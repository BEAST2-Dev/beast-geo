package sphericalGeo.util;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beastfx.app.tools.Application;
import beastfx.app.util.XMLFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.MRCAPrior;
import sphericalGeo.region.KMLRegion;
import sphericalGeo.util.HeatMapMaker.Dot;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;
import sphericalGeo.util.treeset.TreeSet;
import beast.base.inference.Runnable;

@Description("calculates proportion of root locations from posterior tree set fitting in a region")
public class RootDensityFitCalculator extends Runnable  {
	public Input<File> treesetInput = new Input<>("trees","file containing tree set annotated with locations", new File("data.trees"));
	public Input<String> tagInput = new Input<>("tag","tag used in annotated of locations", "location");
	public Input<Integer> burninInput = new Input<>("burnin","burn in percentage, default 10", 10);
	public Input<List<String>> kmlInput = new Input<>("kml","region as stored in KML file to be fitted", new ArrayList<>(), Validate.REQUIRED);
	public Input<Double> lowerAgeInput = new Input<>("lower","lower bound of root age to be taken in account", Double.NEGATIVE_INFINITY);
	public Input<Double> upperAgeInput = new Input<>("upper","upper bound of root age to be taken in account", Double.POSITIVE_INFINITY);

	public Input<XMLFile> xmlInput = new Input<>("clade", "Name of XML file containing a single TaxonSet in BEAST XML format. "
			+ "If specified, the location of the MRCA of a clade is used (and the rootonly flag is ignored).");

	Pattern pattern = Pattern.compile("([0-9\\.Ee-]+),.*=([0-9\\.Ee-]+)");
	Pattern pattern2 = Pattern.compile(".*location=\\{([0-9\\.Ee-]+),([0-9\\.Ee-]+)\\}.*");

	TaxonSet taxonset = null;
	
	
	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		if (xmlInput.get() != null && !xmlInput.get().getName().equals("[[none]]")) {
			taxonset = HeatMapMaker.getTaxonSet(xmlInput.get());
		}

		for (String kmlFile : kmlInput.get()) {
		
		KMLRegion region = new KMLRegion(kmlFile);//.getPath());
		double lowerAge = lowerAgeInput.get();
		double upperAge = upperAgeInput.get();
		
		String tag = tagInput.get();
		pattern = Pattern.compile("([0-9\\.Ee-]+),.*=([0-9\\.Ee-]+)");
		pattern2 = Pattern.compile(".*" + tag + "=\\{([0-9\\.Ee-]+),([0-9\\.Ee-]+)\\}.*");
		int burnIn = burninInput.get();
		if (burnIn < 0) {
			Log.warning.println("Setting burn in to 0");
			burnIn = 0;
		}
		if (burnIn >= 100) {
			Log.err.println("Burnin is a percentage, and must be below 100. Quiting now");
			return;
		}

		// get trees
//		Log.warning.print("Reading trees...");
//		NexusParser parser = new NexusParser();
//		parser.parseFile(treesetInput.get());
//		Log.warning.println("Done");
//		List<Tree> trees = parser.trees;
		
		TreeSet treeSet = new MemoryFriendlyTreeSet(treesetInput.get().getPath(), burnIn);		
    	treeSet.reset();
		double propFit = 0, unFit = 0;
		while (treeSet.hasNext()) {
			Tree tree = treeSet.next();
			Node root = tree.getRoot();
			if (taxonset != null) {
				MRCAPrior p = new MRCAPrior();
				p.initByName("taxonset", taxonset, "tree", tree);
				root = p.getCommonAncestor();
			}

			
			if (root.getHeight() > lowerAge && root.getHeight() <= upperAge) {
				String location = root.metaDataString;
				double [] start = parseLoction(location);
				if (region.isInside(start[0], start[1])) {
					propFit++;
				} else {
					unFit++;
				}
			} else {
				unFit++;
			}
		}
		propFit /= treeSet.size();
		unFit /= treeSet.size();
		Log.info.println("Root location fits " + 100*propFit + "% of the time " + 100*unFit + "% misfit");
		}
	}
	
	private double[] parseLoction(String location) {
		Matcher matcher = pattern2.matcher(location);
		String sMatch0;
		String sMatch1;
		if (matcher.find()) {
			// it is in x=12,y=34 format
			//int nGroups = matcher.groupCount();
			sMatch0 = matcher.group(1);
			sMatch1 = matcher.group(2);
		} else {
			// it is in location={12,34} format
			matcher = pattern.matcher(location);
			if (matcher.find()) {
				// it is in x=12,y=34 format
				sMatch0 = matcher.group(1);
				sMatch1 = matcher.group(2);
			} else {
				throw new RuntimeException("no match");
			}
		}
		double [] coords = new double[2];
		coords[0] = Double.parseDouble(sMatch0);
		coords[1] = Double.parseDouble(sMatch1);
		return coords;
	}

    

	public static void main(String[] args) throws Exception {
		new Application(new RootDensityFitCalculator(), "RootAreaDensityCalculator", args);
	}
}
