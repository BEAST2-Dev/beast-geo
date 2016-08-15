package sphericalGeo.util;

import java.io.File;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.NexusParser;
import sphericalGeo.region.KMLRegion;
import beast.core.Runnable;

@Description("calculates proportion of root locations from posterior tree set fitting in a region")
public class RootDensityFitCalculator extends Runnable  {
	public Input<File> treesetInput = new Input<>("trees","file containing tree set annotated with locations", new File("data.trees"));
	public Input<String> tagInput = new Input<>("tag","tag used in annotated of locations", "location");
	public Input<Integer> burninInput = new Input<>("burnin","burn in percentage, default 10", 10);
	public Input<File> kmlInput = new Input<>("kml","region as stored in KML file to be fitted", Validate.REQUIRED);

	Pattern pattern = Pattern.compile("([0-9\\.Ee-]+),.*=([0-9\\.Ee-]+)");
	Pattern pattern2 = Pattern.compile(".*location=\\{([0-9\\.Ee-]+),([0-9\\.Ee-]+)\\}.*");

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		KMLRegion region = new KMLRegion(kmlInput.get().getPath());
		
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
		Log.warning.print("Reading trees...");
		NexusParser parser = new NexusParser();
		parser.parseFile(treesetInput.get());
		Log.warning.println("Done");
		List<Tree> trees = parser.trees;
		double propFit = 0, unFit = 0;
		for (int i = trees.size() * burnIn / 100; i < trees.size(); i++) {
			Node root = trees.get(i).getRoot();
			String location = (String) root.metaDataString;
			double [] start = parseLoction(location);
			if (region.isInside(start[0], start[1])) {
				propFit++;
			} else {
				unFit++;
			}
		}
		propFit /= (trees.size() - trees.size() * burnIn / 100);
		unFit /= (trees.size() - trees.size() * burnIn / 100);
		Log.info.println("Root location fits " + propFit + "% of the time " + unFit + "% misfit");
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
