package sphericalGeo.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import sphericalGeo.GreatCircleDistance;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;
import sphericalGeo.util.treeset.TreeSet;

@Description("Annotates a set of trees with their speed on a branch -- requires geographical info in the input tree set")
public class SpeedAnnotator extends Runnable  {
	public Input<File> treesetInput = new Input<>("trees","file containing tree set annotated with locations", new File("data.trees"));
	public Input<String> tagInput = new Input<>("tag","tag used in annotated of locations", "location");
	public Input<OutFile> outputInput = new Input<>("output","where to save the file", new OutFile("/tmp/x.trees"));
	public Input<File> leafLocationsInput = new Input<>("leafLocations","file with tip locations (optional). "
			+ "Tab delimited file: first column is taxon name, second the latitude, and third the longitude");

	String tag;
	double interval;
	Map<String, Double[]> tipLocationMap;
	
	Pattern pattern = Pattern.compile("([0-9\\.Ee-]+),.*=([0-9\\.Ee-]+)");
	Pattern pattern2 = Pattern.compile(".*location=\\{([0-9\\.Ee-]+),([0-9\\.Ee-]+)\\}.*");

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {

		tag = tagInput.get();
		pattern = Pattern.compile("([0-9\\.Ee-]+),.*=([0-9\\.Ee-]+)");
		pattern2 = Pattern.compile(".*" + tag + "=\\{([0-9\\.Ee-]+),([0-9\\.Ee-]+)\\}.*");

		if (leafLocationsInput.get() != null && !leafLocationsInput.get().getName().equals("[[none]]")) {
			tipLocationMap = new HashMap<>();
			BufferedReader fin = new BufferedReader(new FileReader(leafLocationsInput.get()));
			while (fin.ready()) {
				String str = fin.readLine();
				String [] strs = str.split("\t");
				if (strs.length >= 3) {
					try {
						double lat = Double.parseDouble(strs[1]);
						double lon = Double.parseDouble(strs[2]);
						tipLocationMap.put(strs[0], new Double[]{lat, lon});
					} catch (NumberFormatException e) {
						// ignore
					}
				}
			}
			fin.close();
		}
		
		// get trees
		TreeSet treeset = new MemoryFriendlyTreeSet(treesetInput.get().getPath(), 0);
		treeset.reset();

		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
		}
		Tree tree = treeset.next();
		tree.init(out);
		treeset.reset();
		out.println();
		int k = 1;
		while (treeset.hasNext()) {
			tree = treeset.next();
			for (Node node : tree.getNodesAsArray()) {
				if (!node.isRoot()) {
					String parentLocation = node.getParent().metaDataString;
					String location = node.metaDataString;
					double [] start = parseLoction(parentLocation);
					double [] end;
					if (node.isLeaf() && node.metaDataString == null) {
						end = getLocation(node);
					} else {
						end = parseLoction(location);
					}
					double distance = GreatCircleDistance.pairwiseDistance(start, end);
					double time = (node.getParent().getHeight() - node.getHeight());
					double speed = distance / time;
					if (distance == 0) {
						speed = 0;
					}
					node.metaDataString = "speed="+speed;				
				} else {
					node.metaDataString = "speed=0.0";				
				}
			}
            out.print("tree TREE" + k++ + " = ");
            int[] dummy = new int[1];
            String newick = tree.getRoot().toSortedNewick(dummy, true);
            out.print(newick);
            out.println(";");
		}
		tree.close(out);
		if (outputInput.get() != null) {
			out.close();
		}
		Log.warning.println("Done.");
	}
	
	protected double[] parseLoction(String location) {
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

	protected double[] getLocation(Node node) {
		if (tipLocationMap != null && tipLocationMap.containsKey(node.getID())) {
			Double [] location = tipLocationMap.get(node.getID());
			double [] l2 = new double[2];
			l2[0] = location[0];
			l2[1] = location[1];
			return l2;
		}
		throw new IllegalArgumentException("could not find location for node " + node.getID() + ". "
				+ "Perhaps leafLocations is not specified, or the node name is misspelled (note "
				+ "names are case sensitive)");
	}

	public static void main(String[] args) throws Exception {
		new Application(new SpeedAnnotator(), "Speed Annotator", args);
	}

}
