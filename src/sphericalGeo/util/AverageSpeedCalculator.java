package sphericalGeo.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import sphericalGeo.GreatCircleDistance;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;
import sphericalGeo.util.treeset.TreeSet;

@Description("calculates average speed on a posterior set of trees from a phylogeographical analysis")
public class AverageSpeedCalculator extends Runnable  {
	public Input<File> treesetInput = new Input<>("trees","file containing tree set annotated with locations", new File("data.trees"));
	public Input<String> tagInput = new Input<>("tag","tag used in annotated of locations", "location");
	public Input<Integer> burninInput = new Input<>("burnin","burn in percentage, default 10", 10);
	public Input<Double> intervalInput = new Input<>("interval","interval size used to measure through time, default 5", 5.0);
	public Input<File> leafLocationsInput = new Input<>("leafLocations","file with tip locations (optional). "
			+ "Tab delimited file: first column is taxon name, second the latitude, and third the longitude");

	double sumOfTime = 0;
	double sumOfDistance = 0;
	double [] sumOfTimes;
	double [] sumOfDistances;
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
		interval = intervalInput.get();

		tag = tagInput.get();
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

		if (leafLocationsInput.get() != null) {
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
		TreeSet treeset = new MemoryFriendlyTreeSet(treesetInput.get().getPath(), burnIn);
		treeset.reset();
		
		double sum = 0;
		
		int k = 0;
		while (treeset.hasNext()) {
			Tree tree = treeset.next();
			sumOfTime = 0;
			sumOfDistance = 0;
			calcSpeed(tree.getRoot());
			sum += sumOfDistance / sumOfTime;
			k++;
			if (k % 100 == 0) {Log.warning.print('.');}
			//System.err.println(sumOfDistance / sumOfTime);
		}
		Log.warning.println();
		sum /= k;
		Log.info.println("Average displacement is " + sum + " km per unit of time in the trees");

		double maxHeight = 0;
		treeset.reset();
		while (treeset.hasNext()) {
			Tree tree = treeset.next();
			maxHeight = Math.max(maxHeight, tree.getRoot().getHeight());
		}
		int n = 1 + (int)(maxHeight/intervalInput.get());
		sumOfTimes = new double[n];
		sumOfDistances = new double[n];
		treeset.reset();
		k = 0;
		while (treeset.hasNext()) {
			Tree tree = treeset.next();
			for (Node node : tree.getNodesAsArray()) {
				if (!node.isRoot()) {
					String parentLocation = (String) node.getParent().metaDataString;
					String location = (String) node.metaDataString;
					double [] start = parseLoction(parentLocation);
					double [] end;
					if (node.isLeaf() && node.metaDataString == null) {
						end = getLocation(node);
					} else {
						end = parseLoction(location);
					}
					double distance = GreatCircleDistance.pairwiseDistance(start, end);
					distribute(distance, node.getHeight(), node.getParent().getHeight());
				}
			}
			k++;
			if (k % 100 == 0) {Log.warning.print('.');}
		}
		Log.warning.println();
		for (int i = 0; i < n; i++) {
			System.out.println(i*interval + " -- " + (i+1) * interval + ": " + sumOfDistances[i]/sumOfTimes[i] + " km/unit of time");
		}
	}

	private void distribute(double distance, double height, double height2) {
		int i = (int)(height / interval);
		int j = (int)(height2 / interval);
		while (i < j) {
			double t = (i+1) * interval;
			sumOfTimes[i] += t - height;
			double delta = distance * (t - height) / (height2 - height);
			sumOfDistances[i] += delta;
			distance -= delta;
			height = t;
			i++;
		}
		sumOfTimes[i] += height2 - height;
		sumOfDistances[i] += distance;
	}

	private double calcSpeed(Node node) {
		if (!node.isRoot()) {
			sumOfTime += node.getLength();
			String parentLocation = (String) node.getParent().metaDataString;
			String location = (String) node.metaDataString;
			double [] start = parseLoction(parentLocation);
			double [] end;
			if (node.isLeaf() && node.metaDataString == null) {
				end = getLocation(node);
			} else {
				end = parseLoction(location);
			}
			double distance = GreatCircleDistance.pairwiseDistance(start, end);
			sumOfDistance += distance;
		}
		for (Node child : node.getChildren()) {
			calcSpeed(child);
		}
		return 0;
	}

	private double[] getLocation(Node node) {
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
		AverageSpeedCalculator averageSpeedCalculator = new AverageSpeedCalculator();
	
		if (args.length == 0) {
			// create BeautiDoc and beauti configuration
			BeautiDoc doc = new BeautiDoc();
			doc.beautiConfig = new BeautiConfig();
			doc.beautiConfig.initAndValidate();
		
			// create panel with entries for the application
			BEASTObjectPanel panel = new BEASTObjectPanel(averageSpeedCalculator, averageSpeedCalculator.getClass(), doc);
			
			// wrap panel in a dialog
			BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
	
			// show the dialog
			if (dialog.showDialog()) {
				dialog.accept(averageSpeedCalculator, doc);
				// create a console to show standard error and standard output
				//ConsoleApp app = new ConsoleApp("PathSampler", "Path Sampler: " + sampler.model1Input.get().getPath());
				averageSpeedCalculator.initAndValidate();
				averageSpeedCalculator.run();
			}
			return;
		}

		Application main = new Application(averageSpeedCalculator);
		main.parseArgs(args, false);
		averageSpeedCalculator.initAndValidate();
		averageSpeedCalculator.run();
	}

}
