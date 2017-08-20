package sphericalGeo.util;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import sphericalGeo.SphericalDiffusionModel;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;
import sphericalGeo.util.treeset.TreeSet;

@Description("Merges two files over independent tree sets by finding most probable location where tree2 fits "
		+ "with tree1. Output written in NEWICK (not NEXUS).")
public class GeoTreeMerger extends SpeedAnnotator {
	public Input<File> treeset2Input = new Input<>("tree2","second file containing tree set annotated with locations");
	public Input<Double> precisionInput = new Input<>("precision","precision use for spherical diffusion process", 100.0);

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		tag = tagInput.get();
		pattern = Pattern.compile("([0-9\\.Ee-]+),.*=([0-9\\.Ee-]+)");
		pattern2 = Pattern.compile(".*" + tag + "=\\{([0-9\\.Ee-]+),([0-9\\.Ee-]+)\\}.*");

		// get trees
		TreeSet treeset = new MemoryFriendlyTreeSet(treesetInput.get().getPath(), 0);
		treeset.reset();

		TreeSet treeset2 = new MemoryFriendlyTreeSet(treeset2Input.get().getPath(), 0);
		treeset2.reset();

		if (treeset.size() != treeset2.size()) {
			throw new IllegalArgumentException("Tree set have different sizes: " + 
					treeset.size() + "!=" + treeset2.size() +". Cannot merge tree sets, giving up!");
		}
		
		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
		}
		int k = 1;
		SphericalDiffusionModel model = new SphericalDiffusionModel(precisionInput.get());
		
		while (treeset.hasNext()) {
			Tree tree = treeset.next();
			Tree tree2 = treeset2.next();
			String rootLocation = tree2.getRoot().metaDataString;
			double [] root = parseLoction(rootLocation);
			double rootAge = tree2.getRoot().getHeight();
			if (rootAge > tree.getRoot().getHeight()) {
				throw new IllegalArgumentException("Cannot merge tree2 into tree1 since tree1 is younger than tree2. "
						+ "Perhaps swap tree1 and tree2?");
			}
			
			List<Node> nodes = new ArrayList<>();
			List<Double> densities = new ArrayList<>();
			List<Double> heights = new ArrayList<>();
			List<Double> lats = new ArrayList<>();
			List<Double> longs = new ArrayList<>();
			
			for (Node node : tree.getNodesAsArray()) {
				if (!node.isRoot()) {
					double tn = node.getHeight();
					double tp = node.getParent().getHeight();

					double [] nodeLocation = parseLoction(node.metaDataString);
					double [] parentLocation = parseLoction(node.getParent().metaDataString);

					double [] end = new double[2];
					for (int i = 0; i < 10; i++) {
						double t = tp + i * (tp - tn) / 10;
						if (t > rootAge) {
							end[0] = nodeLocation[0] + i * (parentLocation[0] - nodeLocation[0]) / 10;
							end[1] = nodeLocation[1] + i * (parentLocation[1] - nodeLocation[1]) / 10;
							double d = model.getLogLikelihood(null, root, end, t);
							nodes.add(node);
							densities.add(d);
							heights.add(t);
							lats.add(end[0]);
							longs.add(end[1]);
						}
					}
				}					
			}
			// normalise densities;
			double sum = 0;
			for (double d : densities) {
				sum += d;
			}
			double r = Randomizer.nextDouble() * sum;
			int i = 0;
			while (r > 0) {
				r -= densities.get(i++);
			}
			Node newNode = new Node();
			newNode.metaDataString= tag+"={" + lats.get(i)+","+longs.get(i)+"}";
			newNode.setHeight(heights.get(i));
			newNode.setLeft(tree2.getRoot());
			newNode.setParent(nodes.get(i).getParent());
			newNode.setRight(nodes.get(i));

			
            out.print("tree TREE" + k++ + " = ");
            String newick = tree.getRoot().toNewick();
            out.print(newick);
            out.println(";");
		}
		if (outputInput.get() != null) {
			out.close();
		}
		Log.warning.println("Done.");
	}

	
	public static void main(String[] args) throws Exception {
		new Application(new GeoTreeMerger(), "Geo Tree Merger", args);
	}
}
