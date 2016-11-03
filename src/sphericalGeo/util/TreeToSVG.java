package sphericalGeo.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.MRCAPrior;
import beast.util.XMLParser;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;
import sphericalGeo.util.treeset.TreeSet;

@Description("Create SVG file for a summary tree -- requires geographical info in the input tree set")
public class TreeToSVG extends SpeedAnnotator  {
	public Input<File> xmlFileInput = new Input<>("xml", "xml file with MRCAPrior constraints. Any branch inside such constraint will be omitted");
	
	
	Node mrca;
	
	public TreeToSVG() {
		outputInput.setValue(new OutFile("/tmp/x.svg"), this);
	}

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		List<Set<String>> mrcaPriors = new ArrayList<>();
		Tree baseTree = null;
		
		if (xmlFileInput.get() != null && xmlFileInput.get().exists()) {
			XMLParser parser = new XMLParser();
			MCMC mcmc = (MCMC) parser.parseFile(xmlFileInput.get());
			State state = mcmc.startStateInput.get();
			for (StateNode n : state.stateNodeInput.get()) {
				if (n instanceof Tree) {
					baseTree = (Tree) n;
					for (BEASTInterface o : n.getOutputs()) {
						if (o instanceof MRCAPrior) {
							Set<String> taxa = ((MRCAPrior) o).taxonsetInput.get().getTaxaNames();
							mrcaPriors.add(taxa);
						}
					}
				}
			}
		}
		

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
		
		out.println("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n" +
"<svg\n" +
"   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n" +
"   xmlns:cc=\"http://creativecommons.org/ns#\"\n" +
"   xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n" +
"   xmlns:svg=\"http://www.w3.org/2000/svg\"\n" +
"   xmlns=\"http://www.w3.org/2000/svg\"\n" +
"   version=\"1.1\"\n" +
"   width=\"360\"\n" +
"   height=\"180\"\n" +
"   viewBox=\"0 0 360 180\"\n" +
"   id=\"svg2\">\n" +
"<style type=\"text/css\">\n" +
".pnline {\n" +
"stroke-width:0.25;stroke-miterlimit:4;stroke-dasharray:none;stroke:#0000e0;stroke-opacity:1;stroke-linejoin:round;stroke-linecap:round\n"+
"}\n" +
".pnline2 {\n"+
"stroke-width:0.2;stroke-miterlimit:4;stroke-dasharray:none;stroke:#ffff00;stroke-opacity:1;stroke-linejoin:round;stroke-linecap:round\n"+
"}\n"+
"</style>\n" +
"<g>\n");
		
		Tree tree = treeset.next();
		treeset.reset();
		out.println();
		while (treeset.hasNext()) {
			tree = treeset.next();
			
			boolean [] suppress = new boolean[tree.getNodeCount()];
			if (mrcaPriors != null) {
				for(Set<String> isInTaxaSet : mrcaPriors) {
					calcMRCAtime(tree.getRoot(), new int[1], isInTaxaSet, isInTaxaSet.size());
					suppress[mrca.getNr()] = true;
				}
				mark(tree.getRoot(), suppress, false);
			}
			
			
			
			
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
					if (!suppress[node.getNr()] || !suppress[node.getParent().getNr()]) {
						out.println("<line x1=\""+ (start[1] + 180) +"\" "+
								"y1=\""+(90-start[0]) +"\" " +
								"x2=\""+(end[1]+180) +"\" " +
								"y2=\""+(90-end[0])+"\" class=\"pnline\"/>\n");
						out.println("<line x1=\""+ (start[1] + 180) +"\" "+
								"y1=\""+(90-start[0]) +"\" " +
								"x2=\""+(end[1]+180) +"\" " +
								"y2=\""+(90-end[0])+"\" class=\"pnline2\"/>\n");
					}
				} else {
					String location = (String) node.metaDataString;
					double [] start = parseLoction(location);
					out.println("<circle cx=\"" + (start[1]+180) +"\" cy=\"" + (90-start[0]) + "\" r=\"0.5\" stroke=\"#000000\" stroke-width=\"0.05\" fill=\"#ffff00d\" />\n");
				}
			}
		}
		out.println("</g>\n" + 
				"</svg>\n");
		out.println("");
		if (outputInput.get() != null) {
			out.close();
			Log.warning("Output written to " + outputInput.get().getPath());
		}
		Log.warning.println("Done.");
	}
	
	private void mark(Node node, boolean[] suppress, boolean b) {
		suppress[node.getNr()] = b;
		for (Node child : node.getChildren()) {
			b = suppress[node.getNr()] || suppress[child.getNr()];
			mark(child, suppress, b);
		}
	}
	
	
    int calcMRCAtime(final Node node, final int[] taxonCount2, Set<String> isInTaxaSet, int nrOfTaxa) {
        if (node.isLeaf()) {
            taxonCount2[0]++;
            if (isInTaxaSet.contains(node.getID())) {
                return 1;
            } else {
                return 0;
            }
        } else {
            int taxonCount = calcMRCAtime(node.getLeft(), taxonCount2, isInTaxaSet, nrOfTaxa);
            final int leftTaxa = taxonCount2[0];
            taxonCount2[0] = 0;
            if (node.getRight() != null) {
                taxonCount += calcMRCAtime(node.getRight(), taxonCount2, isInTaxaSet, nrOfTaxa);
                final int rightTaxa = taxonCount2[0];
                taxonCount2[0] = leftTaxa + rightTaxa;
                if (taxonCount == nrOfTaxa) {
                	mrca = node;
                    return taxonCount + 1;
                }
            }
            return taxonCount;
        }
    }


	public static void main(String[] args) throws Exception {
		new Application(new TreeToSVG(), "Tree to SVG", args);
	}

}
