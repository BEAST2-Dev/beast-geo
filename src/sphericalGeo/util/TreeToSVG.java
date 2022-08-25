package sphericalGeo.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.MCMC;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.util.Randomizer;
import beast.base.parser.XMLParser;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;
import sphericalGeo.util.treeset.TreeSet;

@Description("Create SVG file for a summary tree -- requires geographical info in the input tree set")
public class TreeToSVG extends SpeedAnnotator  {
	public Input<File> xmlFileInput = new Input<>("xml", "xml file with MRCAPrior constraints. Any branch inside such constraint will be omitted");
	public Input<Boolean> suppressCladesInput = new Input<>("suppressClades", "do not export branches in clades that have MRCAPrior", false);
	public Input<Double> legendOffsetXInput = new Input<>("legendOffsetX", "offset of legend (specified by MRCAs, if any)", 0.0);
	public Input<Double> legendOffsetYInput = new Input<>("legendOffsetY", "offset of legend (specified by MRCAs, if any)", 0.0);
	
	Node mrca;
	String [] cssclass;
	
	public TreeToSVG() {
		outputInput.setValue(new OutFile("/tmp/x.svg"), this);
	}

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		boolean suppressClades = suppressCladesInput.get();
		
		List<Set<String>> mrcaPriors = new ArrayList<>();
		List<String> classIDs = new ArrayList<>();
		Set<String> classes = new HashSet<>();
		if (xmlFileInput.get() != null && xmlFileInput.get().exists()) {
			XMLParser parser = new XMLParser();
			MCMC mcmc = (MCMC) parser.parseFile(xmlFileInput.get());
			State state = mcmc.startStateInput.get();
			for (StateNode n : state.stateNodeInput.get()) {
				if (n instanceof Tree) {
					for (BEASTInterface o : n.getOutputs()) {
						if (o instanceof MRCAPrior) {
							Set<String> taxa = ((MRCAPrior) o).taxonsetInput.get().getTaxaNames();
							mrcaPriors.add(taxa);
							classIDs.add(o.getID());
							classes.add(o.getID());
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
"}\n" +
".mrca {\n" +
"stroke-width:0.1;\n" +
"}\n" +
".legend {\n" +
"font-size:1pt;\n" +
"font-family:helvetica;\n" +
"}\n" +
".curly {\n" +
"stroke-width:0.2;\n" +
"stroke:red;\n" +
"fill:none;\n" +
"}\n" +
".curly2 {\n" +
"stroke-width:0.25;\n" +
"stroke:black;\n" +
"fill:none;\n" +
"}\n" +
".label{\n"+
"font-size:0.1;\n" +
"font-family:helvetica;\n" +
"}\n"
);
		for (String c : classes) {
			// Will produce only bright / light colours:
			int r = Randomizer.nextInt(128) + 128;
			int g = Randomizer.nextInt(128) + 128;
			int b = Randomizer.nextInt(128) + 128;
			String color= Integer.toHexString(r) + Integer.toHexString(g) + Integer.toHexString(b);
			
			out.println("." + c.replaceAll("\\..*", "") + "{\nstroke:#" + color + ";\n}");
		}
	
		out.println("</style>\n" +
"<g>\n");
		
		Tree tree = treeset.next();
		treeset.reset();
		out.println();
		while (treeset.hasNext()) {
			tree = treeset.next();
		
			
			boolean [] suppress = new boolean[tree.getNodeCount()];
			if (mrcaPriors != null) {
				cssclass = new String[tree.getNodeCount()];
				Arrays.fill(cssclass, "");
				for (int i = 0; i < mrcaPriors.size(); i++) {
					System.err.println(i);
					Set<String> isInTaxaSet = mrcaPriors.get(i);
					calcMRCAtime(tree.getRoot(), new int[1], isInTaxaSet, isInTaxaSet.size(), classIDs.get(i).replaceAll("\\..*", ""));
					suppress[mrca.getNr()] = true;
				}
				mark(tree.getRoot(), suppress, false);
			}
			
			
			
			Node [] nodes = tree.getNodesAsArray();
			Arrays.sort(nodes, (n1, n2) -> {return Double.compare(n1.getHeight(), n2.getHeight());});
			for (Node node : nodes) {
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
					if (!suppress[node.getNr()] || !suppress[node.getParent().getNr()] || !suppressClades) {
						out.print("<g" + (node.isLeaf() ? " id=\"" + node.getID() + "\"": "") + ">");
						String mrca = "";
						if (cssclass != null && !cssclass[node.getParent().getNr()].equals("")) {
							mrca = " mrca";
						}
						out.println("<line x1=\""+ (start[1] + 180) +"\" "+
								"y1=\""+(90-start[0]) +"\" " +
								"x2=\""+(end[1]+180) +"\" " +
								"y2=\""+(90-end[0])+"\" class=\"pnline " + (cssclass == null ? "" : cssclass[node.getParent().getNr()] + "2" + mrca) + "\"/>\n");
						out.print("<line x1=\""+ (start[1] + 180) +"\" "+
								"y1=\""+(90-start[0]) +"\" " +
								"x2=\""+(end[1]+180) +"\" " +
								"y2=\""+(90-end[0])+"\" class=\"pnline2 " + (cssclass == null ? "" : cssclass[node.getParent().getNr()] + mrca) + "\"/>\n");
						
						if (!suppress[node.getNr()] || !suppress[node.getParent().getNr()]) {
							out.print("<path d=\"M" + (start[1] + 180)+ " " + (90-start[0]) + 
									" c " +(end[1]-start[1])/2.0 + " " + (start[0]-end[0])/2.0 +
									" , "+ (end[1]-start[1])/2.0 + " " + (start[0]-end[0])/2.0 +
									" ," + (end[1]-start[1]) + " " + (start[0]-end[0])+" \" class=\"curly2 \"/>\n"); 
							out.print("<path d=\"M" + (start[1] + 180)+ " " + (90-start[0]) + 
									" c " +(end[1]-start[1])/2.0 + " " + (start[0]-end[0])/2.0 +
									" , "+ (end[1]-start[1])/2.0 + " " + (start[0]-end[0])/2.0 +
									" ," + (end[1]-start[1]) + " " + (start[0]-end[0])+" \" class=\"curly \"/>\n"); 
						}
						if (node.isLeaf()) {
//							out.print("\n<text x=\"" + (end[1] + 180) + "\" y=\"" + (90 - end[0]) + "\" class=\"label\">" + node.getID() + "</text>");
						}
						out.println("</g>");
					}
				} else {
					String location = node.metaDataString;
					double [] start = parseLoction(location);
					out.println("<circle cx=\"" + (start[1]+180) +"\" cy=\"" + (90-start[0]) + "\" r=\"0.5\" stroke=\"#000000\" stroke-width=\"0.05\" fill=\"#ffff00\" />\n");
				}
			}
			out.print("<g>\n");
			for (Node node : nodes) {
				if (!node.isRoot()) {
					boolean isMRCA = false;
					if (cssclass != null && !cssclass[node.getNr()].equals("")) {
						isMRCA = true;;
					}
					boolean isParentMRCA = false;
					if (cssclass != null && !cssclass[node.getParent().getNr()].equals("")) {
						isParentMRCA = true;;
					}
	
					if (isMRCA && !isParentMRCA) {
						String location = node.metaDataString;
						double [] start = parseLoction(location);
						double x = start[1] + 180;
						double y = 90 - start[0]; 
						out.print("<circle cx=\"" + x +"\" cy=\"" + y + "\" r=\"0.5\" class=\"clabel\"/>");
						out.print("<text id=\""+cssclass[node.getNr()]+"\" x=\"" + x +"\" y=\"" + (y + 0.3) +"\" class=\"label\">1</text>\n");
					}
				}

			}
			out.print("</g>\n");
		}
		out.println("</g>\n");
		int k = 1;
		String [] classes2 = classes.toArray(new String[]{});
		Arrays.sort(classes2);
		double offsetX = legendOffsetXInput.get();
		double offsetY = legendOffsetYInput.get();
		double dy = 1;
		double len = 5;
		for (String s : classes2) {
			s = s.replaceAll("\\..*", "");
			out.print("<g class=\"legend\">");
			out.print("<line x1=\"" + offsetX + "\" y1=\"" + (k * dy + offsetY) + "\" x2=\""+( + offsetX + len)+"\" y2=\"" + (k * dy + offsetY) + "\" class=\"" + s +"\"/>");
			out.print("<text x=\""+( + offsetX + len + 1)+"\" y=\"" + (k * dy + offsetY + dy/2) + "\">"+ s + "</text>");
			out.println("</g>");
			k++;
		}
		out.println("</svg>\n");
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
	
	
    int calcMRCAtime(final Node node, final int[] taxonCount2, Set<String> isInTaxaSet, int nrOfTaxa, String css) {
        if (node.isLeaf()) {
            taxonCount2[0]++;
            if (isInTaxaSet.contains(node.getID())) {
            	cssclass[node.getNr()] = css;
                return 1;
            } else {
                return 0;
            }
        } else {
            int taxonCount = calcMRCAtime(node.getLeft(), taxonCount2, isInTaxaSet, nrOfTaxa, css);
            final int leftTaxa = taxonCount2[0];
            taxonCount2[0] = 0;
            if (node.getRight() != null) {
                taxonCount += calcMRCAtime(node.getRight(), taxonCount2, isInTaxaSet, nrOfTaxa, css);
                final int rightTaxa = taxonCount2[0];
                taxonCount2[0] = leftTaxa + rightTaxa;
                if (taxonCount > 0 && taxonCount <= nrOfTaxa) {
                	cssclass[node.getNr()] = css;
                }
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
