package sphericalGeo.util;

import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.List;

import sphericalGeo.GreatCircleDistance;
import sphericalGeo.SphericalDiffusionModel;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.NexusParser;

@Description("Takes a summary tree and insert single child nodes with locations following "
		+ "the shortes path on a sphere.")
public class CirclePathMaker extends BEASTObject {
	// constants for defaults
	final static double MIN_ARC = 5.0;
	final static String LOCATION_NAME = "location";
	
	// members
	PrintStream out = System.out;
	String treeFile = null;
	List<Tree> trees;

	double minarc = MIN_ARC;
	String locationName = LOCATION_NAME;

	public CirclePathMaker(String[] args) {
		try {
			parseArgs(args);
			NexusParser parser = new NexusParser();
			if (treeFile != null) {
				parser.parseFile(new File(treeFile));
			} else {
				parser.parseFile("x", new InputStreamReader(System.in));
			}
			trees = parser.trees;
			if (trees.size() == 0) {
				throw new Exception("no trees found in file");
			}
		} catch (Exception e) {
			e.printStackTrace();
			printUsageAndExit();
		}
	}

	private void parseArgs(String[] args) throws Exception {
		for (int i = 0; i < args.length; i++) {
			String s = args[i];
			switch (s) {
			case "-tree":
				if (i == args.length - 1) {
					throw new Exception("not enough arguments");
				}
				treeFile = args[i+1];
				i++;
			break;
			case "-out":
				if (i == args.length - 1) {
					throw new Exception("not enough arguments");
				}
				out = new PrintStream(new File(args[i+1]));
				i++;
				break;
			case "-minarc":
				if (i == args.length - 1) {
					throw new Exception("not enough arguments");
				}
				minarc = Double.parseDouble(args[i+1]);
				i++;
				break;
			case "-name":
				if (i == args.length - 1) {
					throw new Exception("not enough arguments");
				}
				locationName = args[i+1];
				i++;
				break;
			default:
				throw new Exception("unexpected argument");
			}				
		}
		
	}

	private void printUsageAndExit() {
		System.out.println("CirclePathMaker [args]");
		System.out.println(getDescription()); 
		System.out.println("-tree <nexus file> nexus file containing summary tree to be annotated (default stdin)");
		System.out.println("-out <tree file> where to save results (default to stdout)");
		System.out.println("-minarc <number> minimum arc in degrees to be split (default " + MIN_ARC + ")");
		System.out.println("-name <name> name of the attribute storing the location info (default " + LOCATION_NAME + ")");
		System.exit(0);
	}

	@Override
	public void initAndValidate() {
	}

	private void process() throws Exception {
		for (Tree tree : trees) {
			annotateTree(tree.getRoot());
		}
		trees.get(0).init(out);
		out.println();
		int nSample = 0;
		for (Tree tree : trees) {
			out.print("tree STATE" + nSample++ + " = ");
			printNewick(tree.getRoot()); 
			out.println(";");
		}
		trees.get(0).close(out);
	}
	
    void printNewick(Node node) {
        if (node.getLeft() != null) {
            out.print("(");
            printNewick(node.getLeft());
            if (node.getRight() != null) {
                out.print(',');
                printNewick(node.getRight());
            }
            out.print(")");
        } else {
            out.print(node.getNr()+ 1);
        }
        out.print(node.getNewickMetaData());
        out.print(":");
        out.print(node.getLength());
    }

	
	
	private void annotateTree(Node node) {
		if (node == null) {
			return;
		}
		for (Node child : node.getChildren()) {
			annotateTree(child);
		}
		for (Node child : node.getChildren()) {
			annotateBranch(node, child);
		}
		return;
	}
	
	private void annotateBranch(Node from, Node to) {
		double [] pos1 = nodePosition(from, locationName);
		double [] pos2 = nodePosition(to, locationName);
		if (pos1 == null) {
			Log.err.println("Could not find position in metadata -- perhaps the name need to be specified (with the -name argument)?");
			System.exit(0);
		}
		double angle = GreatCircleDistance.pairwiseDistance(pos1, pos2) / GreatCircleDistance.EARTHRADIUS;
		angle = Math.abs(angle) * 180 / Math.PI;
		if (angle < minarc) {
			return;
		}
		
		double [] cart1 = SphericalDiffusionModel.spherical2Cartesian(pos1[0], pos1[1]);
		double [] cart2 = SphericalDiffusionModel.spherical2Cartesian(pos2[0], pos2[1]);
		// find mid-point
		double [] mid = new double[3];
		mid[0] = (cart1[0] + cart2[0]) / 2.0;
		mid[1] = (cart1[1] + cart2[1]) / 2.0;
		mid[2] = (cart1[2] + cart2[2]) / 2.0;

		double len = Math.sqrt(mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]);
		mid[0] /= len;
		mid[1] /= len;
		mid[2] /= len;
		
		double [] mid0 = SphericalDiffusionModel.cartesian2Sperical(mid, true);
		
		Node node = new Node();
		node.metaDataString = locationName + "={" + mid0[0] + "," + mid0[1] + "}"; 
		node.setHeight((from.getHeight() + to.getHeight()) / 2.0);
		to.setParent(node);
		if (from.getLeft() == to) {
			from.setLeft(node);
		} else {
			from.setRight(node);
		}
		node.setLeft(to);
		node.setParent(from);
		annotateBranch(from, node);
		annotateBranch(node, to);
	}

	private double [] nodePosition(Node node, String sLocationName) {
		if (node == null) {
			return null;
		}
		if (node.metaDataString != null && 
				((node.metaDataString.matches(".*x=[0-9\\.-]+.*") && node.metaDataString.matches(".*y=[0-9\\.-]+.*")) ||
				 (node.metaDataString.matches(".*" + sLocationName + "1=[0-9\\.-]+.*") && node.metaDataString.matches(".*" + sLocationName + "2=[0-9\\.-]+.*")
			    ||(node.metaDataString.matches(".*xy1_median=[0-9\\.-]+.*") && node.metaDataString.matches(".*xy2_median=[0-9\\.-]+.*"))
			    ||(node.metaDataString.matches(".*xy1=[0-9\\.-]+.*") && node.metaDataString.matches(".*xy2=[0-9\\.-]+.*"))
				||(node.metaDataString.matches(".*" + sLocationName + "=\\{.*")))
				||(node.metaDataString.matches(".*xy=\\{.*"))
						)) {
			String [] metaDatas = node.metaDataString.split(",");
			double x = Double.NaN;
			double y = Double.NaN;
			int k = 0;
			for (String str : metaDatas) {
				if (str.startsWith("&")) {
					str = str.substring(1);
				}
				// x, y for landscape aware phylogeography
				if (str.startsWith("x=")) {
					x = (int) (Double.parseDouble(str.substring(2)) + 0.5);
				}
				if (str.startsWith("xy1=") && Double.isNaN(x)) {
					x = (int) (Double.parseDouble(str.substring(4)) + 0.5);
				}
				if (str.startsWith("xy1_median=")) {
					x = (int) (Double.parseDouble(str.substring(11)) + 0.5);
				}
				if (str.startsWith("y=")) {
					y = (int) (Double.parseDouble(str.substring(2)) + 0.5);
				}
				if (str.startsWith("xy2=") && Double.isNaN(y)) {
					y = (int) (Double.parseDouble(str.substring(4)) + 0.5);
				}
				if (str.startsWith("xy2_median=")) {
					y = (int) (Double.parseDouble(str.substring(11)) + 0.5);
				}
				// location.geo1, location.geo2 for continuous phylogeography
				if (str.startsWith(sLocationName + "1=")) {
					y =  Double.parseDouble(str.substring(sLocationName.length() + 2));
				}
				if (str.startsWith(sLocationName + "2=")) {
					x =  Double.parseDouble(str.substring(sLocationName.length() + 2));
				}
				if (str.startsWith(sLocationName + "={")) {
					y = Double.parseDouble(str.substring(sLocationName.length() + 2));
					str = metaDatas[k+1];
					str = str.substring(0, str.indexOf('}'));
					x = Double.parseDouble(str);
				}
				if (str.startsWith("xy={")) {
					x = Double.parseDouble(str.substring(4));
					str = metaDatas[k+1];
					str = str.substring(0, str.indexOf('}'));
					y = Double.parseDouble(str);
				}
				k++;
			}
			return new double[]{y, x};
		}
		return null;
	}	
	
	public static void main(String[] args) throws Exception {
		
		CirclePathMaker app = new CirclePathMaker(args);
		app.process();
	}



}
