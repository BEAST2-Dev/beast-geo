package sphericalGeo.util;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.FunctionEvaluationException;

import sphericalGeo.SphericalDiffusionModel;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

@Description("Simulates spherical diffusion process on a tree")
public class SphericalDiffusionSimulator extends beast.core.Runnable {
    public Input<TreeInterface> treeInput = new Input<TreeInterface>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    public Input<SiteModelInterface> siteModelInput = new Input<SiteModelInterface>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    
    public Input<BranchRateModel.Base> branchRateModelInput = new Input<BranchRateModel.Base>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    public Input<RealParameter> startInput = new Input<RealParameter>("start","start location at the root of the tree");

    public Input<String> fileInput = new Input<String>("file","file where to store the leaf locations (default stdout)");
    public Input<String> treeFileInput = new Input<String>("treefile","file where to store the tree (default stdout)");
    public Input<String> treePNGFileInput = new Input<String>("pngfile","file where to store an image of the tree (default 'tree.png')", "tree.png");
    
    public Input<Boolean> forceAngleInput = new Input<Boolean>("forceAngle","enforce that angle of next point is restricted by its parents direction", true);
	
	SphericalDiffusionModel substModel;
	TreeInterface tree;
	BranchRateModel clockModel;

	double [][] position;
	double precision;
	boolean forceAngle;
	
	@Override
	public void initAndValidate() throws Exception {
		clockModel = branchRateModelInput.get();
		SiteModel siteModel = (SiteModel) siteModelInput.get();
		substModel = (SphericalDiffusionModel) siteModel.substModelInput.get();
		tree = treeInput.get();
		
		position = new double[tree.getNodeCount()][2];
		precision = substModel.precisionInput.get().getValue();
		forceAngle = forceAngleInput.get();
	}

	@Override
	public void run() {
		// starting position
		if (startInput.get() != null) {
			RealParameter start = startInput.get();
			position[tree.getRoot().getNr()] = new double[]{start.getValue(0), start.getValue(1)};
		} else {
			position[tree.getRoot().getNr()] = new double[]{0.0, 0.0};
		}
		
		// propagate along the tree
		if (forceAngle) {
			traverse(tree.getRoot(), new double[]{0,Math.PI*2});
		} else {
			traverse(tree.getRoot());
		}
		
		// output leaf locations
		if (fileInput.get() == null) {
			System.out.println(getLeafLocations());
		} else {
			try {
				FileWriter outfile = new FileWriter(fileInput.get());
				outfile.write(getLeafLocations());
				outfile.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		// output Newick tree
		StringBuilder buf = new StringBuilder();
		printTree(tree.getRoot(), buf);
		if (treeFileInput.get() == null) {
			System.out.println(buf.toString());
		} else {
			try {
				FileWriter outfile = new FileWriter(treeFileInput.get());
				outfile.write(buf.toString());
				outfile.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		// output png tree
		drawTreeToFile(treePNGFileInput.get());
	}
	

	private void printTree(Node node, StringBuilder buf) {
		int i = node.getNr();
		if (node.isLeaf()) {
			buf.append(node.getID());
		} else {
			buf.append("(");
			printTree(node.getLeft(), buf);
			buf.append(",");
			printTree(node.getRight(), buf);
			buf.append(")");
		}
		buf.append("[&location={" + position[i][0] +"," +position[i][1] +"}]");
		buf.append(":");
		buf.append(node.getLength());
		
	}

	String getLeafLocations() {
		StringBuilder buf = new StringBuilder();
		// report leaf locations
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			buf.append(tree.getNode(i).getID() + " = " + position[i][0] + " " + position[i][1] + ",\n");
		}
		// remove trailing comma
		buf.delete(buf.length()-2, buf.length());
		return buf.toString();
	}

	
	private void traverse(Node node, double[] range) {
		if (!node.isRoot()) {
			double t = node.getLength() * clockModel.getRateForBranch(node);
			double [] start = position[node.getParent().getNr()];
			try {
				position[node.getNr()] = substModel.sample(start, t, precision, range);
			} catch (ConvergenceException | FunctionEvaluationException | IllegalArgumentException e) {
				e.printStackTrace();
			}
		}
		if (!node.isLeaf()) {
			traverse(node.getLeft(), new double[]{range[0], (range[0] + range[1])/2.0});
			traverse(node.getRight(), new double[]{(range[0] + range[1])/2.0, range[1]});
		}
	}

	private void traverse(Node node) {
		if (!node.isRoot()) {
			double t = node.getLength() * clockModel.getRateForBranch(node);
			double [] start = position[node.getParent().getNr()];
			try {
				position[node.getNr()] = substModel.sample(start, t, precision);
			} catch (ConvergenceException | FunctionEvaluationException | IllegalArgumentException e) {
				e.printStackTrace();
			}
		}
		if (!node.isLeaf()) {
			traverse(node.getLeft());
			traverse(node.getRight());
		}
	}

	
	void drawTreeToFile(String sFileName) {
		Panel panel = new Panel();
		panel.min0 = -30;
		panel.max0 = 30;
		panel.min1 = -30;
		panel.max1 = 30;
		panel.setSize(1024, 768);
		BufferedImage bi;
		Graphics g;
		bi = new BufferedImage(panel.getWidth(), panel.getHeight(), BufferedImage.TYPE_INT_RGB);
		g = bi.getGraphics();
		g.setPaintMode();
		panel.print(g);
		try {
			ImageIO.write(bi, "png", new File(sFileName));
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	class Panel extends JPanel {
		double min0 = -90, max0 = 90, min1 = -180, max1 = 180;
		
		@Override
		protected void paintComponent(Graphics g) {
			int h = getHeight();
			int w = getWidth();
			g.setColor(Color.WHITE);
			g.fillRect(0,  0, w, h);
			((Graphics2D) g).setStroke(new BasicStroke(2.0f));
			for (Node node : tree.getNodesAsArray()) {
				g.setColor(Color.RED);
				int i = node.getNr();
				if (!node.isRoot()) {
					int y0 = (int)((position[i][0]-min0) * h/(max0-min0));
					int x0 = (int)((position[i][1]-min1) * w/(max1-min1));
					int j = node.getParent().getNr(); 
					int y1 = (int)((position[j][0]-min0) * h/(max0-min0));
					int x1 = (int)((position[j][1]-min1) * w/(max1-min1));
					g.drawLine(x0, y0, x1, y1);
					
				} else {
					g.setColor(Color.BLUE);
				}
				int y = (int)((position[i][0]-min0) * h/(max0-min0));
				int x = (int)((position[i][1]-min1) * w/(max1-min1));
				int width = 10; int height = 10;
				g.fillArc(x-width/2, y-width/2, width, height, 0, 360);
			}
			
			// grid
			g.setColor(Color.GRAY);
			((Graphics2D) g).setStroke(new BasicStroke(2.0f, // line width
				      /* cap style */BasicStroke.CAP_BUTT,
				      /* join style, miter limit */BasicStroke.JOIN_BEVEL, 1.0f,
				      /* the dash pattern */new float[] { 8.0f, 3.0f, 2.0f, 3.0f },
				      /* the dash phase */0.0f));
			for (int i = -360; i < 360; i += 10) {
				int x = (int)((i-min1) * w/(max1-min1));
				g.drawLine(x, 0, x, h);
				g.drawString(i+"", x, 15);
			}
			for (int i = -90; i < 90; i += 10) {
				int y = (int)((i-min0) * h/(max0-min0));
				g.drawLine(0, y, w, y);
				g.drawString(i+"", 5, y);
			}
		}
	}
}
