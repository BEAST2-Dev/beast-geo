package sphericalGeo.util;


import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.WindowConstants;

//import com.itextpdf.awt.PdfGraphics2D;
//import com.itextpdf.text.pdf.PdfContentByte;
//import com.itextpdf.text.pdf.PdfWriter;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.util.Randomizer;
import beast.base.evolution.tree.TreeParser;

public class RandomWalkSimulator extends JPanel implements MouseListener, KeyListener {
	
	private static final long serialVersionUID = 1L;


	enum MODE {simple, levy, biased, correlated};
	class RandomWalk {
		MODE mode;
		int nSteps;
		double [][] point;
		double lastAngle;
		
		RandomWalk(int nSteps, MODE mode, double startX, double startY, double lastAngle) {
			this.nSteps = nSteps;
			this.mode = mode;
			point = new double[nSteps][2]; 
			this.lastAngle = lastAngle;
			walk(startX, startY);
		}

		void walk(double startX, double startY) {
			point[0][0] = startX;
			point[0][1] = startY;
			for (int i = 1; i < nSteps; i++) {
				
				double angle = 0, length = 0;
				switch (mode) {
				case simple:
					angle= Randomizer.nextDouble() * Math.PI * 2.0;
					length = Randomizer.nextExponential(0.125);
					break;
				case levy:
					angle= Randomizer.nextDouble() * Math.PI * 2.0;
					length = Math.tan(Math.PI * (Randomizer.nextDouble() - 0.5)/ 0.1);
					break;
				case biased:
					double angle0 = -Math.atan2(point[0][1] - dp[1], point[0][0] - dp[0]);
					angle = angle0 + Randomizer.nextGaussian() * Math.PI / 2;
					length = Randomizer.nextExponential(0.25);
					break;
				case correlated:
					if (i > 2) {
						angle0 = Math.atan2(point[i-1][1] - point[i-2][1], point[i-1][0] - point[i-2][0]);
					} else {
						angle0 = lastAngle;
					}
					angle = angle0 + Randomizer.nextGaussian() * Math.PI /8.33;
					length = Randomizer.nextExponential(0.25);
					break;
				}
				point[i][0] = startX + Math.cos(angle) * length;
				point[i][1] = startY + Math.sin(angle) * length;
				startX = point[i][0];
				startY = point[i][1];
				lastAngle = angle;
			}
		}
		
	}
	
	
	int H, W;
	double [] dp;
	int nBlues = 1;
	MODE mode = MODE.simple;
	double lastAngle = 0;
	boolean showReconstruction = false;
	
	@Override
	protected void paintComponent(Graphics g) {
		H = getHeight();
		W = getWidth();
		g.setColor(Color.white);
		g.fillRect(0, 0, W, H);
		Graphics2D g2 = (Graphics2D) g;
		
		double [] A;double [] B = null;double [] C;
		double [] end = new double[]{W/2, H/2};
		dp = new double[]{W/2, H/2};
		
		int nSteps = 200;
		
		g.setColor(Color.black);
		g2.drawString("mode=" + mode, 10, 10);
		g.fillOval(W/2-10, H/2-10, 20, 20);
		end = walk(end[0], end[1], g2, nSteps, lastAngle);
		double lastAngle0 = lastAngle;
		
		g.setColor(Color.red);
		A = walk(end[0], end[1], g2, nSteps, lastAngle0);
		
		g.setColor(new Color(0x0066CC));
		for (int i = 0; i < nBlues; i++) {
			B = walk(end[0], end[1], g2, nSteps, lastAngle0);
		}

		dp = end;
		g.setColor(new Color(0x009900));
		C = walk(W/2.0, H/2.0, g2, nSteps * 2, Math.PI);


		if (showReconstruction) {
		try {
			double [][] p = getTree(A, B, C);
			g.setColor(new Color(0x9933FF));
			g.fillOval((int)p[4][0]-10, (int)p[4][1]-10, 20, 20);
			line(g, p[0], p[3]);
			line(g, p[1], p[3]);
			line(g, p[3], p[4]);
			line(g, p[2], p[4]);
		} catch (Exception e) {
			e.printStackTrace();
		}
		}
	}		
		
		
	void exportPDF(String sFileName) {
		try {
//			com.itextpdf.text.Document doc = new com.itextpdf.text.Document();
//			PdfWriter writer = PdfWriter.getInstance(doc, new FileOutputStream(sFileName));
//			doc.setPageSize(new com.itextpdf.text.Rectangle(getWidth(), getHeight()));
//			doc.open();
//			PdfContentByte cb = writer.getDirectContent();
//			Graphics2D g = new PdfGraphics2D(cb, getWidth(), getHeight());
//			 
//			BufferedImage bi;
//			bi = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
//			//g = bi.getGraphics();
//			g.setPaintMode();
//			g.setColor(getBackground());
//			g.fillRect(0, 0, getWidth(), getHeight());
//			paint(g);
//			//m_Panel.printAll(g);
//		
//			g.dispose();
//			doc.close();
		} catch (Exception e) {
			JOptionPane.showMessageDialog(this, "Export may have failed: " + e.getMessage());
		}
	}

	private void line(Graphics g, double[] p0, double[] p1) {
		g.drawLine((int) p0[0], (int) p0[1], (int) p1[0], (int) p1[1]);
		
	}


	public double [][] getTree(double[]pA, double[]pB, double[]pC) throws Exception {
	        Sequence A = new Sequence("A", "?");
	        Sequence B = new Sequence("B", "?");
	        Sequence C = new Sequence("C", "?");

	        Alignment alignment = new Alignment();
	        alignment.initByName("sequence", A, "sequence", B, "sequence", C);

	        TreeParser t = new TreeParser();
	        t.initByName("taxa", alignment,
	                "newick", "((A:1.0,B:1.0):1.0,C:2.0)",
	                "IsLabelledNewick", true);

	        sphericalGeo.SphericalDiffusionModel sModel = new sphericalGeo.SphericalDiffusionModel();
	        sModel.initByName("precision","1.0");
	        SiteModel siteModel = new SiteModel();
	        siteModel.initByName("substModel", sModel);
	        
            RealParameter location = new RealParameter();
            location.initByName("dimension", 2, "minordimension", 2, "value", "0.0 0.0");

            sphericalGeo.TreeTraitMap map = new sphericalGeo.TreeTraitMap();
            double scale = 100;
            map.initByName("tree", t, "traitName", "location", "randomizeupper", "90 180", "randomizelower","-90 -180","initByMean", true,
            		"parameter", location,
            		"value", "A=" + (pA[0] - W/2)/scale + " " + (pA[1] - H/2)/ scale+ 
            		       ", B=" + (pB[0] - W/2)/scale + " " + (pB[1] - H/2)/ scale+
            		       ", C=" + (pC[0] - W/2)/scale + " " + (pC[1] - H/2)/ scale);
            
            sphericalGeo.AlignmentFromTraitMap data = new sphericalGeo.AlignmentFromTraitMap();
            data.initByName("userDataType", new sphericalGeo.LocationDataType(), "traitMap", map);
            
            
            sphericalGeo.ApproxMultivariateTraitLikelihood likelihood = new sphericalGeo.ApproxMultivariateTraitLikelihood();
            likelihood.initByName("tree", t, "scale", false, "logAverage", true, "location", location, "siteModel", siteModel, "data", data);
            double [][] p = likelihood.getPositions();
            for (int i = 0; i < p.length; i++) {
            	p[i][0] = (p[i][0] * scale + W/2);
            	p[i][1] = (p[i][1] * scale + H/2);
            }
            return p;
	}

	private double[] walk(double startX, double startY, Graphics2D g, int nSteps, double lastAngle) {
		RandomWalk walk = new RandomWalk(nSteps, mode, startX, startY, lastAngle);
		
		double [] p0 = walk.point[0];
		g.setStroke(new BasicStroke(1.0f));
		for (int i = 1; i < nSteps; i++) {
			double [] p1 = walk.point[i];
			g.drawLine((int) p0[0], (int) p0[1], (int) p1[0], (int) p1[1]);
			p0 = p1;
		}
		g.setStroke(new BasicStroke(3.0f));
		p0 = walk.point[0];
		double [] p1 = walk.point[nSteps - 1];
		g.drawLine((int) p0[0], (int) p0[1], (int) p1[0], (int) p1[1]);
		
		lastAngle = walk.lastAngle;
		return walk.point[nSteps-1];
	}

	public static void main(String[] args) {
		JFrame frame = new JFrame();
		RandomWalkSimulator panel =  new RandomWalkSimulator();
		frame.add(panel);
		frame.addMouseListener(panel);
		frame.addKeyListener(panel);
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.setSize(1024,768);
		frame.setVisible(true);

	}
	@Override
	public void mouseClicked(MouseEvent e) {
		if (e.isShiftDown() && e.isAltDown()) {
			exportPDF("/tmp/x.pdf");
			return;
		}
		repaint();
		if (e.isShiftDown()) {
			switch (mode) {
			case simple: mode = MODE.biased; break;
			case biased: mode = MODE.correlated; break;
			case correlated: mode = MODE.levy; break;
			case levy: mode = MODE.simple; break;
			}
		}
		if (e.isControlDown()) {
			nBlues = 1000-nBlues;
		}
		if (e.isAltDown()) {
			showReconstruction = !showReconstruction;
		}
	}
	@Override
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void keyTyped(KeyEvent e) {
	}


	@Override
	public void keyPressed(KeyEvent e) {
		if (e.getKeyCode() == KeyEvent.VK_P) {
			exportPDF("/tmp/x.pdf");
		}
	}


	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}

}
