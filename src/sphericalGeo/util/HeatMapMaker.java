package sphericalGeo.util;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.imageio.ImageIO;

import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.app.util.Application;
import beast.app.util.ConsoleApp;
import beast.app.util.OutFile;
import beast.app.util.Utils;
import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.MRCAPrior;
import beast.util.NexusParser;
import beast.util.Randomizer;
import beast.util.XMLParser;
import beast.util.XMLParserException;

@Description("Creates heat map of locations with colour representing time")
public class HeatMapMaker extends Runnable {
	public Input<File> treesetInput = new Input<>("treeSet","file containing tree set annotated with locations");//, new File("IE.trees"));
	public Input<File> bgInput = new Input<>("background","image file with a map for the background");//, new File("/Users/remco/data/map/World98b.png"));
	public Input<String> bboxInput = new Input<>("boundingBox","bounding box of the background image", "-90 -180 90 180");
	public Input<String> tagInput = new Input<>("tag","tag used in annotated of locations", "location");
	public Input<Integer> widthInput = new Input<>("width","width of the heat map", 1024);
	public Input<Integer> heightInput = new Input<>("height","heightof the heat map", 1024);
	public Input<Double> maxTimeInput = new Input<>("maxTime","maximum time (all older nodes will be coloured red)", Double.POSITIVE_INFINITY);
	public Input<OutFile> outputInput = new Input<>("output","where to save the file", new OutFile("heatmap.png"));
	public Input<Integer> discSizeInput = new Input<>("discSize","size of the dots used to draw heat map", 10);
	public Input<Double> translucencyInput = new Input<>("translucency","translucency of the dots used to draw heat map (a number between 0=invisible and 1=solid)", 0.2);
	public Input<Double> saturationInput = new Input<>("saturation","saturation of colour for the dots", 0.9);
	public Input<Double> brightnessInput = new Input<>("brightness","brightnessof colour for the dots ", 0.9);
	public Input<Integer> jitterInput = new Input<>("jitter", "jitter applied to dot locations (in pixels)", 0);
	public Input<File> maskInput = new Input<>("mask", "image file with a mask: dots will not be shown outside mask");//, new File("/Users/remco/data/map/World98b.png"));
	public Input<Integer> repeatsInput = new Input<>("repeats", "number of times a dot should be drasn", 1);
	public Input<Boolean> rootOnlyInput = new Input<>("rootOnly", "only draw root locations and ignore the rest", false);
	
	public Input<XMLFile> xmlInput = new Input<>("clade", "Name of XML file containing a single TaxonSet in BEAST XML format. "
			+ "If specified, the location of the MRCA of a clade is used (and the rootonly flag is ignored).");
	
	BufferedImage image;

	double minLat = -90, maxLat = 90, minLong = -180, maxLong = 180;
	
	int width;
	int height;
	String tag;
	int repeats = 1;
	boolean rootOnly;
	TaxonSet taxonset = null;
	
	final static String DIR_SEPARATOR = (Utils.isWindows() ? "\\\\" : "/");
	@Override
	public void initAndValidate() {
		width = widthInput.get();
		height = heightInput.get();
        image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		tag = tagInput.get();
		rootOnly = rootOnlyInput.get();
	}

	@Override
	public void run() throws Exception {
		long start = System.currentTimeMillis();
		if (xmlInput.get() != null && !xmlInput.get().getName().equals("[[none]]")) {
			taxonset = getTaxonSet(xmlInput.get());
		}
		
		Graphics g = image.getGraphics();
		int jitter = jitterInput.get();
		repeats = repeatsInput.get();
		rootOnly = rootOnlyInput.get();

		File bg = bgInput.get();
		if (bg != null && ! bg.getName().equals("[[none]]")) {
			if (bg.exists()) {
				System.out.println("Loading background image " + bg.getPath());
				BufferedImage bgImage = ImageIO.read(bg);
				parseBBox();
				g.drawImage(bgImage, 0, 0, width, height, 
						(int) (bgImage.getWidth() * (180 + minLong) / 360.0),
						(int) (bgImage.getHeight() * (90 - maxLat) / 180.0), 
						(int) (bgImage.getWidth() * (180 + maxLong) / 360.0), 
						(int) (bgImage.getHeight() * (90 - minLat) / 180.0), null);
			} else {
				System.out.println("Could not find backgroung image " + bg.getPath());
			}
		} else {
			System.out.println("No background image");
		}
		File maskFile = maskInput.get();
		BufferedImage mask = null;
		if (maskFile != null && !maskFile.getName().equals("[[none]]")) {
			if (maskFile.exists()) {
				System.out.println("Loading mask image " + maskFile.getPath());
				BufferedImage maskImage = ImageIO.read(maskFile);
				parseBBox();
		        mask = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		        Graphics g2 = mask.getGraphics();
				g2.drawImage(maskImage, 0, 0, width, height, 
						(int) (maskImage.getWidth() * (180 + minLong) / 360.0),
						(int) (maskImage.getHeight() * (90 - maxLat) / 180.0), 
						(int) (maskImage.getWidth() * (180 + maxLong) / 360.0), 
						(int) (maskImage.getHeight() * (90 - minLat) / 180.0), null);
				try {
					ImageIO.write(mask, "png", new File("/tmp/mask.png"));					
				} catch (Exception e) {
					// ignore -- this is for debugging only
				}

			} else {
				System.out.println("Could not find mask image " + maskFile.getPath());
			}
		} else {
			System.out.println("No mask image");
		}

		// get trees
		System.out.print("Parsing trees " + treesetInput.get().getPath() + " ... ");
		NexusParser parser = new NexusParser();
		parser.parseFile(treesetInput.get());
		List<Tree> trees = parser.trees;
		if (trees == null) {
			System.out.println("File did not contain any trees. Is this a tree set file?");
			return;
		}
		System.out.println(" done");
		
		// get dots
		System.out.println("Drawing dots...");
		System.out.println("|=========|=========|=========|=========|=========|=========|=========|=========|=========|=========|");
		List<Dot> dots = new ArrayList<>();
		for (Tree tree: trees) {
			if (taxonset == null) {
				collectDots(tree.getRoot(),dots);
			} else {
				MRCAPrior p = new MRCAPrior();
				p.initByName("taxonset", taxonset, "tree", tree);
				Node mrca = p.getCommonAncestor();
				String meta =  mrca.metaDataString;
				dots.add(new Dot(mrca.getHeight(), meta));
			}
		}
		Collections.sort(dots, new Comparator<Dot>() {
			@Override
			public int compare(Dot o1, Dot o2) {
				if (o1.age > o2.age) {
					return 1;
				} else if (o1.age < o2.age) {
					return -1;
				}
				return 0;
			}
		});
		
		// draw dots
		double oldest = calcOldest(trees);
		int radius = discSizeInput.get();
		int halfRadius = radius / 2;
		double translucency = translucencyInput.get();
		((Graphics2D) g).setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, (float) translucency));
		float brightness = (float) (double) brightnessInput.get();
		float saturation = (float) (double) saturationInput.get();
		
		int k = 0;
		for (Dot dot : dots) {
			double f = Math.min(dot.age, oldest)/oldest;
			
//			 int red   = (int) (Math.sin(2*3.14*f + 0) * 127 + 128);
//			 int green = (int) (Math.sin(2*3.14*f + 2) * 127 + 128);
//			 int blue  = (int) (Math.sin(2*3.14*f + 4) * 127 + 128);
//			 g.setColor(new Color(red,green,blue));
			
			g.setColor(Color.getHSBColor(0.2f + (float) (f/0.8f), saturation, brightness));
			
			for (int i = 0; i < repeats; i++) {
				int y = height - (int)( (dot.latitude - minLat) * height/(maxLat-minLat));  
				int x = (int)( (dot.longitude - minLong) * width/(maxLong-minLong)); 
				if (jitter > 0) {
					x += Randomizer.nextGaussian() * jitter;
					y += Randomizer.nextGaussian() * jitter;
				}
				if (mask == null) {
					g.fillOval(x-halfRadius, y-halfRadius, radius, radius);
				} else {
					if (x > 0 && x < width && y > 0 && y < height && ((mask.getRGB(x,  y) & 0xff0000) == 0)) {
						g.fillOval(x-halfRadius, y-halfRadius, radius, radius);
					}
				}
			}
			k++;
			if (k % (Math.max(1, dots.size()/100)) == 0) {
				System.err.print('.');
			}
		}
		System.out.println("\ndone");
		
		// write file
		System.out.print("Writing file " + outputInput.get().getPath());
		ImageIO.write(image, "png", outputInput.get());
		System.out.println(" done");
		
		
		// create legend
		System.out.print("Writing legend " + outputInput.get().getParent() + DIR_SEPARATOR + "legend.png ");
        image = new BufferedImage(200, 200, BufferedImage.TYPE_INT_RGB);
        g = image.getGraphics();
        g.setColor(Color.white);
        g.fillRect(0, 0, 200, 200);
        
        for (int i = 0; i < 200; i++) {
        	double f = (200 - i) / 200.0;
			g.setColor(Color.getHSBColor(0.2f + (float) (f/0.8f), 1.0f, 1.0f));
        	g.drawLine(25, i, 75, i);
        }
		g.setColor(Color.black);
        for (int i = 0; i < 10; i++) {
        	g.drawString(""+ (oldest * (10 - i)/ 10), 100, 200*i/10);
        }
		ImageIO.write(image, "png", new File((outputInput.get().getParent() != null ? outputInput.get().getParent() + DIR_SEPARATOR : "") + "legend.png"));
		long end = System.currentTimeMillis();
		System.out.println(" done in " + (end - start) / 1000 + " seconds");
		
	}

	static TaxonSet getTaxonSet(File xmlFile) throws IOException, XMLParserException {
		XMLParser parser = new XMLParser();
		String xml = BeautiDoc.load(xmlFile);
		Object o = parser.parseBareFragment(xml, true);
		if (o instanceof TaxonSet) {
			return (TaxonSet) o;
		}
		return null;
	}

	private void collectDots(Node node, List<Dot> dots) {
		String meta =  node.metaDataString;
		dots.add(new Dot(node.getHeight(), meta));
		if (!rootOnly) {
			for (Node child:  node.getChildren()) {
				collectDots(child, dots);
			}
		}
	}

	class Dot {
		public Dot(double height, String meta) {
			age = height;
			String [] metaDatas = meta.split(",");
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
				if (str.startsWith(tag + "1=")) {
					y =  Double.parseDouble(str.substring(tag.length() + 2));
				}
				if (str.startsWith(tag + "2=")) {
					x =  Double.parseDouble(str.substring(tag.length() + 2));
				}
				if (str.startsWith(tag + "={")) {
					y = Double.parseDouble(str.substring(tag.length() + 2));
					str = metaDatas[k+1];
					str = str.substring(0, str.indexOf('}'));
					x = Double.parseDouble(str);
					//isRawLatLong = true;
				}
				if (str.startsWith("xy={")) {
					x = Double.parseDouble(str.substring(4));
					str = metaDatas[k+1];
					str = str.substring(0, str.indexOf('}'));
					y = Double.parseDouble(str);
				}
				k++;
			}
			latitude = y;
			longitude = x;
		}
		
		double age;
		double latitude;
		double longitude;
	}
	
	private double calcOldest(List<Tree> trees) {
		double h = 0;
		for(Tree tree : trees) {
			h = Math.max(h,  tree.getRoot().getHeight());
		}
		double maxTime = maxTimeInput.get();
		h = Math.min(h, maxTime);
		return h;
	}
	
    void parseBBox() {
        if (bboxInput.get() == null) {
                return;
        }
        String str = bboxInput.get().trim();
        String [] strs = str.split("\\s+");
        if (strs.length != 4) {
                throw new RuntimeException("bbox input must contain 4 numbers");
        }
        minLat = Double.parseDouble(strs[0]);
        minLong = Double.parseDouble(strs[1]);
        maxLat = Double.parseDouble(strs[2]);
        maxLong = Double.parseDouble(strs[3]);
        if (minLat >= maxLat) {
                double tmp = maxLat; maxLat = minLat; minLat = tmp;
        }
        if (minLong >= maxLong) {
            double tmp = maxLong; maxLong = minLong; minLong = tmp;
        }
    }

    static ConsoleApp app;
	public static void main(String[] args) throws Exception {
		HeatMapMaker heatMapMaker = new HeatMapMaker();
	
		if (args.length == 0) {
			// create BeautiDoc and beauti configuration
			BeautiDoc doc = new BeautiDoc();
			doc.beautiConfig = new BeautiConfig();
			doc.beautiConfig.initAndValidate();
			
			// suppress a few inputs that we don't want to expose to the user
			doc.beautiConfig.suppressBEASTObjects.add(heatMapMaker.getClass().getName() + ".mcmc");
			doc.beautiConfig.suppressBEASTObjects.add(heatMapMaker.getClass().getName() + ".value");
			doc.beautiConfig.suppressBEASTObjects.add(heatMapMaker.getClass().getName() + ".hosts");
		
			// create panel with entries for the application
			BEASTObjectPanel panel = new BEASTObjectPanel(heatMapMaker, heatMapMaker.getClass(), doc);
			
			// wrap panel in a dialog
			BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
	
			// show the dialog
			if (dialog.showDialog()) {
				dialog.accept(heatMapMaker, doc);
				// create a console to show standard error and standard output
				//ConsoleApp 
				app = new ConsoleApp("HeatMapMaker", "Heat Map Maker: " + heatMapMaker.treesetInput.get().getPath(), null);
				heatMapMaker.initAndValidate();
				heatMapMaker.run();
			}
			System.out.println("All done");
			return;
		}

		Application main = new Application(heatMapMaker);
		main.parseArgs(args, false);
		heatMapMaker.initAndValidate();
		heatMapMaker.run();
		System.out.println("All done");
	}



}
