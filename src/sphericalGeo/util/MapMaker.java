package sphericalGeo.util;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.core.util.Log;

@Description("Creates a map in mercator projection from a KML file")
public class MapMaker extends Runnable {
	public Input<File> kmlFileInput = new Input<>("kmlfile", "kml file contiaining polygons to draw", Validate.REQUIRED);
	public Input<Integer> widthInput = new Input<>("width","width of the map", 8000);
	public Input<Integer> heightInput = new Input<>("height","heightof the map", 4000);
	public Input<File> outputInput = new Input<>("output","where to save the file", new File("/tmp/map.png"));
	public Input<Color> fillColorInput = new Input<>("fillColour","colour used to fill", Color.white);
	public Input<Color> lineColorInput = new Input<>("lineColour","colour used to draw borders", Color.black);
	public Input<Color> bgColorInput = new Input<>("bgColour","colour used for background", new Color(0xa0a0ff));
	public Input<String> bboxInput = new Input<>("boundingBox","bounding box of the background image", "-90 -180 90 180");
	public Input<Double> lineWidthInput = new Input<>("lineWidth","line width used for drawing borders", 1.0);

	BufferedImage image;
	double minLat = -90, maxLat = 90, minLong = -180, maxLong = 180;

	@Override
	public void initAndValidate() throws Exception {
		if (!kmlFileInput.get().exists()) {
			throw new RuntimeException("kml file " + kmlFileInput.get().getPath() + "does not exist");
		}
		parseBBox();
	}
	
	public void run() throws Exception {
		Log.info.println("Parsing KML file...");
		List<List<Double>> coordinates = parseKML();
		Log.info.println("Creating map...");
		
		Color fillColour = fillColorInput.get();
		Color lineColour = lineColorInput.get();

		int width = widthInput.get();
		int height = heightInput.get();
		
		image = new BufferedImage(widthInput.get(), heightInput.get(), BufferedImage.TYPE_INT_RGB);
		Graphics g = image.getGraphics();
		((Graphics2D)g).setStroke(new BasicStroke((float)(double)lineWidthInput.get()));
		g.setColor(bgColorInput.get());
		g.fillRect(0,  0, width, height);
		
		double w = width / (maxLong - minLong);
		double h = height / (maxLat - minLat);
		
		for (List<Double> coords : coordinates) {
			int nPoints = coords.size() / 2;
			int[] xPoints = new int[nPoints];
			int[] yPoints = new int[nPoints];
			for (int i = 0; i < coords.size(); i += 2) {
				xPoints[i / 2] = (int) ((coords.get(i + 1) - minLong) * w);
				yPoints[i / 2] = height - (int) ((coords.get(i) - minLat) * h);
			}
			g.setColor(fillColour);
			g.fillPolygon(xPoints, yPoints, nPoints);
			g.setColor(lineColour);
			g.drawPolygon(xPoints, yPoints, nPoints);
		}

		
		Log.info.println("Saving KML...");
		ImageIO.write(image, "png", outputInput.get());
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
    
	private List<List<Double>> parseKML() throws Exception {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		factory.setValidating(false);
		org.w3c.dom.Document doc = factory.newDocumentBuilder().parse(kmlFileInput.get());
		doc.normalize();

		List<List<Double>> coordinates = new ArrayList<>();

		// grab 'coordinates' elements out of the KML file
		NodeList oCoordinates = doc.getElementsByTagName("coordinates");
		for (int iNode = 0; iNode < oCoordinates.getLength(); iNode++) {
			Node oCoordinate = oCoordinates.item(iNode);
			String sCoordinates = oCoordinate.getTextContent();
			List<Double> polygon = new ArrayList<>();
			String[] sStrs = sCoordinates.split("\\s+");
			for (String sStr : sStrs) {
				if (sStr.contains(",")) {
					String[] sCoords = sStr.split(",");
					polygon.add(Double.parseDouble(sCoords[1].trim()));
					polygon.add(Double.parseDouble(sCoords[0].trim()));
				}
			}
			coordinates.add(polygon);
		}
		return coordinates;
	}

	public static void main(String[] args) throws Exception {
		MapMaker mm = new MapMaker();
		mm.initByName("kmlfile", "/home/remco/data/map/borders/KML_zip_20141212060828.kml");
		
		if (args.length == 0) {
			// create BeautiDoc and beauti configuration
			BeautiDoc doc = new BeautiDoc();
			doc.beautiConfig = new BeautiConfig();
			doc.beautiConfig.initAndValidate();
		
			// create panel with entries for the application
			BEASTObjectPanel panel = new BEASTObjectPanel(mm, mm.getClass(), doc);
			
			// wrap panel in a dialog
			BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
	
			// show the dialog
			if (dialog.showDialog()) {
				dialog.accept(mm, doc);
				// create a console to show standard error and standard output
				//ConsoleApp app = new ConsoleApp("PathSampler", "Path Sampler: " + sampler.model1Input.get().getPath());
				mm.initAndValidate();
				mm.run();
			}
			Log.info.println("Done");
			return;
		}

		Application main = new Application(mm);
		main.parseArgs(args, false);
		mm.initAndValidate();
		mm.run();
		Log.info.println("Done");

	}
}
