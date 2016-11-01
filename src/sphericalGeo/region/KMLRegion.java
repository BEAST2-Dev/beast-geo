package sphericalGeo.region;

import java.io.File;
import java.io.IOException;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

@Description("Geographical region defined by one or more polygons in a KML file")
public class KMLRegion extends Region {
	public Input<File> kmlFileInput = new Input<>("kml", "kml file with polygons over admissable locations.", Validate.REQUIRED);

	public KMLRegion() {}
	public KMLRegion(String kmlFile) throws Exception {initByName("kml", kmlFile);}
	
	@Override
	public void initAndValidate() {
		Log.info.println("Processing " + kmlFileInput.get());
		List<List<Double>> coordinates;
		try {
			width = 1024;
			height = 1024;
			coordinates = parseKML();
			calcAdmissableNodes(coordinates);
		} catch (Exception e) {
			throw new IllegalArgumentException(e);
		}
	}

	protected void calcAdmissableNodes(List<List<Double>> coordinates) {
		boolean debug = Boolean.valueOf(System.getProperty("beast.debug"));
		boolean traversesMapCenter = false;
		for (List<Double> coords : coordinates) {
			for (int i = 0; i < coords.size() - 2; i += 2) {
				double longitude = coords.get(i + 1);
				double longitude2 = coords.get(i + 3);
				if ((longitude > 150 && longitude2 < 150) ||
					(longitude2 > 150 && longitude < 150)) {
					traversesMapBoundary = true;
				}
				if ((longitude > -30 && longitude < 0 && longitude2 > 0 && longitude2 < 30 ) ||
						(longitude2 > -30 && longitude2 < 0 && longitude > 0 && longitude < 30 )) {
					traversesMapCenter = true;
				}
			}
		}
		if (traversesMapBoundary) {
			if (traversesMapCenter) {
				Log.warning.println("region too large: it crosses both center (long=0) and map boundary (long=180), expect problems");
			} else {
				for (List<Double> coords : coordinates) {
					for (int i = 0; i < coords.size(); i += 2) {
						double longitude = coords.get(i + 1);
						if (longitude < 0) {
							coords.set(i + 1, longitude + 360);
						}
					}
				}
			}
		}

		minLong = 360 * 2; maxLong = -360; maxLat = -90; minLat = 180;
		for (List<Double> coords : coordinates) {
			for (int i = 0; i < coords.size(); i += 2) {
				double latitude = coords.get(i);
				double longitude = coords.get(i + 1);
				minLong = Math.min(minLong, longitude);
				maxLong = Math.max(maxLong, longitude);
				minLat = Math.min(minLat, latitude);
				maxLat = Math.max(maxLat, latitude);
			}
		}

		// create bitmap to draw in
		double w = width / (maxLong - minLong);
		double h = height / (maxLat - minLat);

		// draw polygons in an image
		image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics g = image.getGraphics();

		if (debug) {
			try {
				final BufferedImage worldimage = ImageIO.read(new File("World98b.png"));
				System.err.println((int) (worldimage.getWidth() * (180 + minLong) / 360.0));
				System.err.println((int) (worldimage.getHeight() * (90 + minLat) / 180.0));
				System.err.println((int) (worldimage.getWidth() * (180 + maxLong) / 360.0));
				System.err.println((int) (worldimage.getHeight() * (90 + maxLat) / 180.0));
	
				g.drawImage(worldimage, 0, 0, width, height, 
						(int) (worldimage.getWidth() * (180 + minLong) / 360.0),
						(int) (worldimage.getHeight() * (90 - maxLat) / 180.0), 
						(int) (worldimage.getWidth() * (180 + maxLong) / 360.0), 
						(int) (worldimage.getHeight() * (90 - minLat) / 180.0), null);
				ImageIO.write(image, "png", new File("/tmp/kmlrange" + getID() + ".png"));
			} catch (Exception e) {
				// ignore
			}
		}

		for (List<Double> coords : coordinates) {
			int nPoints = coords.size() / 2;
			int[] xPoints = new int[nPoints];
			int[] yPoints = new int[nPoints];
			for (int i = 0; i < coords.size(); i += 2) {
				xPoints[i / 2] = (int) ((coords.get(i + 1) - minLong) * w);
				yPoints[i / 2] = height - (int) ((coords.get(i) - minLat) * h);
			}
			g.setColor(new Color(regionColor));
			g.fillPolygon(xPoints, yPoints, nPoints);
		}

		if (debug) {
			try {
//				g.setColor(Color.red);
//				double [] centre = centroid();
//				int x = (int)(width * (centre[1] - minLong)/(maxLong-minLong));
//				int y = (int)(height * (centre[0] - minLat)/(maxLat-minLat));
//				if ((image.getRGB(x, y) & 0xFFFFFF) != regionColor) {
//					System.err.println("centroid mismatch: " + getID());
//				}
//				g.fillOval(x - 2, y - 2, 5, 5);
				ImageIO.write(image, "png", new File("/tmp/kmlregion" + getID() + ".png"));
				
//				String s = BeautiDoc.load("/tmp/" + getID() + ".dat");
//				String [] strs = s.split("\n");
//				for (int i = 1; i < strs.length; i++) {
//					String [] x = strs[i].split("\\s+");
//					g.setColor(Color.blue);
//					if (x.length == 2) {
//						int x0 = (int) ((Double.parseDouble(x[1]) - minLong) * w);
//						int y0 = height - (int) ((Double.parseDouble(x[0]) - minLat) * h);
//						g.fillOval(x0-5, y0-5, 10, 10);
//					}
//				}
//				ImageIO.write(image, "png", new File("/tmp/kmlsamples" + getID() + ".png"));
			} catch (Exception e) {
				// ignore
			}
		}

	}

	private List<List<Double>> parseKML() throws SAXException, IOException, ParserConfigurationException {
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
	
}
