package sphericalGeo.region;

import java.io.File;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

@Description("Geographical region defined by one or more polygons in a KML file")
public class KMLRegion extends Region {
	public Input<File> kmlFileInput = new Input<File>("kml", "kml file with polygons over admissable locations.", Validate.REQUIRED);

	@Override
	public void initAndValidate() throws Exception {
		List<List<Double>> coordinates = parseKML();
		calcAdmissableNodes(coordinates);
	}

	private void calcAdmissableNodes(List<List<Double>> coordinates) throws Exception {
		boolean debug = true;// Boolean.valueOf(System.getProperty("beast.debug"));

		minLong = 360; maxLong = -360; maxLat = -90; minLat = 180;
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
		width = 1024;
		height = 1024;
		double w = width / (maxLong - minLong);
		double h = height / (maxLat - minLat);

		// draw polygons in an image
		image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics g = image.getGraphics();

		if (debug) {
			final BufferedImage worldimage = ImageIO.read(new File("World98b.png"));
			System.err.println((int) (worldimage.getWidth() * (180 + minLong) / 360.0));
			System.err.println((int) (worldimage.getHeight() * (90 + minLat) / 180.0));
			System.err.println((int) (worldimage.getWidth() * (180 + maxLong) / 360.0));
			System.err.println((int) (worldimage.getHeight() * (90 + maxLat) / 180.0));

			g.drawImage(worldimage, 0, 0, width, height, (int) (worldimage.getWidth() * (180 + minLong) / 360.0),
					(int) (worldimage.getHeight() * (minLat) / 180.0), (int) (worldimage.getWidth() * (180 + maxLong) / 360.0), (int) (worldimage.getHeight()
							* (maxLat) / 180.0), null);
			ImageIO.write(image, "png", new File("/tmp/kmlrange" + getID() + ".png"));
		}

		for (List<Double> coords : coordinates) {
			int nPoints = coords.size() / 2;
			int[] xPoints = new int[nPoints];
			int[] yPoints = new int[nPoints];
			for (int i = 0; i < coords.size(); i += 2) {
				xPoints[i / 2] = (int) ((coords.get(i + 1) - minLong) * w);
				yPoints[i / 2] = (int) ((coords.get(i) - minLat) * h);
			}
			g.setColor(new Color(regionColor));
			g.fillPolygon(xPoints, yPoints, nPoints);
		}

		if (debug) {
			ImageIO.write(image, "png", new File("/tmp/kmlregion" + getID() + ".png"));
		}

	}

	private List<List<Double>> parseKML() throws Exception {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		factory.setValidating(false);
		org.w3c.dom.Document doc = factory.newDocumentBuilder().parse(kmlFileInput.get());
		doc.normalize();

		List<List<Double>> coordinates = new ArrayList<List<Double>>();

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
					polygon.add(90 - Double.parseDouble(sCoords[1].trim()));
					polygon.add(Double.parseDouble(sCoords[0].trim()));
				}
			}
			coordinates.add(polygon);
		}
		return coordinates;
	}

}