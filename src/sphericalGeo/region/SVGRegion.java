package sphericalGeo.region;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import beast.app.beauti.BeautiDoc;
import beast.core.Input.Validate;

public class SVGRegion extends KMLRegion {
	double lat, lon;
	
	List<List<Double>> coordinates;

	public SVGRegion(String id, double lat, double lon, String path) {
		kmlFileInput.setRule(Validate.OPTIONAL);
		
		width = 512;
		height = 512;
		
		setID(id);
		this.lat = lat;
		this.lon = lon;
		try {
			coordinates = parseSVG(path);
			calcAdmissableNodes(coordinates);
		} catch (Exception e) {
			throw new IllegalArgumentException(e);
		}
		double [] p = centroid();
	}
	
	
	private List<List<Double>> parseSVG(String path) {
		List<List<Double>> coordinates = new ArrayList<>();

		List<Double> polygon = new ArrayList<>();
		String [] paths = path.split("Z");
		for(String p : paths) {
			String [] points = p.split("L");
			for (String sStr : points) {
				if (sStr.contains(",")) {
					String[] sCoords = sStr.split(",");
					polygon.add(Double.parseDouble(sCoords[1].trim()));
					polygon.add(Double.parseDouble(sCoords[0].trim()));
				}
			}
			// close the loop
			polygon.add(polygon.get(0));
			polygon.add(polygon.get(1));
			coordinates.add(polygon);
		}
		return coordinates;
	}

	
	static public List<SVGRegion> parseCFGFile(String fileName) throws IOException {
		String [] strs = BeautiDoc.load(fileName).split("\n");
		List<SVGRegion> list = new ArrayList<>();
		for (String s : strs) {
			String [] x = s.split("=");
			String [] y = x[1].split(",");
			double lat = Double.parseDouble(y[0]);
			double lon = Double.parseDouble(y[1]);
			SVGRegion region = new SVGRegion(x[0],lat,lon,x[2]);
			list.add(region);
		}
		return list;
	}


	public double squaredDistance(double lat2, double lon2) {
		return (lat-lat2)*(lat-lat2) + (lon-lon2)*(lon-lon2);
	}
	
}
