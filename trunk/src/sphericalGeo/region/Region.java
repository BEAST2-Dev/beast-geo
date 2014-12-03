package sphericalGeo.region;

import java.awt.image.BufferedImage;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.util.Randomizer;

@Description("Defines a geographical region")
public class Region extends BEASTObject {

	double minLat = -90, maxLat = 90, minLong = -180, maxLong = 180;
	BufferedImage image;
	int regionColor =  0x00FF00;
    int width = 1024;
    int height = 1024;

    /** flag indicating the region contains points both sides of longitude=180 **/
	boolean traversesMapBoundary = false;

	@Override
	public void initAndValidate() throws Exception {
        image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
	}
	
	public boolean isInside(double latitude, double longitude) {
		if (traversesMapBoundary && longitude < 0) {
			longitude += 360;
		}
		int iLat = (int)(height * (latitude - minLat) / (maxLat - minLat));
		int iLong = (int)(width * (longitude - minLong) / (maxLong - minLong));
		if (iLat < 0 || iLat >= width || iLong < 0 || iLong >= height) {
			// outside bitmap
			return false;
		}
		// check bitmap
		int color = image.getRGB(iLong, iLat) & 0xFFFFFF;
		return (color == regionColor);
	}

	public double [] sample() {
		int i;
		do {
			i = Randomizer.nextInt(width * height);
		} while ((image.getRGB(i % width , i / width) & 0xFFFFFF) != regionColor);

		double [] location = new double[2];
		location[0] = minLat + (maxLat - minLat) * (i + 0.5) /(width * height);
		location[1] = minLong + (maxLong - minLong) * (i % width + 0.5)/width;

		// sanity check
		//if (!isInside(location[0], location[1])) {
		//	System.err.println("Failed to sample region!");
		//}
		
		return location;
	}
}
