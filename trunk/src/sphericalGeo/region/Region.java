package sphericalGeo.region;

import java.awt.image.BufferedImage;

import beast.core.BEASTObject;
import beast.core.Description;

@Description("Defines a geographical region")
public class Region extends BEASTObject {

	double minLat = -90, maxLat = 90, minLong = -180, maxLong = 180;
	BufferedImage image;
	int regionColor =  0x00FF00;
    int width = 1024;
    int height = 1024;
	
	@Override
	public void initAndValidate() throws Exception {
        image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
	}
	
	public boolean isInside(double latitude, double longitude) {
		int iLat = (int)(height * (latitude - minLat) / (maxLat - minLat));
		int iLong = (int)(width * (longitude - minLong) / (maxLong - minLong));
		if (iLat < 0 || iLat >= width || iLong < 0 || iLong >= height) {
			// outside bitmap
			return false;
		}
		// check bitmap
		int color = image.getRGB(iLat, iLong) & 0xFFFFFF;
		return (color == regionColor);
	}

}
