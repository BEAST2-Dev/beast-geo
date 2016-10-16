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
	public void initAndValidate() {
        image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
	}
	
	public boolean isInside(double latitude, double longitude) {
		if (traversesMapBoundary && longitude < 0) {
			longitude += 360;
		}
		int iLat = height - 1 - (int)(height * (latitude - minLat) / (maxLat - minLat));
		int iLong = (int)(width * (longitude - minLong) / (maxLong - minLong));
		if (iLat < 0 || iLat >= width || iLong < 0 || iLong >= height) {
			// outside bitmap
			return false;
		}
		// check bitmap
		int color = image.getRGB(iLong, iLat) & 0xFFFFFF;
		return (color == regionColor);
	}

	public double [] sample(boolean isInside) {
		int i;
		if (isInside) {
			do {
				i = Randomizer.nextInt(width * height);
			} while ((image.getRGB(i % width , i / width) & 0xFFFFFF) != regionColor);
		} else {
			do {
				i = Randomizer.nextInt(width * height);
			} while ((image.getRGB(i % width , i / width) & 0xFFFFFF) == regionColor);
		}

		double [] location = new double[2];
		location[0] = maxLat - (maxLat - minLat) * (i - 0.5) /(width * height);
		location[1] = minLong + (maxLong - minLong) * (i % width + 0.5)/width;

		// sanity check
		if (isInside) {
			if (!isInside(location[0], location[1])) { 
				location = sample(isInside);
			}
		} else {
			if (isInside(location[0], location[1])) { 
				location = sample(isInside);
			}
		}
		
		return location;
	}
	
	/** return lat/long pair that represents the center of the region 
	 * If the centroid is outside the region, return closest point on image
	 * to the centroid. **/
	public double [] centroid() {
		double lat = 0;
		double lon = 0;
		int k = 0;
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				if ((image.getRGB(i, j) & 0xFFFFFF) == regionColor) {
					lat += j;
					lon += i;
					k++;
				}
			}
		}
		lat = lat / k;
		lon = lon / k;
		
		lat = minLat + (lat / height) * (maxLat - minLat);
		lon = minLong + (lon / width) * (maxLong - minLong);
		
		int x = (int)(width * (lon - minLong)/(maxLong-minLong));
		int y = (int)(height * (lat - minLat)/(maxLat-minLat));
		if ((image.getRGB(x, y) & 0xFFFFFF) != regionColor) {
			double minDist = Double.MAX_VALUE;
			int minX = 0;
			int minY = 0;
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					if ((image.getRGB(i, j) & 0xFFFFFF) == regionColor) {
						double d = (x-i)*(x-i)+(y-j)*(y-j);
						if (d < minDist) {
							minDist = d;
							minX = i;
							minY = j;
						}
					}
				}
			}
			
			lat = minLat + ((minY + 0.5) * (maxLat - minLat) / height);
			lon = minLong + ((minX + 0.5) * (maxLong - minLong) / width);
		}

		
		System.out.println(getID() + "=" + lat + " " + lon);
		return new double[]{lat, lon};
	}

}
