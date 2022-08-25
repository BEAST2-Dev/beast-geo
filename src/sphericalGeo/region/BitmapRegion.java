package sphericalGeo.region;

import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;

@Description("Defines geographical region through a bitmap")
public class BitmapRegion extends Region {
	public Input<File> bitmapFileInput = new Input<File>("file", "location of bitmap file on file system", Validate.REQUIRED);
	public Input<String> bboxInput = new Input<>("bbox", "bounding box covered by bitmap file, "
			+ "It is specified as space separated list (in min-latitude min-longitude max-latitude max-longitude), "
            + " e.g. \"-47 120 -5 154\" to cover Australia. "
			+ "default is -90 -180 90 180");
	
	public Input<Integer> regionColorInput = new Input<>("color", "color in bitmap that specifies the region", 0x00FF00);

	@Override
	public void initAndValidate() {
		regionColor = regionColorInput.get();
		parseBBox();
		File file = bitmapFileInput.get();
		if (!file.exists()) {
			throw new RuntimeException("Cannot find file " + file.getName());
		}
		try {
			image = ImageIO.read(file);
		} catch (IOException e) {
			throw new IllegalArgumentException(e);
		}
		width = image.getWidth();
		height = image.getHeight();
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
                throw new RuntimeException("bbox input first latitude must be smaller than second latitude");
        }
        if (minLong >= maxLong) {
                throw new RuntimeException("bbox input first longitude must be smaller than second longitude");
        }
    }

}
