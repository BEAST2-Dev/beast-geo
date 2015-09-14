package sphericalGeo;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

@Description("Transform latitude/longitude through an affine transform. "
		+ "Rotates around a center point, then scales latitude and longitude")
public class AffineTransformer extends BEASTObject implements Transformer {
	public Input<RealParameter> centerInput = new Input<>("center", "center of transformation (default (0,0))");
	public Input<RealParameter> angleInput = new Input<>("angle", "angle used for rotation around center (default 0)");
	public Input<RealParameter> scaleInput = new Input<>("scale", "scale input, multiplies latitude, divides longitude (default 1)");
	
	RealParameter center;
	RealParameter angle;
	RealParameter scale;
	
	public AffineTransformer(String center, String angle, String scale) throws Exception {
		centerInput.setValue(center, this);
		angleInput.setValue(angle, this);
		scaleInput.setValue(scale, this);
		initAndValidate();
	}
	
	
	@Override
	public void initAndValidate() throws Exception {
		center = centerInput.get();
		if (center != null && center.getDimension() != 2) {
			throw new Exception("Center should have a dimension of 2");
		}
		angle = angleInput.get();
		if (angle != null && angle.getDimension() != 1) {
			throw new Exception("angle should have a dimension of 1");
		}
		scale = scaleInput.get();
		if (scale != null && scale.getDimension() != 1) {
			throw new Exception("scale should have a dimension of 1");
		}
	}

	@Override
	public double[] project(double lat0, double long0) {
		if (center != null) {
			lat0 += center.getValue(0);
			long0 += center.getValue(1);
		}
		if (angle != null) {
			double a = Math.cos(angle.getValue());
			double b = Math.sin(angle.getValue());
			double lat1 = a * lat0 - b * long0;
			double long1 = b * lat0 + a * long0;
			lat0 = lat1;
			long0 = long1;
		}
		if (scale != null) {
			lat0 *= scale.getValue();
			long0 /= scale.getValue();
		}
		return new double[]{lat0, long0};
	}

	@Override
	public double[] projectInverse(double lat0, double long0) {
		if (scale != null) {
			lat0 /= scale.getValue();
			long0 *= scale.getValue();
		}
		if (angle != null) {
			double a = Math.cos(-angle.getValue());
			double b = Math.sin(-angle.getValue());
			double lat1 = a * lat0 - b * long0;
			double long1 = b * lat0 + a * long0;
			lat0 = lat1;
			long0 = long1;
		}
		if (center != null) {
			lat0 -= center.getValue(0);
			long0 -= center.getValue(1);
		}
		return new double[]{lat0, long0};
	}
	
}
