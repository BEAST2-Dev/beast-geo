package sphericalGeo;

import beast.base.inference.parameter.RealParameter;

public class LocationParameter extends RealParameter {

	public void setValueSilently(Double[] d) {
		System.arraycopy(d, 0, values, 0, d.length);
	}

}
