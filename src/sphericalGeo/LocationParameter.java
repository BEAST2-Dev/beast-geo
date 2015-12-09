package sphericalGeo;

import beast.core.parameter.RealParameter;

public class LocationParameter extends RealParameter {

	public void setValueSilently(Double[] d) {
		System.arraycopy(d, 0, values, 0, d.length);
	}

}
