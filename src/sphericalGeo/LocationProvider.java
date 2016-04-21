package sphericalGeo;

import beast.core.Function;

public interface LocationProvider extends Function {
	
//    @Override
//	default public int getMinorDimension1() {return 2;}
	
	public double[] getPosition(int iDim);
	
	@Override
	default double getArrayValue() {
		return 0;
	}

	@Override
	default double getArrayValue(int i) {
		return getPosition(i/2)[i%2];
	}
}
