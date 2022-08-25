package sphericalGeo;

import beast.base.core.Description;


@Description("Provides landscape transormation to capture some inhomogenuity in the diffusion process")
public interface Transformer {

	double [] project(double lat0, double long0);

	double [] projectInverse(double lat0, double long0);
	
	@Description("Default transformer leaves the map undistorteds")
	public class Default implements Transformer {

		@Override
		public double[] project(double lat0, double long0) {
			return new double[]{lat0, long0};
		}

		@Override
		public double[] projectInverse(double lat0, double long0) {
			return new double[]{lat0, long0};
		}
		
	}
	
}
