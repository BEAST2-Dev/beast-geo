package test.sphericalGeo;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import sphericalGeo.GreatCircleDistance;
import sphericalGeo.SphericalDiffusionModel;
import junit.framework.TestCase;

public class SphericalDiffusionModelTest extends TestCase {

	@Test
	public void testModelIntegratesTo1() throws Exception {
		for (double time : new double[]{0.0001, 0.001, 0.01, 0.1, 1.0}) {
			for (int i = 1; i < 1000; i+= 100) {
				//double sum = integrateAlpha(time, i * 10, 100000, 10e-4);
				double sum = integrateAlpha(time, i, 1000000, 10e-6);
				System.err.println("t=" + time + " p=" + i + " " + sum);
			}
		}
	}
		
	/** integrate over angle from 0 to PI **/
	double integrateAlpha(double time0, double precision0, int NR_OF_STEPS, double accuracy) throws Exception {	
		SphericalDiffusionModel model = new SphericalDiffusionModel();
		RealParameter precision = new RealParameter(precision0 +"");
		model.initByName("precision", precision, "fast", false);
		
		// simple integration by stepping  
		//final int NR_OF_STEPS = 1000000;
		double [] start = new double[2];
		double [] end = new double[2];
		
		double time = time0;
		double sum = 0;
		for (int i = 1; i < NR_OF_STEPS; i++) {
			end[1] = (180.0 * i/NR_OF_STEPS)/180;
			sum += Math.exp(model.getLogLikelihood(null, start, end, time))/NR_OF_STEPS;
		}
		
		// multiply by total length of interval = PI
		sum *= Math.PI/180;
//		assertEquals(1.0, sum, accuracy);
		return sum;
	}
	
	/** integrate over angle from 0 to PI **/
	double expectedDistanceInKiloMeters(double time0, double precision0, int NR_OF_STEPS, double accuracy) throws Exception {	
		SphericalDiffusionModel model = new SphericalDiffusionModel();
		RealParameter precision = new RealParameter(precision0 +"");
		model.initByName("precision", precision, "fast", false);
		
		// simple integration by stepping  
		//final int NR_OF_STEPS = 1000000;
		double [] start = new double[2];
		double [] end = new double[2];
		
		double time = time0;
		double sum = 0;
		for (int i = 1; i < NR_OF_STEPS; i++) {
			end[1] = 180.0 * i/NR_OF_STEPS;
			sum += GreatCircleDistance.EARTHRADIUS * Math.PI * i/NR_OF_STEPS * Math.exp(model.getLogLikelihood(null, start, end, time));
		}
		
		// multiply by total length of interval = PI
		sum *= Math.PI/NR_OF_STEPS;
		//assertEquals(1.0, sum, accuracy);
		return sum;
	}
	
	@Test
	public void testExpectedDistance() throws Exception {
		for (double time : new double[]{0.1}) { //, 2, 20, 50}) {
			for (int i = 1; i < 20; i+= 1) {
				double precision = 50 + i * 50;//Math.pow(2, i - 4);
				//double sum = integrateAlpha(time, i * 10, 100000, 10e-4);
				double sum = expectedDistanceInKiloMeters(time, precision, 1000000, 10e-6);
				System.err.println("t=" + time + " p=" + precision + " " + sum);
			}
		}
	}

}
