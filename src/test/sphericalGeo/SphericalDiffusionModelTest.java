package test.sphericalGeo;

import org.junit.Test;

import beast.core.parameter.RealParameter;
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
			end[1] = 180.0 * i/NR_OF_STEPS;
			sum += Math.exp(model.getLogLikelihood(start, end, time))/NR_OF_STEPS;
		}
		
		// multiply by total length of interval = PI
		sum *= Math.PI;
		assertEquals(1.0, sum, accuracy);
		return sum;
	}
	
}
