package test.sphericalGeo;

import org.junit.Test;

import beast.util.Randomizer;
import sphericalGeo.AffineTransformer;
import junit.framework.TestCase;

public class AffineTransformerTest extends TestCase {
	
	@Test
	public void testAffineTransformer() throws Exception {
		for (int i = 0; i < 10; i++) {
			tryAT();
		}
	}
	
	private void tryAT() throws Exception {
		// null transform
		AffineTransformer t = new AffineTransformer(null, null, null);
		singleRun(t);
		
		// translate
		double c0 = Randomizer.nextDouble() * 180 - 90;
		double l0 = Randomizer.nextDouble() * 360 - 180;
		t = new AffineTransformer(c0 + " " + l0, null, null);
		singleRun(t);

		// rotate
		double angle = Randomizer.nextDouble() * Math.PI * 2;
		t = new AffineTransformer(null, angle+"", null);
		singleRun(t);
		
		// translate + rotate
		t = new AffineTransformer(c0 + " " + l0, angle+"", null);
		singleRun(t);

		// scale
		double scale = Randomizer.nextDouble() * 5;
		t = new AffineTransformer(null, null, scale + "");
		singleRun(t);

		// translate + scale
		t = new AffineTransformer(c0 + " " + l0, null, scale + "");
		singleRun(t);

		// translate + rotate + scale
		t = new AffineTransformer(c0 + " " + l0, angle + "", scale + "");
		singleRun(t);
	}

	private void singleRun(AffineTransformer t) {
		for  (int i = 0; i < 100; i++) {
			double lat0 = Randomizer.nextDouble() * 180 - 90;
			double long0 = Randomizer.nextDouble() * 360 - 180;
			double [] proj = t.project(lat0, long0);
			double [] revers = t.projectInverse(proj[0], proj[1]);
			assertEquals(lat0, revers[0], 1e-10);
			assertEquals(long0, revers[1], 1e-10);
		}
	}

}
