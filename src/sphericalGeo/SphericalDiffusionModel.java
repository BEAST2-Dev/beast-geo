/*
 * GreatCircleDiffusionModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package sphericalGeo;



import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

// while FastMath is supposed to be faster, this is not evident in a small timing test.
//import static org.apache.commons.math3.util.FastMath;
import static java.lang.Math.*;

@Description("Diffusion model that assumes a normal diffusion process on a sphere")
public class SphericalDiffusionModel extends SubstitutionModel.Base {

    private static final double RAD2DEG = 180 / Math.PI;
    private static final double DEG2RAD = Math.PI / 180;
    public Input<RealParameter> precisionInput = new Input<RealParameter>("precision", "precision of diffusion process", Validate.REQUIRED);
    public Input<Boolean> m_fast = new Input<>("fast", "Use an approximation for arccos for angles close to 0 "
    + "(|cos(x) > 0.9). In this range the approximation has an error of at most 1e-10, " +
            "and is faster than the Java version.", false);

    public SphericalDiffusionModel() {
		frequenciesInput.setRule(Validate.OPTIONAL);
	}

    RealParameter precision;
    boolean fast = false;

    @Override
    public void initAndValidate() throws Exception {
        precision = precisionInput.get();
        fast = m_fast.get();

        super.initAndValidate();
    }

    public double getLogLikelihood(Node node, double [][] position, double [] branchLengths) {
		double [] start = position[node.getParent().getNr()];
		double [] stop = position[node.getNr()];
		double time = branchLengths[node.getNr()];
		return getLogLikelihood(node, start, stop, time);
    }
    
    // assumes start = {latitude, longitude}
    //         stop  = {latitude, longitude}
    // and -90 < latitude < 90, -180 < longitude < 180
    public double getLogLikelihood(Node node, double [] start, double [] stop, double time) {

    	if (fast) {
    		return getLogLikelihood2(node, start, stop, time);
    	}
           
            if (time <= 1e-20) {
                    return -1e100;
            }
            if (start[0] == stop[0] && start[1] == stop[1]) {
                    return -1e100;
            }
           
           
            double latitude1 = start[0];
            double longitude1 = start[1];
            double theta1 = (latitude1) * Math.PI/180.0;
            if (longitude1 < 0) longitude1 += 360;
            double phi1 = longitude1 * Math.PI/180;

            double latitude2 = stop[0];
            double longitude2 = stop[1];
            double theta2 = (latitude2) * Math.PI/180.0;
            if (longitude2 < 0) longitude2 += 360;
            double phi2 = longitude2 * Math.PI/180;
           
            double Deltalambda = phi2 - phi1;
            
            // See http://en.wikipedia.org/wiki/Great-circle_distance#Formulas
            double angle = Math.acos(Math.sin(theta1) * Math.sin(theta2) + Math.cos(theta1) * Math.cos(theta2) * Math.cos(Deltalambda));

            final double tau = time / precision.getValue(0);
            final double logN = calcLogN(tau, node);
            final double logP = 0.5 * Math.log(angle * sin(angle)) - Math.log(tau) + -angle * angle / (tau * 2.0);

//            double inverseVariance = precision.getValue(0) / time;
//            // See Equation (8) from http://arxiv.org/pdf/1303.1278v1.pdf
//            // logN normalising 'constant'
//            double logP = 0.5 * Math.log(angle * Math.sin(angle)) + 0.5 * Math.log(inverseVariance) -0.5 * angle*angle * inverseVariance;
//            //double logP = - 0.5 * angle*angle * inverseVariance;
//            System.err.println(start[0] + " " + start[1] + " -> " + stop[0] + " " + stop[1] + " => " + logP + " " + angle + " " + 
//            		(0.5 * Math.log(angle * Math.sin(angle))) + " " + ( -0.5 * angle*angle * inverseVariance) + " " + logN);
            return logP - logN;
    }

    
	static double [] in = new double[]{0.0000006743496, 0.0000009536743, 0.000001348699, 0.000001907349, 0.000002697398, 0.000003814697, 0.000005394797, 0.000007629395, 0.00001078959, 0.00001525879, 0.00002157919, 0.00003051758, 0.00004315837, 0.00006103516, 0.00008631675, 0.0001220703, 0.0001726335, 0.0002441406, 0.000345267, 0.0004882812, 0.000690534, 0.0009765625, 0.001381068, 0.001953125, 0.002762136, 0.00390625, 0.005524272, 0.0078125, 0.01104854, 0.015625, 0.02209709, 0.03125, 0.04419417, 0.0625, 0.08838835, 0.125, 0.1767767, 0.25, 0.3535534, 0.5, 0.7071068, 1, 1.414214, 2, 2.828427, 4, 5.656854, 8, 11.31371, 16, 22.62742, 32, 45.25483, 64, 90.50967, 128, 181.0193, 256, 362.0387, 512, 724.0773, 1024, 1448.155, 2048, 2896.309, 4096, 5792.619, 8192, 11585.24, 16384, 23170.48, 32768, 46340.95, 65536, 92681.9, 131072, 185363.8, 262144, 370727.6, 524288, 741455.2, 1048576, 1482910, 2097152, 2965821, 4194304, 5931642, 8388608, 11863280, 16777220, 23726570, 33554430, 47453130, 67108860, 94906270, 134217700, 189812500, 268435500, 379625100, 536870900, 759250100};
	static double [] out = new double[]{0.9999999, 0.9999998, 0.9999998, 0.9999997, 0.9999996, 0.9999994, 0.9999991, 0.9999987, 0.9999982, 0.9999975, 0.9999964, 0.9999949, 0.9999928, 0.9999898, 0.9999856, 0.9999797, 0.9999712, 0.9999593, 0.9999425, 0.9999186, 0.9998849, 0.9998372, 0.9997698, 0.9996745, 0.9995397, 0.999349, 0.9990795, 0.9986983, 0.9981593, 0.9973972, 0.9963198, 0.994797, 0.992645, 0.9896045, 0.9853106, 0.9792494, 0.9706989, 0.9586452, 0.9416613, 0.9177062, 0.8837365, 0.8353602, 0.7681435, 0.6813698, 0.5805998, 0.4757169, 0.3765316, 0.2896687, 0.2179294, 0.1612042, 0.1177539, 0.08522718, 0.06127549, 0.0438445, 0.0312647, 0.02223985, 0.01579263, 0.01120063, 0.007936867, 0.005620646, 0.003978664, 0.002815455, 0.001991886, 0.001409006, 0.0009965826, 0.0007048228, 0.0004984513, 0.0003524914, 0.0002492657, 0.0001762657, 0.0001246428, 0.00008813786, 0.00006232392, 0.00004407018, 0.00003116258, 0.0000220354, 0.00001558145, 0.00001101778, 0.000007790763, 0.00000550891, 0.000003895391, 0.00000275446, 0.000001947698, 0.000001377231, 0.0000009738497, 0.0000006886158, 0.000000486925, 0.000000344308, 0.0000002434625, 0.000000172154, 0.0000001217313, 0.00000008607701, 0.00000006086564, 0.00000004303851, 0.00000003043282, 0.00000002151925, 0.00000001521641, 0.00000001075963, 0.000000007608205, 0.000000005379814, 0.000000003804103};
	static double [] logout;
	static double minin, maxin, logminin;
	static {
		logout = new double[out.length];
		for (int i = 0; i < out.length; i++) {
			logout[i] = Math.log(out[i]);
		}
		minin = in[0];
		logminin = FastMath.log(minin);
		maxin = in[in.length - 1];
	}
	
    private double [] lastTau;// = -1;
    private double [] lastN;// = 0;
    
    double calcLogN(double tau, Node node) {
    	if (node != null) {
    		// check the per node cache
	    	if (lastTau == null) {
	    		lastTau = new double[node.getTree().getNodeCount()];
	    		lastN = new double[lastTau.length];
	    	}
	    	if (tau == lastTau[node.getNr()]) {
	    		return lastN[node.getNr()];
	    	}
			lastTau[node.getNr()] = tau;
    	}

    	// try approximation of pieces of the curve by polynomials
    	double logN;
        if( 1e-6 <= tau && tau < 1e-1 ) {
            final double a = -0.0087466200605258813, b = -0.16664983515519927, c = -1.4116091739448584e-07;
            logN = (a * tau + b) * tau + c;
        } else if( 1e-1 <= tau && tau <= 1 ) {
        	final double a = -0.017361715075749726, b = -0.1617395191006793, c = -0.00064481343766429537;
            logN = (a * tau + b) * tau + c;
        } else if( tau < 1e-6 ) {
        	logN = 0;
        } else if (tau < 10) {
	        // this is terrible for tau>10 but those kinds of values should not matter
        	final double a = -0.14247308312783535, b = 1.2366027188629642, c= -1.7402299884011929;
        	final double lt = FastMath.log(tau);
	        logN = -FastMath.exp(lt * (a*lt + b) +c);
        } else { // tau >= 10
        	// fall back on original approximation
        	logN  = calcLogN2(tau);
        }

        if (node != null) {
        	lastN[node.getNr()] = logN;
        }
        return logN;
    }

    
    /** approximation described in the GEO_SPHERE paper **/
    double calcLogN2(double tau) {

    	if (tau < minin) {
    		return 0;
    	}
    	if (tau > maxin) {
    		return FastMath.log(2.888266/tau);
    	}
    	
    	int i = (int) ((FastMath.log(tau)-logminin)/0.3465735);
    	if (i == in.length-1) {
    		return logout[i];
    	}
//      int i = Arrays.binarySearch(in, tau);
//    	if (i >= 0) {
//    		// found exact match
//    		return logout[i];
//    	}
//    	i = - i - 2;
    	// interpolation needed to reconstruct trans prob
    	//double fiTime = CACHE_SIZE * fTime / m_fMaxTime; 
    	//int iTime =  (int) fiTime;
    	double logout1 = logout[i];
    	double logout2 = logout[i + 1];
    	double fWeight2 = (tau - in[i])/(in[i+1] - in[i]);
    	double fWeight1 = 1.0 - fWeight2;
		// TODO Auto-generated method stub
		double logN = logout1 * fWeight1 + logout2 * fWeight2;
    	return logN;
	}

    
//        static double log0(double x) {
//            double y = (x-1.)/(1+x) ;
//            // making this 2*y reduces accuracy to 1e-8 over [1e-8,2]
//            return 2*(y*(1 + y*y/3));
//        }
//        static int N = 100;
//
//        static public double log(double x) {
//            int v = 0;
//            while( x < 1 ) {
//                x *= Math.E;
//                v -= 1;
//            }
//            while ( x > Math.E ) {
//                x /= Math.E;
//                v += 1;
//            }
//            // assert 1 <= x <= math.e
//            int i = (int)(N*(x - 1)/(Math.E-1));
//            return v + tab1[i] + log0(x * tab0[i]);
//        }
//
//        static double[] tab0 = new double[N+1];
//        static double[] tab1 = new double[N+1];
//
//        static {
//            for(int i = 0; i < 101; ++i) {
//                tab0[i] = 1/(1 + (1./N) * i * (Math.E - 1));
//                tab1[i] = -Math.log(tab0[i]);
//            }
//        }
//    
//    public static double logBadApprox(double x) {
//    	System.err.println(x);
//        return 6 * (x - 1) / (x + 1 + 4 * (Math.sqrt(x)));
//    }
    
	//@Override
    public double getLogLikelihood2(Node node, double[] start, double[] stop, double time) {

        if( time <= 1e-20 ) {
            return -1e100;
        }
        if( start[0] == stop[0] && start[1] == stop[1] ) {
            return -1e100;
        }

        // assumes start = {latitude, longitude}
        // assumes stop = {latitude, longitude}
        // and -90 < latitude < 90, -180 < longitude < 180

        double latitude1 = start[0];
        double longitude1 = start[1];
        //final double DEG2RAD = Math.PI / 180.0;
        final double theta1 = latitude1 * DEG2RAD;
        if( longitude1 < 0 ) longitude1 += 360;
        //final double phi1 = longitude1 * DEG2RAD;

        double latitude2 = stop[0];
        double longitude2 = stop[1];
        final double theta2 = latitude2 * DEG2RAD;
        if( longitude2 < 0 ) longitude2 += 360;
        //final double phi2 = longitude2 * DEG2RAD;

        final double deltaLambda = (longitude2 - longitude1) * DEG2RAD; //phi2 - phi1, in radians;

        // Use trigonometric equalities to reduce cost of computing both sin(x)*sin(y) and cos(x)*cos(y)
        // to two cos() calls and 4 +/- and one '/2'.
        final double cosplus = cos(theta1 + theta2);
        final double cosminus = cos(theta1 - theta2);
        final double twicecoscos = (cosminus + cosplus);
        final double twicesinsin = (cosminus - cosplus);
        final double x = (twicesinsin + twicecoscos * cos(deltaLambda)) / 2;
//        final double x = FastMath.sin(theta1) * FastMath.sin(theta2) +
//                FastMath.cos(theta1) * FastMath.cos(theta2) * FastMath.cos(deltaLambda);
        final double angle;
        if( fast ) {
          angle = (abs(x) > .9 ? acos_parts_fast7(x) : FastMath.acos(x));
        } else {
          angle = acos(x);
        }

//        final double inverseVariance = precision.getValue(0) / time;
        final double tau = time/precision.getValue(0);
        final double logN = calcLogN(tau, node);
        //final double logP = -angle * angle * inverseVariance / 2.0 + 0.5 * Math.log(angle * sin(angle)) + Math.log(inverseVariance);

        //final double logP = Math.log(Math.sqrt(angle * sin(angle)) / tau * exp(-angle * angle / (tau * 2.0));
        //                  = Math.log(Math.sqrt(angle * sin(angle)) / tau) + log(exp(-angle * angle / (tau * 2.0));
        //                  = Math.log(Math.sqrt(angle * sin(angle)) - Math.log(tau) + -angle * angle / (tau * 2.0);
        //                  = 0.5 * Math.log(angle * sin(angle)) - Math.log(tau) + -angle * angle / (tau * 2.0);
        //                  = 0.5 * Math.log(angle * sin(angle)) + Math.log(inverseVariance) + -angle * angle * inverseVariance / 2.0;
        final double logP = 0.5 * Math.log(angle * sin(angle)) - Math.log(tau) + -angle * angle / (tau * 2.0);

        
//		System.err.println(start[0] + " " + start[1] + " -> " + stop[0] + " " + stop[1] + " => " + logP);
        return logP - logN;

    }

	final static int NR_OF_STEPS = 1000;
	
	public double[] sample(double[] start, double time, double precision) throws ConvergenceException, FunctionEvaluationException, IllegalArgumentException {
		return sample(start, time, precision, new double[]{0.0, Math.PI * 2.0});
	}
	
	/** anglerange[0] = lower bound, anglerange[1] = upper bound of range for angle2
	 * **/
	public double[] sample(double[] start, double time, double precision, double [] angleRange)
            throws ConvergenceException, FunctionEvaluationException, IllegalArgumentException {
		
		// first, sample an angle from the spherical diffusion density
		final double inverseVariance = precision / time;

		UnivariateRealFunction function = new UnivariateRealFunction() {
			@Override
			public double value(double x) throws FunctionEvaluationException {
				double logR = -x*x * inverseVariance /2.0 + 0.5 * Math.log(x * Math.sin(x) * inverseVariance);

				return Math.exp(logR);
			}
		};
		
        UnivariateRealIntegrator integrator = new TrapezoidIntegrator();
        integrator.setAbsoluteAccuracy(1.0e-10);
        integrator.setRelativeAccuracy(1.0e-14);
        integrator.setMinimalIterationCount(2);
        integrator.setMaximalIterationCount(25);
		
		double [] cumulative = new double[NR_OF_STEPS];
		for (int i = 1; i < NR_OF_STEPS; i++) {
			cumulative[i] = cumulative[i - 1] + integrator.integrate(function, (i - 1) * Math.PI / NR_OF_STEPS, i * Math.PI / NR_OF_STEPS);
		}
		
		// normalise 
		double sum = cumulative[NR_OF_STEPS - 1];
		for (int i = 0; i < NR_OF_STEPS; i++) {
			cumulative[i] /= sum;
		}

		int i = Randomizer.randomChoice(cumulative);
		double angle = i* Math.PI / NR_OF_STEPS;
		
		
		// now we have an angle, use this to rotate the point [0,0] over
		// this angle in a random direction angle2
		double angle2 = angleRange[0] + Randomizer.nextDouble() * (angleRange[1] - angleRange[0]);
		//angleRange[0] = angle2;

	    double [] xC = new double[] {Math.cos(angle), Math.sin(angle)*Math.cos(angle2), Math.sin(angle)*Math.sin(angle2)};

	    // convert back to latitude, longitude relative to (lat=0, long=0)
		double [] xL = cartesian2Sperical(xC, true);
	    //double [] sC = spherical2Cartesian(start[0], start[1]);
		double [] position = reverseMap(xL[0], xL[1], start[0], start[1]);
		return position;
	}

    public static double getAngle(double[] start, double[] stop) {
        double latitude1 = start[0];
        double longitude1 = start[1];
        double theta1 = (latitude1) * Math.PI/180.0;
        if (longitude1 < 0) longitude1 += 360;
        double phi1 = longitude1 * Math.PI/180;

        double latitude2 = stop[0];
        double longitude2 = stop[1];
        double theta2 = (latitude2) * Math.PI/180.0;
        if (longitude2 < 0) longitude2 += 360;
        double phi2 = longitude2 * Math.PI/180;

        double Deltalambda = phi2 - phi1;

        // See http://en.wikipedia.org/wiki/Great-circle_distance#Formulas
        double angle = Math.acos(Math.sin(theta1) * Math.sin(theta2) + Math.cos(theta1) * Math.cos(theta2) * Math.cos(Deltalambda));
        return angle;
    }

	/** Convert spherical coordinates (latitude,longitude) in degrees on unit sphere 
	 * to Cartesian (x,y,z) coordinates **/
	public static double [] spherical2Cartesian(double fLat, double fLong) {
		double fPhi = fLong * DEG2RAD;
		double fTheta = (90 - fLat) * DEG2RAD;
		//double fTheta = (fLat) * Math.PI / 180.0;
	    //{x}=\rho \, \sin\theta \, \cos\phi  
	    //{y}=\rho \, \sin\theta \, \sin\phi  
	    //{z}=\rho \, \cos\theta 
		double [] fNorm = new double[3];
        final double sinTheta = FastMath.sin(fTheta);
        fNorm[0] = sinTheta * FastMath.cos(fPhi);
		fNorm[1] = sinTheta * FastMath.sin(fPhi);
		fNorm[2] = FastMath.cos(fTheta);
//		fNorm[0] = Math.sin(fTheta) * Math.cos(fPhi);
//		fNorm[1] = Math.sin(fTheta) * Math.sin(fPhi);
//		fNorm[2] = Math.cos(fTheta);
		return fNorm;
	} // spherical2Cartesian

	/** inverse of spherical2Cartesian **/
	public static double [] cartesian2Sperical(double[] f3dRotated2, boolean fastApprox) {
        final double v = usemyat ? at2(f3dRotated2[1], f3dRotated2[0], 30) : fast_atan2(f3dRotated2[1], f3dRotated2[0]);
        //assert abs(v - v1) < 1e-5;
        return 	new double[]{
				//Math.acos(-f3dRotated2[2]) * RAD2DEG - 90,
				!fastApprox ? FastMath.acos(-f3dRotated2[2]) * RAD2DEG - 90 :
				             myacos(-f3dRotated2[2]) * RAD2DEG - 90, // <- faster but considerably less accurate
				//Math.atan2(f3dRotated2[1], f3dRotated2[0]) * 180.0/Math.PI
				//FastMath.atan2(f3dRotated2[1], f3dRotated2[0]) * 180.0/Math.PI
				v * 180.0/Math.PI
				
		};
	}

	// from http://stackoverflow.com/questions/523531/fast-transcendent-trigonometric-functions-for-java
	/** Fast approximation of 1.0 / sqrt(x).
	   * See <a href="http://www.beyond3d.com/content/articles/8/">http://www.beyond3d.com/content/articles/8/</a>
	   * @param x Positive value to estimate inverse of square root of
	   * @return Approximately 1.0 / sqrt(x)
	   **/
	  public static double
	  invSqrt(double x)
	  {
	    double xhalf = 0.5 * x; 
	    long i = Double.doubleToRawLongBits(x);
	    i = 0x5FE6EB50C7B537AAL - (i>>1); 
	    x = Double.longBitsToDouble(i);
	    x = x * (1.5 - xhalf*x*x); 
	    return x; 
	  }

	  /** Approximation of arctangent.
	   *  Slightly faster and substantially less accurate than
	   *  {@link Math#atan2(double, double)}.
	   **/
	  public static double fast_atan2(double y, double x)
	  {
	    double d2 = x*x + y*y;

	    // Bail out if d2 is NaN, zero or subnormal
	    if (Double.isNaN(d2) ||
	        (Double.doubleToRawLongBits(d2) < 0x10000000000000L))
	    {
	      return Double.NaN;
	    }

	    // Normalise such that 0.0 <= y <= x
	    boolean negY = y < 0.0;
	    if (negY) {y = -y;}
	    boolean negX = x < 0.0;
	    if (negX) {x = -x;}
	    boolean steep = y > x;
	    if (steep)
	    {
	      double t = x;
	      x = y;
	      y = t;
	    }

	    // Scale to unit circle (0.0 <= y <= x <= 1.0)
	    double rinv = invSqrt(d2); // rinv ≅ 1.0 / hypot(x, y)
	    x *= rinv; // x ≅ cos θ
	    y *= rinv; // y ≅ sin θ, hence θ ≅ asin y

	    // Hack: we want: ind = floor(y * 256)
	    // We deliberately force truncation by adding floating-point numbers whose
	    // exponents differ greatly.  The FPU will right-shift y to match exponents,
	    // dropping all but the first 9 significant bits, which become the 9 LSBs
	    // of the resulting mantissa.
	    // Inspired by a similar piece of C code at
	    // http://www.shellandslate.com/computermath101.html
	    double yp = FRAC_BIAS + y;
	    int ind = (int) Double.doubleToRawLongBits(yp);

	    // Find φ (a first approximation of θ) from the LUT
	    double φ = ASIN_TAB[ind];
	    double cφ = COS_TAB[ind]; // cos(φ)

	    // sin(φ) == ind / 256.0
	    // Note that sφ is truncated, hence not identical to y.
	    double sφ = yp - FRAC_BIAS;
	    double sd = y * cφ - x * sφ; // sin(θ-φ) ≡ sinθ cosφ - cosθ sinφ

	    // asin(sd) ≅ sd + ⅙sd³ (from first 2 terms of Maclaurin series)
	    double d = (6.0 + sd * sd) * sd * ONE_SIXTH;
	    double θ = φ + d;

	    // Translate back to correct octant
	    if (steep) { θ = Math.PI * 0.5 - θ; }
	    if (negX) { θ = Math.PI - θ; }
	    if (negY) { θ = -θ; }

	    return θ;
	  }

	  private static final double ONE_SIXTH = 1.0 / 6.0;
	  private static final int FRAC_EXP = 8; // LUT precision == 2 ** -8 == 1/256
	  private static final int LUT_SIZE = (1 << FRAC_EXP) + 1;
	  private static final double FRAC_BIAS =
	    Double.longBitsToDouble((0x433L - FRAC_EXP) << 52);
	  private static final double[] ASIN_TAB = new double[LUT_SIZE];
	  private static final double[] COS_TAB = new double[LUT_SIZE];

	  static
	  {
	    /* Populate trig tables */
	    for (int ind = 0; ind < LUT_SIZE; ++ ind)
	    {
	      double v = ind / (double) (1 << FRAC_EXP);
	      double asinv = Math.asin(v);
	      COS_TAB[ind] = Math.cos(asinv);
	      ASIN_TAB[ind] = asinv;
	    }
	  }
	
	
	
	
	
	public static double [] reverseMap(double fLat, double fLong, double fLatT, double fLongT) {
		// from spherical to Cartesian coordinates
		double [] f3DPoint = spherical2Cartesian(fLat, fLong);

		// double [] p = cartesian2Sperical(f3DPoint);
		// rotate, first latitude, then longitude
		double [] f3DRotated = new double[3];
		double fC = Math.cos(fLongT * DEG2RAD);
		double fS = Math.sin(fLongT * DEG2RAD);
		double [] f3DRotated2 = new double[3];
		double fC2 = Math.cos(-fLatT * DEG2RAD);
		double fS2 = Math.sin(-fLatT * DEG2RAD);

		// rotate over latitude
		f3DRotated[0] = f3DPoint[0] * fC2 + f3DPoint[2] * fS2;
		f3DRotated[1] = f3DPoint[1];
		f3DRotated[2] = -f3DPoint[0] * fS2 + f3DPoint[2] * fC2;

		// rotate over longitude
		f3DRotated2[0] = f3DRotated[0] * fC - f3DRotated[1] * fS; 
		f3DRotated2[1] = f3DRotated[0] * fS + f3DRotated[1] * fC; 
		f3DRotated2[2] = f3DRotated[2]; 

		double [] point = cartesian2Sperical(f3DRotated2, true);
		//System.err.println(Arrays.toString(point) + " " + Arrays.toString(f3DRotated2) + " " + Arrays.toString(f3DRotated));
		return point;
	} // map
	
	/** convert spherical coordinates (latitude,longitude) to
	 * 2D point in interval [-180, -90] x [180, 90] 
	 * wrt plane defined by fNorm
	 * 
	 * http://en.wikipedia.org/wiki/Sinusoidal_projection
	 *
	 * Alternatives:
	 * http://en.wikipedia.org/wiki/Hammer_projection <= looks very promissing, comes with inverse
	 * http://en.wikipedia.org/wiki/Category:Equal-area_projections
	 */
	public static double [] map(double fLat, double fLong, double fLatT, double fLongT) {
		// from spherical to Cartesian coordinates
		double [] f3DPoint = spherical2Cartesian(fLat, fLong);
		// rotate, first longitude, then latitude
		double [] f3DRotated = new double[3];
		double fC = Math.cos(-fLongT * DEG2RAD);
		double fS = Math.sin(-fLongT * DEG2RAD);
		double [] f3DRotated2 = new double[3];
		double fC2 = Math.cos(fLatT * DEG2RAD);
		double fS2 = Math.sin(fLatT * DEG2RAD);

		// rotate over longitude
		f3DRotated[0] = f3DPoint[0] * fC - f3DPoint[1] * fS; 
		f3DRotated[1] = f3DPoint[0] * fS + f3DPoint[1] * fC; 
		f3DRotated[2] = f3DPoint[2]; 

		// rotate over latitude
//		f3DRotated2 = f3DRotated;
		f3DRotated2[0] = f3DRotated[0] * fC2 + f3DRotated[2] * fS2; 
		f3DRotated2[1] = f3DRotated[1];
		f3DRotated2[2] = -f3DRotated[0] * fS2 + f3DRotated[2] * fC2;
		
		System.err.println(Arrays.toString(f3DPoint) + " " + Arrays.toString(f3DRotated) + " " + Arrays.toString(f3DRotated2));
		// translate back to (longitude, latitude)
		double [] point = cartesian2Sperical(f3DRotated2, true);
		return point;
	} // map
	
	
    //  first 7 terms from the classic expansion
    // coefficients for series expansion
    static final double c_3 = 1.0 / 6.0;
    static final double c_5 = 3.0 / 40.0;
    static final double c_7 = 5.0 / 112.0;
    static final double c_9 = 35.0 / 1152.0;
    static final double c_11 = 63.0 / 2816.0;
    static final double c_13 = 231.0 / 13312;

    static final double halfpi = Math.PI / 2;

    private static double th = sqrt(2.0)/2;

    static double acos_fast7(double x) {
        assert -1 <= x && x <= 1;
        final double u2 = x * x;
        final double u3 = u2 * x;
        final double u5 = u3 * u2;
        final double u7 = u5 * u2;
        final double u9 = u7 * u2;
        final double u11 = u9 * u2;
        final double u13 = u11 * u2;

        return halfpi - x - c_3 * u3 - c_5 * u5 - c_7 * u7 - c_9 * u9 - c_11 * u11 - c_13 * u13;
    }
    
    static public double acos_parts_fast7(double x) {
        if ( x > th ) {
            return (halfpi - acos_fast7(sqrt(1 - x * x)));
        }
        else if ( x < -th ) {
            return halfpi +  acos_fast7(sqrt(1 - x * x));
        }
        else {
            return acos_fast7(x);
        }
    }
    static private double myacos(double x) {
		// This can happen if calling code is not carefull. Hack at the moment
		if( Double.isInfinite(x) || Double.isNaN(x) ) {
			return x;
		}
        assert -1 <= x && x <= 1;
		return acos_parts_fast7(x);
    }

   static final int NN = 100;
   static private double[] phi = new double[NN]; // [2.**(-k) for k in range(1,n+1)]
   static private double [] ap = new double[NN]; // [math.atan(phi[n]) for n in range(len(phi))]

   static double at2(double y, double x, int N) {
       //#ap = [math.atan(x) for x in phi]
       //#v = [cos(0.25), sin(0.25)]'; % Starting vector (input)
       boolean r1 = false;
       if( x < 0 ) {
           x = -x;
           r1 = true;
       }
       boolean r2 = false;
       if( y < 0 ) {
           y = -y;
           r2 = true;
       }
       boolean r = false;
       if( y > x ) {
           final double t = x;
           x = y;
           y = t;
           r = true;
       }
       double a = 0;
       for (int n = 0; n < N; ++n) {
           //final double s = y > 0 ? 1 : -1;
           final double p = y > 0 ? phi[n] : -phi[n];
           a += y > 0 ? ap[n] : -ap[n];

           final double tx = x + p * y;
           y -= p * x;
           x = tx;
           //System.err.println("a " + a + " " + x + " " + y);
       }
       if( r ) {
           a = Math.PI / 2 - a;
       }
       if( r1 ) {
           a = Math.PI - a;
       }
       if( r2 ) {
           a = -a;
       }
       return a;
   }

	@Override
	public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
		// TODO Auto-generated method stub
	}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		return null;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return false;
	}
	
	static boolean usemyacos = false;
    static boolean usemyat = false;

    // Args - N seed , 0/1 (0 fastmath/ 1 jhacos), 0/1 (0 fastmath/ 1 jhatn) 0/1 (0 gen rands inline / 1 pre allocate)
	public static void main(String[] args) {
        for(int k = 1; k < NN+1; ++k) {
            phi[k-1] = Math.pow(2.0, -k);
            ap[k-1] = Math.atan(phi[k-1]);
        }

		// speed test
        final int N = args.length > 0 ? Integer.parseInt(args[0]) : 10000000;
        final int seed = args.length > 1 ? Integer.parseInt(args[1]) : 123;
		Randomizer.setSeed(seed);
        if( args.length >  2 ) {
            usemyacos = (Integer.parseInt(args[2]) != 0);
        }
        if( args.length >  3 ) {
            usemyat = (Integer.parseInt(args[3]) != 0);
        }
        boolean pre = true;
        if( args.length >  4 ) {
            pre = (Integer.parseInt(args[4]) != 0);
        }

        long start ;
        double a0 = 0, a1 = 1;
        double x0 = 0, x1 = 1;
        if( pre ) {
            double[][] point = new double[N][2];
            for (int i = 0; i < N; i++) {
                point[i][0] = Randomizer.nextDouble() * 180 - 90;
                point[i][1] = Randomizer.nextDouble() * 360 - 180;
            }

            System.err.println("Starting...");

            start = System.currentTimeMillis();

            for (int i = 0; i < N; i++) {
                double[] cart = SphericalDiffusionModel.spherical2Cartesian(point[i][0], point[i][1]);
                double[] sper = SphericalDiffusionModel.cartesian2Sperical(cart, usemyacos);
                x0 += sper[0];
                x1 += sper[1];
                a0 += point[i][0];
                a1 += point[i][1];
            }
        } else {
            double[] point = new double[2];
            System.err.println("Starting...");

            start = System.currentTimeMillis();

            for (int i = 0; i < N; i++) {
                point[0] = Randomizer.nextDouble() * 180 - 90;
                point[1] = Randomizer.nextDouble() * 360 - 180;
                double[] cart = SphericalDiffusionModel.spherical2Cartesian(point[0], point[1]);
                double[] sper = SphericalDiffusionModel.cartesian2Sperical(cart, usemyacos);
                x0 += sper[0];
                x1 += sper[1];
                a0 += point[0];
                a1 += point[1];
            }
        }

        long end = System.currentTimeMillis();
        System.err.println("N : " + N + " seed: " + seed + (usemyacos ? " jhacos" : " fastmath") + (usemyat ? " jhatan2" : " atanapprox"));
		System.err.println("Expeted sum: " + a0 + " " + a1);
		System.err.println("Calculated : " + x0 + " " + x1);
        System.err.println("Diff : " + Math.abs(a0-x0)/Math.max(a0,x0) + " " + Math.abs(a1-x1)/Math.max(a1,x1));

		System.err.println("Runtime: " + ((end - start)/1000.0) + " seconds");
	}
}
