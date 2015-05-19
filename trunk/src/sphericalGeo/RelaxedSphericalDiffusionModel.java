/*
 * GreatCircleDiffusionModel.java
 *
 * Copyright (C) Remco Bouckaert
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




import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import org.apache.commons.math.MathException;
import org.apache.commons.math3.util.FastMath;

import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.math.distributions.ParametricDistribution;

@Description("Diffusion model that assumes a normal diffusion process on a sphere and integrates out relaxed clock")
public class RelaxedSphericalDiffusionModel extends SphericalDiffusionModel {
    public Input<ParametricDistribution> rateDistInput = new Input<ParametricDistribution>("distr", "the distribution governing the rates among branches. Must have mean of 1. The precission parameter can be used to change the mean rate.", Input.Validate.REQUIRED);
    public Input<Integer> numberOfDiscreteRatesInput = new Input<Integer>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by", 10);

    int nrOfDiscreteRates;
    private double[] rates;
    private double[] storedRates;
    boolean recompute = true;
    ParametricDistribution distribution;
    
    final static int CACHESIZE = 10000;
    double [][] cache;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        nrOfDiscreteRates = numberOfDiscreteRatesInput.get();
        rates = new double[nrOfDiscreteRates];
        storedRates = new double[nrOfDiscreteRates];
        distribution = rateDistInput.get();
        
        prepare();
        
        Log.warning.println("Setting up cache...");
        cache = new double[180][CACHESIZE];
        for (int i = 0; i < 180; i++) {
        	for (int j = 1; j < CACHESIZE; j++) {
        		cache[i][j] = calc(i * Math.PI / 180.0, j);
        	}
        }
        Log.warning.println("Done with cache");
    }

    // assumes start = {latitude, longitude}
    //         stop  = {latitude, longitude}
    // and -90 < latitude < 90, -180 < longitude < 180
    public double getLogLikelihood(double[] start, double[] stop, double time) {
        if( start[0] == stop[0] && start[1] == stop[1] ) {
            return -1e100;
        }

        if (recompute) {
   	        // this must be synchronized to avoid being called simultaneously by
   	        // two different likelihood threads
   	    	synchronized (this) {
    			prepare();
    			recompute = false;
    		}
    	}
        // assumes start = {latitude, longitude}
        // assumes stop = {latitude, longitude}
        // and -90 < latitude < 90, -180 < longitude < 180

        double latitude1 = start[0];
        double longitude1 = start[1];
        final double DEG2RAD = Math.PI / 180.0;
        final double theta1 = (latitude1) * DEG2RAD;
        if( longitude1 < 0 ) longitude1 += 360;
        //final double phi1 = longitude1 * DEG2RAD;

        double latitude2 = stop[0];
        double longitude2 = stop[1];
        final double theta2 = (latitude2) * DEG2RAD;
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
          angle = (abs(x) > .9 ? acos_parts_fast7(x) : acos(x));
        } else {
           angle = acos(x);
        }

    	double precision0 = precision.getValue(0);
    	double invVariance = precision0 / time;
        double cachedResult = getCachedResult(angle, invVariance);
        if (Double.isFinite(cachedResult)) {
        	return cachedResult;
        }
        return calc(angle, invVariance);
        
    }
    
    double calc(double angle, double invVariance) {
   		// calc contribution for each rate
    	double [] logP = new double[nrOfDiscreteRates];

    	for (int i = 0; i < nrOfDiscreteRates; i++) {
    		logP[i] = getLogLikelihood2(angle, invVariance/rates[i]); // precision0/(time * rates[i]));
    	}
    	
    	// find maximum contribution
    	double max = logP[0];
    	for (double d : logP) {
    		max = Math.max(max, d);
    	}
    	
    	// calc sum of contributions in a numerically stable way
    	double sum = 0;
    	for (double d : logP) {
    		sum += FastMath.exp(d - max);
    	}
    	double logLikelihood = max + FastMath.log(sum)/nrOfDiscreteRates;
    	
    	return logLikelihood;
	}
	
    private double getCachedResult(double angle, double invVariance) {
    	angle = angle * 180.0 / Math.PI;
    	if (angle <= 1.0) {
    		return Double.POSITIVE_INFINITY;
    	}
    	if (invVariance >= CACHESIZE-1) {
    		return Double.POSITIVE_INFINITY;
    	}
    	int i = (int) (angle);
    	double w = angle - i;
    	int j = (int) invVariance;
    	double w2 = invVariance - j;
    	double A = cache[i][j];
    	double B = cache[i+1][j];
    	double C = cache[i][j+1];
    	double D = cache[i+1][j+1];
    	
    	double X = (1-w) * A + w * B; 
    	double Y = (1-w) * C + w * D;
    	
    	double R = (1-w2) * X + w2 * Y;
    	
		return R;
	}

	public double getLogLikelihood2(double angle, double inverseVariance) {
//        if( inverseVariance > 1e20 ) {
//            return -1e100;
//        }
        final double tau = 1.0/inverseVariance;
        double logN = calcLogN(tau);
        final double logP = 0.5 * Math.log(angle * sin(angle)) - Math.log(tau) + -angle * angle / (tau * 2.0);

        //double logN = calcLogN(inverseVariance);
        //final double logP = -angle * angle * inverseVariance / 2.0 + 0.5 * Math.log(angle * sin(angle)) + Math.log(inverseVariance);
        
        return logP - logN;

    }

    
    private void prepare() {
    	for (int i = 0; i < nrOfDiscreteRates; i++) {
            try {
        		rates[i] = distribution.inverseCumulativeProbability((i + 0.5) / nrOfDiscreteRates);
            } catch (MathException e) {
                throw new RuntimeException("Failed to compute inverse cumulative probability!");
            }
    	}
	}

    @Override
    protected boolean requiresRecalculation() {
        recompute = false;
        if (rateDistInput.get().isDirtyCalculation()) {
            recompute = true;
            return true;
        }
        return recompute;
    }

    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;
        super.restore();
    }
}
