package sphericalGeo;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

public class LocationOperator2 extends Operator {
	public Input<RealParameter> locationInput = new Input<>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree", Validate.REQUIRED);
	public Input<ApproxMultivariateTraitLikelihoodF> likelihoodInput = new Input<>("likelihood", 
			"likelihood over the locations", Validate.REQUIRED);
//	public Input<Operator> operatorInput = new Input<>("operator", "(tree) operator to user before changing the locations");

	RealParameter sampledLocations;
	int [] taxonNrs;
	ApproxMultivariateTraitLikelihoodF likelihood;
	TreeInterface tree;
	
	double [][] position;
	double [][] sphereposition;
	boolean [] isSampled;
	double precision;
	double [] parentweight;
	double [] branchLengths;
//	Operator operator;
	
	double stdDevLat;
	double stdDevLong;

	@Override
	public void initAndValidate() throws Exception {
		sampledLocations = locationInput.get();
		likelihood = likelihoodInput.get();
		tree = likelihood.tree;
		position = new double[tree.getNodeCount()][2];
		sphereposition = new double[position.length][3];
		isSampled = new boolean[tree.getNodeCount()];
		branchLengths = new double[tree.getNodeCount()];
		parentweight = new double[tree.getNodeCount()];
//		operator = operatorInput.get();
	}
	
	@Override
	public double proposal() {
		double logHR = 0;
//		if (operator != null) {
//			logHR = operator.proposal();
//			if (logHR == Double.NEGATIVE_INFINITY) {
//				return logHR;
//			}
//		}
		
		taxonNrs = likelihood.taxonNrs;
		if (taxonNrs == null || taxonNrs.length == 0) {
			return logHR;
		}
		
//		if (true || Randomizer.nextDouble() < 0.9) {
//			// use random walk
//			MultiGeoPrior multiGeoPrior = likelihood.multiGeopriorsInput.get();
//			List<GeoPrior> geopriors = multiGeoPrior.geopriorsInput.get();
//			int i = Randomizer.nextInt(geopriors.size());
//
//			int k = multiGeoPrior.getCladeTopNodeNr(i);
//			double window = 2;
//			sampledLocations.setValue(k * 2, sampledLocations.getValue(k*2) + window * (Randomizer.nextDouble() - 0.5));
//			sampledLocations.setValue(k * 2 + 1, sampledLocations.getValue(k*2+1) + window * (Randomizer.nextDouble() - 0.5));
//			return logHR;
//		}
		
		
		
		MultiGeoPrior multiGeoPrior = likelihood.multiGeopriorsInput.get();
		List<GeoPrior> geopriors = multiGeoPrior.geopriorsInput.get();
		int i = Randomizer.nextInt(geopriors.size());
		GeoPrior prior = geopriors.get(i);
		if (prior.allInternalNodes) {
			return Double.NEGATIVE_INFINITY;
		}
		if (!prior.isMonophyletic) {
			prior.initialise();
			if (!prior.isMonophyletic) {
				// fail -- we should never move a non-monophyletic clade
				return Double.NEGATIVE_INFINITY;
			}
		}
		
		List<Integer> sampleNumber = new ArrayList<>();
		for (GeoPrior geoprior : geopriors) {
			if (geoprior.region != null || Randomizer.nextDouble() < 0.25) {
				int taxonNr = geoprior.getTaxonNr();
				sampleNumber.add(taxonNr);
				isSampled[geoprior.getTaxonNr()] = true;
			}
		}
		precision = likelihood.substModel.precisionInput.get().getValue();
		calcBranchLengths();
		calcPositions(sampleNumber, precision);
		

		calcStdDevs();
		double forwardStdDevLat = stdDevLat;
		double forwardStdDevLong = stdDevLong;
		Double [] orgValues = sampledLocations.getValues();
		
		for (GeoPrior geoprior : geopriors) {
			if (geoprior.region == null) {
				int taxonNr = geoprior.getTaxonNr();
				double lat = position[taxonNr][0] + Randomizer.nextGaussian() * stdDevLat;
				sampledLocations.setValue(taxonNr * 2, lat);
				
				double long_ = position[taxonNr][1] + Randomizer.nextGaussian() * stdDevLong;
				sampledLocations.setValue(taxonNr * 2 + 1, long_);
			}
		}
		
		calcStdDevs();
		double backwardStdDevLat = stdDevLat;
		double backwardStdDevLong = stdDevLong;

		NormalDistribution normalLat1 = new NormalDistributionImpl(0, forwardStdDevLat);
		NormalDistribution normalLong1 = new NormalDistributionImpl(0, forwardStdDevLong);
		NormalDistribution normalLat2 = new NormalDistributionImpl(0, backwardStdDevLat);
		NormalDistribution normalLong2 = new NormalDistributionImpl(0, backwardStdDevLong);
		
		for (GeoPrior geoprior : geopriors) {
			if (geoprior.region == null) {
				int taxonNr = geoprior.getTaxonNr();
				double lat = position[taxonNr][0] - sampledLocations.getValue(taxonNr * 2);
				double long_ = position[taxonNr][1] - sampledLocations.getValue(taxonNr * 2 + 1);
				logHR -= normalLat1.logDensity(lat) + normalLong1.logDensity(long_);
				
				lat = position[taxonNr][0] - orgValues[taxonNr * 2];
				long_ = position[taxonNr][1] - orgValues[taxonNr * 2 + 1];
				logHR += normalLat2.logDensity(lat) + normalLong2.logDensity(long_);
			}
		}
		
		
		return logHR;
	}

	void calcStdDevs() {
		MultiGeoPrior multiGeoPrior = likelihood.multiGeopriorsInput.get();
		List<GeoPrior> geopriors = multiGeoPrior.geopriorsInput.get();

		double meanLat = 0, meanLong = 0, squaredLat = 0, squaredLong = 0;
		int k = 0;
		for (GeoPrior geoprior : geopriors) {
			if (geoprior.region == null) {
				int taxonNr = geoprior.getTaxonNr();
				double lat = sampledLocations.getValue(taxonNr * 2);
				double diff = (lat - position[taxonNr][0]);
				meanLat += diff;
				squaredLat += diff * diff;
				double long_ = sampledLocations.getValue(taxonNr * 2 + 1);
				diff = (long_ - position[taxonNr][1]);
				meanLong += diff;
				squaredLong += diff * diff;
				
				k++;
			}
		}
		
		
		final double  epsilon = 1e-10;
		meanLat /= k;
		squaredLat /= k;
		stdDevLat = Math.sqrt(squaredLat - meanLat * meanLat) + epsilon;
		meanLong /= k;
		squaredLong /= k;
		stdDevLong = Math.sqrt(squaredLong - meanLong * meanLong) + epsilon;

	}
	
	void calcBranchLengths() {
		BranchRateModel clockModel = likelihood.branchRateModelInput.get();
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				branchLengths[node.getNr()] = node.getLength() * clockModel.getRateForBranch(node);
//				if (!node.isLeaf()) {
//					Node child1 = node.getLeft();
//					Node child2 = node.getRight();
//					sumLengths[node.getNr()] = branchLengths[node.getNr()] +
//							child1.getLength() * clockModel.getRateForBranch(child1) +
//							child2.getLength() * clockModel.getRateForBranch(child2);
//				}
			}
		}
	}

	
	void calcPositions(List<Integer> sampleNumber, double precision) {
		final double EPSILON = 1e-8;
		
		// process sampled locations
		for (int i : sampleNumber) {
			double lat1 = sampledLocations.getMatrixValue(i, 0);
			double long1 = sampledLocations.getMatrixValue(i, 1);
			if (position[i][0] != lat1 || position[i][1] != long1) {
				position[i][0] = lat1;
				position[i][1] = long1;
				sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
			}
		}
		
		initByMean(tree.getRoot());
		resetMeanDown(tree.getRoot());
		
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			if (!isSampled[i]) {
				position[i] = SphericalDiffusionModel.cartesian2Sperical(sphereposition[i], true);
			} else {
				//System.err.print("[skip " + i +"]");
			}
		}
		//System.err.print("pos153B=["+position[153][0] + "," + position[153][1] + ":" + sphereposition[153][0]+","+sphereposition[153][1]+","+sphereposition[153][2]+"]");
		//System.err.println();
		
	}

	void initByMean(Node node) {
		if (!node.isLeaf()) {
			initByMean(node.getLeft());
			initByMean(node.getRight());
			
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (!isSampled[nodeNr]) {
				setHalfWayPosition(nodeNr, child1, child2);
			} else {
				parentweight[nodeNr] = 0;
				//System.err.print("skip" + nodeNr);
			}
		}
	}		

	/** top down recalculation **/
	void resetMeanDown(Node node) {
		if (!node.isLeaf()) {
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (node.isRoot() || isSampled[node.getNr()]) {
				//setHalfWayPosition(nodeNr, child1, child2);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent);
			}

			resetMeanDown(node.getLeft());
			resetMeanDown(node.getRight());
		}
	}

	void setHalfWayPosition(int nodeNr, int child1, int child2) {
		double b1 = 1.0/Math.sqrt(branchLengths[child1]/precision);
		double b2 = 1.0/Math.sqrt(branchLengths[child2]/precision);
		if (tree.getNode(nodeNr).isRoot()) {
			double len = b1 + b2;
			b1 = b1 / len;
			b2 = b2 / len;
			double w = (1.0 + b1 * parentweight[child1] + b2 * parentweight[child2]);
			b1 /= w;
			b2 /= w;
		} else {
			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
			final double len = b1 + b2 + p;
            if( true ) {
                //same as other branch
                final double s = 1.0 / (len - b1 * parentweight[child1] - b2 * parentweight[child2]);
                b1 *= s;
                b2 *= s;
                parentweight[nodeNr] = p * s;
            } else {
                b1 /= len;
                b2 /= len;
                p /= len;
                double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
                b1 /= w;
                b2 /= w;
                p /= w;
                parentweight[nodeNr] = p;
            }
		}
		sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2);
		sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2);
		sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2);
		normalise(sphereposition[nodeNr]);
		//System.err.print("["+sphereposition[nodeNr][0]+","+sphereposition[nodeNr][1]+","+sphereposition[nodeNr][2]+"]");
	}

	void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
        final double sp = 1.0/precision;
		double b1 = 1.0/Math.sqrt(branchLengths[child1] * sp);
		double b2 = 1.0/Math.sqrt(branchLengths[child2] * sp);
		double p = 1.0/Math.sqrt(branchLengths[nodeNr] * sp);
		final double len = b1 + b2 + p;
        if( true ) {
            final double s = (len - b1 * parentweight[child1] - b2 * parentweight[child2]);
            p /= s;
        } else {
            b1 /= len;
            b2 /= len;
            p /= len;
            double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
            p /= w;
        }
		sphereposition[nodeNr][0] += sphereposition[parent][0] * p;
		sphereposition[nodeNr][1] += sphereposition[parent][1] * p;
		sphereposition[nodeNr][2] += sphereposition[parent][2] * p;
		normalise(sphereposition[nodeNr]);
	}

	private void normalise(double[] position) {
		double len = Math.sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
		position[0] /= len;
		position[1] /= len;
		position[2] /= len;
	}
	
	
}
