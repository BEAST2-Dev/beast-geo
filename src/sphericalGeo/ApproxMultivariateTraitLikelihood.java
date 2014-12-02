package sphericalGeo;



import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.AlignmentFromTraitMap;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeTraitMap;

@Description("Approximate likelihood by MAP approximation of internal states")
public class ApproxMultivariateTraitLikelihood extends GenericTreeLikelihood {
	public Input<Boolean> scaleByBranchLengthInput = new Input<Boolean>("scale", "scale by branch lengths for initial position", true);
	public Input<List<GeoPrior>> geopriorsInput = new Input<List<GeoPrior>>("geoprior", "geographical priors on tips, root or clades restricting these nodes to a region", new ArrayList<>());
	public Input<RealParameter> locationInput = new Input<RealParameter>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree");

	
	
	SphericalDiffusionModel substModel;
	TreeInterface tree;
	BranchRateModel clockModel;
	
	double [][] position;
	double [][] sphereposition;
	double [] branchLengths;
	double [] sumLengths;
	double [] parentweight;

	boolean needsUpdate = true;
	boolean scaleByBranchLength;
	
	PFApproxMultivariateTraitLikelihood loggerLikelihood;
	
	RealParameter sampledLocations;
	boolean [] isSampled;
	List<Integer> sampleNumber;
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		clockModel = branchRateModelInput.get();
		SiteModel siteModel = (SiteModel) siteModelInput.get();
		substModel = (SphericalDiffusionModel) siteModel.substModelInput.get();
		tree = treeInput.get();
		scaleByBranchLength = scaleByBranchLengthInput.get();
		
		// initialise leaf positions
		position = new double[tree.getNodeCount()][2];
		AlignmentFromTraitMap data = (AlignmentFromTraitMap) dataInput.get();
		TreeTraitMap traitMap = data.getTraitMap(); 
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			position[i] = traitMap.getTrait(tree, nodes[i]);
		}
		
		sphereposition = new double[position.length][3];
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
		}

		branchLengths = new double[tree.getNodeCount()];
		sumLengths = new double[tree.getNodeCount()];
		parentweight = new double[tree.getNodeCount()];
	
		List<GeoPrior> geopriors = geopriorsInput.get();
		isSampled = new boolean[tree.getNodeCount()];
		sampleNumber = new ArrayList<Integer>();
		if (geopriors.size() > 0) {
			sampledLocations = locationInput.get();
			for (GeoPrior prior : geopriors) {
				isSampled[prior.taxonNr] = true;
				sampleNumber.add(prior.taxonNr);
				double [] location = prior.region.sample();
				sampledLocations.setValue(prior.taxonNr * 2, location[0]);
				sampledLocations.setValue(prior.taxonNr * 2 + 1, location[1]);
			}
		}
		
		loggerLikelihood = new PFApproxMultivariateTraitLikelihood();
		if (geopriors.size() > 0) {
			loggerLikelihood.initByName("scale", scaleByBranchLength, "tree", tree, "siteModel", siteModel, 
				"branchRateModel", clockModel, "data", data,
				"geoprior", geopriors, "location", sampledLocations);
		} else {
			loggerLikelihood.initByName("scale", scaleByBranchLength, "tree", tree, "siteModel", siteModel, 
					"branchRateModel", clockModel, "data", data);
		}
	}
	
	@Override
    public double getCurrentLogP() {
        double logP = Double.NaN;
		try {
			// check prior
			if (sampledLocations != null) {
				for (GeoPrior prior : geopriorsInput.get()) {
					if (prior.calculateLogP() != 0) {
						logP = Double.NEGATIVE_INFINITY;
						return logP;
					}
				}
			}
			
			
			// calc likelihood
			logP = calculateLogP();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return logP;
    }

	
	@Override
	public double calculateLogP() throws Exception {
		logP = 0.0;
		calcBranchLengths();
		caclPositions();
		logP = calcLogP();
		needsUpdate = false;
		return logP;
	}

	void calcBranchLengths() {
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				branchLengths[node.getNr()] = node.getLength() * clockModel.getRateForBranch(node);
				if (!node.isLeaf()) {
					Node child1 = node.getLeft();
					Node child2 = node.getRight();
					sumLengths[node.getNr()] = branchLengths[node.getNr()] +
							child1.getLength() * clockModel.getRateForBranch(child1) +
							child2.getLength() * clockModel.getRateForBranch(child2);
				}
			}
		}
	}

	/** traverse tree **/
	double calcLogP() {
		double logP = 0;
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				logP += substModel.getLogLikelihood(
						position[node.getParent().getNr()], 
						position[node.getNr()], 
						branchLengths[node.getNr()]);
			}
		}
		return logP;
	}
	
//	void caclPositions() {
//		final double EPSILON = 1e-8;
//		
//		initByMean(tree.getRoot());
//		resetMeanDown(tree.getRoot());
//		
//		if (scaleByBranchLength) {
//			for (int i = 0; i < 50; i++) {
//				double [][] oldPosition = new double[tree.getNodeCount()][2];
//				for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
//					oldPosition[j][0] = position[j][0];
//					oldPosition[j][1] = position[j][1];
//				}
//				resetMeanDown(tree.getRoot());
//				resetMeanUp(tree.getRoot());
//				
//				double max = 0;
//				for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
//					double delta0 = oldPosition[j][0] - position[j][0];
//					double delta1 = oldPosition[j][1] - position[j][1];
//					max = Math.max(max, Math.max(Math.abs(delta0), Math.abs(delta1)));
//				}
//				if (max < EPSILON) {
//					// System.err.println(i + " maxdelta = " + max);
//					break;
//				}
//			}
//		}
//	}
	
	
	void caclPositions() {
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
		
			
		if (scaleByBranchLength) {
		double [][] oldPosition = new double[tree.getNodeCount()][3];
		for (int i = 0; i < 50; i++) {
			if (i % 5 == 0)
			for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
				oldPosition[j][0] = sphereposition[j][0];
				oldPosition[j][1] = sphereposition[j][1];
				oldPosition[j][2] = sphereposition[j][2];
			}

				resetMeanDown(tree.getRoot());
				resetMeanUp(tree.getRoot());

				if (i % 5 == 0) {
				double max = 0;
				for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
					double delta0 = oldPosition[j][0] - sphereposition[j][0];
					double delta1 = oldPosition[j][1] - sphereposition[j][1];
					double delta2 = oldPosition[j][2] - sphereposition[j][2];
					max = Math.max(max, Math.max(Math.abs(delta0), Math.max(Math.abs(delta1), Math.abs(delta2))));
				}
				if (max < EPSILON) {
					break;
				}
				}
					
			}
		}
		
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			position[i] = SphericalDiffusionModel.cartesian2Sperical(sphereposition[i]);
		}

//		System.err.println("maxdelta2 = " + max);
		
	}
	
	void setHalfWayPosition(int nodeNr, int child1, int child2) {
		// start in weighted middle of the children
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double len = b1 + b2;
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2) / len;
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2) / len;
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2) / len;
//			double len = 1.0/branchLengths[child1] + 1.0/branchLengths[child2];
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] / branchLengths[child1] + sphereposition[child2][0] / branchLengths[child2]) / len;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] / branchLengths[child1] + sphereposition[child2][1] / branchLengths[child2]) / len;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] / branchLengths[child1] + sphereposition[child2][2] / branchLengths[child2]) / len;
		} else{
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			if (tree.getNode(nodeNr).isRoot()) {
				double len = b1 + b2;
				b1 = b1 / len;
				b2 = b2 / len;
				double w = (1.0 + b1 * parentweight[child1] + b2 * parentweight[child2]);
				b1 /= w;
				b2 /= w;
			} else {
				double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
				double len = b1 + b2 + p;
				b1 /= len;
				b2 /= len;
				p /= len;
				double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
				b1 /= w;
				b2 /= w;
				p /= w;
				parentweight[nodeNr] = p;
			}
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2);
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2);
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2);
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] + sphereposition[child2][0]) / 2.0;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] + sphereposition[child2][1]) / 2.0;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] + sphereposition[child2][2]) / 2.0;
		}
		normalise(sphereposition[nodeNr]);
	}

	void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
			double len = b1 + b2 + p;
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2 + sphereposition[parent][0] * p) / len;
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2 + sphereposition[parent][1] * p) / len;
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2 + sphereposition[parent][2] * p) / len;

//			double len = 1.0/branchLengths[child1] + 1.0/branchLengths[child2] + 1.0/branchLengths[nodeNr];
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] / branchLengths[child1] + sphereposition[child2][0] / branchLengths[child2] + sphereposition[parent][0] / branchLengths[nodeNr]) / len;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] / branchLengths[child1] + sphereposition[child2][1] / branchLengths[child2] + sphereposition[parent][1] / branchLengths[nodeNr]) / len;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] / branchLengths[child1] + sphereposition[child2][2] / branchLengths[child2] + sphereposition[parent][2] / branchLengths[nodeNr]) / len;
		} else {
			double b1 = branchLengths[child1];
			double b2 = branchLengths[child2];
			double p = branchLengths[nodeNr];
			double len = b1 + b2 + p;
			b1 /= len;
			b2 /= len;
			p /= len;
			double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
			p /= w;
			sphereposition[nodeNr][0] += sphereposition[parent][0] * p;
			sphereposition[nodeNr][1] += sphereposition[parent][1] * p;
			sphereposition[nodeNr][2] += sphereposition[parent][2] * p;
		}
		normalise(sphereposition[nodeNr]);
	}


	private void normalise(double[] position) {
		double len = Math.sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
		position[0] /= len;
		position[1] /= len;
		position[2] /= len;
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
			}
		}
	}		
	
	/** bottom up recalculation **/
	void resetMeanUp(Node node) {
		if (!node.isLeaf()) {
			resetMeanUp(node.getLeft());
			resetMeanUp(node.getRight());
			
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (node.isRoot()) {
				if (!isSampled[nodeNr]) {
					setHalfWayPosition(nodeNr, child1, child2);
				}
			} else {
				if (!isSampled[nodeNr]) {
					int parent = node.getParent().getNr();
					setHalfWayPosition(nodeNr, child1, child2, parent);
				}
			}
		}
	}

	/** top down recalculation **/
	void resetMeanDown(Node node) {
		if (!node.isLeaf()) {
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();
			if (node.isRoot()) {
				//setHalfWayPosition(nodeNr, child1, child2);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent);
			}

			resetMeanDown(node.getLeft());
			resetMeanDown(node.getRight());
		}
	}

//	void setHalfWayPosition(int nodeNr, int child1, int child2) {
//		// start in weighted middle of the children
//		
//
//		if (scaleByBranchLength) {
//			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
//			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
//			double len = b1 + b2;
//			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2) / len;
//			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2) / len;
//			//position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
//			//position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
//		} else{
//			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
//			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
//			if (tree.getNode(nodeNr).isRoot()) {
//				double len = b1 + b2;
//				b1 = b1 / len;
//				b2 = b2 / len;
//				double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
//				b1 /= w;
//				b2 /= w;
//			} else {
//				double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
//				double len = b1 + b2 + p;
//				b1 /= len;
//				b2 /= len;
//				p /= len;
//				double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
//				b1 /= w;
//				b2 /= w;
//				p /= w;
//				parentweight[nodeNr] = p;
//			}
//			position[nodeNr][0] = position[child1][0] * b1 + position[child2][0] * b2;
//			position[nodeNr][1] = position[child1][1] * b1 + position[child2][1] * b2;
//		}
//		if (MAX_ITER <= 0) {return;}
//
//		// optimise for subst model
//		double [] newPos = new double[2];
//		newPos[0] = position[nodeNr][0]; 
//		newPos[1] = position[nodeNr][1]; 
//
//		double logP = 
//				substModel.getLogLikelihood(position[nodeNr], position[child1], branchLengths[child1]) +
//				substModel.getLogLikelihood(position[nodeNr], position[child2], branchLengths[child2]);
//		
//		// optimise by random walk
//		int i = 0;
//		double epsilon = 0.5;
//		while (i < MAX_ITER && epsilon > MIN_EPSILON) {
////			boolean progress = false;
////			newPos[0] += epsilon;
////			double newLogP = 
////					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
////					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]);
////			double dX = newLogP - logP;
////			if (dX > 0) {
////				logP = newLogP;
////				progress = true;
////			} else {
////				newPos[0] -= epsilon;
////			}
////			newPos[1] += epsilon;
////			newLogP = 
////					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
////					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]);
////			double dY= newLogP - logP;
////			if (dY > 0) {
////				logP = newLogP;
////				progress = true;
////			} else {
////				newPos[1] -= epsilon;
////			}
////			if (progress) {
////				position[nodeNr][0] = newPos[0]; 
////				position[nodeNr][1] = newPos[1];
////				epsilon *= 1.5;
////			} else {
////				epsilon /= 1.5;
////			}
//			
//			newPos[0] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
//			newPos[1] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
//			double newLogP = 
//					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
//					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]);
//			if (newLogP > logP) {
//				logP = newLogP;
//				position[nodeNr][0] = newPos[0]; 
//				position[nodeNr][1] = newPos[1]; 
//			} else {
//				newPos[0] = position[nodeNr][0]; 
//				newPos[1] = position[nodeNr][1]; 
//			}
//			i++;
//		}
//	}
//
//	void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
//		// start in weighted middle of the children and parent location
//		if (scaleByBranchLength) {
//			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
//			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
//			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
//			double len = b1 + b2 + p;
//			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2 + position[parent][0] * p) / len;
//			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2 + position[parent][1] * p) / len;
//
//			//position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2] + position[parent][0] * branchLengths[nodeNr]) / sumLengths[nodeNr];
//			//position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2] + position[parent][1] * branchLengths[nodeNr]) / sumLengths[nodeNr];
//		} else {
//			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
//			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
//			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
//			double len = b1 + b2 + p;
//			b1 /= len;
//			b2 /= len;
//			p /= len;
//			double w = (1.0 - b1 * parentweight[child1] - b2 * parentweight[child2]);
//			p /= w;
//			position[nodeNr][0] += position[parent][0] * p;
//			position[nodeNr][1] += position[parent][1] * p;
//
//		}
//
//		if (MAX_ITER <= 0) {return;}
//		
//		// optimise for subst model
//		double [] newPos = new double[2];
//		newPos[0] = position[nodeNr][0]; 
//		newPos[1] = position[nodeNr][1]; 
//
//		double logP = 
//				substModel.getLogLikelihood(position[nodeNr], position[child1], branchLengths[child1]) +
//				substModel.getLogLikelihood(position[nodeNr], position[child2], branchLengths[child2]) +
//				substModel.getLogLikelihood(position[parent], position[nodeNr], branchLengths[nodeNr]);
//		
//		// optimise by random walk
//		int i = 0;
//		double epsilon = 0.5;
//		while (i < MAX_ITER && epsilon > MIN_EPSILON) {
////			boolean progress = false;
////			newPos[0] += epsilon;
////			double newLogP = 
////					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
////					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]) +
////					substModel.getLogLikelihood(position[parent], newPos, branchLengths[nodeNr]);
////			double dX = newLogP - logP;
////			if (dX > 0) {
////				logP = newLogP;
////				progress = true;
////			} else {
////				newPos[0] -= epsilon;
////			}
////			newPos[1] += epsilon;
////			newLogP = 
////					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
////					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]) +
////					substModel.getLogLikelihood(position[parent], newPos, branchLengths[nodeNr]);
////			double dY= newLogP - logP;
////			if (dY > 0) {
////				logP = newLogP;
////				progress = true;
////			} else {
////				newPos[1] -= epsilon;
////			}
////			if (progress) {
////				position[nodeNr][0] = newPos[0]; 
////				position[nodeNr][1] = newPos[1];
////				epsilon *= 1.5;
////			} else {
////				epsilon /= 1.5;
////			}
//
//			newPos[0] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
//			newPos[1] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
//			double newLogP = 
//					substModel.getLogLikelihood(newPos, position[child1], branchLengths[child1]) +
//					substModel.getLogLikelihood(newPos, position[child2], branchLengths[child2]) +
//					substModel.getLogLikelihood(position[parent], newPos, branchLengths[nodeNr]);
//			if (newLogP > logP) {
//				logP = newLogP;
//				position[nodeNr][0] = newPos[0]; 
//				position[nodeNr][1] = newPos[1]; 
//			} else {
//				newPos[0] = position[nodeNr][0]; 
//				newPos[1] = position[nodeNr][1]; 
//			}
//			i++;
//		}
//	}
	final static int MAX_ITER = 0;
	final static double MIN_EPSILON = 0.001;
	
	@Override
	public boolean isStochastic() {
		return true;
	}

	public double[] getPostion(int iDim) {
		return loggerLikelihood.getPostion(iDim);
//		if (needsUpdate) {
//			try {
//				calculateLogP();
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
//		}
//		return position[iDim];
	}
	
	@Override
	public void store() {
		needsUpdate = true;
		loggerLikelihood.needsUpdate = true;
		super.store();
	}
	@Override
	public void restore() {
		needsUpdate = true;
		loggerLikelihood.needsUpdate = true;
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		loggerLikelihood.needsUpdate = true;
		return super.requiresRecalculation();
	}
	
}
