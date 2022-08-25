package sphericalGeo;




import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
//import beast.base.evolution.alignment.AlignmentFromTraitMap;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
//import beast.base.evolution.tree.TreeTraitMap;
import beast.base.util.Randomizer;

@Description("Approximate likelihood by particle filter approximation")
public class PFApproxMultivariateTraitLikelihood extends GenericTreeLikelihood implements LocationProvider {
	public Input<Integer> nrOfParticlesInput = new Input<>("nrOfParticles", "number of particles to use", 1);//100
	public Input<Integer> nrOfIterationsInput = new Input<>("nrOfIterations", "number of iterations to run the particle filter", 100);
	public Input<Integer> rangeSizeInput = new Input<>("nrrange", "number of random samples for placing a node", 20);//10
	
	public Input<Boolean> scaleByBranchLengthInput = new Input<>("scale", "scale by branch lengths for initial position", true);

	public Input<List<GeoPrior>> geopriorsInput = new Input<>("geoprior", "geographical priors on tips, root or clades restricting these nodes to a region", new ArrayList<>());
	public Input<RealParameter> locationInput = new Input<>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree");

	public Input<Transformer> transformerInput = new Input<>("transformer","landscape transformer to capture some inheterogenuity in the diffusion process");
	public Input<Double> epsilonInput = new Input<>("epsilon", "size of interval to use when randomly choosing new latitude/longitude", 3.0);


	double epsilon = 2.0;
	boolean scaleByBranchLength;

	RealParameter sampledLocations;
	boolean [] isSampled;
	List<Integer> sampleNumber;

	double [][] sphereposition;
	double [] parentweight;

	class LeafParticleSet {
		int particleCount;
		int nodeNr;
		double [] initialPosition;
		
		LeafParticleSet(int nodeNr, int particleCount, double [] initialPosition) {
			this.nodeNr = nodeNr;
			this.particleCount = particleCount;
			this.initialPosition = initialPosition;
		}
		double [] getPosition(int iParticle) {
			return initialPosition;
		}
		void init(double [] position) {
			throw new RuntimeException("thou shalt not call init on LeafParticleSet");
		}
		
		void resample(SphericalDiffusionModel substModel) {}
		public void randomize() {}
	} // class LeafParticleSet

	
	class ParticleSet extends LeafParticleSet {
		double [][] pposition;
		double [] logP;
		
		ParticleSet(int nodeNr, int particleCount, double [] initialPosition) {
			super(nodeNr, particleCount, initialPosition);
			pposition = new double[particleCount][2];
			logP = new double[particleCount];
		}
		
		@Override
		double [] getPosition(int iParticle) {
			return pposition[iParticle];
		}

		@Override
		void init(double [] position) {
			for (double [] pos : pposition) {
				pos[0] = position[0];
				pos[1] = position[1];
			}
		}

		@Override
		void resample(SphericalDiffusionModel substModel) {
			Node node = tree.getNode(nodeNr);
			int iChild1 = node.getLeft().getNr();
			LeafParticleSet child1 = particleSets[iChild1];
			int iChild2 = node.getRight().getNr();
			LeafParticleSet child2 = particleSets[iChild2];
			
			if (node.isRoot()) {
				for (int i = 0; i < particleCount; i++) {
					double [] nodePosition = getPosition(i);
					logP[i] = substModel.getLogLikelihood(null, child1.getPosition(i), nodePosition, branchLengths[iChild1]) +
							substModel.getLogLikelihood(null, child2.getPosition(i), nodePosition, branchLengths[iChild2]);
				}
			} else {
				int iParent = node.getParent().getNr();
				LeafParticleSet parent = particleSets[iParent];
				for (int i = 0; i < particleCount; i++) {
					double [] nodePosition = getPosition(i);
					logP[i] = substModel.getLogLikelihood(null, child1.getPosition(i), nodePosition, branchLengths[iChild1]) +
							substModel.getLogLikelihood(null, child2.getPosition(i), nodePosition, branchLengths[iChild2]) +
							substModel.getLogLikelihood(null, nodePosition, parent.getPosition(i), branchLengths[node.getNr()]);
				}
			}
			double max = logP[0];
			for (double d : logP) {
				max = Math.max(d,  max);
			}
			for (int i = 0; i < particleCount; i++) {
				logP[i] = Math.exp(logP[i] - max);
			}
			double [][] newposition = new double[particleCount][2];
			for (int i = 0; i < particleCount; i++) {
				int k = Randomizer.randomChoicePDF(logP);
				newposition[i][0] = pposition[k][0];
				newposition[i][1] = pposition[k][1];
			}
			pposition = newposition;
		}
		
		@Override
		public void randomize() {
			for (int i = 0; i < particleCount; i++) {
				pposition[i][0] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
				pposition[i][1] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
			}			
		}

	} // class ParticleSet
	
	
	
	
	SphericalDiffusionModel substModel;
	TreeInterface tree;
	BranchRateModel clockModel;
	Transformer transformer;
	int particleCount;
	int iterationCount;
	int rangeSize;
	
	double [][] position;
	double [] branchLengths;
	double [] sumLengths;
	public boolean needsUpdate = true;
	
	LeafParticleSet [] particleSets;
	
	/** particles x positions x 2 **/
	double [][][] particlePosition;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		scaleByBranchLength = scaleByBranchLengthInput.get();
		clockModel = branchRateModelInput.get();
		SiteModel siteModel = (SiteModel) siteModelInput.get();
		substModel = (SphericalDiffusionModel) siteModel.substModelInput.get();
		tree = treeInput.get();
		
		transformer = transformerInput.get();
		
		// initialise leaf positions
		position = new double[tree.getNodeCount()][2];
		AlignmentFromTraitMap data = (AlignmentFromTraitMap) dataInput.get();
		TreeTraitMap traitMap = data.getTraitMap(); 
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			position[i] = traitMap.getTrait(tree, nodes[i]);
			if (transformer != null) {
				position[i] = transformer.project(position[i][0], position[i][1]);
			}
			
		}
		branchLengths = new double[tree.getNodeCount()];
		sumLengths = new double[tree.getNodeCount()];
		parentweight = new double[tree.getNodeCount()];

		
		particleCount = nrOfParticlesInput.get();
		iterationCount = nrOfIterationsInput.get();
		rangeSize = rangeSizeInput.get();

		// set up particles
//		particleSets = new LeafParticleSet[tree.getNodeCount()];
//		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
//			particleSets[i] = new LeafParticleSet(i, particleCount, position[i]);
//		}
//		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
//			particleSets[i] = new ParticleSet(i, particleCount, position[i]);
//		}
		
		
		// set up particles
		particlePosition = new double[particleCount][tree.getNodeCount()][2];
		for (double [][] pposition : particlePosition) {
			for (int i = 0; i < tree.getLeafNodeCount(); i++) {
				pposition[i] = position[i];
			}
		}
		
		List<GeoPrior> geopriors = geopriorsInput.get();
		isSampled = new boolean[tree.getNodeCount()];
		sampleNumber = new ArrayList<>();
		if (geopriors.size() > 0) {
			sampledLocations = locationInput.get();
			for (GeoPrior prior : geopriors) {
				if (prior.allInternalNodes) {
					sampleNumber.add(-1);
					for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
						isSampled[i] = true;
					}
				} else {
					isSampled[prior.getTaxonNr()] = true;
					sampleNumber.add(prior.getTaxonNr());
				}
			}
		}
		
		sphereposition = new double[position.length][3];
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
		}

		epsilon = epsilonInput.get();
	}

	@Override
    public double getCurrentLogP() {
        double logP = Double.NaN;
		try {
			logP = calculateLogP();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return logP;
    }

	@Override
	public double calculateLogP() {
		logP = 0.0;
		calcBranchLengths();
		setUpInitialPositions();
		for (int i = 0; i < iterationCount; i++) {
			sampleUp(tree.getRoot());
			resample();
			sampleDown(tree.getRoot());
			resample();
		}

		logP = calcLogP();
		needsUpdate = false;
		return logP;
	}

	void sampleUp(Node node) {
		if (!node.isLeaf()) {
			for (double [][] pposition : particlePosition) {
				choosePosition(node.getNr(), pposition);
			}			

			//particleSets[node.getNr()].randomize();
			//particleSets[node.getNr()].resample(substModel);

			sampleUp(node.getLeft());
			sampleUp(node.getRight());
		}
	}

	void sampleDown(Node node) {
		if (!node.isLeaf()) {
			sampleDown(node.getLeft());
			sampleDown(node.getRight());
			
			for (double [][] pposition : particlePosition) {
				choosePosition(node.getNr(), pposition);
			}
			
			//particleSets[node.getNr()].randomize();
			//particleSets[node.getNr()].resample(substModel);

		}
	}

	private void normalise(double[] position) {
		double len = Math.sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
		position[0] /= len;
		position[1] /= len;
		position[2] /= len;
	}

	private void choosePosition(int nodeNr, double[][] pposition) {
		// randomize
		double [][] newPosition = new double[rangeSize][2];
		for (int i = 0; i < rangeSize; i++) {
			if (!isSampled[nodeNr]) {
//				double [] spherePosition = SphericalDiffusionModel.spherical2Cartesian(pposition[i][0], pposition[i][1]);
//				spherePosition[0] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
//				spherePosition[1] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
//				spherePosition[2] += Randomizer.nextDouble() * epsilon - epsilon / 2.0;
//				normalise(spherePosition);
//				newPosition[i] = SphericalDiffusionModel.cartesian2Sperical(spherePosition, true);
				newPosition[i][0] = pposition[nodeNr][0] + Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
				newPosition[i][1] = pposition[nodeNr][1] + Randomizer.nextDouble() * epsilon - epsilon / 2.0;
			} else {
				newPosition[i][0] = pposition[nodeNr][0];
				newPosition[i][1] = pposition[nodeNr][1];
			}
		}
		
		// resample
		Node node = tree.getNode(nodeNr);
		int iChild1 = node.getLeft().getNr();
		int iChild2 = node.getRight().getNr();
		
		
		double [] logP = new double[rangeSize];
		if (node.isRoot()) {
			for (int i = 0; i < rangeSize; i++) {
				double [] nodePosition = newPosition[i];
				logP[i] = substModel.getLogLikelihood(null, pposition[iChild1], nodePosition, branchLengths[iChild1]) +
						substModel.getLogLikelihood(null, pposition[iChild2], nodePosition, branchLengths[iChild2]);
			}
		} else {
			int iParent = node.getParent().getNr();
			for (int i = 0; i < rangeSize; i++) {
				double [] nodePosition = newPosition[i];
				logP[i] = substModel.getLogLikelihood(null, pposition[iChild1], nodePosition, branchLengths[iChild1]) +
						substModel.getLogLikelihood(null, pposition[iChild2], nodePosition, branchLengths[iChild2]) +
						substModel.getLogLikelihood(null, nodePosition, pposition[iParent], branchLengths[node.getNr()]);
			}
		}
		double max = logP[0];
		for (double d : logP) {
			max = Math.max(d,  max);
		}
		for (int i = 0; i < rangeSize; i++) {
			logP[i] = Math.exp(logP[i] - max);
		}
		int k = Randomizer.randomChoicePDF(logP);
		pposition[nodeNr][0] = newPosition[k][0];
		pposition[nodeNr][1] = newPosition[k][1];
	}

	
	void resample() {
		double [] logP = new double[particleCount];
		// logP is mean of logP's of the particles(?)
		for (int i = 0; i < particleCount;  i++) {
			for (Node node : tree.getNodesAsArray()) {
				if (!node.isRoot()) {
//					logP += substModel.getLogLikelihood(
//							particleSets[node.getParent().getNr()].getPosition(i), 
//							particleSets[node.getNr()].getPosition(i), 
//							branchLengths[node.getNr()]);
					logP[i] += substModel.getLogLikelihood(null, 
							particlePosition[i][node.getParent().getNr()], 
							particlePosition[i][node.getNr()], 
							branchLengths[node.getNr()]);
				}
			}
		}
		double max = logP[0];
		for (double d : logP) {
			max = Math.max(d,  max);
		}
		for (int i = 0; i < particleCount; i++) {
			logP[i] = Math.exp(logP[i] - max);
		}
		
		double [][][] newposition = new double[particleCount][tree.getNodeCount()][2];
		for (int i = 0; i < particleCount; i++) {
			int k = Randomizer.randomChoicePDF(logP);
			copy(newposition[i], particlePosition[k]);
		}
		particlePosition = newposition;
	}

	
	private void copy(double[][] dest, double[][] src) {
		for (int i = 0; i < dest.length; i++) {
			dest[i][0] = src[i][0];
			dest[i][1] = src[i][1];
		}
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
		// logP is mean of logP's of the particles(?)
		for (int i = 0; i < particleCount;  i++) {
			for (Node node : tree.getNodesAsArray()) {
				if (!node.isRoot()) {
//					logP += substModel.getLogLikelihood(
//							particleSets[node.getParent().getNr()].getPosition(i), 
//							particleSets[node.getNr()].getPosition(i), 
//							branchLengths[node.getNr()]);
					logP += substModel.getLogLikelihood(null, 
							particlePosition[i][node.getParent().getNr()], 
							particlePosition[i][node.getNr()], 
							branchLengths[node.getNr()]);
				}
			}
		}
		logP /= particleCount;
		return logP;
	}
	
	void setUpInitialPositions() {
		// process sampled locations
//		for (int i : sampleNumber) {
//			position[i][0] = sampledLocations.getMatrixValue(i, 0);
//			position[i][1] = sampledLocations.getMatrixValue(i, 1);
//			if (transformer != null) {
//				double [] t = transformer.project(position[i][0], position[i][1]);
//				position[i][0] = t[0];
//				position[i][1] = t[1];
//			}
//		}
		
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
		}

		// process sampled locations
		for (int i : sampleNumber) {
			double lat1 = sampledLocations.getMatrixValue(i, 0);
			double long1 = sampledLocations.getMatrixValue(i, 1);
			if (transformer != null) {
				double [] t = transformer.project(lat1, long1);
				lat1 = t[0];
				long1 = t[1];
			}
			if (position[i][0] != lat1 || position[i][1] != long1) {
				position[i][0] = lat1;
				position[i][1] = long1;
				sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
			}
		}

		
		initByMean(tree.getRoot());
		resetMean2(tree.getRoot());

		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			if (!isSampled[i]) {
				position[i] = SphericalDiffusionModel.cartesian2Sperical(sphereposition[i], true);
			} else {
				//System.err.print("[skip " + i +"]");
			}
		}
		
		for (double [][] pposition : particlePosition) {
			for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
				pposition[i][0] = position[i][0];
				pposition[i][1] = position[i][1];
			}
			for (int i : sampleNumber) {
				pposition[i][0] = position[i][0];
				pposition[i][1] = position[i][1];
				
			}
		}
		
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
	void resetMean(Node node) {
		if (!node.isLeaf()) {
			resetMean(node.getLeft());
			resetMean(node.getRight());
			
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
	void resetMean2(Node node) {
		if (!node.isLeaf()) {
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

			resetMean2(node.getLeft());
			resetMean2(node.getRight());
		}
	}

	void setHalfWayPosition(int nodeNr, int child1, int child2) {
		// start in weighted middle of the children
		double precision = 1.0;
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]/precision);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]/precision);
			double len = b1 + b2;
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2) / len;
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2) / len;
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2) / len;
//			double len = 1.0/branchLengths[child1] + 1.0/branchLengths[child2];
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] / branchLengths[child1] + sphereposition[child2][0] / branchLengths[child2]) / len;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] / branchLengths[child1] + sphereposition[child2][1] / branchLengths[child2]) / len;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] / branchLengths[child1] + sphereposition[child2][2] / branchLengths[child2]) / len;
		} else{
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
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] + sphereposition[child2][0]) / 2.0;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] + sphereposition[child2][1]) / 2.0;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] + sphereposition[child2][2]) / 2.0;
		}
		normalise(sphereposition[nodeNr]);
		//System.err.print("["+sphereposition[nodeNr][0]+","+sphereposition[nodeNr][1]+","+sphereposition[nodeNr][2]+"]");
//		double [] pos = SphericalDiffusionModel.cartesian2Sperical(sphereposition[nodeNr], true);
//		System.err.print("["+pos[0]+","+pos[1]+"]");
//		pos = SphericalDiffusionModel.map(pos[0], pos[1], -80, 0);
//		System.err.println("["+pos[0]+","+pos[1]+"]");
	}

//	private void setHalfWayPosition(int nodeNr, int child1, int child2) {
//		// start in weighted middle of the children
//		if (scaleByBranchLength) {
//			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
//			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
//			double len = b1 + b2;
//			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2) / len;
//			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2) / len;
//			// position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
//			// position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
//		} else {
//			position[nodeNr][0] = (position[child1][0] + position[child2][0]) / 2.0;
//			position[nodeNr][1] = (position[child1][1] + position[child2][1]) / 2.0;
//		}
//		
//	}

	void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
		final double precision = 1.0;
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]/precision);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]/precision);
			double p = 1.0/Math.sqrt(branchLengths[nodeNr]/precision);
			double len = b1 + b2 + p;
			sphereposition[nodeNr][0] = (sphereposition[child1][0] * b1 + sphereposition[child2][0] * b2 + sphereposition[parent][0] * p) / len;
			sphereposition[nodeNr][1] = (sphereposition[child1][1] * b1 + sphereposition[child2][1] * b2 + sphereposition[parent][1] * p) / len;
			sphereposition[nodeNr][2] = (sphereposition[child1][2] * b1 + sphereposition[child2][2] * b2 + sphereposition[parent][2] * p) / len;

//			double len = 1.0/branchLengths[child1] + 1.0/branchLengths[child2] + 1.0/branchLengths[nodeNr];
//			sphereposition[nodeNr][0] = (sphereposition[child1][0] / branchLengths[child1] + sphereposition[child2][0] / branchLengths[child2] + sphereposition[parent][0] / branchLengths[nodeNr]) / len;
//			sphereposition[nodeNr][1] = (sphereposition[child1][1] / branchLengths[child1] + sphereposition[child2][1] / branchLengths[child2] + sphereposition[parent][1] / branchLengths[nodeNr]) / len;
//			sphereposition[nodeNr][2] = (sphereposition[child1][2] / branchLengths[child1] + sphereposition[child2][2] / branchLengths[child2] + sphereposition[parent][2] / branchLengths[nodeNr]) / len;
		} else {
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
		}
		normalise(sphereposition[nodeNr]);

//		double [] pos = SphericalDiffusionModel.cartesian2Sperical(sphereposition[nodeNr], true);
//		System.err.print("["+pos[0]+","+pos[1]+"]");
//		pos = SphericalDiffusionModel.map(pos[0], pos[1], -80, 0);
//		System.err.println("["+pos[0]+","+pos[1]+"]");
	}
//	private void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
//		// start in weighted middle of the children and parent location
//		if (scaleByBranchLength) {
//			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
//			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
//			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
//			double len = b1 + b2 + p;
//			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2 + position[parent][0] * p) / len;
//			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2 + position[parent][1] * p) / len;
//			// position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2] + position[parent][0] * branchLengths[nodeNr]) / sumLengths[nodeNr];
//			// position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2] + position[parent][1] * branchLengths[nodeNr]) / sumLengths[nodeNr];
//		} else {
//			position[nodeNr][0] = (position[child1][0] + position[child2][0] + position[parent][0]) / 3.0;
//			position[nodeNr][1] = (position[child1][1] + position[child2][1] + position[parent][1]) / 3.0;
//		}
//
//	}

	@Override
	public boolean isStochastic() {
		return true;
	}

	@Override
	public double[] getPosition(int iDim) {
		if (needsUpdate) {
			try {
				calculateLogP();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		if (transformer != null) {
			double lat = particlePosition[0][iDim][0];
			double lon = particlePosition[0][iDim][1];
			while (lon < -180) {
				lon += 360;
			}
			while (lon > 180) {
				lon -= 360;
			}
			double [] pos = transformer.projectInverse(lat, lon);
			if (pos == null) {
				pos = transformer.projectInverse(lat, lon);
			}
			return pos;
		}
		return particlePosition[0][iDim];
	}
	
	@Override
	public void store() {
		needsUpdate = true;
		super.store();
	}
	@Override
	public void restore() {
		needsUpdate = true;
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		return super.requiresRecalculation();
	}
	
}
