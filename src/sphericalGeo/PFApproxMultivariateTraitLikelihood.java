package sphericalGeo;



import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.AlignmentFromTraitMap;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeTraitMap;
import beast.util.Randomizer;

@Description("Approximate likelihood by particle filter approximation")
public class PFApproxMultivariateTraitLikelihood extends GenericTreeLikelihood {
	public Input<Integer> nrOfParticlesInput = new Input<Integer>("nrOfParticles", "number of particles to use", 25);//100
	public Input<Integer> nrOfIterationsInput = new Input<Integer>("nrOfIterations", "number of iterations to run the particle filter", 10);
	public Input<Integer> rangeSizeInput = new Input<Integer>("nrrange", "number of random samples for placing a node", 10);//10
	
	public Input<Boolean> scaleByBranchLengthInput = new Input<Boolean>("scale", "scale by branch lengths for initial position", true);

	double epsilon = 2.0;
	boolean scaleByBranchLength;

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
					logP[i] = substModel.getLogLikelihood(child1.getPosition(i), nodePosition, branchLengths[iChild1]) +
							substModel.getLogLikelihood(child2.getPosition(i), nodePosition, branchLengths[iChild2]);
				}
			} else {
				int iParent = node.getParent().getNr();
				LeafParticleSet parent = particleSets[iParent];
				for (int i = 0; i < particleCount; i++) {
					double [] nodePosition = getPosition(i);
					logP[i] = substModel.getLogLikelihood(child1.getPosition(i), nodePosition, branchLengths[iChild1]) +
							substModel.getLogLikelihood(child2.getPosition(i), nodePosition, branchLengths[iChild2]) +
							substModel.getLogLikelihood(nodePosition, parent.getPosition(i), branchLengths[node.getNr()]);
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
	int particleCount;
	int iterationCount;
	int rangeSize;
	
	double [][] position;
	double [] branchLengths;
	double [] sumLengths;
	boolean needsUpdate = true;
	
	LeafParticleSet [] particleSets;
	
	/** particles x positions x 2 **/
	double [][][] particlePosition;
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		scaleByBranchLength = scaleByBranchLengthInput.get();
		clockModel = branchRateModelInput.get();
		SiteModel siteModel = (SiteModel) siteModelInput.get();
		substModel = (SphericalDiffusionModel) siteModel.substModelInput.get();
		tree = treeInput.get();
		
		// initialise leaf positions
		position = new double[tree.getNodeCount()][2];
		AlignmentFromTraitMap data = (AlignmentFromTraitMap) dataInput.get();
		TreeTraitMap traitMap = data.getTraitMap(); 
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			position[i] = traitMap.getTrait(tree, nodes[i]);
		}
		branchLengths = new double[tree.getNodeCount()];
		sumLengths = new double[tree.getNodeCount()];
		
		
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
	public double calculateLogP() throws Exception {
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

	private void choosePosition(int nodeNr, double[][] pposition) {
		// randomize
		double [][] newPosition = new double[rangeSize][2];
		for (int i = 0; i < rangeSize; i++) {
			newPosition[i][0] = pposition[nodeNr][0] + Randomizer.nextDouble() * epsilon - epsilon / 2.0;			
			newPosition[i][1] = pposition[nodeNr][1] + Randomizer.nextDouble() * epsilon - epsilon / 2.0;
		}
		
		// resample
		Node node = tree.getNode(nodeNr);
		int iChild1 = node.getLeft().getNr();
		int iChild2 = node.getRight().getNr();
		
		
		double [] logP = new double[rangeSize];
		if (node.isRoot()) {
			for (int i = 0; i < rangeSize; i++) {
				double [] nodePosition = newPosition[i];
				logP[i] = substModel.getLogLikelihood(pposition[iChild1], nodePosition, branchLengths[iChild1]) +
						substModel.getLogLikelihood(pposition[iChild2], nodePosition, branchLengths[iChild2]);
			}
		} else {
			int iParent = node.getParent().getNr();
			for (int i = 0; i < rangeSize; i++) {
				double [] nodePosition = newPosition[i];
				logP[i] = substModel.getLogLikelihood(pposition[iChild1], nodePosition, branchLengths[iChild1]) +
						substModel.getLogLikelihood(pposition[iChild2], nodePosition, branchLengths[iChild2]) +
						substModel.getLogLikelihood(nodePosition, pposition[iParent], branchLengths[node.getNr()]);
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
					logP[i] += substModel.getLogLikelihood(
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
					logP += substModel.getLogLikelihood(
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
		initByMean(tree.getRoot());
		resetMean2(tree.getRoot());
		
		for (double [][] pposition : particlePosition) {
			for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
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
			setHalfWayPosition(nodeNr, child1, child2);
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
				setHalfWayPosition(nodeNr, child1, child2);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent);
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
				setHalfWayPosition(nodeNr, child1, child2);
			} else {
				int parent = node.getParent().getNr();
				setHalfWayPosition(nodeNr, child1, child2, parent);
			}

			resetMean2(node.getLeft());
			resetMean2(node.getRight());
		}
	}

	private void setHalfWayPosition(int nodeNr, int child1, int child2) {
		// start in weighted middle of the children
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double len = b1 + b2;
			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2) / len;
			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2) / len;
			// position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
			// position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2]) / (branchLengths[child1] + branchLengths[child2]);
		} else {
			position[nodeNr][0] = (position[child1][0] + position[child2][0]) / 2.0;
			position[nodeNr][1] = (position[child1][1] + position[child2][1]) / 2.0;
		}
		
	}

	private void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
		if (scaleByBranchLength) {
			double b1 = 1.0/Math.sqrt(branchLengths[child1]);
			double b2 = 1.0/Math.sqrt(branchLengths[child2]);
			double p = 1.0/Math.sqrt(branchLengths[nodeNr]);
			double len = b1 + b2 + p;
			position[nodeNr][0] = (position[child1][0] * b1 + position[child2][0] * b2 + position[parent][0] * p) / len;
			position[nodeNr][1] = (position[child1][1] * b1 + position[child2][1] * b2 + position[parent][1] * p) / len;
			// position[nodeNr][0] = (position[child1][0] * branchLengths[child1] + position[child2][0] * branchLengths[child2] + position[parent][0] * branchLengths[nodeNr]) / sumLengths[nodeNr];
			// position[nodeNr][1] = (position[child1][1] * branchLengths[child1] + position[child2][1] * branchLengths[child2] + position[parent][1] * branchLengths[nodeNr]) / sumLengths[nodeNr];
		} else {
			position[nodeNr][0] = (position[child1][0] + position[child2][0] + position[parent][0]) / 3.0;
			position[nodeNr][1] = (position[child1][1] + position[child2][1] + position[parent][1]) / 3.0;
		}

	}
	
	@Override
	public boolean isStochastic() {
		return true;
	}

	public double[] getPostion(int iDim) {
		if (needsUpdate) {
			try {
				calculateLogP();
			} catch (Exception e) {
				e.printStackTrace();
			}
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
