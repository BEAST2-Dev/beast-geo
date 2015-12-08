package sphericalGeo;



import java.util.*;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;


@Description(value="Approximate likelihood by MAP approximation of internal states optimised to take sampled nodes in account", 
			isInheritable=false)
@Citation("Remco R. Bouckaert. Phylogeography by diffusion on a sphere. bioRxiv, BIORXIV/2015/016311, 2015.")
public class ApproxMultivariateTraitLikelihood2 extends ApproxMultivariateTraitLikelihood { 
	
	/** maps each node to the partition containing the node **/
	int [] nodeToPartitionMap;
	
	/** flag to indicate for each node whether it is a 'root' for a partition **/
	boolean [] isTopOfPartition;
	/** array of root nodes, one for each partition **/
	int [] rootNode;
	
	/** partition root nodes can affect 3 (or at least 2 for global root) **/
	int [][] rootNodeToPartitionMap;
	int partitionCount;
	
	/** contribution to logP by each of the partitions **/
	double [] logPContributions;
	double [] storedLogPContributions;
	
	/** list of partitions that need recalculating **/
	List<Integer> dirtyPartitionList = new ArrayList<>();
	
	
	

	
    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] storedBranchLengths;
	double [][] storedPosition;
	double [][] storedSphereposition;
	//double [] storedSumLengths;
	//double [] storedParentweight;
	
	boolean wasInitialised;
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();

		branchLengths = new double[tree.getNodeCount()];
		storedBranchLengths = new double[tree.getNodeCount()];
		
		storedPosition = new double[tree.getNodeCount()][2];
		storedSphereposition = new double[tree.getNodeCount()][3];
	}
	
	@Override
	void initialiseSampledStates() {
		wasInitialised = true;
		super.initialiseSampledStates();
		if (!isMonoPhyletic) {
			return;
		}

		// initialise partition information
		isTopOfPartition = new boolean[tree.getNodeCount()];
		for (int i : sampleNumber) {
			isTopOfPartition[i] = true;
		}
		nodeToPartitionMap = new int[tree.getNodeCount()];
		Arrays.fill(nodeToPartitionMap, -1);
		int [] nextParitionNr = new int[1];
		nextParitionNr[0] = isSampled[tree.getRoot().getNr()] ? -1 : 0;
		initNodeToPartitionMap(tree.getRoot(), nextParitionNr, 0);
		partitionCount = nextParitionNr[0] + 1;

		rootNodeToPartitionMap = new int[tree.getNodeCount()][3];
		rootNode = new int[partitionCount];
		for (int i: sampleNumber) {
			Node node = tree.getNode(i);
			int left = node.getLeft().getNr();
			int right = node.getRight().getNr();
			if (!node.isRoot()) {
				// normal behaviour for internal node
				rootNodeToPartitionMap[i][0] = nodeToPartitionMap[i];
				rootNode[nodeToPartitionMap[i]] = firstSampledAncestor(node.getParent());
			} else {
				// hack to deal with root node
				rootNodeToPartitionMap[i][0] = nodeToPartitionMap[left];
			}
			rootNodeToPartitionMap[i][1] = nodeToPartitionMap[left];
			rootNodeToPartitionMap[i][2] = nodeToPartitionMap[right];
			rootNode[nodeToPartitionMap[left]] = i;
			rootNode[nodeToPartitionMap[right]] = i;
		}
		
		if (logPContributions == null) {
			logPContributions = new double[partitionCount];
			storedLogPContributions = new double[partitionCount];
		}
	}
	
	private int firstSampledAncestor(Node node) {
		if (node.isRoot()) {
			return node.getNr(); 
		}
		if (isSampled[node.getNr()]) {
			return node.getNr();
		}
		return firstSampledAncestor(node.getParent());
	}

	private void initNodeToPartitionMap(Node node, int [] nextParitionNr, int currentPartition) {
		int nodeNr = node.getNr();
		nodeToPartitionMap[nodeNr] = currentPartition;
		if (isSampled[nodeNr]) {
			for (Node child : node.getChildren()) {
				nextParitionNr[0]++;
				initNodeToPartitionMap(child, nextParitionNr, nextParitionNr[0]);
			}
		} else {
			for (Node child : node.getChildren()) {
				initNodeToPartitionMap(child, nextParitionNr, currentPartition);
			}
		}	
	}

	private boolean nodeToPartitionMapChanged(Node node, int [] nextParitionNr, int currentPartition) {
		int nodeNr = node.getNr();
		if (nodeToPartitionMap[nodeNr] != currentPartition) {
			return true;
		}
		if (isSampled[nodeNr]) {
			for (Node child : node.getChildren()) {
				nextParitionNr[0]++;
				if (nodeToPartitionMapChanged(child, nextParitionNr, nextParitionNr[0])) {
					return true;
				}
			}
		} else {
			for (Node child : node.getChildren()) {
				if (nodeToPartitionMapChanged(child, nextParitionNr, currentPartition)) {
					return true;
				}
			}
		}
		return false;
	}


	@Override
	public double calculateLogP() throws Exception {
		if (!initialised) {
			initialiseSampledStates();
			initialised = isMonoPhyletic;
		}
		
		if (!isMonoPhyletic) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}
		
        logP = Double.NaN;
		try {
			// check prior
			if (sampledLocations != null) {
				for (GeoPrior prior : geopriorsInput.get()) {
					if (Double.isInfinite(prior.calculateLogP())) {
						logP = Double.NEGATIVE_INFINITY;
						return logP;
					}
				}
			}
			
			// calc likelihood
			logP = 0.0;
			// calcBranchLengths(); already done in requiresRecalculation
			
			
			// only update dirty partitions
			for (int i : dirtyPartitionList) {
				Node node = tree.getNode(rootNode[i]);
				if (node.isRoot() && !isSampled[node.getNr()]) {
					calcPositions(node);
				} else if (nodeToPartitionMap[node.getLeft().getNr()] == i) {
					calcPositions(node.getLeft());
				} else {
					if (nodeToPartitionMap[node.getRight().getNr()] != i) {
						throw new RuntimeException("Programmer error: one of nodeToPartitionMaps should be equal to " + i);
					}
					calcPositions(node.getRight());
				}
			}

			for (int i : dirtyPartitionList) {
				Node node = tree.getNode(rootNode[i]);
				if (node.isRoot() && !isSampled[node.getNr()]) {
					logPContributions[i] = calcLogPContribution(node);
				} else if (nodeToPartitionMap[node.getLeft().getNr()] == i) {
					logPContributions[i] = calcLogPContribution(node.getLeft());
				} else {
					logPContributions[i] = calcLogPContribution(node.getRight());					
				}
			}

			// logP = sum of contributions of all partitions
			for (double c : logPContributions) {
				logP += c;
				//System.err.println(logP);
			}
			
// TODO: slightly more efficient way to sum instead of the loop just above		
// what if everything is set to dirty?
//			logP = storedLogP;
//			for (int i : dirtyPartitionList) {
//				logP += logPContributions[i] - storedLogPContributions[i];
//			}
			needsUpdate = false;
		} catch (Exception e) {
			e.printStackTrace();
		}
		//System.err.print("locP2(" + logP +")\n");
		sanitycheck();
		return logP;
	}

	private void calcPositions(Node node) {		
		precision = substModel.precisionInput.get().getValue();
		initByMean(node);
		resetMeanDown(node);
		
			
		if (scaleByBranchLength) {
			// not implemented yet
			throw new RuntimeException("scaleByBrancheLength=true is not implemented yet");
		}
		
		// TODO: loop over nodes in partition stored in List<Integer>[partitionCount] instead of over all internal nodes
		int currentPartition = nodeToPartitionMap[node.getNr()];
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			if (nodeToPartitionMap[i] == currentPartition && !isSampled[i]) {
				position[i] = SphericalDiffusionModel.cartesian2Sperical(sphereposition[i], true);
			}
		}
		
	}

	@Override
	void initByMean(Node node) {
		if (!node.isLeaf()) {
			int nodeNr = node.getNr();
			int child1 = node.getLeft().getNr();
			int child2 = node.getRight().getNr();

			if (!isSampled[child1]) {
				initByMean(node.getLeft());
			}
			if (!isSampled[child2]) {
				initByMean(node.getRight());
			}
			
			if (!isSampled[nodeNr]) {
				setHalfWayPosition(nodeNr, child1, child2);
			}
		}
	}		

	@Override
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

			if (!isSampled[child1]) {
				resetMeanDown(node.getLeft());
			}
			if (!isSampled[child2]) {
				resetMeanDown(node.getRight());
			}
		}
	}


	/** traverse partition of a tree **/
	private double calcLogPContribution(Node node) {
		double logP = 0.0;
		if (!node.isRoot()) {
			logP += substModel.getLogLikelihood(node, position, branchLengths);
		}
		if (!isSampled[node.getNr()]) {
			for (Node child : node.getChildren()) {
				logP += calcLogPContribution(child);
			}
		}
		return logP;
	}

	
	@Override
	public void store() {
        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);

		System.arraycopy(logPContributions, 0, storedLogPContributions, 0, partitionCount);
		
		// TODO: instead of copying all positions, only copy those of the partitions that changed in the last update
		double [] p, sp;
		for (int i = 0; i < position.length; i++) {
			p = position[i];
			sp = storedPosition[i];
			sp[0] = p[0];
			sp[1] = p[1];
		}
		for (int i = 0; i < position.length; i++) {
			p = sphereposition[i];
			sp = storedSphereposition[i];
			sp[0] = p[0];
			sp[1] = p[1];
			sp[2] = p[2];
		}

		super.store();
		
		wasInitialised = false;
	}
	
	@Override
	public void restore() {
		super.restore();
        
		double[] tmp = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
        
        tmp = logPContributions;
        logPContributions = storedLogPContributions;
        storedLogPContributions = tmp;
        
        double [][]tmp2 = position;
        position = storedPosition;
        storedPosition = tmp2;
        
        tmp2 = sphereposition;
        sphereposition = storedSphereposition;
        storedSphereposition = tmp2;
        
        if (wasInitialised) {
			initialiseSampledStates();
        }
		wasInitialised = false;
	}
	
	@Override
	protected boolean requiresRecalculation() {
		if (!initialised) {
			initialiseSampledStates();
		}
		boolean [] dirtyPartitions = new boolean[partitionCount];
		
		needsUpdate = true;
		if (loggerLikelihood != null) {
			loggerLikelihood.needsUpdate = true;
		}
		
		if (geoPriorChanged()) {
			initialiseSampledStates();
			Arrays.fill(dirtyPartitions, true);
		}
		
		if (tree.somethingIsDirty()){
			int [] nextParitionNr = new int[1];
			nextParitionNr[0] = isSampled[tree.getRoot().getNr()] ? -1 : 0;
			if (nodeToPartitionMapChanged(tree.getRoot(), nextParitionNr, 0)) {
				initialiseSampledStates();
				Arrays.fill(dirtyPartitions, true);
			}
		}

		if (substModel.isDirtyCalculation()) {
			Arrays.fill(dirtyPartitions, true);
		}
		
		if (((CalculationNode) clockModel).isDirtyCalculation() || tree.somethingIsDirty()) {
			calcBranchLengths();
			for (int i = 0; i < branchLengths.length; i++) {
				if (Math.abs(branchLengths[i] - storedBranchLengths[i]) > 1e-10) {
					if (isTopOfPartition[i]) {
						dirtyPartitions[rootNodeToPartitionMap[i][0]] = true;
						dirtyPartitions[rootNodeToPartitionMap[i][1]] = true;
						dirtyPartitions[rootNodeToPartitionMap[i][2]] = true;
					} else {
						dirtyPartitions[nodeToPartitionMap[i]] = true;
					}
				}
			}
		}
		
		if (sampledLocations.somethingIsDirty()) {
			for (int i : sampleNumber) {
				double lat1  = sampledLocations.getMatrixValue(i, 0);
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
					if (isTopOfPartition[i]) {
						dirtyPartitions[rootNodeToPartitionMap[i][0]] = true;
						dirtyPartitions[rootNodeToPartitionMap[i][1]] = true;
						dirtyPartitions[rootNodeToPartitionMap[i][2]] = true;
					} else {
						throw new RuntimeException("Programmer error: Sampled nodes should be root of partitions");
					}
				}
			}
		}
		
		dirtyPartitionList.clear();
		for (int i = 0; i < partitionCount; i++) {
			if (dirtyPartitions[i]) {
				dirtyPartitionList.add(i);
			}
		}
		return true;
	}

}
