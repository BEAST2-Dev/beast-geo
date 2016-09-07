package sphericalGeo;



import java.util.*;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;


@Description(value="Approximate likelihood by MAP approximation of internal states optimised to take sampled nodes in account. "
		+ "Efficiently calculates geo priors as well.", 
			isInheritable=false)
@Citation(value="Remco R. Bouckaert. Phylogeography by diffusion on a sphere: whole world phylogeography. 2016, PeerJ 4:e2406 https://doi.org/10.7717/peerj.2406", DOI="10.7717/peerj.2406", year=2016)
public class ApproxMultivariateTraitLikelihoodF2 extends ApproxMultivariateTraitLikelihoodF { 
	
	/** maps each node to the partition containing the node **/
	int [] nodeToPartitionMap;
	int [] storedNodeToPartitonMap;
	
	/** flag to indicate for each node whether it is a 'root' for a partition **/
	boolean [] isTopOfPartition;
	boolean [] storedIsTopOfPartition;
	/** array of root nodes, one for each partition **/
	int [] rootNode;
	int [] storedRootNode;
	
	/** partition root nodes can affect 3 (or at least 2 for global root) **/
	int [][] rootNodeToPartitionMap;
	int [][] storedRootNodeToPartitionMap;
	int partitionCount;
	//int storedPartitionCount;
	
	/** contribution to logP by each of the partitions **/
	double [] logPContributions;
	double [] storedLogPContributions;
	
	
	int [/*partitionCount*/][/*2*/] partitionNrCandidates;
	
	/** list of partitions that need recalculating **/
	List<Integer> dirtyPartitionList = new ArrayList<>();
	
	List<Integer> storedDirtyPartitionList = new ArrayList<>();
	List<Integer> storedDirtyPartitionList2 = new ArrayList<>();
	
	
	boolean wasInitialised;
	

	
    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] storedBranchLengths;
	//double [] storedSumLengths;
	double [] storedParentweight;
	
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		branchLengths = new double[tree.getNodeCount()];
		storedBranchLengths = new double[tree.getNodeCount()];
		storedParentweight = new double[tree.getNodeCount()];
		
	}
	
	@Override
	void initialiseSampledStates() {
		wasInitialised = true;
		//System.err.print("M");
		
		super.initialiseSampledStates();
		if (!isMonoPhyletic) {
			return;
		}

		
		// initialise partition information
		if (isTopOfPartition == null) {
			isTopOfPartition = new boolean[tree.getNodeCount()];
			storedIsTopOfPartition = new boolean[tree.getNodeCount()];
			nodeToPartitionMap = new int[tree.getNodeCount()];
			storedNodeToPartitonMap = new int[tree.getNodeCount()];
			rootNodeToPartitionMap = new int[tree.getNodeCount()][3];
			storedRootNodeToPartitionMap = new int[tree.getNodeCount()][3];
		}
		Arrays.fill(isTopOfPartition, false);
		for (int i : taxonNrs) {
			isTopOfPartition[i] = true;
		}
		Arrays.fill(nodeToPartitionMap, -1);
		
		int [] nextParitionNr = new int[1];
		nextParitionNr[0] = isSampled[tree.getRoot().getNr()] ? -1 : 0;
		if (partitionNrCandidates == null) {
			partitionNrCandidates = new int[taxonNrs.length][];
		}
		initNodeToPartitionMap(tree.getRoot(), nextParitionNr, 0);
		partitionCount = Math.max(partitionCount, nextParitionNr[0] + 1);

		if (rootNode == null) {
			rootNode = new int[partitionCount];
			storedRootNode = new int[partitionCount];
		}
		
		for (int i: taxonNrs) {
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
		List<Integer> sampleNumber = asList(taxonNrs);
		if (isSampled[nodeNr]) {
			int k = sampleNumber.indexOf(nodeNr);
			if (partitionNrCandidates[k] == null) {
				partitionNrCandidates[k] = new int[node.getChildCount()];
				int i = 0;
				for (Node child : node.getChildren()) {
					nextParitionNr[0]++;
					partitionNrCandidates[k][i] = nextParitionNr[0];
					initNodeToPartitionMap(child, nextParitionNr, nextParitionNr[0]);
					i++;
				} 
			} else {
				int i = 0;
				for (Node child : node.getChildren()) {
					nextParitionNr[0] = partitionNrCandidates[k][i];
					initNodeToPartitionMap(child, nextParitionNr, nextParitionNr[0]);
					i++;
				} 
			}
		} else {
			for (Node child : node.getChildren()) {
				initNodeToPartitionMap(child, nextParitionNr, currentPartition);
			}
		}	
	}
	  private List<Integer> asList(final int[] is)
	    {
	            return new AbstractList<Integer>() {
	                    public Integer get(int i) { return is[i]; }
	                    public int size() { return is.length; }
	            };
	    }
	
	private boolean nodeToPartitionMapChanged(Node node, int [] nextParitionNr, int currentPartition, boolean[] dirtyPartitions) {
		int nodeNr = node.getNr();
		boolean changed = false;
		if (nodeToPartitionMap[nodeNr] != currentPartition) {
			dirtyPartitions[currentPartition] = true;
			dirtyPartitions[nodeToPartitionMap[nodeNr]] = true;
			changed = true;
		}
		List<Integer> sampleNumber = asList(taxonNrs);
		if (isSampled[nodeNr]) {
			int k = sampleNumber.indexOf(nodeNr);			
			int i = 0;
			for (Node child : node.getChildren()) {
				nextParitionNr[0] = partitionNrCandidates[k][i];
				if (nodeToPartitionMapChanged(child, nextParitionNr, nextParitionNr[0], dirtyPartitions)) {
					changed = true;
				}
				i++;
			}
		} else {
			for (Node child : node.getChildren()) {
				if (nodeToPartitionMapChanged(child, nextParitionNr, currentPartition, dirtyPartitions)) {
					changed = true;
				}
			}
		}
		return changed;
	}


	@Override
	public double calculateLogP() {
		if (!initialised) {
			initialiseSampledStates();
			initialised = isMonoPhyletic;
		}
		
		preRecalculation();
		
		if (!isMonoPhyletic) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}

		
		
		logP = 0.0;
		try {
			// calc prior
			if (sampledLocations != null) {
				logP = multiGeopriorsInput.get().calculateLogP();
				if (Double.isInfinite(logP)) {
					return logP;
				}
			}
			
			// calc likelihood
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
//for (int i = 0; i < position.length; i++) {
//	System.err.print("["+position[i][0] + "," + position[i][1] + "]");
//}

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
			double[] l = logPContributions.clone();
            Arrays.sort(l);
            for(int k = l.length-1; k >= 0; --k) {
                logP += l[k];
            }
						
//			for (double c : logPContributions) {
//				logP += c;
//				//System.err.println(logP);
//			}
			
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
		
		storedDirtyPartitionList2.clear();
		storedDirtyPartitionList2.addAll(storedDirtyPartitionList);
		storedDirtyPartitionList.clear();
		storedDirtyPartitionList.addAll(dirtyPartitionList);
		
		
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
			} else {
				parentweight[child1] = 0;
			}
			if (!isSampled[child2]) {
				initByMean(node.getRight());
			} else {
				parentweight[child2] = 0;
			}
			
			if (!isSampled[nodeNr]) {
				setHalfWayPosition(nodeNr, child1, child2);
			} else {
				parentweight[nodeNr] = 0;
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
		//System.err.print("store");
		super.store();

		System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
		System.arraycopy(parentweight, 0, storedParentweight, 0, parentweight.length);

		//if (storedLogPContributions.length != logPContributions.length) {
		//	storedLogPContributions = new double[logPContributions.length];
		//}

		System.arraycopy(logPContributions, 0, storedLogPContributions, 0, logPContributions.length);
		
		System.arraycopy(isTopOfPartition, 0, storedIsTopOfPartition, 0, isTopOfPartition.length);
		System.arraycopy(nodeToPartitionMap, 0, storedNodeToPartitonMap, 0, nodeToPartitionMap.length);
		System.arraycopy(rootNode, 0, storedRootNode, 0, rootNode.length);
		
		int [] p, sp;
		for (int i = 0; i < rootNodeToPartitionMap.length; i++) {
			p = rootNodeToPartitionMap[i];
			sp= storedRootNodeToPartitionMap[i];
			sp[0] = p[0];
			sp[1] = p[1];
			sp[2] = p[2];
		}

//		// TODO: instead of copying all positions, only copy those of the partitions that changed in the last update
//		double [] p, sp;
//		for (int i = 0; i < position.length; i++) {
//			p = position[i];
//			sp = storedPosition[i];
//			sp[0] = p[0];
//			sp[1] = p[1];
//		}
//		for (int i = 0; i < position.length; i++) {
//			p = sphereposition[i];
//			sp = storedSphereposition[i];
//			sp[0] = p[0];
//			sp[1] = p[1];
//			sp[2] = p[2];
//		}
		wasInitialised = false;
	}
	
	@Override
	public void restore() {
		//System.err.print("restore");
		super.restore();
        
		double[] tmp = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
        
		tmp = parentweight;
		parentweight = storedParentweight;
		storedParentweight = tmp;

        tmp = logPContributions;
        logPContributions = storedLogPContributions;
        storedLogPContributions = tmp;
        
        int [][]tmp3 = rootNodeToPartitionMap;
        rootNodeToPartitionMap = storedRootNodeToPartitionMap;
        storedRootNodeToPartitionMap = tmp3;
    
		boolean [] tmp4 = isTopOfPartition;
		isTopOfPartition = storedIsTopOfPartition;
		storedIsTopOfPartition = tmp4;
		
		int [] tmp5 = nodeToPartitionMap;
		nodeToPartitionMap = storedNodeToPartitonMap;
		storedNodeToPartitonMap = tmp5;
		
		tmp5 = rootNode;
		rootNode = storedRootNode;
		storedRootNode = tmp5;
		
//		  double [][]tmp2 = position;
//        position = storedPosition;
//        storedPosition = tmp2;
//        
//        tmp2 = sphereposition;
//        sphereposition = storedSphereposition;
//        storedSphereposition = tmp2;
        if (wasInitialised) {
			//initialiseSampledStates();
        }
		wasInitialised = false;
    }
	
	protected boolean geoPriorChanged(boolean[] dirtyPartitions) {
		MultiGeoPrior multiPrior = multiGeopriorsInput.get();

		boolean changed = false;
		for (int k = 0; k < multiPrior.size(); k++) {
			int taxonNr = multiPrior.getCladeTopNodeNr(k);
			if (taxonNr != taxonNrs[k]) {
				int partition = nodeToPartitionMap[taxonNrs[k]];
				int rootNode = this.rootNode[partition];
				dirtyPartitions[rootNodeToPartitionMap[rootNode][0]] = true;
				dirtyPartitions[rootNodeToPartitionMap[rootNode][1]] = true;
				dirtyPartitions[rootNodeToPartitionMap[rootNode][2]] = true;

				// need updated nodeToPartitionMap, rootNode and rootNodeToPartitionMap for this
				priorChanged.add(k);
				changed = true;
			}
		}
		return changed;
	}

	List<Integer> priorChanged = new ArrayList<>();
	
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}
	
	protected boolean preRecalculation() {
		boolean [] dirtyPartitions = new boolean[partitionCount];
		if (!initialised) {
			initialiseSampledStates();
			dirtyPartitions = new boolean[partitionCount];
			Arrays.fill(dirtyPartitions, true);
		}
		
		needsUpdate = true;
		if (loggerLikelihood != null) {
			loggerLikelihood.needsUpdate = true;
		}
		
		boolean needsReInit = false;
		priorChanged.clear();
		if (geoPriorChanged(dirtyPartitions)) {
			needsReInit = true;
			//dirtyPartitions = new boolean[partitionCount];
			//Arrays.fill(dirtyPartitions, true);
			//System.err.print("c");
		} else {
			//System.err.print("u");
		}
		
		if (tree.somethingIsDirty() && needsReInit == false){
			// make sure something is filthy
			boolean isFilthy = false;
			for (Node node : tree.getNodesAsArray()) {
				if (node.isDirty() == Tree.IS_FILTHY) {
					isFilthy = true;
					break;
				}
			}
			if (isFilthy) {
				// now do the more expensive check that partition numbers changed
				int [] nextPartitionNr = new int[1];
				nextPartitionNr[0] = isSampled[tree.getRoot().getNr()] ? -1 : 0;
				if (nodeToPartitionMapChanged(tree.getRoot(), nextPartitionNr, 0, dirtyPartitions)) {
					needsReInit = true;
					//Arrays.fill(dirtyPartitions, true);
				}
			}
		}
		if (needsReInit) {
			initialiseSampledStates();
			
			for (int k : priorChanged) {
				int partition = nodeToPartitionMap[taxonNrs[k]];
				int rootNode = this.rootNode[partition];
				dirtyPartitions[rootNodeToPartitionMap[rootNode][0]] = true;
				dirtyPartitions[rootNodeToPartitionMap[rootNode][1]] = true;
				dirtyPartitions[rootNodeToPartitionMap[rootNode][2]] = true;
			}
		}
		
		if (substModel.isDirtyCalculation()) {
			Arrays.fill(dirtyPartitions, true);
		}

		if (((CalculationNode) clockModel).isDirtyCalculation() || tree.somethingIsDirty()) {
			calcBranchLengths();
			//System.err.print("t");
			for (int i = 0; i < branchLengths.length; i++) {
				if (branchLengths[i] != storedBranchLengths[i]) {
					//System.err.print(" " + i);
					if (isTopOfPartition[i]) {
						dirtyPartitions[rootNodeToPartitionMap[i][0]] = true;
						dirtyPartitions[rootNodeToPartitionMap[i][1]] = true;
						dirtyPartitions[rootNodeToPartitionMap[i][2]] = true;
						//System.err.print("["+rootNodeToPartitionMap[i][0] + " " + rootNodeToPartitionMap[i][1] + " " + rootNodeToPartitionMap[i][2] + "]");
					} else {
						dirtyPartitions[nodeToPartitionMap[i]] = true;
						//System.err.print("["+nodeToPartitionMap[i] + "]");
					}
					if (storedIsTopOfPartition[i]) {
						dirtyPartitions[storedRootNodeToPartitionMap[i][0]] = true;
						dirtyPartitions[storedRootNodeToPartitionMap[i][1]] = true;
						dirtyPartitions[storedRootNodeToPartitionMap[i][2]] = true;
						//System.err.print("["+storedRootNodeToPartitionMap[i][0] + " " + storedRootNodeToPartitionMap[i][1] + " " + storedRootNodeToPartitionMap[i][2] + "]");
					} else {
						dirtyPartitions[storedNodeToPartitonMap[i]] = true;						
						//System.err.print("["+storedNodeToPartitonMap[i] + "]");
					}
				}
			}
		}
		
//		System.err.println(dirtyPartitions[7] + " " + dirtyPartitions[8]);
		
		if (sampledLocations.somethingIsDirty()) {
			for (int i : taxonNrs) {
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
		
		//System.err.println(Arrays.toString(dirtyPartitions));
		//System.err.println(dirtyPartitions);
		
		dirtyPartitionList.clear();
		for (int i = 0; i < partitionCount; i++) {
			if (dirtyPartitions[i]) {
				dirtyPartitionList.add(i);
			}
		}
//		System.err.print(partitionCount + ">>>" + dirtyPartitionList);
		return true;
	}

	public double getNodeLogP(Node node) {
		if (node.isRoot()) {
			return 0.0;
		}
		return substModel.getLogLikelihood(node, position, branchLengths);
	}

	
	
	
//	@Override
//	public String getID() {
//		if (logPContributions == null) {
//			return super.getID();
//		}
//		StringBuilder buf = new StringBuilder();
//		for (int i = 0; i < logPContributions.length; i++) {
//			buf.append(storedLogPContributions[i] + " " + logPContributions[i] + " " + (storedLogPContributions[i] - logPContributions[i]) + "\n");
//		}
//		
//		buf.append(dirtyPartitionList.toString().replaceAll("[\\[\\]]",""));
//		buf.append("<=\n");
//		buf.append(storedDirtyPartitionList.toString().replaceAll("[\\[\\]]",""));
//		buf.append("<=\n");
//		buf.append(storedDirtyPartitionList2.toString().replaceAll("[\\[\\]]",""));
//		buf.append("<=\n");
//		return buf.toString();
//	}
	
}
