package sphericalGeo;




import java.util.*;

import beast.core.BEASTInterface;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;


@Description("Approximate likelihood by MAP approximation of internal states")
@Citation("Remco R. Bouckaert. Phylogeography by diffusion on a sphere. bioRxiv, BIORXIV/2015/016311, 2015.")
public class ApproxMultivariateTraitLikelihood extends GenericTreeLikelihood implements StateNodeInitialiser { 
	public Input<Boolean> scaleByBranchLengthInput = new Input<>("scale", "scale by branch lengths for initial position", false);
	public Input<List<GeoPrior>> geopriorsInput = new Input<>("geoprior", "geographical priors on tips, root or clades restricting these nodes to a region", new ArrayList<GeoPrior>());
	public Input<LocationParameter> locationInput = new Input<>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree");

	public Input<Transformer> transformerInput = new Input<>("transformer","landscape transformer to capture some inheterogenuity in the diffusion process");
	public Input<Boolean> logAverageInput = new Input<>("logAverage", "when logging, use average position instead of sample from particle filter. "
			+ "This is faster, but also artificially reduces uncertainty in locations. ", false);


	SphericalDiffusionModel substModel;
	TreeInterface tree;
	BranchRateModel clockModel;
	Transformer transformer;
	
	double [][] position;
	double [][] storedPosition;

	double [][] sphereposition;
	double [][] storedSphereposition;
	
	double [] branchLengths;
	//double [] sumLengths;
	double [] parentweight;

	public boolean needsUpdate = true;
	boolean scaleByBranchLength;
	boolean logAverage;
	
	PFApproxMultivariateTraitLikelihood loggerLikelihood;
	
	LocationParameter sampledLocations;
	boolean [] isSampled;
	boolean [] storedIsSampled;
	List<Integer> sampleNumber;
	List<Integer> storedSampleNumber;
	
	double precision;
	int [] taxonNrs;
	int [] storedTaxonNrs;
	
	boolean isMonoPhyletic = false;
	boolean storedIsMonoPhyletic = false;
	
	@Override
	public void initAndValidate() throws Exception {
		super.initAndValidate();
		clockModel = branchRateModelInput.get();
		if (clockModel == null) {
			clockModel = new StrictClockModel();
		}
		SiteModel siteModel = (SiteModel) siteModelInput.get();
		substModel = (SphericalDiffusionModel) siteModel.substModelInput.get();
		tree = treeInput.get();
		scaleByBranchLength = scaleByBranchLengthInput.get();
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
		
		sphereposition = new double[position.length][3];
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
		}
		
		storedPosition = new double[tree.getNodeCount()][2];
		storedSphereposition = new double[tree.getNodeCount()][3];


		branchLengths = new double[tree.getNodeCount()];
		//sumLengths = new double[tree.getNodeCount()];
		parentweight = new double[tree.getNodeCount()];
	
		
		List<GeoPrior> geopriors = geopriorsInput.get();
		{
			// geopriors.size: sanity check to see whether there really are no geo-priors
			String ids ="";
			for (Object o : ((BEASTInterface) tree).getOutputs()) {
				if (o instanceof GeoPrior &&  !geopriors.contains(o)) {
					ids += ((BEASTInterface)o).getID() + ", ";
				}
			}
			if (ids.length() > 1) {
				ids = ids.substring(0, ids.length() - 2);
				Log.warning.println("\nWARNING: this analysis contains GeoPriors (" + ids + "), but these are not "
						+ "connected to the ApproxMultivariateTraitLikelihood (" + getID() + ").");
				Log.warning.println("For every GeoPrior, there should be a geoprior entry.");
				Log.warning.println("Expect this analysis to fail.\n");
			}
		}

//		initialiseSampledStates();

		logAverage = logAverageInput.get();
	}
	
	boolean initialised = false;
	
	
	void initialiseSampledStates() {
		
		List<GeoPrior> geopriors = geopriorsInput.get();
		if (isSampled == null) {
			isSampled = new boolean[tree.getNodeCount()];
			storedIsSampled = new boolean[tree.getNodeCount()];
			sampleNumber = new ArrayList<>();
			storedSampleNumber = new ArrayList<>();
			taxonNrs = new int[geopriors.size()];
			storedTaxonNrs = new int[geopriors.size()];
		}
		Arrays.fill(isSampled, false);
		sampleNumber.clear();
		Arrays.fill(taxonNrs, 0);
		
		if (geopriors.size() > 0) {
			sampledLocations = locationInput.get();
			if (sampledLocations == null) {
				Log.warning.println("\nWARNING: 'location' needs to be specified when geopriors are defined, but location=null\n");
			}
			Double [] d = new Double[sampledLocations.getDimension()];
			for (int i = 0; i < d.length; i++) {
				d[i] = sampledLocations.getValue(i);
			}
			Double [] d2 = d.clone();
			
			int k = 0;
			for (GeoPrior prior : geopriors) {
				prior.initialise();
				if (prior.allInternalNodes) {
					for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
						isSampled[i] = true;
						sampleNumber.add(i);
						double [] location = prior.sample();
						// check if the location is already initialised (e.g. through resuming a chain)
						if (Math.abs(d[i * 2]) < 1e-10 && Math.abs(d[i * 2 + 1]) < 1e-10) {
							d[i * 2] = location[0];
							d[i * 2 + 1] = location[1];
							setPosition(i, location[0], location[1]);
						} else {
							if (!initialised) {
								Log.warning.println("location of " + prior.getID() + " already set at " + d[i*2]+","+d[i*2+1]);
							}
						}
					}
					taxonNrs[k] = -1;
				} else {
					if (prior.isMonoPhyletic()) {
						int taxonNr = prior.getTaxonNr();
						isSampled[taxonNr] = true;
						sampleNumber.add(taxonNr);
						int storedTaxonNr = prior.getStoredTaxonNr();
						if (storedTaxonNr >= 0) {
							d[taxonNr * 2]     = d2[storedTaxonNr * 2];
							d[taxonNr * 2 + 1] = d2[storedTaxonNr * 2 + 1];						
							setPosition(taxonNr, d2[storedTaxonNr * 2], d2[storedTaxonNr * 2 + 1]);
						} else {
							double [] location = prior.sample();
							// check if the location is already initialised (e.g. through resuming a chain)
							if (Math.abs(d[taxonNr * 2]) < 1e-10 && Math.abs(d[taxonNr * 2 + 1]) < 1e-10) {
								d[taxonNr * 2] = location[0];
								d[taxonNr * 2 + 1] = location[1];
								setPosition(taxonNr, location[0], location[1]);
							} else {
								if (!initialised) {
									Log.warning.println("location of " + prior.getID() + " already set at " + d[taxonNr*2]+","+d[taxonNr*2+1]);
								}
							}
						}
						taxonNrs[k] = taxonNr;
					} else {
						// only works when clades are monophyletic, so we give up here
						isMonoPhyletic = false;
						taxonNrs[0] = - (1+storedTaxonNrs[0]);
						return;
					}
				}
				k++;
			}
			isMonoPhyletic = true;
			// RRB: assignFromWithoutID overrides sampledLocations.storedValues, so sampledLocations is invalid after restore
			// RealParameter tmp = new RealParameter(d);
			// sampledLocations.assignFromWithoutID(tmp);
			sampledLocations.setValueSilently(d);
			//System.err.println("Initialised");
		} else {
			isMonoPhyletic = true;
		}
		

		if (initialised || logAverage) {
			return;
		}

		SiteModel siteModel = (SiteModel) siteModelInput.get();
		AlignmentFromTraitMap data = (AlignmentFromTraitMap) dataInput.get();

		loggerLikelihood = new PFApproxMultivariateTraitLikelihood();
		try {
			if (geopriors.size() > 0) {
				loggerLikelihood.initByName("scale", scaleByBranchLength, "tree", tree, "siteModel", siteModel, 
					"branchRateModel", clockModel, "data", data,
					"geoprior", geopriors, "location", sampledLocations, "transformer", transformer);
			} else {
				loggerLikelihood.initByName("scale", scaleByBranchLength, "tree", tree, "siteModel", siteModel, 
						"branchRateModel", clockModel, "data", data, "transformer", transformer);
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	
	
	private void setPosition(int i, double lat, double long_) {
		//System.err.print(i+"(" + lat + "," + long_ + ")");
		this.position[i][0] = lat;
		this.position[i][1] = long_;
		sphereposition[i] = SphericalDiffusionModel.spherical2Cartesian(lat, long_);		
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

        // check prior
		if (sampledLocations != null) {
			for (GeoPrior prior : geopriorsInput.get()) {
				if (Double.isInfinite(prior.calculateLogP())) {
					prior.calculateLogP();
					logP = Double.NEGATIVE_INFINITY;
					return logP;
				}
			}
		}
		// calc likelihood
		logP = 0.0;
		calcBranchLengths();
		calcPositions();
		logP = calcLogP();
		needsUpdate = false;

		//System.err.println("\nlocP(" + logP +") ");
		sanitycheck();
		return logP;
	}

	void sanitycheck() {
		for (int i : sampleNumber) {
			if (position[i][0] != sampledLocations.getValue(i*2) ||
				position[i][1] != sampledLocations.getValue(i*2+1)) {
				System.err.println(position[i][0] +"!="+ sampledLocations.getValue(i*2));
				System.err.println(position[i][1] +"!="+ sampledLocations.getValue(i*2+1));
				System.exit(0);
			}
		}
//		for (int i = 0; i < position.length; i++) {
//			System.err.print("["+position[i][0] + "," + position[i][1] + "]");
//		}
//		System.err.println();
//		System.err.println(tree.getRoot().toNewick());
//		System.err.println(sampleNumber.toString());
//		System.err.println(Arrays.toString(branchLengths));
//		System.err.println(Arrays.toString(isSampled));
//		System.err.println(Arrays.toString(parentweight));

		
	}

	void calcBranchLengths() {
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

	/** traverse tree **/
	double calcLogP() {
		double logP = 0;
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				double d = substModel.getLogLikelihood(node, position,
						branchLengths);
				//System.err.print(d+" ");
				logP += d;
			}
		}
		//System.err.println();
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
	
	



	void calcPositions() {
		final double EPSILON = 1e-8;
		
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
		
		precision = substModel.precisionInput.get().getValue();
		initByMean(tree.getRoot());
		//for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			//System.err.print("[" + sphereposition[i][0]+","+sphereposition[i][1]+","+sphereposition[i][2]+"]");
			//System.err.println("pos153A=["+position[153][0] + "," + position[153][1] + ":" + sphereposition[153][0]+","+sphereposition[153][1]+","+sphereposition[153][2]+"]");
		//}
		//System.err.println();
		//System.err.print(Arrays.toString(parentweight));
		resetMeanDown(tree.getRoot());
		
			
		if (scaleByBranchLength) {
			double [][] oldPosition = new double[tree.getNodeCount()][3];
			for (int i = 0; i < 50; i++) {
				if (i % 5 == 0) {
					for (int j = tree.getLeafNodeCount(); j < oldPosition.length; j++) {
						oldPosition[j][0] = sphereposition[j][0];
						oldPosition[j][1] = sphereposition[j][1];
						oldPosition[j][2] = sphereposition[j][2];
					}
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
			if (!isSampled[i]) {
				position[i] = SphericalDiffusionModel.cartesian2Sperical(sphereposition[i], true);
			} else {
				//System.err.print("[skip " + i +"]");
			}
		}
		//System.err.print("pos153B=["+position[153][0] + "," + position[153][1] + ":" + sphereposition[153][0]+","+sphereposition[153][1]+","+sphereposition[153][2]+"]");
		//System.err.println();
		
	}
	
	void setHalfWayPosition(int nodeNr, int child1, int child2) {
		// start in weighted middle of the children
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
	}

	void setHalfWayPosition(int nodeNr, int child1, int child2, int parent) {
		// start in weighted middle of the children and parent location
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
			} else {
				parentweight[nodeNr] = 0;
				//System.err.print("skip" + nodeNr);
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


	@Override
	public boolean isStochastic() {
		return !logAverage;
	}

	/** return randomized position **/
	public double[] getPosition(int iDim) {
		if (!logAverage) {
			return loggerLikelihood.getPostion(iDim);
		} else {
			if (needsUpdate) {
				try {
					if (geoPriorChanged()) {
						initialiseSampledStates();
					}
					calculateLogP();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			if (transformer != null) {
				double [] p = transformer.projectInverse(position[iDim][0], position[iDim][1]);
				return p;
			}
			sanitycheck();
			return position[iDim];
		}
	}

	/** return non-randomized positions **/
	public double [][] getPositions() {
		if (needsUpdate) {
			try {
				calculateLogP();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		if (transformer != null) {
			double [][] position0 = new double[position.length][];
			for (int i = 0; i < position.length; i++) {
				position0[i] = transformer.projectInverse(position[i][0], position[i][1]);
			}
			return position0;
		}
		return position;
	}

	
	@Override
	public void store() {
		super.store();

		needsUpdate = true;
		if (loggerLikelihood != null) {
			loggerLikelihood.needsUpdate = true;
		}
		
		System.arraycopy(isSampled, 0, storedIsSampled, 0, isSampled.length);
		System.arraycopy(taxonNrs, 0, storedTaxonNrs, 0, taxonNrs.length);
		
		storedSampleNumber.clear();
		storedSampleNumber.addAll(sampleNumber);

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

		storedIsMonoPhyletic = isMonoPhyletic;
	}
	
	@Override
	public void restore() {
		super.restore();

		needsUpdate = true;
		if (loggerLikelihood != null) {
			loggerLikelihood.needsUpdate = true;
		}
		
		boolean [] tmp = storedIsSampled;
		storedIsSampled = isSampled;
		isSampled = tmp;

		int [] tmp3 = storedTaxonNrs;
		storedTaxonNrs = taxonNrs;
		taxonNrs = tmp3;

		
		List<Integer> tmp4 = storedSampleNumber;
		storedSampleNumber = sampleNumber;
		sampleNumber = tmp4;

        double [][]tmp2 = position;
        position = storedPosition;
        storedPosition = tmp2;
        
        tmp2 = sphereposition;
        sphereposition = storedSphereposition;
        storedSphereposition = tmp2;

		isMonoPhyletic = storedIsMonoPhyletic;

	}
	
	
	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		if (loggerLikelihood != null) {
			loggerLikelihood.needsUpdate = true;
		}
		
		if (initialised && geoPriorChanged()) {
			initialiseSampledStates();
		}
		
		super.requiresRecalculation();
		return true;
	}



	protected boolean geoPriorChanged() {
		int k = 0;
		for (GeoPrior prior : geopriorsInput.get()) {
			int taxonNr = prior.getTaxonNr();
			if (taxonNr != taxonNrs[k]) {
				return true;
			}
			k++;
		}
		return false;
	}



	@Override
	public void initStateNodes() throws Exception {
		initialiseSampledStates();
		initialised = false;
	}



	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(locationInput.get());
	}
}
