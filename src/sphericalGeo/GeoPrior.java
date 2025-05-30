package sphericalGeo;

import java.io.PrintStream;
import java.util.List;
import java.util.Random;

import sphericalGeo.region.Region;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

@Description("Flat prior over a region -- emforces monophyly of the clade (if any)")
public class GeoPrior extends Distribution {
	public Input<Region> regionInput = new Input<>("region", "region to be in (or not, depending on 'isInside' flag). "
			+ "If not specified the MRCA node will be sampled (which may or may not be more efficient), but no restriction "
			+ "is placed on the node.");//, Validate.REQUIRED);
	public Input<Boolean> isInsideInput = new Input<>("isInside", "whether the prior is for being inside the region, instead of outside", true);
	public Input<RealParameter> locationInput = new Input<>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree", Validate.REQUIRED);
	public Input<Tree> treeInput = new Input<>("tree", "beast tree (from which to get the taxon set)", Validate.REQUIRED);

	public Input<Boolean> allInternalNodesInput = new Input<>("allInternalNodes", "if true, apply prior to all internal nodes", false);
	public Input<Taxon> taxonInput = new Input<>("taxon", "taxon associated with this region, if only a tip is restricted. Otherwise use 'taxonset'");
	public Input<TaxonSet> taxonSetInput = new Input<>("taxonset",
			"specify the prior over an internal node that is the MRCA of set of taxa. Select all taxa for the root");//, Validate.XOR, taxonInput);

	Region region;
	RealParameter location;
	Tree tree;
	//TaxonSet taxonSet;
	List<String> names;
	private int storedTaxonNr = -1;
	private int taxonNr = -1;
	
	boolean isMonophyletic = true;
	boolean storedIsMonophyletic = true;
	
	boolean isRoot;
	boolean isTip = false;
	boolean allInternalNodes = false;
	
    // number of taxa in taxon set
    int nrOfTaxa = -1;
    //Set<Integer> cladeSet = null;
    //Set<Integer> storedCladeSet = null;
    
	boolean initialised = false;
	
    // array of indices of taxa
    int[] taxonIndex;
    int[] storedTaxonIndex;
	boolean [] nodesTraversed;
	boolean [] storedNodesTraversed;
    int nseen;

	@Override
	public void initAndValidate() {
		if (taxonInput.get() != null && taxonSetInput.get() != null) {
			throw new RuntimeException("At most one of taxon and taxonset should be specified");
			// if none specified, inputs can still be valid if allinternalnodes = true
		}
		if (taxonInput.get() == null && taxonSetInput.get() == null && !allInternalNodesInput.get()) {
			throw new RuntimeException("Specify one of taxon and taxonset, or set allInternalNodes=true");
		}
		
		
		region = regionInput.get();
		location = locationInput.get();
		if (location.getMinorDimension1() != 2) {
			throw new RuntimeException("expected that location parameter to have minor dimension 2");
		}
		tree = treeInput.get();
		TaxonSet taxonSet = tree.getTaxonset();
		if (location.getDimension() != taxonSet.getTaxonCount() * 4 - 2) {
			Log.warning.println("Setting dimension of location parameter to have dimension 2 time number of taxa - 1 = " + (taxonSet.getTaxonCount() * 4 - 2)
					+ " (from " + location.getDimension() +")");
			location.setDimension(taxonSet.getTaxonCount() * 4 - 2);
		}
		names = taxonSet.asStringList();
		allInternalNodes = allInternalNodesInput.get();
		
		if (taxonSetInput.get() != null) {
			nrOfTaxa = taxonSetInput.get().asStringList().size();
		} else if (taxonInput.get() != null) {
			nrOfTaxa = 1;
		} else {
			// allInternalNodes = true when we got here
			nrOfTaxa = names.size();
		}
		
        //storedCladeSet = new HashSet<>();
        storedNodesTraversed = new boolean[tree.getNodeCount()];

		super.initAndValidate();
		//initialise();
	}

	public int getTaxonNr() {
		if (taxonNr == -1 || nodesTraversed == null) {
			return taxonNr;
		}
		Node node = tree.getNode(taxonNr);
		
//		for (int i : cladeSet) {
//			if (!nodesTraversed[i]) {
//				int h = 3;
//				h++;
//			}
//		}
		
		while (!node.isRoot() && nodesTraversed[node.getParent().getNr()]) {		
			node = node.getParent();
		}
		taxonNr = node.getNr();
		return taxonNr;
	}


	/** 
    * Need delayed initialisation in order for the tree to get set up.
	* If this happens through a StateNodeInitialiser, node numbering can change.
	**/
	public void initialise() {
		if (allInternalNodes) {
		} else if (taxonInput.get() != null) {
			// prior over a single tip
			isTip = true;
			String taxonName = taxonInput.get().getID();
			taxonNr = names.indexOf(taxonName);
			if (taxonNr < 0) {
				throw new RuntimeException("Could not find taxon " + taxonName + ". Typo perhaps?");
			}
		} else {
			// prior over a clade
			TaxonSet taxonset2 = taxonSetInput.get();
			if (taxonset2.getTaxonCount() == names.size()) { //taxonSet.getTaxonCount()) {
				isRoot = true;
				taxonNr = tree.getRoot().getNr();
			} else {

				//isInTaxaSet = new boolean[taxonSet.getTaxonCount()];
				//List<String> names = taxonSet.asStringList();
				if (!initialised) {
					int k = 0;
		            taxonIndex = new int[nrOfTaxa];
		            if (storedTaxonIndex == null) {
		            	storedTaxonIndex = new int[nrOfTaxa];
		            }
		            for (final String sTaxon : taxonset2.asStringList()) {
		                final int iTaxon = names.indexOf(sTaxon);
		                if (iTaxon < 0) {
		                    throw new RuntimeException("Cannot find taxon " + sTaxon + " in data");
		                }
		                //if (isInTaxaSet[iTaxon]) {
		                //    throw new RuntimeException("Taxon " + sTaxon + " is defined multiple times, while they should be unique");
		                //}
		                //isInTaxaSet[iTaxon] = true;
		                taxonIndex[k++] = iTaxon;
		            }
				}
				
				// set up taxonNr
	            isMonophyletic = false;
                Node m;
	            //if(false) {
	            //    calcMRCAtime(tree.getRoot(), new int[1]);
	            //    m = tree.getNode(taxonNr);
	            //} else {

	                nodesTraversed = new boolean[tree.getNodeCount()];
	                nseen = 0;
	                m = getCommonAncestor();
	            	taxonNr = m.getNr();
	                isMonophyletic = nseen == 2 * taxonIndex.length - 1;
	            //}
				
	           //cladeSet = new HashSet<>();
	           //setUpCladeSet(m);
			}
		}
		initialised = true;
	}

    // would be nice to use nodeRef's, but they are not preserved :(
    public Node getCommonAncestor() {
    	if (isTip) {
    		return tree.getNode(taxonNr);
    	}
        Node cur = tree.getNode(taxonIndex[0]);

        for (int k = 1; k < taxonIndex.length; ++k) {
            cur = getCommonAncestor(cur, tree.getNode(taxonIndex[k]));
        }
        return cur;
    }

    private Node getCommonAncestor(Node n1, Node n2) {
        // assert n1.getTree() == n2.getTree();
        if( ! nodesTraversed[n1.getNr()] ) {
           nodesTraversed[n1.getNr()] = true;
            nseen += 1;
        }
        if( ! nodesTraversed[n2.getNr()] ) {
            nodesTraversed[n2.getNr()] = true;
            nseen += 1;
        }
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            if( ! nodesTraversed[n1.getNr()] ) {
	                nodesTraversed[n1.getNr()] = true;
	                nseen += 1;
	            }
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            if( ! nodesTraversed[n2.getNr()] ) {
	                nodesTraversed[n2.getNr()] = true;
	                nseen += 1;
	            }
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
                    n = n1 = n.getParent();
                } else {
                    n = n2 = n.getParent();  
                }
	            if( ! nodesTraversed[n.getNr()] ) {
	                nodesTraversed[n.getNr()] = true;
	                nseen += 1;
	            } 
	        }
        }
        return n1;
    }

	//private void setUpCladeSet(Node node) {
	//	cladeSet.add(node.getNr());
	//	for (Node child : node.getChildren()) {
	//		setUpCladeSet(child);
	//	}
	//}

	@Override
	public double calculateLogP() {
		if (!initialised || region == null) {
			logP = 0;
			return logP;
			//initialise();
		}

		logP = Double.NEGATIVE_INFINITY;
		boolean isInside = isInsideInput.get();

		double[] location = new double[2];
		if (allInternalNodes) {
			logP = 0;
			for (int i = 0; i < tree.getNodeCount(); i++) {
				this.location.getMatrixValues1(i, location);
				if (region.isInside(location[0], location[1])) {
					if (!isInside) {
						logP -= 1e-20;
					}
				} else {
					if (isInside) {
						logP -= 1e-20;
					}
				}
			}
			return logP;
		}
		
		
		if (isRoot) {
			taxonNr = tree.getRoot().getNr();
		} else {
			if (!isTip) {
				//calcMRCAtime(tree.getRoot(), new int[1]);
                nodesTraversed = new boolean[tree.getNodeCount()];
                nseen = 0;
                Node m = getCommonAncestor();
            	taxonNr = m.getNr();
                isMonophyletic = nseen == 2 * taxonIndex.length - 1;
			}
		}
		this.location.getMatrixValues1(taxonNr, location);
		if (region.isInside(location[0], location[1])) {
			if (isInside) {
				logP = 0;
			}
		} else {
			if (!isInside) {
				logP = 0;
			}
		}
		return logP;
	}
	
	
    /**
     * Recursively visit all leaf nodes, and collect number of taxa in the taxon
     * set. When all taxa in the set are visited, record the time.
     * *
     * @param node
     * @param nTaxonCount
     */
//    int calcMRCAtime(final Node node, final int[] nTaxonCount) {
//        if (node.isLeaf()) {
//            nTaxonCount[0]++;
//            if (isInTaxaSet[node.getNr()]) {
//                return 1;
//            } else {
//                return 0;
//            }
//        } else {
//            int iTaxons = calcMRCAtime(node.getLeft(), nTaxonCount);
//            int nTaxa = nTaxonCount[0];
//            nTaxonCount[0] = 0;
//            for (int i = 1; i < node.getChildCount(); i++) {
//            	Node child = node.getChild(i);
//                iTaxons += calcMRCAtime(child, nTaxonCount);
//                nTaxa += nTaxonCount[0];
//            }
//            nTaxonCount[0] = nTaxa;
//            if (iTaxons == nrOfTaxa) {
//            	taxonNr = node.getNr();
//            	isMonophyletic = (nTaxonCount[0] == nrOfTaxa);
//                return iTaxons + 1;                	
//           }
//            return iTaxons;
//        }
//    }

    
    
    @Override
    public void init(PrintStream out) {
    	if (allInternalNodes) {
    		out.print(getID() + "\t");
    	} else {
    		out.print(getID() + ".latitude\t");
    		out.print(getID() + ".longitude\t");
    	}
    }
    
    @Override
    public void log(long nSample, PrintStream out) {
		if (!initialised) {
			initialise();
		}
    	if (allInternalNodes) {
    		out.print(logP + "\t");
    		return;
    	}
		double[] location = new double[2];
		this.location.getMatrixValues1(taxonNr, location);
        out.print(location[0] + "\t");
        out.print(location[1] + "\t");
    }
    

    @Override
    public void store() {
    	storedIsMonophyletic = isMonophyletic;
    	if (taxonIndex != null) {
    		System.arraycopy(taxonIndex, 0, storedTaxonIndex, 0, taxonIndex.length);
    		System.arraycopy(nodesTraversed, 0, storedNodesTraversed, 0, nodesTraversed.length);
    		//storedCladeSet.clear();
            //storedCladeSet.addAll(cladeSet);
    	}
    	
    	
    	storedTaxonNr = taxonNr;
    	super.store();
    }
    
    @Override
    public void restore() {
    	isMonophyletic = storedIsMonophyletic;
    	
    	int [] tmp = taxonIndex;
    	taxonIndex = storedTaxonIndex;
    	storedTaxonIndex = tmp;
    	
    	//Set<Integer> tmp3 = cladeSet;
    	//cladeSet = storedCladeSet;
    	//storedCladeSet = tmp3;
    	
    	boolean [] tmp2 = nodesTraversed;
    	nodesTraversed = storedNodesTraversed;
    	storedNodesTraversed = tmp2;
    	
    	taxonNr = storedTaxonNr;
    	super.restore();
    }
    
	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}

	public double[] sample() {
		if (region != null) {
			return region.sample(isInsideInput.get());
		}
		// no region available to sample from, so uniformly sample from the whole planet.
		double [] location = new double[2];
		location[0] = -90 + Randomizer.nextDouble() * 180;
		location[1] = -180 + Randomizer.nextDouble() * 360;
		return location;
	}

	public int getStoredTaxonNr() {
		return storedTaxonNr;
	}

	public boolean isMonoPhyletic() {
		return isMonophyletic;
	}

}
