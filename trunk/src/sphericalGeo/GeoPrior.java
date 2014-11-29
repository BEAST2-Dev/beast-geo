package sphericalGeo;

import java.util.List;
import java.util.Random;

import sphericalGeo.region.Region;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Flat prior over a region")
public class GeoPrior extends Distribution {
	public Input<Region> regionInput = new Input<Region>("region", "region to be in (or not depending on 'isinside' flag)", Validate.REQUIRED);
	public Input<Boolean> isInsideInput = new Input<Boolean>("isInside", "whether the prior is for being inside the region, instead of outside", true);
	public Input<RealParameter> locationInput = new Input<RealParameter>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree", Validate.REQUIRED);
	public Input<Tree> treeInput = new Input<Tree>("tree", "beast tree (from which to get the taxon set)");

	public Input<Taxon> taxonInput = new Input<Taxon>("taxon", "taxon associated with this region, if only a tip is restricted. Otherwise use 'taxonset'");
	public Input<TaxonSet> taxonSetInput = new Input<TaxonSet>("taxonset",
			"specify the prior over an internal node that is the MRCA of set of taxa. Select all taxa for the root", Validate.XOR, taxonInput);

	Region region;
	RealParameter location;
	Tree tree;
	TaxonSet taxonSet;
	int taxonNr;
	boolean isRoot;
	boolean isTip = false;
	
    // number of taxa in taxon set
    int nrOfTaxa = -1;
    // array of flags to indicate which taxa are in the set
    boolean[] isInTaxaSet;


	@Override
	public void initAndValidate() throws Exception {
		region = regionInput.get();
		location = locationInput.get();
		if (location.getMinorDimension1() != 2) {
			throw new RuntimeException("expected that location parameter to have minor dimension 2");
		}
		tree = treeInput.get();
		taxonSet = tree.getTaxonset();
		if (location.getDimension() != taxonSet.getTaxonCount() * 4 - 2) {
			throw new RuntimeException("expected that location parameter to have dimension 2 time number of taxa - 1 = " + (taxonSet.getTaxonCount() * 4 - 2)
					+ " not " + location.getDimension());
		}

		if (taxonInput.get() != null) {
			isTip = true;
			String taxonName = taxonInput.get().getID();
			List<String> names = taxonSet.asStringList();
			taxonNr = names.indexOf(taxonName);
			if (taxonNr < 0) {
				throw new RuntimeException("Could not find taxon " + taxonName + ". Typo perhaps?");
			}
		} else {
			TaxonSet taxonset2 = taxonSetInput.get();
			if (taxonset2.getTaxonCount() == taxonSet.getTaxonCount()) {
				isRoot = true;
				taxonNr = tree.getRoot().getNr();
			} else {
				isInTaxaSet = new boolean[taxonSet.getTaxonCount()];
				List<String> names = taxonSet.asStringList();
				//int k = 0;
	            for (final String sTaxon : taxonset2.asStringList()) {
	                final int iTaxon = names.indexOf(sTaxon);
	                if (iTaxon < 0) {
	                    throw new RuntimeException("Cannot find taxon " + sTaxon + " in data");
	                }
	                if (isInTaxaSet[iTaxon]) {
	                    throw new RuntimeException("Taxon " + sTaxon + " is defined multiple times, while they should be unique");
	                }
	                isInTaxaSet[iTaxon] = true;
	                //taxonIndex[k++] = iTaxon;
	            }
			}
		}

		super.initAndValidate();
	}

	@Override
	public double calculateLogP() throws Exception {
		logP = Double.NEGATIVE_INFINITY;
		boolean isInside = isInsideInput.get();

		double[] location = new double[2];
		if (isRoot) {
			taxonNr = tree.getRoot().getNr();
		} else {
			if (!isTip) {
				calcMRCAtime(tree.getRoot(), new int[1]);
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
    int calcMRCAtime(final Node node, final int[] nTaxonCount) {
        if (node.isLeaf()) {
            nTaxonCount[0]++;
            if (isInTaxaSet[node.getNr()]) {
                return 1;
            } else {
                return 0;
            }
        } else {
            int iTaxons = calcMRCAtime(node.getLeft(), nTaxonCount);
            final int nLeftTaxa = nTaxonCount[0];
            nTaxonCount[0] = 0;
            if (node.getRight() != null) {
                iTaxons += calcMRCAtime(node.getRight(), nTaxonCount);
                final int nRightTaxa = nTaxonCount[0];
                nTaxonCount[0] = nLeftTaxa + nRightTaxa;
                if (iTaxons == nrOfTaxa) {
                	taxonNr = node.getNr();
                    return iTaxons + 1;                	
               }
            }
            return iTaxons;
        }
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

}
