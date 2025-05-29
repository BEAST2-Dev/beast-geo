package sphericalGeo;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.RealParameter;
import sphericalGeo.region.Region;

@Description("Initialise location so it conforms to geoprior")
public class LocationInitialiser extends BEASTObject implements StateNodeInitialiser {
	public Input<RealParameter> locationInput = new Input<>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree", Validate.REQUIRED);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
	public Input<List<GeoPrior>> geopriorsInput = new Input<>("geoprior", "geographical priors on tips, root or clades restricting these nodes to a region", new ArrayList<>(), Validate.REQUIRED);
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void initStateNodes() {
		RealParameter locationParam = locationInput.get();
		TreeInterface tree = treeInput.get();
		TaxonSet taxonSet = tree.getTaxonset();
		if (locationParam.getDimension() != taxonSet.getTaxonCount() * 4 - 2) {
			Log.warning.println("Setting dimension of location parameter to have dimension 2 time number of taxa - 1 = " + (taxonSet.getTaxonCount() * 4 - 2)
					+ " (from " + locationParam.getDimension() +")");
			locationParam.setDimension(taxonSet.getTaxonCount() * 4 - 2);
		}
		
		for (GeoPrior prior : geopriorsInput.get()) {
			if (prior.allInternalNodesInput.get()) {
				Region region = prior.region;
				double[] location = new double[2];
				for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
					locationParam.getMatrixValues1(i, location);
					if (!region.isInside(location[0], location[1])) {
						double [] newLocation = region.sample(true);
						locationParam.setMatrixValue(i, 0, newLocation[0]);
						locationParam.setMatrixValue(i, 1, newLocation[1]);
					}
				}
				prior.initialised = true;
			}
		}
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
	}

}
