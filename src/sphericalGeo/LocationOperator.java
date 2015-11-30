package sphericalGeo;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

@Description("operator to samples locations from geo-priors")
public class LocationOperator extends Operator {
	public Input<RealParameter> locationInput = new Input<RealParameter>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree", Validate.REQUIRED);
	public Input<ApproxMultivariateTraitLikelihood> likelihoodInput = new Input<ApproxMultivariateTraitLikelihood>("likelihood", 
			"likelihood over the locations", Validate.REQUIRED);
	
	RealParameter sampledLocations;
	List<Integer> sampleNumber;
	ApproxMultivariateTraitLikelihood likelihood;
	
	@Override
	public void initAndValidate() throws Exception {
		sampledLocations = locationInput.get();

		likelihood = likelihoodInput.get();
		
	}

	@Override
	public double proposal() {
		sampleNumber = likelihood.sampleNumber;
		if (sampleNumber.size() == 0) {
			return 0;
		}
		List<GeoPrior> geopriors = likelihood.geopriorsInput.get();
		int i = Randomizer.nextInt(geopriors.size());
		GeoPrior prior = geopriors.get(i);
		if (prior.allInternalNodes) {
			TreeInterface tree = likelihood.tree;
			int k = tree.getLeafNodeCount() + Randomizer.nextInt(tree.getInternalNodeCount());
			double [] location = prior.sample();
			sampledLocations.setValue(k * 2, location[0]);
			sampledLocations.setValue(k * 2 + 1, location[1]);
			return 0;
		}
		double [] location = prior.sample();
		int k = prior.getTaxonNr();
		sampledLocations.setValue(k * 2, location[0]);
		sampledLocations.setValue(k * 2 + 1, location[1]);
		
		return 0;
	}

}
