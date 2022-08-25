package sphericalGeo;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

@Description("operator to samples locations from geo-priors")
public class LocationOperatorF extends Operator {
	public Input<RealParameter> locationInput = new Input<>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree", Validate.REQUIRED);
	public Input<ApproxMultivariateTraitLikelihoodF> likelihoodInput = new Input<>("likelihood", 
			"likelihood over the locations", Validate.REQUIRED);
	
	RealParameter sampledLocations;
	//List<Integer> sampleNumber;
	int [] taxonNrs;
	ApproxMultivariateTraitLikelihoodF likelihood;
	
	@Override
	public void initAndValidate() {
		sampledLocations = locationInput.get();

		likelihood = likelihoodInput.get();
		
	}

	@Override
	public double proposal() {
		taxonNrs = likelihood.taxonNrs;
		if (taxonNrs == null || taxonNrs.length == 0) {
			return 0;
		}
		MultiGeoPrior multiGeoPrior = likelihood.multiGeopriorsInput.get();
		List<GeoPrior> geopriors = multiGeoPrior.geopriorsInput.get();
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
		if (!prior.isMonophyletic) {
			prior.initialise();
			if (!prior.isMonophyletic) {
				// fail -- we should never move a non-monophyletic clade
				return Double.NEGATIVE_INFINITY;
			}
		}
		double [] location = prior.sample();
		int k = multiGeoPrior.getCladeTopNodeNr(i);
		//System.out.println(prior.getID() + ":" + k + " ("+sampledLocations.getValue(k*2)+","+sampledLocations.getValue(k*2+1) + ") => (" + location[0] + "," + location[1]+")");
		sampledLocations.setValue(k * 2, location[0]);
		sampledLocations.setValue(k * 2 + 1, location[1]);
		
		return 0;
	}

}
