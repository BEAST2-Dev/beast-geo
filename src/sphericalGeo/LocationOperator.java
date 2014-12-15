package sphericalGeo;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

@Description("operator to samples locations from geo-priors")
public class LocationOperator extends Operator {
	public Input<RealParameter> locationInput = new Input<RealParameter>("location",
			"2 dimensional parameter representing locations (in latitude, longitude) of nodes in a tree", Validate.REQUIRED);
	public Input<ApproxMultivariateTraitLikelihood> likelihoodInput = new Input<ApproxMultivariateTraitLikelihood>("likelihood", 
			"likelihood over the locations", Validate.REQUIRED);
	
	RealParameter sampledLocations;
	boolean [] isSampled;
	List<Integer> sampleNumber;
	List<GeoPrior> geopriors;
	
	@Override
	public void initAndValidate() throws Exception {
		sampledLocations = locationInput.get();

		ApproxMultivariateTraitLikelihood likelihood = likelihoodInput.get();
		isSampled = likelihood.isSampled;
		sampleNumber = likelihood.sampleNumber;
		geopriors = likelihood.geopriorsInput.get();
		
	}

	@Override
	public double proposal() {
		if (sampleNumber.size() == 0) {
			return 0;
		}
		int i = Randomizer.nextInt(sampleNumber.size());
		GeoPrior prior = geopriors.get(i);
		double [] location = prior.sample();
		i = sampleNumber.get(i);
		sampledLocations.setValue(i * 2, location[0]);
		sampledLocations.setValue(i * 2 + 1, location[1]);
		
		return 0;
	}

}
