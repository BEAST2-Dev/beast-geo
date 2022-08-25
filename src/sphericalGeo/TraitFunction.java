package sphericalGeo;




import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.Tree;

@Description("Helper class for logging locations from ApproxMultivariateTraitLikelihood")
public class TraitFunction extends RealParameter {
	public Input<LocationProvider> likelihoodInput = new Input<>("likelihood", "trait likelihood to be logged", Validate.REQUIRED);

	LocationProvider likelihood;
	Tree tree;

	@Override
	public void initAndValidate() {
		likelihood = likelihoodInput.get();
		if (likelihood instanceof GenericTreeLikelihood) {
	        tree = ((Tree) ((GenericTreeLikelihood) likelihood).treeInput.get());
		} else {
			throw new RuntimeException("likelihood should have derived from GenericTreeLikelihood");
		}
	}

	@Override
	public int getDimension() {
		return tree.getNodeCount() * 2;
	}

	@Override
	public int getMinorDimension1() {
		return 2;
	}

	@Override
	public double getArrayValue() {
		return 0;
	}

	@Override
	public Double getMatrixValue(int i, int j) {
		return likelihood.getPosition(i)[j];
	}

}
