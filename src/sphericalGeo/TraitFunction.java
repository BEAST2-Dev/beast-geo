package sphericalGeo;




import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Tree;

@Description("Helper class for logging locations from ApproxMultivariateTraitLikelihood")
public class TraitFunction extends RealParameter {
	public Input<LocationProvider> likelihoodInput = new Input<>("likelihood", "trait likelihood to be logged", Validate.REQUIRED);

	LocationProvider likelihood;
	Tree tree;

	@Override
	public void initAndValidate() {
		if (likelihoodInput.get() instanceof GenericTreeLikelihood) {
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
