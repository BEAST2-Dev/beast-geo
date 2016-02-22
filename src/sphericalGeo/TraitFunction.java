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
	public Input<GenericTreeLikelihood> likelihoodInput = new Input<>("likelihood", "trait likelihood to be logged", Validate.REQUIRED);

	GenericTreeLikelihood likelihood;
	Tree tree;

	Method getPosition;
	
	@Override
	public void initAndValidate() {
		if (likelihoodInput.get() instanceof ApproxMultivariateTraitLikelihood) {
			likelihood = (ApproxMultivariateTraitLikelihood)likelihoodInput.get();
		} else if (likelihoodInput.get() instanceof ApproxMultivariateTraitLikelihoodF) {
			likelihood = (ApproxMultivariateTraitLikelihoodF)likelihoodInput.get();
		} else if (likelihoodInput.get() instanceof PFApproxMultivariateTraitLikelihood) {
			likelihood= (PFApproxMultivariateTraitLikelihood)likelihoodInput.get();
		} else {
			throw new RuntimeException("likelihood should be one of ApproxMultivariateTraitLikelihood or PFApproxMultivariateTraitLikelihood");
		}
        tree = ((Tree) (likelihood.treeInput.get()));
		try {
			getPosition = likelihood.getClass().getMethod("getPosition", int.class);
		} catch (NoSuchMethodException | SecurityException e) {
			throw new IllegalArgumentException(e);
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
		try {
			return ((double[]) getPosition.invoke(likelihood, i))[j];
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 0.0;
	}

}
