package sphericalGeo;




import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Tree;

@Description("Helper class for logging locations from ApproxMultivariateTraitLikelihood")
public class TraitFunction extends RealParameter {
	public Input<GenericTreeLikelihood> likelihoodInput = new Input<GenericTreeLikelihood>("likelihood", "trait likelihood to be logged", Validate.REQUIRED);

	ApproxMultivariateTraitLikelihood likelihood;
	PFApproxMultivariateTraitLikelihood likelihood2;
	Tree tree;
	
	@Override
	public void initAndValidate() throws Exception {
		if (likelihoodInput.get() instanceof ApproxMultivariateTraitLikelihood) {
			likelihood = (ApproxMultivariateTraitLikelihood)likelihoodInput.get();
			tree = ((Tree) (likelihood.treeInput.get()));
		} else if (likelihoodInput.get() instanceof PFApproxMultivariateTraitLikelihood) {
			likelihood2 = (PFApproxMultivariateTraitLikelihood)likelihoodInput.get();
	        tree = ((Tree) (likelihood2.treeInput.get()));
		} else {
			throw new RuntimeException("likelihood should be one of ApproxMultivariateTraitLikelihood or PFApproxMultivariateTraitLikelihood");
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
		if (likelihood != null)
		return likelihood.getPostion(i)[j];
		return likelihood2.getPostion(i)[j];
	}

}
