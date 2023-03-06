package sphericalGeo.util;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.tree.Node;

@Description("Branch rate model that reads branch rates from the tree -- to be used in combination with PostHocLoationSampler")
public class PostHocBranchRateModel extends StrictClockModel {
	final public Input<String> tagInput = new Input<>("tag","meta data tag used in node for the branch rate","rate");

	String tag;
	
	@Override
	public void initAndValidate() {
		tag = tagInput.get();
		tag = tag + "=";
	}

	@Override
	public double getRateForBranch(Node node) {
		String str = node.metaDataString;
		int i = str.indexOf(tag);
		if (i < 0) {
			final double meanRate = meanRateInput.get().getArrayValue();
			return meanRate;
		}
		
		str = str.substring(i + tag.length());
		double rate = Double.parseDouble(str);
		return rate;
	}

}
