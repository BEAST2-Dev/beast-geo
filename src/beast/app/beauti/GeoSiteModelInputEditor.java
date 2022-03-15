package beast.app.beauti;

import javax.swing.JComponent;

import beast.app.draw.BEASTObjectInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import sphericalGeo.sitemodel.GeoSiteModel;

public class GeoSiteModelInputEditor extends BEASTObjectInputEditor {
	
	private static final long serialVersionUID = 1L;

	@Override
	public Class<?> type() {
		return GeoSiteModel.class;
	}
	
	public GeoSiteModelInputEditor(BeautiDoc doc) {
		super(doc);
	}
	
	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
			boolean addButtons) {
		GeoSiteModel siteModel = (GeoSiteModel)input.get(); 
		input = siteModel.substModelInput;
		beastObject = (BEASTInterface) input.get();
		super.init(input, beastObject, itemNr, ExpandOption.TRUE, addButtons);
	}

	@Override
	protected void addComboBox(JComponent box, Input<?> input, BEASTInterface beastObject0) {
		// suppress addition of combobox
	}
}
