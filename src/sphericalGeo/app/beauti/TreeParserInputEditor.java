package sphericalGeo.app.beauti;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor.Base;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.tree.TreeParser;

public class TreeParserInputEditor extends Base {

	private static final long serialVersionUID = 1L;

	public TreeParserInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return TreeParser.class;
	}
	
	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, itemNr, isExpandOption, addButtons);
	}

}
