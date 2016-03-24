package beast.app.beauti;

import beast.app.draw.InputEditor.Base;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.util.TreeParser;

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
