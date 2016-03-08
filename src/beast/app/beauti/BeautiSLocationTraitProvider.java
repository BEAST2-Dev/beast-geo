package beast.app.beauti;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.border.EmptyBorder;

import sphericalGeo.ApproxMultivariateTraitLikelihood;
import beast.app.beauti.BeautiAlignmentProvider;
import beast.app.beauti.BeautiDoc;
import beast.app.beauti.PartitionContext;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.State;
import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;


/** trait provided for spherical diffusion model **/
@Description("Location provider for BEAUti template to set up spherical diffusion models")
public class BeautiSLocationTraitProvider extends BeautiAlignmentProvider {

	@Override
	protected List<BEASTInterface> getAlignments(BeautiDoc doc) {
		try {
            List<String> trees = new ArrayList<>();
            doc.scrubAll(true, false);
            State state = (State) doc.pluginmap.get("state");
            for (StateNode node : state.stateNodeInput.get()) {
                if (node instanceof Tree) { // && ((Tree) node).m_initial.get() != null) {
                    trees.add(BeautiDoc.parsePartition(((Tree) node).getID()));
                }
            }
            TraitDialog2 dlg = new TraitDialog2(doc, trees);
            if (dlg.showDialog("Create new location")) {
            	String tree = dlg.tree;
            	String name = dlg.name;
            	PartitionContext context = new PartitionContext(name, name, name, tree);

            	Alignment alignment = (Alignment) doc.addAlignmentWithSubnet(context, template.get());
            	List<BEASTInterface> list = new ArrayList<>();
            	list.add(alignment);
            	editAlignment(alignment, doc);
            	return list;
            }
		} catch (Exception e) {
			e.printStackTrace();
		}
        return null;
	}
	
	@Override
	protected int matches(Alignment alignment) {
		for (Object o : alignment.getOutputs()) {
			BEASTInterface  output = (BEASTInterface) o;
			if (output instanceof sphericalGeo.ApproxMultivariateTraitLikelihood) {
				return 10;
			}
		}
		return 0;
	}
	
	
	@Override
	void editAlignment(Alignment alignment, BeautiDoc doc) {
		SLocationInputEditor editor = new SLocationInputEditor(doc);
		ApproxMultivariateTraitLikelihood likelihood = null;
		for (Object o : alignment.getOutputs()) {
			BEASTInterface  output = (BEASTInterface) o;
			if (output instanceof ApproxMultivariateTraitLikelihood) {
				likelihood = (ApproxMultivariateTraitLikelihood) output;
				editor.initPanel(likelihood);
		        JOptionPane optionPane = new JOptionPane(editor, JOptionPane.PLAIN_MESSAGE,
		                JOptionPane.CLOSED_OPTION, null, new String[]{"Close"}, "Close");
		        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

		        final JDialog dialog = optionPane.createDialog(Frame.getFrames()[0], "Location trait editor");
		    	dialog.setName("LocationTraitEditor");
		        // dialog.setResizable(true);
		        dialog.pack();

		        dialog.setVisible(true);
		        editor.convertTableDataToTrait();
		        try {
			        // TODO: any post-processing...
			        // AlignmentFromTraitMap traitData = (AlignmentFromTraitMap) likelihood.m_data.get();
			        
		        } catch (Exception e) {
					e.printStackTrace();
				}

				return;
			}
		}
	}

}
