package beast.app.beauti;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.border.EmptyBorder;

import sphericalGeo.AlignmentFromTraitMap;
import sphericalGeo.ApproxMultivariateTraitLikelihood;
import sphericalGeo.TraitFunction;
import sphericalGeo.TreeTraitMap;
import beast.app.beauti.BeautiAlignmentProvider;
import beast.app.beauti.BeautiDoc;
import beast.app.beauti.PartitionContext;
import beast.app.util.Utils;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeWithMetaDataLogger;
import beast.util.NexusParser;
import beast.util.TreeParser;


/** trait provided for spherical diffusion model **/
@Description("Location provider for BEAUti template to set up spherical diffusion models")
public class BeautiSLocationTraitProvider extends BeautiAlignmentProvider {
	final static String FIXED_TREE = "[[Fixed tree]]";

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc) {
		try {
            List<String> trees = new ArrayList<>();
            doc.scrubAll(true, false);
            State state = (State) doc.pluginmap.get("state");
            for (StateNode node : state.stateNodeInput.get()) {
                if (node instanceof Tree) { // && ((Tree) node).m_initial.get() != null) {
                    trees.add(BeautiDoc.parsePartition(((Tree) node).getID()));
                }
            }
            trees.add(FIXED_TREE);
            TraitDialog2 dlg = new TraitDialog2(doc, trees);
            if (dlg.showDialog("Create new location")) {
            	String tree = dlg.tree;
            	if (tree.equals(FIXED_TREE)) {
            		return getTree(doc);
            	}
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
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files, String[] args) {
		String name = args[0];
		String tree = args[1];
    	PartitionContext context = new PartitionContext(name, name, name, tree);

    	Alignment alignment = (Alignment) doc.addAlignmentWithSubnet(context, template.get());
    	List<BEASTInterface> list = new ArrayList<>();
    	list.add(alignment);

		AlignmentFromTraitMap traitData = (AlignmentFromTraitMap) alignment;
		try {
	        BufferedReader fin = new BufferedReader(new FileReader(files[0]));
	        StringBuffer buf = new StringBuffer();
	        // do not eat up header -- it might contain a useful entry, 
	        // but if not, it will not hurt
	        // fin.readLine();
	        // process data
	        while (fin.ready()) {
	            String str = fin.readLine();
	            str = str.replaceFirst("\t", "=") + ",";
	            // only add entries that are non-empty
	            if (str.indexOf("=") > 0 && !str.matches("^\\s+=.*$")) {
	                buf.append(str);
	            }
	        }
	        fin.close();
	        TreeTraitMap trait = traitData.traitInput.get();
	        trait.value.setValue(buf.toString().trim(), trait);
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}

    	return list;
	}
	
	protected List<BEASTInterface> getTree(BeautiDoc doc) {
		try {
            File file = Utils.getLoadFile("Open tree file with fixed tree", new File(Beauti.g_sDir), "Tree file", "tree","tre","txt","nxs");
            if (file != null) {
            	NexusParser parser = new NexusParser();
            	parser.parseFile(file);
            	if (parser.trees == null || parser.trees.size() == 0) {
            		JOptionPane.showMessageDialog(null, "Did not find any tree in the file -- giving up.");
            		return null;
            	}
            	if (parser.trees.size() > 1) {
            		JOptionPane.showMessageDialog(null, "Found more than one tree in the file -- expected only 1!");
            		return null;
            	}
            	Tree tree = parser.trees.get(0);
            	
            	// create dummy alignment
            	List<Taxon> seqs = new ArrayList<>();
            	for (String name : tree.getTaxaNames()) {
            		seqs.add(new Taxon(name));
            	}
            	TaxonSet taxa = new TaxonSet();
            	taxa.initByName("taxon", seqs);

            	TreeParser treeParser = new TreeParser();
            	treeParser.initByName("newick", tree.getRoot().toNewick(),
            			"IsLabelledNewick", true,
            			"taxonset", taxa,
            			"adjustTipHeights", false,
            			"estimate", false
            			);
            	String id = file.getName();
            	if (id.lastIndexOf('.') > -1) {
            		id = id.substring(0, id.lastIndexOf('.'));
            	}

            	treeParser.setID("Tree.t:" + id);
            	doc.registerPlugin(treeParser);
            	
            	PartitionContext context = new PartitionContext(id, id, id, id);

            	Alignment alignment = (Alignment) doc.addAlignmentWithSubnet(context, template.get());
            	List<BEASTInterface> list = new ArrayList<>();
            	list.add(alignment);
            	editAlignment(alignment, doc);
            	
            	// add tree logger
            	GenericTreeLikelihood likelihood = null;
            	for (BEASTInterface o : alignment.getOutputs()) {
            		if (o instanceof GenericTreeLikelihood) {
            			likelihood = (GenericTreeLikelihood) o;
            		}
            	}
            	
            	TraitFunction locationTrait = new TraitFunction();
            	locationTrait.initByName("likelihood", likelihood, "value", "0.0");
            	locationTrait.setID("location");
            	TreeWithMetaDataLogger treelogger = new TreeWithMetaDataLogger();
            	treelogger.initByName("tree", treeParser, "metadata", locationTrait);
            	treelogger.setID("TreeWithMetaDataLogger");

            	Logger logger = new Logger();
            	logger.initByName("log", treelogger, "logEvery", 1000, "mode", "tree", "fileName", id + ".trees");
            	logger.setID("TreeLogger");
            	
            	MCMC mcmc = (MCMC) doc.pluginmap.get("mcmc");
            	mcmc.loggersInput.setValue(logger, mcmc);
            	return list;
            }
        } catch (Exception e) {
        	e.printStackTrace();
        	JOptionPane.showMessageDialog(null, "Something went wrong: " + e.getMessage());
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
