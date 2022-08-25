package sphericalGeo.app.beauti;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.beauti.PriorProvider;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BEASTObjectPanel;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.inputeditor.SmallButton;
import beastfx.app.inputeditor.StringInputEditor;
import beastfx.app.inputeditor.TaxonSetDialog;
import beastfx.app.util.FXUtils;
import beastfx.app.util.Utils;
import javafx.scene.Node;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.CompoundDistribution;
import beast.base.core.Log;
import beast.base.core.ProgramStatus;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.MRCAPrior;
import sphericalGeo.ApproxMultivariateTraitLikelihood;
import sphericalGeo.GeoPrior;
import sphericalGeo.LocationOperator;
import sphericalGeo.region.KMLRegion;;

public class GeoPriorProvider extends BEASTObjectInputEditor implements PriorProvider {
	private static final long serialVersionUID = 1L;

	public GeoPriorProvider() {
		super(null);
	}
	
	public GeoPriorProvider(BeautiDoc doc) {
		super(doc);
	}
	
	@Override
	public Class<?> type() {
		return GeoPrior.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
			boolean addButtons) {
        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr= itemNr;
		
        HBox itemBox = FXUtils.newHBox();

        GeoPrior prior = (GeoPrior) beastObject;
        String text = prior.getID();

        Button taxonButton = new Button(text);
//        taxonButton.setMinimumSize(Base.PREFERRED_SIZE);
//        taxonButton.setPreferredSize(Base.PREFERRED_SIZE);
        itemBox.getChildren().add(taxonButton);
        taxonButton.setOnAction(e -> {
                List<?> list = (List<?>) m_input.get();
                GeoPrior prior2 = (GeoPrior) list.get(itemNr);
                try {
                    TaxonSet taxonset = prior2.taxonSetInput.get();
                    List<Taxon> originalTaxa = new ArrayList<>();
                    originalTaxa.addAll(taxonset.taxonsetInput.get());
                    Set<Taxon> candidates = getTaxonCandidates(prior2);
                    TaxonSetDialog dlg = new TaxonSetDialog(taxonset, candidates, doc);
                    if (dlg.showDialog()) {
        	            if (dlg.taxonSet.taxonsetInput.get().size() == 0) {
        	            	JOptionPane.showMessageDialog(doc.beauti, "At least one taxon should be included in the taxon set",
        	            			"Error specifying taxon set", JOptionPane.ERROR_MESSAGE);
        	            	taxonset.taxonsetInput.get().addAll(originalTaxa);
        	            	return;
        	            }

                        prior2.taxonSetInput.setValue(dlg.taxonSet, prior2);
                        int i = 1;
                        String id = dlg.taxonSet.getID();
                        while (doc.pluginmap.containsKey(dlg.taxonSet.getID()) && doc.pluginmap.get(dlg.taxonSet.getID()) != dlg.taxonSet) {
                        	dlg.taxonSet.setID(id + i);
                        	i++;
                        }
                        BEASTObjectPanel.addPluginToMap(dlg.taxonSet, doc);
                        prior2.setID(dlg.taxonSet.getID() + ".prior");

                    }
                } catch (Exception e1) {
                    e1.printStackTrace();
                }
                refreshPanel();
            });


        itemBox.getChildren().add((Node) createRegionEditor());

        CheckBox isInsidedBox = new CheckBox(doc.beautiConfig.getInputLabel(prior, prior.isInsideInput.getName()));
        isInsidedBox.setId(text+".isInside");
        isInsidedBox.setSelected(prior.isInsideInput.get());
        isInsidedBox.setTooltip(new Tooltip(prior.isInsideInput.getHTMLTipText()));
        isInsidedBox.setOnAction(e->new GeoPriorActionListener(prior));
        itemBox.getChildren().add(isInsidedBox);
        
        CheckBox isAllInternaldBox = new CheckBox(doc.beautiConfig.getInputLabel(prior, prior.allInternalNodesInput.getName()));
        isAllInternaldBox.setId(text+".allInternalNodes");
        isAllInternaldBox.setSelected(prior.allInternalNodesInput.get());
        isAllInternaldBox.setTooltip(new Tooltip(prior.allInternalNodesInput.getHTMLTipText()));
        isAllInternaldBox.setOnAction(e->new GeoPriorActionListener2(prior));
        itemBox.getChildren().add(isAllInternaldBox);
        

        
        Button deleteButton = new SmallButton("-", true);
        deleteButton.setTooltip(new Tooltip("Delete this geo-prior"));
        deleteButton.setOnAction(e -> {
				Log.warning.println("Trying to delete a geoprior");
				List<?> list = (List<?>) m_input.get();
				GeoPrior prior0 = (GeoPrior) list.get(itemNr);
				doc.disconnect(prior0, "prior", "distribution");
				doc.disconnect(prior0, "tracelog", "log");
				doc.unregisterPlugin(prior0);
				// remove operator and location parameter from state, iff there are no other geo-priors left
				boolean hasMoreGeoPriors = false;
				CompoundDistribution priors = (CompoundDistribution) doc.pluginmap.get("prior");
				for (BEASTInterface o : priors.pDistributions.get()) {
					if (o instanceof GeoPrior) {
						hasMoreGeoPriors = true;
						break;
					}
				}
				if (!hasMoreGeoPriors) {
					BEASTInterface operator = doc.pluginmap.get("location.sampler");
					doc.disconnect(operator, "mcmc", "operator");
					doc.disconnect(prior0.locationInput.get(), "state", "stateNode");
				}
				refreshPanel();
			}        	
        );
        //itemBox.getChildren().add(Box.createGlue());
        itemBox.getChildren().add(deleteButton);

        getChildren().add(itemBox);
	}
	
	
    /**
     * class to deal with toggling isInside flag on an GeoPrior *
     */
    class GeoPriorActionListener implements ActionListener {
        GeoPrior m_prior;

        GeoPriorActionListener(GeoPrior prior) {
            m_prior = prior;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            try {
                m_prior.isInsideInput.setValue(((JCheckBox) e.getSource()).isSelected(), m_prior);
                refreshPanel();
            } catch (Exception ex) {
            	Log.warning.println("PriorListInputEditor " + ex.getMessage());
            }
        }
    }

	
    /**
     * class to deal with toggling isInside flag on an GeoPrior *
     */
    class GeoPriorActionListener2 implements ActionListener {
        GeoPrior m_prior;

        GeoPriorActionListener2(GeoPrior prior) {
            m_prior = prior;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            try {
                m_prior.allInternalNodesInput.setValue(((JCheckBox) e.getSource()).isSelected(), m_prior);
                refreshPanel();
            } catch (Exception ex) {
            	Log.warning.println("PriorListInputEditor " + ex.getMessage());
            }
        }
    }

    Set<Taxon> getTaxonCandidates(GeoPrior prior) {
        Set<Taxon> candidates = new HashSet<>();
        Tree tree = prior.treeInput.get();
        String [] taxa = null;
        if (tree.m_taxonset.get() != null) {
        	try {
            	TaxonSet set = tree.m_taxonset.get();
        		set.initAndValidate();
            	taxa = set.asStringList().toArray(new String[0]);
        	} catch (Exception e) {
            	taxa = prior.treeInput.get().getTaxaNames();
			}
        } else {
        	taxa = prior.treeInput.get().getTaxaNames();
        }
        
        for (String taxon : taxa) {
            candidates.add(doc.getTaxon(taxon));
        }
        return candidates;
	}
	
	@Override
	public List<Distribution> createDistribution(BeautiDoc doc) {
		this.doc = doc;
    	GeoPrior prior = new GeoPrior();
        try {

            List<Tree> trees = new ArrayList<>();
            getDoc().scrubAll(true, false);
            
            for (BEASTInterface partition : doc.getPartitions("SiteModel")) {
            	if (partition instanceof ApproxMultivariateTraitLikelihood) {
            		ApproxMultivariateTraitLikelihood likelihood = (ApproxMultivariateTraitLikelihood) partition;
            		trees.add((Tree) likelihood.treeInput.get());
            	}
            }
                        
            if (trees.size() == 0) {
            	JOptionPane.showMessageDialog(this, "Could not find a spherical geography partition\nNo prior will be added.");
            	return null;
            }

            int treeIndex = 0;
            if (trees.size() > 1) {
                String[] treeIDs = new String[trees.size()];
                for (int j = 0; j < treeIDs.length; j++) {
                    treeIDs[j] = trees.get(j).getID();
                }
                treeIndex = JOptionPane.showOptionDialog(null, "Select a tree", "GeoPrior selector",
                        JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null,
                        treeIDs, trees.get(0));
            }
            if (treeIndex < 0) {
                return null;
            }
            prior.treeInput.setValue(trees.get(treeIndex), prior);
            TaxonSet taxonSet = new TaxonSet();

            TaxonSetDialog dlg = new TaxonSetDialog(taxonSet, getTaxonCandidates(prior), doc);
            if (!dlg.showDialog() || dlg.taxonSet.getID() == null || dlg.taxonSet.getID().trim().equals("")) {
                return null;
            }
            taxonSet = dlg.taxonSet;
            if (taxonSet.taxonsetInput.get().size() == 0) {
            	JOptionPane.showMessageDialog(doc.beauti, "At least one taxon should be included in the taxon set",
            			"Error specifying taxon set", JOptionPane.ERROR_MESSAGE);
            	return null;
            }
            int i = 1;
            String id = taxonSet.getID();
            while (doc.pluginmap.containsKey(taxonSet.getID()) && doc.pluginmap.get(taxonSet.getID()) != taxonSet) {
            	taxonSet.setID(id + i);
            	i++;
            }
            BEASTObjectPanel.addPluginToMap(taxonSet, doc);
            prior.taxonSetInput.setValue(taxonSet, prior);
            prior.setID(taxonSet.getID() + ".prior");
            
            // set up a region
            File kmlFile = Utils.getLoadFile("KML file specifying a region", new File(ProgramStatus.g_sDir), "KML file", "kml");
            if (kmlFile == null) {
            	return null;
            }
            KMLRegion region = new KMLRegion();
            region.setID(taxonSet.getID() + ".region");
            region.kmlFileInput.setValue(kmlFile.getPath(), region);
            prior.regionInput.setValue(region, prior);
            
            // set up location parameter for prior
    		BEASTInterface location = null;
            for (BEASTInterface partition : doc.getPartitions("SiteModel")) {
            	if (partition instanceof ApproxMultivariateTraitLikelihood) {
	            	String partitionName = BeautiDoc.parsePartition(partition.getID());
	            	if (doc.pluginmap.containsKey("location.s:" + partitionName)) {
	            		location = doc.pluginmap.get("location.s:" + partitionName);
	            		prior.locationInput.setValue(location, prior);
	            		break;
	            	}
            	}
            }
	
	
    		// add location parameter to the state
	        State state = (State) doc.pluginmap.get("state");
	        boolean alreadyInState = false;
	        for (StateNode node : state.stateNodeInput.get()) {
	        	if (node == location) {
	        		alreadyInState = true;
	        		break;
	        	}
	        }
	        if (!alreadyInState) {
	        	state.stateNodeInput.setValue(location, state);
	        }

	        // add location operator, if not already present
	        boolean alreadyAdded = false;
	        for(Operator operator : ((MCMC) doc.mcmc.get()).operatorsInput.get()) {
	        	if (operator.getID().equals("location.sampler")) {
	        		alreadyAdded = true;
	        		break;
	        	}
	        }
	        if (!alreadyAdded) {
		        Operator operator = new LocationOperator();
		        operator.setID("location.sampler");
				CompoundDistribution likelihoods = (CompoundDistribution) doc.pluginmap.get("likelihood");
		        BEASTObject locationtreeLikelihood = null;
				for (Distribution likelihood : likelihoods.pDistributions.get()) {
					if (likelihood instanceof sphericalGeo.ApproxMultivariateTraitLikelihood) {
						locationtreeLikelihood = likelihood;
					}
				}
				operator.initByName("location", location, "likelihood", locationtreeLikelihood, "weight", 30.0);
		    	doc.mcmc.get().setInputValue("operator", operator);
		    	doc.registerPlugin(operator);
	        }

	        
            // add to trace log
            Logger logger = (Logger) doc.pluginmap.get("tracelog");
            logger.loggersInput.setValue(prior, logger);
        } catch (Exception e) {
        	return new ArrayList<>();
        }
        List<Distribution> selectedPlugins = new ArrayList<>();
        selectedPlugins.add(prior);
        //g_collapsedIDs.add(prior.getID());
        return selectedPlugins;
	}

    Set<Taxon> getTaxonCandidates(MRCAPrior prior) {
        Set<Taxon> candidates = new HashSet<>();
        Tree tree = prior.treeInput.get();
        String [] taxa = null;
        if (tree.m_taxonset.get() != null) {
        	try {
            	TaxonSet set = tree.m_taxonset.get();
        		set.initAndValidate();
            	taxa = set.asStringList().toArray(new String[0]);
        	} catch (Exception e) {
            	taxa = prior.treeInput.get().getTaxaNames();
			}
        } else {
        	taxa = prior.treeInput.get().getTaxaNames();
        }
        
        for (String taxon : taxa) {
            candidates.add(doc.getTaxon(taxon));
        }
        return candidates;
    }

	@Override
	public String getDescription() {
		return "Geographical prior";
	}
	
	@Override
	public boolean canProvidePrior(BeautiDoc doc) {
        for (BEASTInterface siteModel : doc.getPartitions("SiteModel")) {
        	String partition = BeautiDoc.parsePartition(siteModel.getID());
        	if (doc.pluginmap.containsKey("location.s:" + partition)) {
        		return true;
        	}
        }
        return false;
	}

    public InputEditor createRegionEditor() {
    	Input<?> input = ((GeoPrior) m_beastObject).regionInput;
    	KMLRegion region = (KMLRegion) input.get();
    	
    	StringInputEditor inputEditor = new StringInputEditor(doc);
    	//FileInputEditor inputEditor = new FileInputEditor(doc);
    	inputEditor.init(region.kmlFileInput, region, itemNr, ExpandOption.FALSE, true);
    	
        // increase size of newick text editor
        for (Node node : inputEditor.getChildren()) {
        	if (node instanceof TextField) {
        		((TextField)node).setPrefColumnCount(30);
        	}
//        	if (inputEditor.getComponent(i) instanceof Label) {
//        		int fontSize = ((Label)inputEditor.getComponent(i)).getFont().getSize();
//        		((Label)inputEditor.getComponent(i)).setPreferredSize(new Dimension(40 * fontSize / 13, 20 * fontSize / 13));        		
//        	}
        }

    	return inputEditor;
    }

	// suppress tree input editing 
	public InputEditor createTreeEditor() {return new DummyEditor(doc);}
	public InputEditor createTaxonEditor() {return new DummyEditor(doc);}
	public InputEditor createTaxonsetEditor() {return new DummyEditor(doc);}
	public InputEditor createLocationEditor() {return new DummyEditor(doc);}
	//public InputEditor createAllInternalNodesEditor() {return new DummyEditor(doc);}
	

	class DummyEditor extends StringInputEditor {
		public DummyEditor(BeautiDoc doc) {
			super(doc);
		}

		@Override
		public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
				boolean addButtons) {
			addInputLabel();
		}
		
	};
}
