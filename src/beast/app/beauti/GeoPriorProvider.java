package beast.app.beauti;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
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

import beast.app.draw.BEASTObjectInputEditor;
import beast.app.draw.BEASTObjectPanel;
import beast.app.draw.InputEditor;
import beast.app.draw.SmallButton;
import beast.app.draw.StringInputEditor;
import beast.app.draw.InputEditor.ExpandOption;
import beast.app.util.Utils;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Tree;
import beast.math.distributions.MRCAPrior;
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
		
        Box itemBox = Box.createHorizontalBox();

        GeoPrior prior = (GeoPrior) beastObject;
        String text = prior.getID();

        JButton taxonButton = new JButton(text);
//        taxonButton.setMinimumSize(Base.PREFERRED_SIZE);
//        taxonButton.setPreferredSize(Base.PREFERRED_SIZE);
        itemBox.add(taxonButton);
        taxonButton.addActionListener(e -> {
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


        itemBox.add((Component) createRegionEditor());

        JCheckBox isInsidedBox = new JCheckBox(doc.beautiConfig.getInputLabel(prior, prior.isInsideInput.getName()));
        isInsidedBox.setName(text+".isInside");
        isInsidedBox.setSelected(prior.isInsideInput.get());
        isInsidedBox.setToolTipText(prior.isInsideInput.getHTMLTipText());
        isInsidedBox.addActionListener(new GeoPriorActionListener(prior));
        itemBox.add(isInsidedBox);
        
        JCheckBox isAllInternaldBox = new JCheckBox(doc.beautiConfig.getInputLabel(prior, prior.allInternalNodesInput.getName()));
        isAllInternaldBox.setName(text+".allInternalNodes");
        isAllInternaldBox.setSelected(prior.allInternalNodesInput.get());
        isAllInternaldBox.setToolTipText(prior.allInternalNodesInput.getHTMLTipText());
        isAllInternaldBox.addActionListener(new GeoPriorActionListener2(prior));
        itemBox.add(isAllInternaldBox);
        

        
        JButton deleteButton = new SmallButton("-", true);
        deleteButton.setToolTipText("Delete this geo-prior");
        deleteButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Log.warning.println("Trying to delete a geoprior");
				List<?> list = (List<?>) m_input.get();
				GeoPrior prior = (GeoPrior) list.get(itemNr);
				doc.disconnect(prior, "prior", "distribution");
				doc.disconnect(prior, "tracelog", "log");
				doc.unregisterPlugin(prior);
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
					doc.disconnect(prior.locationInput.get(), "state", "stateNode");
				}
				refreshPanel();
			}        	
        });
        itemBox.add(Box.createGlue());
        itemBox.add(deleteButton);

        add(itemBox);
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
            File kmlFile = Utils.getLoadFile("KML file specifying a region", new File(Beauti.g_sDir), "KML file", "kml");
            if (kmlFile == null) {
            	return null;
            }
            KMLRegion region = new KMLRegion();
            region.setID(taxonSet.getID() + ".region");
            region.kmlFileInput.setValue(kmlFile, region);
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
	
	//@Override
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
    	FileInputEditor inputEditor = new FileInputEditor(doc);
    	inputEditor.init(region.kmlFileInput, region, itemNr, ExpandOption.FALSE, true);
    	
        // increase size of newick text editor
        for (int i = 0; i < inputEditor.getComponentCount(); i++) {
        	if (inputEditor.getComponent(i) instanceof JTextField) {
        		((JTextField)inputEditor.getComponent(i)).setColumns(30);
        	}
        	if (inputEditor.getComponent(i) instanceof JLabel) {
        		int fontSize = ((JLabel)inputEditor.getComponent(i)).getFont().getSize();
        		((JLabel)inputEditor.getComponent(i)).setPreferredSize(new Dimension(40 * fontSize / 13, 20 * fontSize / 13));        		
        	}
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
