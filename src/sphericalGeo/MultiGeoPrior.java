package sphericalGeo;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.MonoCladesMapping;
import beast.math.distributions.MultiMonophyleticConstraint;

/** adapted from MultiMRCAPrors by Joseph Heled **/
@Description("Set of GeoPriors, efficiently managed")
public class MultiGeoPrior extends MultiMonophyleticConstraint {

	final public Input<List<GeoPrior>> geopriorsInput = new Input<>("geoprior", "Set of geographical priors",
			new ArrayList<>());

	private boolean initialised = false;
	private boolean needsUpdate = true;
	private int[] nodeToCladeGroup = null;
	private int[] cladeTopNodeNrs = null;
	private RealParameter location;
	private List<GeoPrior> geoPriors;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		location = geopriorsInput.get().get(0).locationInput.get();

		// sanity check
		geoPriors = geopriorsInput.get();
		for (GeoPrior prior : geoPriors) {
			if (prior.allInternalNodesInput.get()) {
				throw new IllegalArgumentException("GeoPrior " + prior.getID()
						+ " has allInternalNodes=\"true\", which should be false for use with MultiGeoPrior");
			}
			if (prior.taxonSetInput.get() == null && prior.taxonInput.get() == null) {
				throw new IllegalArgumentException(
						"GeoPrior " + prior.getID() + " should have a \"taxonset\" or \"taxon\" specified");
			}
		}

	}

	@Override
	protected void parse(String newick) {
		super.parse(newick);
		for (GeoPrior m : geopriorsInput.get()) {
			if (!m.tree.equals(this.tree)) {
				throw new IllegalArgumentException("All constraints must be on the same tree");
			}
			List<Integer> list = new ArrayList<Integer>();
			if (m.taxonSetInput.get() !=null) {
				for (String taxon : m.taxonSetInput.get().getTaxaNames()) {
					list.add(indexOf(taxon));
				}
			} else {
				String taxon = m.taxonInput.get().getID();
				list.add(indexOf(taxon));
			}
			boolean add = true;
			Set<Integer> slist = new HashSet<>(list);

			for (List<?> l : taxonIDList) {
				if (l.size() == slist.size() && slist.containsAll(l)) {
					add = false;
					break;
				}
			}
			if (add) {
				taxonIDList.add(list);
			}
		}
	}

	@Override
	public double calculateLogP() {
		if (!needsUpdate) {
			return logP;
		}
		logP = super.calculateLogP();
		if (logP != Double.NEGATIVE_INFINITY) {
			if (!initialised) {
				initialise();
				initialised = true;
			}

			for (int k = 0; k < cladeTopNodeNrs.length; ++k) {
				GeoPrior geoPrior = geoPriors.get(k);
				cladeTopNodeNrs[k] = geoPrior.getTaxonNr();
				
//				int nr = cladeTopNodeNrs[k];
//				Node node = tree.getNode(nr);
//				int nodeGroup = nodeToCladeGroup[nr];
//				while (!node.isRoot() && nodeGroup == nodeToCladeGroup[node.getParent().getNr()]) {
//					node = node.getParent();
//				}
//				cladeTopNodeNrs[k] = node.getNr();
//				{
//					boolean ICC = false;
//					if (ICC) {
//						geoPrior.calculateLogP();
//						assert geoPrior.getCommonAncestor().equals(node);
//					}
//				}

				double[] location = new double[2];
				if (geoPrior.region != null) {
					this.location.getMatrixValues1(cladeTopNodeNrs[k], location);
					if (geoPrior.region.isInside(location[0], location[1])) {
						if (!geoPrior.isInsideInput.get()) {
							logP = Double.NEGATIVE_INFINITY;
						}
					} else {
						if (geoPrior.isInsideInput.get()) {
							logP = Double.NEGATIVE_INFINITY;
						}
					}
				}

				if (logP == Double.NEGATIVE_INFINITY) {
					break;
				}
			}
		}
		needsUpdate = false;
		return logP;
	}

	private void initialise() {
		nodeToCladeGroup = MonoCladesMapping.setupNodeGroup(tree, this);

		int geoPriorCount = geoPriors.size();
		cladeTopNodeNrs = new int[geoPriorCount];
		int i = 0;
		for (GeoPrior geoPrior : geoPriors) {
			geoPrior.initialise();
			if (geoPrior.isRoot) {
				cladeTopNodeNrs[i] = tree.getNodeCount() - 1;
			} else {
				geoPrior.calculateLogP(); // init
				cladeTopNodeNrs[i] = geoPrior.getCommonAncestor().getNr();
			}
			i += 1;
		}
	}
	
	
	protected int getCladeTopNodeNr(int k) {
		if (needsUpdate) {
			try {
				calculateLogP();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return cladeTopNodeNrs[k];
	}

	int size() {
		return geoPriors.size();
	}

	@Override
	public void store() {
		super.store();
	}

	@Override
	public void restore() {
		needsUpdate = true;
		super.restore();
	}

	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		return super.requiresRecalculation();
	}

	@Override
	public void sample(State state, Random random) {
	}

	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	public GeoPrior getPrior(int k) {
		return geoPriors.get(k);
	}

	public boolean isMonoPhyletic() {
		try {
			logP = super.calculateLogP();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return Double.isFinite(logP);
	}
	
	
	
	@Override
	public void init(PrintStream out) {
		for (GeoPrior prior : geoPriors) {
			out.append("Latt(" + prior.getID()+")\t");
			out.append("Long(" + prior.getID()+")\t");
		}
		super.init(out);
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		for (GeoPrior prior : geoPriors) {
			int taxonNr = prior.getTaxonNr();
			out.append(location.getValue(taxonNr * 2)+"\t");
			out.append(location.getValue(taxonNr * 2 + 1)+"\t");
		}
		super.log(sample, out);
	}
	
	@Override
	public void close(PrintStream out) {
		super.close(out);
	}
}
