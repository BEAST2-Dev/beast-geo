package sphericalGeo.util;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.TreeWithMetaDataLogger;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.Distribution;
import beast.base.inference.Runnable;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beastfx.app.tools.Application;
import beastfx.app.tools.LogAnalyser;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beastfx.app.util.XMLFile;
import sphericalGeo.LocationProvider;
import sphericalGeo.SphericalDiffusionModel;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;

@Description("Post Hoc Location Sampler produces location sample from posterior tree and trace sample")
public class PostHocLocationSampler extends Runnable {
	final public Input<XMLFile> xmlInput = new Input<>("xml", "BEAST XML file containing a spherical diffusion analysis",
			new XMLFile("[[none]]"));
	final public Input<TreeFile> treeInput = new Input<>("tree", "tree file containing posterior sample of above XML",
			new TreeFile("[[none]]"));
	final public Input<OutFile> traceInput = new Input<>("trace", "trace file containing posterior sample of above XML",
			new OutFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	final public Input<OutFile> outputInput = new Input<>("out", "output file where trees annotated with locatiosn will be stored."
			+ " Stdout is used if not specified",
			new OutFile("[[none]]"));
	
	@Override
	public void initAndValidate() {
	}

	
	private LocationProvider locationProvider;
	private RealParameter precision, clockRate;
	private int precisionIndex, clockRateIndex;
	private LogAnalyser trace;
	private Tree tree;
	private State state;
	private Distribution posterior;
	private TreeWithMetaDataLogger treeLogger;
	
	@Override
	public void run() throws Exception {
		if (xmlInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("XML file must be specified");
		}
		if (treeInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("Tree file must be specified");
		}
		if (traceInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("Trace file must be specified");
		}
		
		
		parseXML();
		parseTrace();
		
		PrintStream out = System.out;
		if (!outputInput.get().getName().equals("[[none]]")) {
			out = new PrintStream(outputInput.get());
			Log.warning("Output to " + outputInput.get().getPath());
		}
		
		MemoryFriendlyTreeSet srcTreeSet = new MemoryFriendlyTreeSet(treeInput.get().getPath(), burnInPercentageInput.get());
		srcTreeSet.reset();
		int k = 0;
		while (srcTreeSet.hasNext()) {
	        state.store(k);
			Tree tree2 = srcTreeSet.next();
			tree.assignFrom(tree2);
			if (precisionIndex >= 0) {
				precision.setValue(trace.getTrace(precisionIndex)[k]);
			} else {
				clockRate.setValue(trace.getTrace(clockRateIndex)[k]);
			}
			state.robustlyCalcPosterior(posterior);
			state.acceptCalculationNodes();
            state.setEverythingDirty(false);

			treeLogger.log(k, out);
			k++;
		}


		
		if (!outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}
	}

	private void parseTrace() throws IOException {
		trace = new LogAnalyser(traceInput.get().getPath(), burnInPercentageInput.get(), true, false);
		List<String> labels = trace.getLabels();
		precisionIndex = -1;
		clockRateIndex = -1;
		for (int i = 0; i < labels.size(); i++) {
			String label = labels.get(i);
			if (precision.getID().startsWith(label)) {
				precisionIndex = i;
			}
			if (clockRate.getID().startsWith(label)) {
				clockRateIndex = i;
			}
		}
		if (precisionIndex < 0 && clockRateIndex < 0) {
			throw new IllegalArgumentException("Trace file should contain one of presicion and clock rate, but could not find either");
		}
		if (precisionIndex >= 0 && clockRateIndex >= 0) {
			throw new IllegalArgumentException("Trace file should contain one of presicion and clock rate, but found both");
		}
	}

	private void parseXML() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		XMLParser parser = new XMLParser();
		BEASTInterface o = parser.parseFile(xmlInput.get());
		posterior = (Distribution) o.getInput("distribution").get();
		traverse(o, new HashSet<BEASTInterface>());
		
		if (locationProvider == null) {
			throw new IllegalArgumentException("Could not find locationProvider in XML."
					+ " The XML must contain a spherical diffusion analysis.");
		}
		if (precision == null) {
			throw new IllegalArgumentException("Could not find precision parameter in XML."
					+ " The XML must contain a SphericalDiffusionModel.");
		}
		if (tree == null) {
			throw new IllegalArgumentException("Could not find tree in XML."
					+ " The XML must contain a tree.");
		}
		if (state == null) {
			throw new IllegalArgumentException("Could not find state in XML.");
		}
		
		Input<?> logAverageInput = ((BEASTInterface)locationProvider).getInput("logAverage");
		logAverageInput.set(false);
	}

	private void traverse(BEASTInterface o, HashSet<BEASTInterface> done) {
		for (BEASTInterface o2 : o.listActiveBEASTObjects()) {
			if (o2 instanceof LocationProvider) {
				locationProvider = (LocationProvider) o2;
				if (o2 instanceof GenericTreeLikelihood) {
					GenericTreeLikelihood likelihood = (GenericTreeLikelihood) o2;
					BranchRateModel clockModel = likelihood.branchRateModelInput.get();
					Input<?> input = ((BEASTInterface)clockModel).getInput("clock.rate");
					clockRate = (RealParameter) input.get();
				}
			}
			if (o2 instanceof SphericalDiffusionModel) {
				SphericalDiffusionModel m = (SphericalDiffusionModel) o2;
				precision = m.precisionInput.get();
			}
			if (o2 instanceof Tree) {
				tree = (Tree) o2;
			}
			if (o2 instanceof State) {
				state = (State) o2;
			}
			if (o2 instanceof TreeWithMetaDataLogger) {
				treeLogger = (TreeWithMetaDataLogger) o2;
			}
			if (!done.contains(o2)) {
				traverse(o2, done);
			}
		}
		done.add(o);
	}

	public static void main(String[] args) throws Exception {
		new Application(new PostHocLocationSampler(), "Post Hoc Location Sampler", args);
	}
}
