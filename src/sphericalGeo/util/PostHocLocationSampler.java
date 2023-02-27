package sphericalGeo.util;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.*;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeWithMetaDataLogger;
import beast.util.LogAnalyser;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import sphericalGeo.LocationProvider;
import sphericalGeo.SphericalDiffusionModel;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;

@Description("Post Hoc Location Sampler produces location sample from posterior tree and trace sample")
public class PostHocLocationSampler extends beast.core.Runnable {
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
	private TreeWithMetaDataLogger
	treeLogger;
	private PostHocBranchRateModel clockModel;
	
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
			}
			if (clockRateIndex >= 0){
				clockRate.setValue(trace.getTrace(clockRateIndex)[k]);
			}
			state.robustlyCalcPosterior(posterior);
			state.acceptCalculationNodes();
            state.setEverythingDirty(false);

            if (k == 0) {
            	tree2.init(out);
    			out.println();
            }
			treeLogger.log((long)k, out);
			out.println();
			k++;
			if (k % 10 != 0) {
				Log.warning.print('.');
			} else {
				Log.warning.print('|');
				if (k % 100 == 0) {
					Log.warning.println(k);
				}
			}
		}
		tree.close(out);


		
		if (!outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}
		Log.warning("Done");
	}

	private void parseTrace() throws IOException {
		trace = new LogAnalyser(traceInput.get().getPath(), burnInPercentageInput.get(), true, false);
		List<String> labels = trace.getLabels();
		precisionIndex = -1;
		clockRateIndex = -1;
		String precision2 = sanitise(precision.getID());
		String precision3 = precision2.replaceAll("\\.geo", "");
		String clockRate2 = sanitise(clockRate.getID());
		String clockRate3 = clockRate2.replaceAll("\\.geo", "");
		for (int i = 0; i < labels.size(); i++) {
			String label = labels.get(i);
			if (precision.getID().equals(label) || precision2.equals(label)|| precision3.equals(label)) {
				precisionIndex = i + 1;
			}
			if (clockRate.getID().equals(label) || clockRate2.equals(label)|| clockRate3.equals(label)) {
				clockRateIndex = i + 1;
			}
		}
		if (precisionIndex < 0 && clockRateIndex < 0) {
			throw new IllegalArgumentException("Trace file should contain one of presicion and clock rate, but could not find either");
		}
//		if (precisionIndex >= 0 && clockRateIndex >= 0) {
//			throw new IllegalArgumentException("Trace file should contain one of presicion and clock rate, but found both");
//		}
	}

	private String sanitise(String id) {
    	// remove clock/site/tree info
    	id = id.replaceAll("\\.c:", ".");
    	id = id.replaceAll("\\.t:", ".");
    	id = id.replaceAll("\\.s:", ".");
    	// remove trailing dots on labels
    	id = id.replaceAll("\\.\\.", ".");
    	id = id.replaceAll("\\.\t", "\t");
    	id = id.replaceAll("\\.$", "");
		return id;
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
		
		if (clockModel == null) {
			throw new IllegalArgumentException("Could not find PostHocBranchRateModel in XML. "
					+ "Perhaps you need to replace the relaxed clock with PostHocBranchRateModel.");
		}
		
		Input<?> logAverageInput = ((BEASTInterface)locationProvider).getInput("logAverage");
		if ((Boolean)logAverageInput.get()) {
			logAverageInput.set(false);
			((BEASTInterface)locationProvider).initAndValidate();
		}
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
			if (o2 instanceof PostHocBranchRateModel) {
				clockModel = (PostHocBranchRateModel) o2;
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
