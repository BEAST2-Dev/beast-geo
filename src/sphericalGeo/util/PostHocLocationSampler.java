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
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeWithMetaDataLogger;
import beast.util.LogAnalyser;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import sphericalGeo.AlignmentFromTraitMap;
import sphericalGeo.LocationProvider;
import sphericalGeo.PFApproxMultivariateTraitLikelihood;
import sphericalGeo.SphericalDiffusionModel;
import sphericalGeo.TraitFunction;
import sphericalGeo.TreeTraitMap;
import sphericalGeo.util.treeset.MemoryFriendlyTreeSet;

@Description("Post Hoc Location Sampler produces location sample from posterior tree and trace sample")
public class PostHocLocationSampler extends beast.core.Runnable {
//	final public Input<XMLFile> xmlInput = new Input<>("xml", "BEAST XML file containing a spherical diffusion analysis",
//			new XMLFile("[[none]]"));
	final public Input<TreeFile> treeInput = new Input<>("tree", "tree file containing posterior sample of above XML",
			new TreeFile("[[none]]"));
	final public Input<OutFile> traceInput = new Input<>("trace", "trace file containing posterior sample of above XML",
			new OutFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	final public Input<OutFile> outputInput = new Input<>("out", "output file where trees annotated with locatiosn will be stored."
			+ " Stdout is used if not specified",
			new OutFile("[[none]]"));
	final public Input<String> locationtagInput = new Input<>("locationtag", "meta data tag used in node for the location", "location");
	final public Input<String> precisionTagInput = new Input<>("precisionTag", "trace tag used for the precision", "precision");
	final public Input<String> clockRateTagInput = new Input<>("clockRateTag", "trace tag used for the mean clock rate", "ucld.mean");

	@Override
	public void initAndValidate() {
	}

	
//	private LocationProvider locationProvider;
	private RealParameter precision, clockRate;
	private int precisionIndex, clockRateIndex;
	private LogAnalyser trace;
	// private Tree tree;
//	private State state;
//	private Distribution posterior;
	// private TreeWithMetaDataLogger treeLogger;
	// private PostHocBranchRateModel clockModel;
	private String locationtag;
	
	@Override
	public void run() throws Exception {
		locationtag = locationtagInput.get();
		precision = new RealParameter("1000.0");
		precision.setID(precisionTagInput.get());
		clockRate = new RealParameter("1.0");
		clockRate.setID(clockRateTagInput.get());
		
//		if (xmlInput.get().getName().equals("[[none]]")) {
//			throw new IllegalArgumentException("XML file must be specified");
//		}
		if (treeInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("Tree file must be specified");
		}
		if (traceInput.get().getName().equals("[[none]]")) {
			throw new IllegalArgumentException("Trace file must be specified");
		}
		
		
		// parseXML();
		parseTrace();
		
		PrintStream out = System.out;
		if (!outputInput.get().getName().equals("[[none]]")) {
			out = new PrintStream(outputInput.get());
			Log.warning("Output to " + outputInput.get().getPath());
		}
		
		
		MemoryFriendlyTreeSet srcTreeSet = new MemoryFriendlyTreeSet(treeInput.get().getPath(), burnInPercentageInput.get());
		srcTreeSet.reset();
		Tree tree = srcTreeSet.next();


	    SphericalDiffusionModel substModel = new SphericalDiffusionModel();
	    substModel.initByName("fast", true, "precision", precision, "threshold", 1.0);
	    SiteModel siteModel = new SiteModel();
	    siteModel.initByName("substModel", substModel);
	    RealParameter location = new RealParameter();
	    location.initByName("minordimension", 2, "value", "0.0 0.0");
	    

	    StringBuilder traitvalue = new StringBuilder();
	    int n = tree.getLeafNodeCount();
	    for (int i = 0; i < n; i++) {
	    	Node leaf = tree.getNode(i);
	    	String metaData = leaf.metaDataString;
	    	traitvalue.append(parseLocation(leaf.getID(), metaData));
	    	if (i < n - 1) {
	    		traitvalue.append(",\n");
	    	}
	    }
	    
	    
		
	    TreeTraitMap traitMap = new TreeTraitMap();
	    traitMap.initByName("parameter", location, 
	    		"initByMean", true, 
	    		"randomizelower", "-90 -180", 
	    		"randomizeupper", "90 180", 
	    		"traitName", "location", 
	    		"tree", tree,
	    		"value", traitvalue.toString());
	    AlignmentFromTraitMap data = new AlignmentFromTraitMap();
	    data.initByName("userDataType", new sphericalGeo.LocationDataType(),
	    				"traitMap", traitMap);
		
		PostHocBranchRateModel clockModel = new PostHocBranchRateModel();
		clockModel.initByName("tag", "rate");

		PFApproxMultivariateTraitLikelihood loggerLikelihood = new PFApproxMultivariateTraitLikelihood();
		loggerLikelihood.initByName("scale", true, 
				"tree", tree, 
				"siteModel", siteModel, 
				"branchRateModel", clockModel,
				"location", location,
				"data", data, 
				"transformer", null);

		TreeWithMetaDataLogger treeLogger = new TreeWithMetaDataLogger();
		TraitFunction traitFunction = new TraitFunction();
		traitFunction.setID("location");
		traitFunction.initByName("likelihood", loggerLikelihood, "value", "0.0");
		treeLogger.initByName("tree", tree, "branchratemodel", clockModel, "metadata", traitFunction);
		
		srcTreeSet.reset();
		int k = 0;
		while (srcTreeSet.hasNext()) {
	        // state.store(k);
			Tree tree2 = srcTreeSet.next();
			tree.assignFrom(tree2);
			if (precisionIndex >= 0) {
				precision.setValue(trace.getTrace(precisionIndex)[k]);
			}
			if (clockRateIndex >= 0){
				clockRate.setValue(trace.getTrace(clockRateIndex)[k]);
			}
			
			loggerLikelihood.calculateLogP();
//			state.robustlyCalcPosterior((Distribution) locationProvider);
//			state.acceptCalculationNodes();
//            state.setEverythingDirty(false);

            if (k == 0) {
            	tree.init(out);
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

	private String parseLocation(String id, String str) {
		int i = str.indexOf(locationtag);
		if (i < 0) {
			throw new IllegalArgumentException("Expected metadata to contain tag " + locationtag);
		}
		
		str = str.substring(i + locationtag.length() + 1);
		str = str.substring(1, str.indexOf('}'));
		str = str.replaceAll(",", " ");
		String location = id + "=" + str;
		return location;
	}

	private void parseTrace() throws IOException {
		trace = new LogAnalyser(traceInput.get().getPath(), burnInPercentageInput.get(), false, false);
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

////	private void parseXML() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
////		XMLParser parser = new XMLParser();
////		BEASTInterface o = parser.parseFile(xmlInput.get());
////		posterior = (Distribution) o.getInput("distribution").get();
////		traverse(o, new HashSet<BEASTInterface>());
////		
////		if (locationProvider == null) {
////			throw new IllegalArgumentException("Could not find locationProvider in XML."
////					+ " The XML must contain a spherical diffusion analysis.");
////		}
////		if (precision == null) {
////			throw new IllegalArgumentException("Could not find precision parameter in XML."
////					+ " The XML must contain a SphericalDiffusionModel.");
////		}
////		if (tree == null) {
////			throw new IllegalArgumentException("Could not find tree in XML."
////					+ " The XML must contain a tree.");
////		}
////		if (state == null) {
////			throw new IllegalArgumentException("Could not find state in XML.");
////		}
////		
////		if (clockModel == null) {
////			throw new IllegalArgumentException("Could not find PostHocBranchRateModel in XML. "
////					+ "Perhaps you need to replace the relaxed clock with PostHocBranchRateModel.");
////		}
////		
////		Input<?> logAverageInput = ((BEASTInterface)locationProvider).getInput("logAverage");
////		if ((Boolean)logAverageInput.get()) {
////			logAverageInput.set(false);
////			((BEASTInterface)locationProvider).initAndValidate();
////		}
////	}
//
//	private void traverse(BEASTInterface o, HashSet<BEASTInterface> done) {
//		for (BEASTInterface o2 : o.listActiveBEASTObjects()) {
//			if (o2 instanceof LocationProvider) {
//				locationProvider = (LocationProvider) o2;
//				if (o2 instanceof GenericTreeLikelihood) {
//					GenericTreeLikelihood likelihood = (GenericTreeLikelihood) o2;
//					BranchRateModel clockModel = likelihood.branchRateModelInput.get();
//					Input<?> input = ((BEASTInterface)clockModel).getInput("clock.rate");
//					clockRate = (RealParameter) input.get();
//				}
//			}
//			if (o2 instanceof SphericalDiffusionModel) {
//				SphericalDiffusionModel m = (SphericalDiffusionModel) o2;
//				precision = m.precisionInput.get();
//			}
//			if (o2 instanceof Tree) {
//				tree = (Tree) o2;
//			}
//			if (o2 instanceof State) {
//				state = (State) o2;
//			}
//			if (o2 instanceof TreeWithMetaDataLogger) {
//				treeLogger = (TreeWithMetaDataLogger) o2;
//			}
//			if (o2 instanceof PostHocBranchRateModel) {
//				clockModel = (PostHocBranchRateModel) o2;
//			}
//			if (!done.contains(o2)) {
//				traverse(o2, done);
//			}
//		}
//		done.add(o);
//	}

	public static void main(String[] args) throws Exception {
		new Application(new PostHocLocationSampler(), "Post Hoc Location Sampler", args);
	}
}
