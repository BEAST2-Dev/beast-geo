package beast.core;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;

import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.ThreadedTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;

@Description("Logs data for pseudo marginal likelihood estimate (see CPOAnalyser for details)")
public class CPOLogger extends Logger {
	public Input<CompoundDistribution> likelihoodInput = new Input<>("likelihood", "likelihood distribution, containing one or more treelikilihoods", Validate.REQUIRED);
	
	
	public CPOLogger() {
		loggersInput.setRule(Validate.OPTIONAL);
	}
	
	List<GenericTreeLikelihood> likelihoods;
	
	@Override
	public void initAndValidate() {
		RealParameter p = new RealParameter("1.0");
		p.setID("dummy");
		loggersInput.setValue(p, this);
		super.initAndValidate();
		loggersInput.get().clear();
		
		likelihoods = new ArrayList<>();
		CompoundDistribution likelihood = likelihoodInput.get();
		for (Distribution d : likelihood.pDistributions.get()) {
			if (d instanceof TreeLikelihood || d instanceof ThreadedTreeLikelihood) {
				likelihoods.add((GenericTreeLikelihood) d);
			}
		}
	}

	@Override
    public void init() throws IOException {
        final boolean needsHeader = openLogFile();
        if (needsHeader) {
            final ByteArrayOutputStream rawbaos = new ByteArrayOutputStream();
            final PrintStream out = new PrintStream(rawbaos);

            // list of pattern weights
            out.print("#patternweights\t");
            for (GenericTreeLikelihood d : likelihoods) {
            	int [] weights = d.dataInput.get().getWeights();
            	for (int w : weights) {
            		out.print(w + "\t");
            	}
            }
            out.print("\n");             
            
            // header line
            out.print("Sample\t");
            int k = 0;
            for (GenericTreeLikelihood d : likelihoods) {
            	int [] weights = d.dataInput.get().getWeights();
            	for (int w : weights) {
            		out.print("logP" + k + "\t");
            		k++;
            	}
            }

            // Remove trailing tab from header
            String header = rawbaos.toString().trim();

           	m_out.print(header);
            
            m_out.println();
        }
    } // init
	
	
	@Override
	public void log(int sampleNr) {
	        if ((sampleNr < 0) || (sampleNr % every > 0)) {
	            return;
	        }
	        if (sampleOffset >= 0) {
	            if (sampleNr == 0) {
	                // don't need to duplicate the last line in the log
	                return;
	            }
	            sampleNr += sampleOffset;
	        }

	        ByteArrayOutputStream baos = new ByteArrayOutputStream();
	        PrintStream out = new PrintStream(baos);

	        out.print((sampleNr) + "\t");

            for (GenericTreeLikelihood d : likelihoods) {
        		double [] logPs = null;
            	if (d instanceof TreeLikelihood) {
            		TreeLikelihood tl = (TreeLikelihood) d;
            		tl.calculateLogP();
            		logPs = tl.getPatternLogLikelihoods();
            	} else {
            		ThreadedTreeLikelihood tl = (ThreadedTreeLikelihood) d;
            		tl.calculateLogP();
                	logPs = tl.getPatternLogLikelihoods();
            	}
        		for (double f : logPs) {
        			out.print(f + "\t");
        		}
            }
            
	        // Acquire log string and trim excess tab
	        String logContent;
	        try {
	            logContent = baos.toString("ASCII").trim();
	        } catch (UnsupportedEncodingException e) {
	            throw new RuntimeException("ASCII string encoding not supported: required for logging!");
	        }

	        m_out.println(logContent);
	    } // log
	
	@Override
	public void close() {
		super.close();

//		CPOAnalyser analyser = new CPOAnalyser();
//		analyser.cpoLogFileInput.setValue(fileNameInput.get(), analyser);
//		try {
//			analyser.run();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
	}
} // CPOLogger
