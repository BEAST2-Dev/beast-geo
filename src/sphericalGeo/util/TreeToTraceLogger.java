package sphericalGeo.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beastfx.app.util.XMLFile;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Param;
import beast.base.inference.Runnable;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.evolution.tree.TreeParser;
import beast.base.parser.NexusParser;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;

@Description("Creates trace log from a tree log and XML file containing MRCAPriors. "
		+ "For each MRCAPrior, height and location are logged.")
public class TreeToTraceLogger extends Runnable {
	XMLFile xmlFile;
	TreeFile treeFile;
	OutFile outFile;
	Integer burnin = 10;
	
	List<MRCAPrior> mrcas;
	   // array of flags to indicate which taxa are in the set
    List<Set<String>> isInTaxaSet;
    
	Tree tree;
	
	public TreeToTraceLogger() {}	
	public TreeToTraceLogger(@Param(name="xmlFile", description="XML file containing MRCAPriors (and possibly many more items. but these will be ignored)") XMLFile xmlFile,
			@Param(name="treeFile", description="file containing tree log") TreeFile treeFile,
			@Param(name="out", description="file to store trace log, if not specified, results are printed on stdout") OutFile out,
			@Param(name="burnin", description="burn-in percentage, defaults to 10", defaultValue="10") Integer burnin
			) {
		this.xmlFile = xmlFile;
		this.treeFile = treeFile;		
		this.outFile = out;
		this.burnin = burnin;
	}
	

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		grabMRCAs();
		
		PrintStream out = (outFile == null? System.out : new PrintStream(outFile));
		// print header
		out.print("Sample\t");
		for (MRCAPrior p : mrcas) {
			//sout.print(p.getID() + ".node\t");
			out.print(p.getID() + ".height\t");
		}
		for (MRCAPrior p : mrcas) {
			out.print(p.getID() + ".lat\t");
			out.print(p.getID() + ".long\t");
		}
		out.println();
		
		processTrees(out);
	}

	
	private void processTrees(PrintStream out) throws IOException {
		Log.warning.println("loading trees");
		MemoryFriendlyTreeSet parser = new MemoryFriendlyTreeSet(treeFile.getPath(), burnin);
		parser.reset();
		Log.warning.println("processing trees");
		int sample = 0;
		while (parser.hasNext()) {
			Tree  tree = parser.next();
			this.tree.assignFrom(tree);
			out.print(sample + "\t");
			locations = new StringBuilder();
			for (int k = 0; k < mrcas.size(); k++) {
				printData(isInTaxaSet.get(k), mrcas.get(k).useOriginateInput.get(), out);
			}
			out.print(locations.toString());
			sample++;
			out.println();
		}
		Log.warning.println("Done");
	}
	
	Node MRCATime;
	StringBuilder locations;
	class FastTreeSet  {
    	int current = 0;
    	Tree [] trees;
    	
    	public FastTreeSet(String inputFileName, int burninPercentage) throws IOException  {
            Log.warning.println("0              25             50             75            100");
            Log.warning.println("|--------------|--------------|--------------|--------------|");
            NexusParser parser = new NexusParser();
	      	parser.parseFile(new File(inputFileName));
	      	List<Tree> treesInFile = parser.trees;
	      	
	      	int n0 = treesInFile.size();
	      	int start = n0 * burninPercentage / 100;
	      	int n = n0 - start;
	      	trees = new Tree[n];
	      	for (int i = start; i < n0; i++) {
	      		trees[i-start] = treesInFile.get(i);
	      	}
		}

		boolean hasNext() {
			return current < trees.length;
		}

		Tree next()  {
			return trees[current++];
		}

		void reset()  {
			current = 0;
		}
    }	

	private void printData(Set<String> taxonset, boolean useOriginate, PrintStream out) {
		calcMRCAtime(tree.getRoot(), taxonset, useOriginate, new int[1]);
		Node node = MRCATime;
		//out.print(node.getNr() + "\t");
		out.print(node.getHeight() + "\t");
		String str = node.metaDataString;
		int i = str.indexOf("location={");
		if (i >= 0) {
			int j = str.indexOf("}", i);
			String substr = str.substring(i + 10, j-1);
			String [] strs = substr.split(",");
			locations.append(strs[0] + "\t");
			locations.append(strs[1] + "\t");
		} else {
			locations.append("-1000\t-1000\t");
		}
	}
	
	
	/**
     * Recursively visit all leaf nodes, and collect number of taxa in the taxon
     * set. When all taxa in the set are visited, record the time.
     * *
     * @param node
     * @param taxonCount2
     */
    int calcMRCAtime(final Node node, Set<String> isInTaxaSet, boolean useOriginate, final int[] taxonCount2) {
    	final int nrOfTaxa = isInTaxaSet.size();
        if (node.isLeaf()) {
            taxonCount2[0]++;
            if (isInTaxaSet.contains(node.getID())) {
                return 1;
            } else {
                return 0;
            }
        } else {
            int taxonCount = calcMRCAtime(node.getLeft(), isInTaxaSet, useOriginate, taxonCount2);
            final int leftTaxa = taxonCount2[0];
            taxonCount2[0] = 0;
            if (node.getRight() != null) {
                taxonCount += calcMRCAtime(node.getRight(), isInTaxaSet, useOriginate, taxonCount2);
                final int rightTaxa = taxonCount2[0];
                taxonCount2[0] = leftTaxa + rightTaxa;
                if (taxonCount == nrOfTaxa) {
                	if (nrOfTaxa == 1 && useOriginate) {
            			MRCATime = node;
                        return taxonCount + 1;
                	}
                    // we are at the MRCA, so record the height
                	if (useOriginate) {
                		Node parent = node.getParent();
                		if (parent != null) {
                			MRCATime = parent;
                		} else {
                			MRCATime = node;
                		}
                	} else {
                		MRCATime = node;
                	}
                    return taxonCount + 1;
                }
            }
            return taxonCount;
        }
    }

	private void grabMRCAs() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		XMLParser parser = new XMLParser();
		BEASTInterface o = parser.parseFile(xmlFile);
		Set<MRCAPrior> set = new HashSet<>();
		scanForMRCAPriors(o, set);
		mrcas = new ArrayList<>();
		mrcas.addAll(set);
		
		// sanity check
		if (mrcas.size() == 0) {
			throw new IllegalArgumentException("No MRCAPriors found in XML file. There should be at least one MRCAPrior");
		}
		
		// sort alphabetically
		Collections.sort(mrcas, (o1, o2) -> {
			if (o1.getID() == null || o2.getID() == null) {
				return 0;
			}
			return o1.getID().compareTo(o2.getID());
		});
		
		tree  = mrcas.get(0).treeInput.get();
		// sanity check: make sure we are only dealing with 1 tree
		for (MRCAPrior prior : mrcas) {
			if (prior.treeInput.get() != tree) {
				throw new IllegalArgumentException("Multiple trees found. There should be only one tree for all MRCAPriors");
			}
		}
		
		isInTaxaSet = new ArrayList<>();
		for (MRCAPrior prior : mrcas) {
			Set<String> taxonset = new LinkedHashSet<>();
			isInTaxaSet.add(taxonset);	        
	        taxonset.addAll(prior.taxonsetInput.get().asStringList());
		}
	}

	private void scanForMRCAPriors(BEASTInterface o, Set<MRCAPrior> set) {
		for (Object o2 : o.listActiveBEASTObjects()) {
			if (o2 instanceof MRCAPrior) {
				set.add((MRCAPrior) o2);
			}
			if (o2 instanceof BEASTInterface) {
				scanForMRCAPriors((BEASTInterface) o2, set);
			}
		}
	}
	public XMLFile getXmlFile() {
		return xmlFile;
	}
	public void setXmlFile(XMLFile xmlFile) {
		this.xmlFile = xmlFile;
	}
	public TreeFile getTreeFile() {
		return treeFile;
	}
	public void setTreeFile(TreeFile treeFile) {
		this.treeFile = treeFile;
	}
	public OutFile getOut() {
		return outFile;
	}
	public void setOut(OutFile out) {
		this.outFile = out;
	}
	public Integer getBurnin() {
		return burnin;
	}
	public void setBurnin(Integer burnin) {
		this.burnin = burnin;
	}

	public static void main(String[] args) throws Exception {
		new Application(new TreeToTraceLogger(), "TreeToTraceLogger", args);
	}
	
    class MemoryFriendlyTreeSet {
//    	Tree [] trees;
    	int current = 0;
    	int lineNr;
        public Map<String, String> translationMap = null;
        public List<String> taxa;
    	
        int burninCount = 0;
        int totalTrees = 0;
        boolean isNexus = true;
        BufferedReader fin;
        String inputFileName;
        // label count origin for NEXUS trees
        int origin = -1;
       
        MemoryFriendlyTreeSet(String inputFileName, int burninPercentage) throws IOException  {
    		this.inputFileName = inputFileName;
    		init(burninPercentage);
        	Log.warning.println("Processing " + (totalTrees - burninCount) + " trees from file" +
                    (burninPercentage > 0 ? " after ignoring first " + burninPercentage + "% = " + burninCount + " trees." : "."));
    		
    		
    	}

    	/** determine number of trees in the file,
    	 * and number of trees to skip as burnin 
    	 * @throws IOException 
    	 * @throws FileNotFoundException **/
    	private void init(int burninPercentage) throws IOException  {
            fin = new BufferedReader(new FileReader(new File(inputFileName)));
            if (!fin.ready()) {
            	throw new IOException("File appears empty");
            }
        	String str = nextLine();
            if (!str.toUpperCase().trim().startsWith("#NEXUS")) {
            	// the file contains a list of Newick trees instead of a list in Nexus format
            	isNexus = false;
            	if (str.trim().length() > 0) {
            		totalTrees = 1;
            	}
            }
            while (fin.ready()) {
            	str = nextLine();
                if (isNexus) {
                    if (str.trim().toLowerCase().startsWith("tree ")) {
                    	totalTrees++;
                    }
                } else if (str.trim().length() > 0) {
            		totalTrees++;
                }            	
            }
            fin.close();
            
            burninCount = Math.max(0, (burninPercentage * totalTrees)/100);
		}

		void reset() throws FileNotFoundException  {
    		current = 0;
            fin = new BufferedReader(new FileReader(new File(inputFileName)));
            lineNr = 0;
            try {
                while (fin.ready()) {
                    final String str = nextLine();
                    if (str == null) {
                        return;
                    }
                    final String lower = str.toLowerCase();
                    if (lower.matches("^\\s*begin\\s+trees;\\s*$")) {
                        parseTreesBlock();
                        return;
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException("Around line " + lineNr + "\n" + e.getMessage());
            }
        } // parseFile

        /**
         * read next line from Nexus file that is not a comment and not empty 
         * @throws IOException *
         */
        String nextLine() throws IOException  {
            String str = readLine();
            if (str == null) {
                return null;
            }
            if (str.contains("[")) {
                final int start = str.indexOf('[');
                int end = str.indexOf(']', start);
                while (end < 0) {
                    str += readLine();
                    end = str.indexOf(']', start);
                }
                str = str.substring(0, start) + str.substring(end + 1);
                if (str.matches("^\\s*$")) {
                    return nextLine();
                }
            }
            if (str.matches("^\\s*$")) {
                return nextLine();
            }
            return str;
        }

        /**
         * read line from nexus file *
         */
        String readLine() throws IOException {
            if (!fin.ready()) {
                return null;
            }
            lineNr++;
            return fin.readLine();
        }

        private void parseTreesBlock() throws IOException  {
            // read to first non-empty line within trees block
            String str = fin.readLine().trim();
            while (str.equals("")) {
                str = fin.readLine().trim();
            }

            // if first non-empty line is "translate" then parse translate block
            if (str.toLowerCase().contains("translate")) {
                translationMap = parseTranslateBlock();
                origin = getIndexedTranslationMapOrigin(translationMap);
                if (origin != -1) {
                    taxa = getIndexedTranslationMap(translationMap, origin);
                }
            }
            // we got to the end of the translate block
            // read bunrinCount trees
            current = 0;
            while (current < burninCount && fin.ready()) {
    			str = nextLine();
                if (str.toLowerCase().startsWith("tree ")) {
                	current++;
                }
            }
        }

        private List<String> getIndexedTranslationMap(final Map<String, String> translationMap, final int origin) {

            //System.out.println("translation map size = " + translationMap.size());

            final String[] taxa = new String[translationMap.size()];

            for (final String key : translationMap.keySet()) {
                taxa[Integer.parseInt(key) - origin] = translationMap.get(key);
            }
            return Arrays.asList(taxa);
        }

        /**
         * @param translationMap
         * @return minimum key value if keys are a contiguous set of integers starting from zero or one, -1 otherwise
         */
        private int getIndexedTranslationMapOrigin(final Map<String, String> translationMap) {

            final SortedSet<Integer> indices = new java.util.TreeSet<>();

            int count = 0;
            for (final String key : translationMap.keySet()) {
                final int index = Integer.parseInt(key);
                indices.add(index);
                count += 1;
            }
            if ((indices.last() - indices.first() == count - 1) && (indices.first() == 0 || indices.first() == 1)) {
                return indices.first();
            }
            return -1;
        }

        /**
         * @return a map of taxa translations, keys are generally integer node number starting from 1
         *         whereas values are generally descriptive strings.
         * @throws IOException
         */
        private Map<String, String> parseTranslateBlock() throws IOException {

            final Map<String, String> translationMap = new HashMap<>();

            String line = readLine();
            final StringBuilder translateBlock = new StringBuilder();
            while (line != null && !line.trim().toLowerCase().equals(";")) {
                translateBlock.append(line.trim());
                line = readLine();
            }
            final String[] taxaTranslations = translateBlock.toString().split(",");
            for (final String taxaTranslation : taxaTranslations) {
                final String[] translation = taxaTranslation.split("[\t ]+");
                if (translation.length == 2) {
                    translationMap.put(translation[0], translation[1]);
//                    System.out.println(translation[0] + " -> " + translation[1]);
                } else {
                    Log.err.println("Ignoring translation:" + Arrays.toString(translation));
                }
            }
            return translationMap;
        }

    	
    	
		boolean hasNext() {
    		return current < totalTrees;
    	}
    	
		Tree next() throws IOException {
			String str = nextLine();
    		if (!isNexus) {
                TreeParser treeParser;

                if (origin != -1) {
                    treeParser = new TreeParser(taxa, str, origin, false);
                } else {
                    try {
                        treeParser = new TreeParser(taxa, str, 0, false);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        treeParser = new TreeParser(taxa, str, 1, false);
                    }
                }
                return treeParser;
    		}
    		
            // read trees from NEXUS file
            if (str.toLowerCase().startsWith("tree ")) {
            	current++;
                final int i = str.indexOf('(');
                if (i > 0) {
                    str = str.substring(i);
                }
                TreeParser treeParser;

                if (origin != -1) {
                    treeParser = new TreeParser(taxa, str, origin, false);
                } else {
                    try {
                        treeParser = new TreeParser(taxa, str, 0, false);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        treeParser = new TreeParser(taxa, str, 1, false);
                    }
                }

                if (translationMap != null) treeParser.translateLeafIds(translationMap);

                return treeParser;
            }
    		return null;
    	}
    }
}
