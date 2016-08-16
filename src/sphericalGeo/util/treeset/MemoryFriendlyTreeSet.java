package sphericalGeo.util.treeset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

public class MemoryFriendlyTreeSet extends TreeSet {
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
       
        
        @Override
        public int size() {
        	return totalTrees - burninCount;
        };
        
        public MemoryFriendlyTreeSet(String inputFileName, int burninPercentage) throws IOException  {
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

    	@Override
    	public void reset() throws FileNotFoundException  {
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

    	
    	
    	@Override
    	public boolean hasNext() {
    		return current < totalTrees;
    	}
    	
    	@Override
    	public Tree next() throws IOException {
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
