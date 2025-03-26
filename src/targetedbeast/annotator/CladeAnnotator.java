package targetedbeast.annotator;


import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import beastfx.app.treeannotator.CladeSystem;
import beastfx.app.treeannotator.CladeSystem.Clade;
import beastfx.app.treeannotator.services.NodeHeightSettingService;
import beastfx.app.treeannotator.services.TopologySettingService;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.parser.NexusParser;
import beast.base.util.DiscreteStatistics;
import beast.base.util.HeapSort;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * 
 * TreeAnnotator ported from BEAST 1
 */
public class CladeAnnotator {
    
	
	private CladeSystem cladeSystem = null;
    
    private static class CladeAnnotatorOptions {
        String inFile;
        String logFile;
        String outFile = "summary.tree";
        int burninPercentage = 0;        
        boolean lowMemory = false;
        double cladeCutOff = 0.1;
        double ESSCutOff = 50;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
            		"Log file: " + logFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "Low memory mode: " + lowMemory + "\n" +
                    "Clade credibility cutoff: " + cladeCutOff + "\n";
        }
    }

    
    public CladeAnnotator(CladeAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");
        
        // get the tree likelihood from the logfile, a tsv file
        BufferedReader br = new BufferedReader(new FileReader(options.logFile));
        // skip every line starting with #
        String line = br.readLine();
		while (line.startsWith("#")) {
			line = br.readLine();
		}
		// get the header
		String[] header = line.split("\t");
		// get the index of the likelihood column
		int likelihoodIndex = -1;
		for (int i = 0; i < header.length; i++) {
			if (header[i].equals("likelihood")) {
				likelihoodIndex = i;
				break;
			}
		}
		if (likelihoodIndex == -1) {
			System.out.println("Likelihood column not found in the log file");
			return;
		}

		
        try {
        	if (options.lowMemory) {
        		treeSet = new MemoryFriendlyTreeSet(options.inFile, options.burninPercentage);
        	} else {
        		treeSet = new FastTreeSet(options.inFile, options.burninPercentage);
        	}
        } catch (Exception e) {
        	e.printStackTrace();
        	Log.err.println("Error Parsing Input Tree: " + e.getMessage());
        	return;
        }
        
      		// remove the burnin
        
        
        cladeSystem = new CladeSystem();
        
        treeSet.reset();
        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();
            tree.getLeafNodeCount();
        	cladeSystem.add(tree, false);
            totalTreesUsed++;
        }
        
		// get the likelihood values
		List<Double> likelihoods = new ArrayList<>();
		while ((line = br.readLine()) != null) {
			String[] values = line.split("\t");
			likelihoods.add(Double.parseDouble(values[likelihoodIndex]));
		}
		br.close();

        
        int burninCount = Math.max(0, (options.burninPercentage * totalTreesUsed)/100);
		likelihoods = likelihoods.subList(burninCount, totalTreesUsed+burninCount);

        
		// check that the size of likelihoods is the same as the number of trees
        if (likelihoods.size() != totalTreesUsed) {
        	System.out.println(likelihoods.size() + " " + totalTreesUsed);
        	throw new RuntimeException("Number of likelihoods does not match number of trees");
        }

        

        cladeSystem.calculateCladeCredibilities(totalTreesUsed);
        
        System.out.println("Annotating target tree...");
        System.out.println(cladeSystem.getCladeMap().size());
        List<Clade> clades = new ArrayList<>();
        // remove a clade if it has below 0.1 credibility
        for (Clade clade : cladeSystem.getCladeMap().values()) {
        	if (clade.getCredibility() < options.cladeCutOff || clade.getCredibility() > (1-options.cladeCutOff)) {
        	}else {
        		clades.add(clade);
        	}    
        }
        System.out.println(clades.size());
        
        // for each clade, annotate the tree
        treeSet.reset();
        
        boolean[][] cladePresence = new boolean[clades.size()][totalTreesUsed];
        
        int c = 0;
        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();
            CladeSystem treeCladeSystem = new CladeSystem();
            treeCladeSystem.add(tree, false);
            // check the presence absence of each clade in clades
            for (int i = 0; i < clades.size(); i++) {
	        	Clade clade = clades.get(i);
	        	
	        	if (treeCladeSystem.getCladeMap().values().contains(clade)) {
	        		cladePresence[i][c] = true;
	        	}            	
            }
            c++;
        }
        
        List<Double> smoothLikelihoods = new ArrayList<>();
        int smooth = 5;
		for (int i = 0; i < smooth; i++) {
			smoothLikelihoods.add(likelihoods.get(i));
		}
        
        
		for (int i = smooth; i < likelihoods.size()-smooth; i++) {
			double sum=0.0;
			for (int j = i - smooth; j <= i + smooth; j++) {
				if (j!=i)
					sum += likelihoods.get(j);
			}
			smoothLikelihoods.add(sum / (2*smooth));
		}
        
		for (int i = likelihoods.size()-smooth; i < likelihoods.size(); i++) {
			smoothLikelihoods.add(likelihoods.get(i));
		}
		
        
        // for each cladePresence, calculate the ESS value
        
        double[] ESS = new double[clades.size()];
		for (int i = 0; i < clades.size(); i++) {
            ESS[i] = essTracerStyle(cladePresence[i]);
		}
		System.out.println(Arrays.toString(ESS));
		// print the minimum ESS value
		System.out.println("Minimum ESS value: " + Arrays.stream(ESS).min().getAsDouble());
		
		// print the results for each clade with an ESS below 50
		// print to a log file
		PrintStream ps = new PrintStream(new File(options.outFile));
		ps.print("Sample");
		ps.print("\tLikelihood");
		ps.print("\tSmoothLikelihood");
		for (int i = 0; i < clades.size(); i++) {
			if (ESS[i] < options.ESSCutOff) {
				ps.print("\tClade" + i);
			}
//			System.out.println(clades.get(i).toString());
//			ps.print("\tClade"+ i);
		}
		ps.println();
		for (int j = 0; j < totalTreesUsed; j++) {
            ps.print(j);
            ps.print("\t" + likelihoods.get(j));
            ps.print("\t" + smoothLikelihoods.get(j));
            for (int i = 0; i < clades.size(); i++) {
            	if (ESS[i] < options.ESSCutOff) {
            		ps.print("\t"+ (cladePresence[i][j] ? 1 : 0));
            	}
            }
            ps.println();
        }
		// for each clade with ESS below 50, check how often any other clade with ESS below 50 is non present
		// in the same tree
		
		boolean[][] nonOverlap = new boolean[clades.size()][clades.size()];
		// set all to tru
		for (int i = 0; i < clades.size(); i++) {
			for (int j = 0; j < clades.size(); j++) {
				nonOverlap[i][j] = true;
			}
		}
		
		for (int i = 0; i < clades.size(); i++) {
			if (ESS[i] < options.ESSCutOff) {
				for (int j = 0; j < clades.size(); j++) {
					if (ESS[j] < options.ESSCutOff && i != j) {
						for (int k = 0; k < totalTreesUsed; k++) {
							if (cladePresence[i][k] && cladePresence[j][k]) {
								nonOverlap[i][j] = false;
								break;
							}
						}
					}
				}
			}
		}
		
		
		
		
		
		// print the actual tip names involved in the clades
        treeSet.reset();
        Tree tree = treeSet.next();
        
        List<String>[] tipNames = new List[clades.size()];        
		for (int i = 0; i < clades.size(); i++) {
			if (ESS[i] < options.ESSCutOff) {
				System.out.println("\nClade" + i + ": ");
				// reverse from this to get the bits
//		        public String toString() {
//		            return "clade " + bits.toString() + " #" + count;
//		        }
		        String cladeStr = clades.get(i).toString();
		        cladeStr = cladeStr.substring(cladeStr.indexOf("{")+1, cladeStr.indexOf("}"));
		        // replace all " "
		        cladeStr = cladeStr.replaceAll(" ", "");
//		        System.out.println(cladeStr);
		        // get all integer in cladeStr
		        tipNames[i] = new ArrayList<>();
		        
		        String[] bits = cladeStr.split(",");
		        for (String bit : bits) {		        	
                	int bitInt = Integer.parseInt(bit);
                	System.out.print(tree.getNode(bitInt/2).getID() + " ");
                	tipNames[i].add(tree.getNode(bitInt/2).getID());
		        }
			}
		}
		
		// for each combination o
		for (int i = 0; i < clades.size(); i++) {
			if (ESS[i] < options.ESSCutOff) {
				for (int j = i+1; j < clades.size(); j++) {
					if (ESS[j] < options.ESSCutOff && i != j) {	
						if (nonOverlap[i][j]) {
							// get the differences in the tipNames
							List<String> diff = new ArrayList<>(tipNames[i]);
							diff.removeAll(tipNames[j]);
							List<String> diff2 = new ArrayList<>(tipNames[j]);
							diff2.removeAll(tipNames[i]);
							System.out.println("\n\nClade" + i + " and Clade " + j + " do not overlap");
//							System.out.println("Clade" + i + " has the following tips: ");
							System.out.println(diff);
//							System.out.println("Clade" + j + " has the following tips: ");
							System.out.println(diff2);
							
						}						
					}
				}
			}
		}

		// among the clades with the lowest ESS, get the average difference between them being on or off in the likelihoods
		double[] diffs = new double[clades.size()];
		for (int i = 0; i < clades.size(); i++) {
			if (ESS[i] < options.ESSCutOff) {
				List<Double> with = new ArrayList<>();
				List<Double> without = new ArrayList<>();
				for (int j = 0; j < totalTreesUsed; j++) {
					if (cladePresence[i][j]) {
						with.add(likelihoods.get(j));
					} else {
						without.add(likelihoods.get(j));
					}
				}
				double avgWith = with.stream().mapToDouble(a -> a).average().getAsDouble();
				double avgWithout = without.stream().mapToDouble(a -> a).average().getAsDouble();
				diffs[i] = Math.abs(avgWith - avgWithout);
			}
		}
		
		// print the clades by the largest difference
		List<Double> corrCoeffs = new ArrayList<>();
		for (int i = 0; i < diffs.length; i++) {
//			if (ESS[i] < options.ESSCutOff) {
				// calculate the correlation between the likelihoods and the presence of the clade
				corrCoeffs.add(Math.abs(correlation(cladePresence[i], likelihoods)));
	//			System.out.println("Clade" + i + " has a difference of " + diffs[i] + " and a correlation of " + corrCoeff);
//			}
		}
		// print the 10 clades with the highest correlation
		List<Clade> maxClades = new ArrayList<>();
		List<Integer> cladeNumbers = new ArrayList<>();
		for (int i = 0; i < clades.size(); i++) {
			int maxIndex = corrCoeffs.indexOf(Collections.max(corrCoeffs));
//			System.out.println("Clade" + maxIndex + " has a correlation of " + corrCoeffs.get(maxIndex));
			maxClades.add(clades.get(maxIndex));
			cladeNumbers.add(maxIndex);
			corrCoeffs.set(maxIndex, Double.MIN_VALUE);		
		}
		
		// print the clades with the highest correlation
		List<String>[] cladeTipNames = new List[maxClades.size()];
		for (int i = 0; i < maxClades.size(); i++) {
			if (ESS[cladeNumbers.get(i)]<options.ESSCutOff) {
				System.out.println("\nClade" + cladeNumbers.get(i) + ": " + correlation(cladePresence[cladeNumbers.get(i)], likelihoods) + " " + ESS[cladeNumbers.get(i)]);
		        String cladeStr = maxClades.get(i).toString();
		        cladeStr = cladeStr.substring(cladeStr.indexOf("{")+1, cladeStr.indexOf("}"));
		        cladeStr = cladeStr.replaceAll(" ", "");
		        // get all integer in cladeStr
		        cladeTipNames[i] = new ArrayList<>();
		        
		        String[] bits = cladeStr.split(",");
		        for (String bit : bits) {		        	
	            	int bitInt = Integer.parseInt(bit);
	            	System.out.print(tree.getNode(bitInt/2).getID() + " ");
	            	System.out.print(bitInt/2 + " ");
//	            	cladeTipNames[i].add(tree.getNode(bitInt/2).getID());
	            	cladeTipNames[i].add(bitInt/2 +"");
		        }	
			}
		}
		System.out.println();

		List<String> uniqueTipNames = new ArrayList<>();
		for (int i = 0; i < cladeTipNames.length; i++) {
			if (cladeTipNames[i] != null) {
	            for (String tipName : cladeTipNames[i]) {
	                if (!uniqueTipNames.contains(tipName)) {
	                    uniqueTipNames.add(tipName);
	                }
	            }
			}
        }
				
		List<String> uniqueTipNamesWithClades = new ArrayList<>();
		// for each member in unqiqueTipNames, make a new string with[&cladex=1,cladey=1...] if it's present in clade
		for (String tipName : uniqueTipNames) {
            StringBuilder b = new StringBuilder(tipName);
            b.append("[&show=1,");
            for (int i = 0; i < cladeTipNames.length; i++) {
            	if (cladeTipNames[i] != null) {
	                if (cladeTipNames[i].contains(tipName)) {
	                    b.append("clade" + cladeNumbers.get(i) + "=1,");
	                }
            	}
            }
            b.deleteCharAt(b.length()-1);
//            b.append("]");
            uniqueTipNamesWithClades.add(b.toString());
        }
		// compute how often each uniqueTipName is present in cladeTipNames
		int[] counts = new int[uniqueTipNames.size()];
		for (int i = 0; i < uniqueTipNames.size(); i++) {
			for (int j = 0; j < cladeTipNames.length; j++) {
				if (cladeTipNames[j] != null) {
					if (cladeTipNames[j].contains(uniqueTipNames.get(i))) {
						counts[i]++;
					}
				}
			}
		}
//		System.out.println("\n\n\n");
//		// print the 10 most common tip names
//		for (int i = 0; i < 10; i++) {
//			int maxIndex = 0;
//			for (int j = 0; j < counts.length; j++) {
//				if (counts[j] > counts[maxIndex]) {
//					maxIndex = j;
//				}
//			}
//			System.out.println(uniqueTipNames.get(maxIndex) + " " + counts[maxIndex]);
//			counts[maxIndex] = 0;
//		}
//		
//		
//		System.out.println(uniqueTipNames);
//		System.out.println(uniqueTipNamesWithClades);
//		System.exit(0);
		
		// read options.inFile line by line and open a second file options.inFile + "_annotated.trees" and write the newick tree with the annotations
		
		PrintStream ps2 = new PrintStream(new File(options.inFile + "_annotated.trees"));
        treeSet.reset();
        while (treeSet.hasNext()) {
            String treeStr = treeSet.next().toString();
//            treeStr = treeStr.replaceAll("\\[&dirty=3\\]", "\\[&show=0]");
            treeStr = treeStr.replaceAll("\\]\\[&dirty=3\\]", "\\]");
            treeStr = treeStr.replaceAll("&sum", "&show=0,sum");
			for (int i = 0; i < uniqueTipNames.size(); i++) {
				treeStr = treeStr.replace(","+uniqueTipNames.get(i) +"[&show=0", ","+uniqueTipNamesWithClades.get(i));
				treeStr = treeStr.replace("("+uniqueTipNames.get(i) +"[&show=0", "("+uniqueTipNamesWithClades.get(i));
//				treeStr = treeStr.replace(","+uniqueTipNames.get(i) +"[&show=0]", ","+uniqueTipNamesWithClades.get(i));
//				treeStr = treeStr.replace("("+uniqueTipNames.get(i) +"[&show=0]", "("+uniqueTipNamesWithClades.get(i));
			}
			ps2.println(treeStr +";");

//            System.out.println(treeStr);
//            tree.getLeafNodeCount();
//        	cladeSystem.add(tree, false);
//            totalTreesUsed++;
        }

		
		
//		while ((line = br2.readLine()) != null) {
//			if (line.trim().length() > 0) {
//				for (int i =0; i < uniqueTipNames.size(); i++) {
//					line = line.replace("(" + uniqueTipNames.get(i) +":", "("+uniqueTipNamesWithClades.get(i)+":");
//					line = line.replace("," + uniqueTipNames.get(i) +":", ","+uniqueTipNamesWithClades.get(i)+":");
//				}
//				ps2.println(line);
//			}
//		}
//		br2.close();
		ps2.close();
		
    }    
    
    
    public static String helpMessage =
			"CladeAnnotator v" + "\n"
					+ "Usage: java -jar cladeannotator.jar [options] <input file> [<output file>]\n" + "Options:\n"
					+ "  -help               Print this message and exit.\n"
					+ "  -burnin <percent>   Percentage of trees to ignore as burnin (default 0).\n"
					+ "  -lowMemory          Use low memory mode (default false).\n"
					+ "  -cladeCutoff <cutoff>  Clade credibility cutoff (default 0.1).\n";
					
    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve ACGAnnotator options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, CladeAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Integer.parseInt(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;
                case "-lowMemory":
					if (args.length <= i + 1)
						printUsageAndError("-lowMemory must be followed by a boolean value (true or false)");
					try {
						options.lowMemory = Boolean.parseBoolean(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing lowMemory value.");
					}
                    
                    break;
                case "-cladeCutoff":
					if (args.length <= i + 1)
						printUsageAndError("-cladeCutoff must be followed by a number (percent)");

					try {
						options.cladeCutOff = Double.parseDouble(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing clade cutoff.");
					}

					if (options.cladeCutOff < 0 || options.cladeCutOff > 1) {
						printUsageAndError("Clade cutoff must be >= 0 and <= 1.");
					}

					i += 1;
					break;   
                case "-ESSCutoff":
                    if (args.length <= i + 1)
                    	printUsageAndError("-ESSCutoff must be followed by a number (percent)");
                    
                    try {
                    	options.ESSCutOff = Double.parseDouble(args[i + 1]);
                    } catch (NumberFormatException e) {
                    	printUsageAndError("Error parsing ESS cutoff.");                    
                    }
					i += 1;
					break;   

				case "-log":
					if (args.length <= i + 1)
						printUsageAndError("-log must be followed by a filename");
					options.logFile = args[i + 1];
					i += 1;
					break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = args[i];

        if (i+1<args.length)
            options.outFile = args[i+1];
    }
    
    
    
    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	CladeAnnotatorOptions options = new CladeAnnotatorOptions();

//        if (args.length == 0) {
//            // Retrieve options from GUI:
//
//            try {
//                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
//            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
//                Log.warning.println("Error setting cross-platform look and feel.");
//            }
//
//            try {
//                SwingUtilities.invokeAndWait(() -> {
//                    if (!getOptionsGUI(options))
//                        System.exit(0);
//
//                    setupGUIOutput();
//                });
//            } catch (InterruptedException | InvocationTargetException e) {
//                e.printStackTrace();
//            }
//
//
//        } else {
            getCLIOptions(args, options);
//        }

        // Run ACGAnnotator
        try {
            new CladeAnnotator(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
    
    
    
    
    
    
    
    
    
    
    public NodeHeightSettingService nodeHeightSettingService;
    public TopologySettingService topologySettingService;
    private int burninPercentage;
    // arguments that do not set any input option
	public Input<List<File>> filesInput = new  Input<> ("file", "Specify the input filename and (optional) output file name" , new ArrayList<>());

    public abstract class TreeSet {
    	public abstract boolean hasNext();
    	public abstract Tree next() throws IOException;
    	public abstract void reset() throws IOException;


    	public String inputFileName;
        public int burninCount = 0;
        public int totalTrees = 0;
        public boolean isNexus = true;

        /** determine number of trees in the file,
    	 * and number of trees to skip as burnin
    	 * @throws IOException
    	 * @throws FileNotFoundException **/
    	void countTrees(int burninPercentage) throws IOException  {
            BufferedReader fin = new BufferedReader(new FileReader(new File(inputFileName)));
            if (!fin.ready()) {
            	throw new IOException("File appears empty");
            }
        	String str = fin.readLine();
            if (!str.toUpperCase().trim().startsWith("#NEXUS")) {
            	// the file contains a list of Newick trees instead of a list in Nexus format
            	isNexus = false;
            	if (str.trim().length() > 0) {
            		totalTrees = 1;
            	}
            }
            while (fin.ready()) {
            	str = fin.readLine();
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

            System.out.println("Processing " + (totalTrees - burninCount) + " trees from file" +
                    (burninPercentage > 0 ? " after ignoring first " + burninPercentage + "% = " + burninCount + " trees." : "."));
		}

    }    
    
    public class FastTreeSet extends TreeSet {
    	int current = 0;
    	Tree [] trees;

    	public FastTreeSet(String inputFileName, int burninPercentage) throws IOException  {
            this.inputFileName = inputFileName;
            countTrees(burninPercentage);

            List<Tree> parsedTrees;
            if (isNexus) {
                NexusParser nexusParser = new NexusParser();
                nexusParser.parseFile(new File(inputFileName));
                parsedTrees = nexusParser.trees;
            } else {
                BufferedReader fin = new BufferedReader(new FileReader(inputFileName));
                parsedTrees = new ArrayList<>();
                current = 0;
                while (fin.ready()) {
                    String line = fin.readLine().trim();

                    String id = "" + current++;
                    try {
                    	int i = line.indexOf("(");
    	                id = line.substring(5, i).split("=")[0].trim();
                    } catch (Exception e) {
                    	// ignore
                    }

                    Tree thisTree;
                    try {
                        thisTree = new TreeParser(null, line, 0, false);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        thisTree = new TreeParser(null, line, 1, false);
                    }
                    thisTree.setID(id);
                    
                    parsedTrees.add(thisTree);
                }
                fin.close();
                current = 0;
            }

            int treesToUse = parsedTrees.size() - burninCount;
	      	trees = new Tree[treesToUse];
            for (int i=burninCount; i<parsedTrees.size(); i++)
                trees[i-burninCount] = parsedTrees.get(i);
		}

		@Override
		public boolean hasNext() {
			return current < trees.length;
		}

		@Override
		public Tree next()  {
			return trees[current++];
		}

		@Override
		public void reset()  {
			current = 0;
		}
    }
    
    public class MemoryFriendlyTreeSet extends TreeSet {
//    	Tree [] trees;
    	int current = 0;
    	int lineNr;
        public Map<String, String> translationMap = null;
        public List<String> taxa;

        // label count origin for NEXUS trees
        int origin = -1;

        BufferedReader fin;

        public MemoryFriendlyTreeSet(String inputFileName, int burninPercentage) throws IOException  {
    		this.inputFileName = inputFileName;
    		countTrees(burninPercentage);

            fin = new BufferedReader(new FileReader(inputFileName));
    	}


    	@Override
    	public void reset() throws FileNotFoundException  {
    		current = 0;
            fin = new BufferedReader(new FileReader(new File(inputFileName)));
            lineNr = 0;
            try {
                if (isNexus) {
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
                } else {
                    while (fin.ready() && lineNr < burninCount) {
                        final String str = nextLine();
                        if (str == null) {
                            return;
                        }
                        if (str.trim().length() > 2 && !str.trim().startsWith("#")) {
                        	lineNr++;
                        }
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
            if (str.matches("^\\s*\\[.*")) {
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
        	fin.mark(1024*1024);
        	int lineNr = this.lineNr;
            String str = readLine().trim();
            while (str.equals("")) {
            	fin.mark(1024*1024);
            	lineNr = this.lineNr;
                str = readLine().trim();
            }

            // if first non-empty line is "translate" then parse translate block
            if (str.toLowerCase().contains("translate")) {
                translationMap = parseTranslateBlock();
                origin = getIndexedTranslationMapOrigin(translationMap);
                if (origin != -1) {
                    taxa = getIndexedTranslationMap(translationMap, origin);
                }
            } else {
            	this.lineNr = lineNr;
            	fin.reset();
            }
            // we got to the end of the translate block
            // read burninCount trees
            current = 0;
            while (current < burninCount && fin.ready()) {
    			str = nextLine();
                if (str.trim().toLowerCase().startsWith("tree ")) {
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
                if (taxa == null) {
                	collectTaxaNames(str);
                }
            	current++;

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
            if (str.trim().toLowerCase().startsWith("tree ")) {
            	current++;
                final int i = str.indexOf('(');

                String id = "" + current;
                try {
	                 id = str.substring(5, i).split("=")[0].trim();
                } catch (Exception e) {
                	// ignore
                }
                
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

                //if (translationMap != null) treeParser.translateLeafIds(translationMap);
                treeParser.setID(id);
                
                return treeParser;
            }
    		return null;
    	}
    	
    	private void collectTaxaNames(String str) {
    		taxa = new ArrayList<>();
    		int i = 0;
    		while (i < str.length()) {
    			char c = str.charAt(i);
    			switch (c) {
    			case '(':
    			case ')':
    			case ',':
    				// ignore
    				i++;
    				break;
    			case '[':
    				// eat up meta data
    				while (i < str.length() && str.charAt(i) != ']') {
    					i++;
    				}
    				break;
    			case ':':
    				// eat up length
    				while (i < str.length() && !(str.charAt(i) == ')'|| str.charAt(i) == ',')) {
    					i++;
    				}
    				break;
    			default:
    				StringBuilder b = new StringBuilder();
    				boolean done = false;
    				while (i < str.length() && !done) {
    					c = str.charAt(i);
    					done =  c == ')' || c == ':' || c == ',' || c == '(' || c == '[';
    					if (!done) {
    						if (c != '\'' && c != '"') {
    							b.append(c);
    						}
    						i++;
    					} else {
    						taxa.add(b.toString());
    					}
    				}
    				
    			}
    		}		
    	}
    }
    TreeSet treeSet;


//    public void run(final int burninPercentage,
//    					 boolean lowMemory, // allowSingleChild was defunct (always set to false), now replaced by flag to say how much 
//                         // HeightsSummary heightsOption,
//                         double posteriorLimit,
//                         double hpd2D,
//                         // Target targetOption,
//                         String targetTreeFileName,
//                         String inputFileName,
//                         String outputFileName
//    ) throws IOException  {
//
//        topologySettingService = getTopologySettingService();
//        nodeHeightSettingService = getNodeHeightSettingService();
//
//        this.posteriorLimit = posteriorLimit;
//        this.hpd2D = hpd2D;
//
//        attributeNames.add("height");
//        attributeNames.add("length");
//
//        totalTrees = 10000;
//        totalTreesUsed = 0;
//
//        try {
//        	if (lowMemory) {
//        		treeSet = new MemoryFriendlyTreeSet(inputFileName, burninPercentage);
//        	} else {
//        		treeSet = new FastTreeSet(inputFileName, burninPercentage);
//        	}
//        } catch (Exception e) {
//        	e.printStackTrace();
//        	Log.err.println("Error Parsing Input Tree: " + e.getMessage());
//        	return;
//        }
//        
//        if (!topologySettingService.getServiceName().equals(UserTargetTreeTopologyService.SERVICE_NAME)) {
//            // even when a user specified target tree is provided we still need to count the totalTreesUsed for subsequent steps.
//            treeSet.reset();
//            while (treeSet.hasNext()) {
//                Tree tree = treeSet.next();
//                tree.getLeafNodeCount();
//                if (tree.getDirectAncestorNodeCount() > 0 && !SAmode && processSA) {
//                    SAmode = true;
//                    Log.err.println("A tree with a sampled ancestor is found. Turning on\n the sampled ancestor " +
//                            "summary analysis.");
//                    if (nodeHeightSettingService.getServiceName().equals("CA")) {
//                        throw new RuntimeException("The common ancestor height is not \n available for trees with sampled " +
//                                "ancestors. Please choose \n another height summary option");
//                    }
//                }
//                totalTreesUsed++;
//            }
//        }
//
//        Tree targetTree = topologySettingService.setTopology(treeSet, progressStream, this);
//
//     
//        cladeSystem = getCladeSystem(targetTree);
//        
//
//        progressStream.println("Annotating target tree...");
//
//        try {
//            annotateTree(cladeSystem, targetTree.getRoot(), null);//, heightsOption);
//
//            nodeHeightSettingService.setNodeHeights(targetTree, progressStream, this);
////            if( heightsOption == HeightsSummary.CA_HEIGHTS ) {
////                setTreeHeightsByCA(targetTree, targetOption);
////            }
//        } catch (Exception e) {
//        	e.printStackTrace();
//            Log.err.println("Error to annotate tree: " + e.getMessage() + "\nPlease check the tree log file format.");
//            return;
//        }
//
//        progressStream.println("Writing annotated tree....");
//
//        
//        processMetaData(targetTree.getRoot());
//        try {
//            final PrintStream stream = outputFileName != null ?
//                    new PrintStream(new FileOutputStream(outputFileName)) :
//                    System.out;
//            targetTree.init(stream);
//            stream.println();
//            
//            stream.print("tree TREE_" + 
//            		topologySettingService.getServiceName() + "_" + 
//            		nodeHeightSettingService.getServiceName() + " = ");
//            int[] dummy = new int[1];
//            String newick = targetTree.getRoot().toSortedNewick(dummy, true);
//            stream.print(newick);
//            stream.println(";");
////            stream.println(targetTree.getRoot().toShortNewick(false));
////            stream.println();
//            targetTree.close(stream);
//            stream.println();
//        } catch (Exception e) {
//            Log.err.println("Error to write annotated tree file: " + e.getMessage());
//            return;
//        }
//
//    }
//
    private void processMetaData(Node node) {
		for (Node child : node.getChildren()) {
			processMetaData(child);
		}
		Set<String> metaDataNames = node.getMetaDataNames(); 
		if (metaDataNames != null && !metaDataNames.isEmpty()) {
			String metadata = "";
			for (String name : metaDataNames) {
				Object value = node.getMetaData(name);
				metadata += name + "=";
				if (value instanceof Object[]) {
					Object [] values = (Object[]) value;
					metadata += "{";
					for (int i = 0; i < values.length; i++) {
						metadata += values[i].toString();
						if (i < values.length - 1) {
							metadata += ",";
						}
					}
					metadata += "}";
				} else {
					 metadata += value.toString();
				}
				metadata += ",";
			}
			metadata = metadata.substring(0, metadata.length() - 1);
			node.metaDataString = metadata;
		}		
	}


    public static void annotateMeanAttribute(Node node, String label, double[] values) {
        double mean = DiscreteStatistics.mean(values);
        node.setMetaData(label, mean);
    }

    public static void annotateMedianAttribute(Node node, String label, double[] values) {
        double median = DiscreteStatistics.median(values);
        node.setMetaData(label, median);

    }

    public static void annotateModeAttribute(Node node, String label, HashMap<Object, Integer> values) {
        Object mode = null;
        int maxCount = 0;
        int totalCount = 0;
        int countInMode = 1;

        for (Object key : values.keySet()) {
            int thisCount = values.get(key);
            if (thisCount == maxCount) {
                // I hope this is the intention
                mode = mode.toString().concat("+" + key);
                countInMode++;
            } else if (thisCount > maxCount) {
                mode = key;
                maxCount = thisCount;
                countInMode = 1;
            }
            totalCount += thisCount;
        }
        double freq = (double) maxCount / (double) totalCount * countInMode;
        node.setMetaData(label, mode);
        node.setMetaData(label + ".prob", freq);
    }

    public static void annotateFrequencyAttribute(Node node, String label, HashMap<Object, Integer> values) {
        double totalCount = 0;
        Set<?> keySet = values.keySet();
        int length = keySet.size();
        String[] name = new String[length];
        Double[] freq = new Double[length];
        int index = 0;
        for (Object key : values.keySet()) {
            name[index] = key.toString();
            freq[index] = Double.valueOf(values.get(key));
            totalCount += freq[index];
            index++;
        }
        for (int i = 0; i < length; i++)
            freq[i] /= totalCount;

        node.setMetaData(label + ".set", name);
        node.setMetaData(label + ".set.prob", freq);
    }

    public static void annotateRangeAttribute(Node node, String label, double[] values) {
        double min = DiscreteStatistics.min(values);
        double max = DiscreteStatistics.max(values);
        node.setMetaData(label, new Object[]{min, max});
    }

    public static void annotateHPDAttribute(Node node, String label, double hpd, double[] values) {
        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(hpd * values.length);
        for (int i = 0; i <= (values.length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = values[indices[hpdIndex]];
        double upper = values[indices[hpdIndex + diff - 1]];
        node.setMetaData(label, new Object[]{lower, upper});
    }



    public static final String CORDINATE = "cordinates";

    private String formattedLocation(double x) {
        return String.format("%5.2f", x);
    }


    int totalTrees = 0;
    int totalTreesUsed = 0;
    double posteriorLimit = 0.0;
    double hpd2D = 0.80;

// 
//
//	public CladeSystem getCladeSystem(Tree targetTree) {
//    	if (cladeSystem != null) {
//    		return cladeSystem;
//    	}
//        System.out.println("Collecting node information...");
//        System.out.println("0              25             50             75            100");
//        System.out.println("|--------------|--------------|--------------|--------------|");
//
//        int stepSize = Math.max(totalTreesUsed / 60, 1);
//        int reported = 0;
//
//        // this call increments the clade counts and it shouldn't
//        // this is remedied with removeClades call after while loop below
//        cladeSystem = new CladeSystem();
//        cladeSystem.setProcessSA(processSA);
//        cladeSystem.add(targetTree, true);
//        int totalTreesUsedNew = 0;
//        try {
//            int counter = 0;
//            treeSet.reset();
//            while (treeSet.hasNext()) {
//            	Tree tree = treeSet.next();
//            	if (counter == 0) {
//                    setupAttributes(tree);
//            	}
//                cladeSystem.collectAttributes(tree, attributeNames);
//    			while (reported < 61 && 1000.0*reported < 61000.0 * (counter + 1) / this.totalTreesUsed) {
//    				System.out.print("*");
//  		          	reported++;
//  		          System.out.flush();
//    			}
//                totalTreesUsedNew++;
//                counter++;
//        	}
//        	
//            cladeSystem.removeClades(targetTree.getRoot(), true);
//            this.totalTreesUsed = totalTreesUsedNew;
//            cladeSystem.calculateCladeCredibilities(totalTreesUsedNew);
//        } catch (Exception e) {
//            Log.err.println("Error Parsing Input Tree: " + e.getMessage());
//            return null;
//        }
//        System.out.println();
//        System.out.println();
//        return cladeSystem;
//    }
//    
//    public int getTotalTreesUsed() {return totalTreesUsed;}
//    public TreeSet getTreeSet() {return treeSet;}
//    public boolean isProcessSA() {return processSA;}
//    public int getTotalTrees() {return totalTrees;}
//    public int getBurninCount() {return burnInPercentageInput.get();}

    public static double essTracerStyle(boolean[] x) {
        // 1. Convert booleans to doubles (1.0 if true, 0.0 if false)
        int N = x.length;
        double[] data = new double[N];
        for (int i = 0; i < N; i++) {
            data[i] = x[i] ? 1.0 : 0.0;
        }

        // 2. Compute the mean
        double sum = 0.0;
        for (double val : data) {
            sum += val;
        }
        double mean = sum / N;

        // 3. Demean the data
        double[] xDemeaned = new double[N];
        for (int i = 0; i < N; i++) {
            xDemeaned[i] = data[i] - mean;
        }

        // 4. Decide on a maximum lag (similar to the R code)
        int maxLag = Math.min(N - 1, 10000);

        // 5. Compute variance of xDemeaned
        double sumSq = 0.0;
        for (double val : xDemeaned) {
            sumSq += (val * val);
        }
        double var = sumSq / N;
        if (var == 0.0) {
            // If all values are the same, variance is zero => no autocorrelation => ESS = N
            return N;
        }

        // 6. Compute autocorrelations for lags 1..maxLag
        //    Then sum the positive autocorrelations
        double posSum = 0.0;
        for (int lag = 1; lag <= maxLag; lag++) {
            double covSum = 0.0;
            // Naive approach to covariance at this lag
            for (int i = lag; i < N; i++) {
                covSum += xDemeaned[i] * xDemeaned[i - lag];
            }
            // Average covariance at this lag
            double cov = covSum / (N - lag);
            // Correlation = covariance / variance
            double r = cov / var;

            if (r <= 0.0) {
                break;
            }
            posSum += r;
        }

        // 7. Compute ESS
        double ess = N / (1.0 + 2.0 * posSum);
        return ess;
    }
    
    public static double correlation(boolean[] x, List<Double> y) {
        if (x.length != y.size()) {
            throw new IllegalArgumentException("Lengths of boolean[] and List<Double> must match.");
        }
        
        int n = x.length;
        if (n < 2) {
            throw new IllegalArgumentException("Need at least two data points to compute correlation.");
        }
        
        // 1. Compute means
        double sumX = 0.0;
        double sumY = 0.0;
        for (int i = 0; i < n; i++) {
            double xVal = x[i] ? 1.0 : 0.0;
            double yVal = y.get(i);
            sumX += xVal;
            sumY += yVal;
        }
        double meanX = sumX / n;
        double meanY = sumY / n;
        
        // 2. Accumulate covariances and variances
        double covXY = 0.0;
        double varX = 0.0;
        double varY = 0.0;
        for (int i = 0; i < n; i++) {
            double xVal = x[i] ? 1.0 : 0.0;
            double yVal = y.get(i);
            double dx = xVal - meanX;
            double dy = yVal - meanY;
            covXY += dx * dy;
            varX += dx * dx;
            varY += dy * dy;
        }
        
        // 3. If one variable has zero variance, correlation is undefined (NaN) or could be 0.0
        if (varX == 0.0 || varY == 0.0) {
            // Typically, correlation is undefined if there is no variation in X or Y.
            // You could throw an exception, return 0.0, or return Double.NaN.
            return Double.NaN;
        }
        
        // 4. Compute the sample correlation coefficient
        // For sample correlation we divide by (n - 1) in both numerator and denominator,
        // which cancels out, so equivalently:
        double correlation = covXY / Math.sqrt(varX * varY);
        
        return correlation;
    }


}

