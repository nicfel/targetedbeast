package targetedbeast.alignment;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.inference.parameter.IntegerParameter;



@Description("Alignment based on a filter operation on another alignment")
public class ConsensusAlignment extends Alignment {
    final public Input<Alignment> alignmentInput = new Input<>("data", "alignment to be filtered", Validate.REQUIRED);

    public ConsensusAlignment() {
        sequenceInput.setRule(Validate.OPTIONAL);
        // it does not make sense to set weights on sites, since they can be scrambled by the filter
        siteWeightsInput.setRule(Validate.FORBIDDEN);
    }

    @Override
    public void initAndValidate() {
        Alignment data = alignmentInput.get();
        m_dataType = data.getDataType();

        
   		counts = data.getCounts();
        taxaNames = data.getTaxaNames();
        stateCounts = data.getStateCounts();

        calcPatterns();
        setupAscertainment();
        
//        // print out all the data to sysout by taxa
//		System.out.print("[");
//		for (int i = 0; i < taxaNames.size(); i++) {
//			for (int j = 0; j < getPatternCount()-1; j++) {
//				System.out.print(getPattern(i, j) + ",");
//			}
//			System.out.print(getPattern(i, getPatternCount()-1) + ";\n");		
//			
//		}
//		System.out.println("]");
//		System.exit(0);
        
        // get all   sites   public int getPatternIndex(int site)    for which pattern is 0 and print
        for (int i = 0; i < getTaxonCount(); i++) {
        	System.out.println(getPattern(i, 65) + " " + taxaNames.get(i));
        }
		for (int i = 0; i < getSiteCount(); i++) {
			if (getPatternIndex(i) == 65) {
				System.out.println(i + " " + getPatternIndex(i));
			}
		}
        // print out the data for the 
        
//        System.exit(0);

        
        
	
    }
    
    @Override
    public List<List<Integer>> getCounts() {
    	if (counts == null) {
			counts = new ArrayList<>();
			for (int j = 0; j < sitePatterns[0].length; j++) {
				counts.add(new ArrayList<>());
			}
			for (int i = 0; i < getSiteCount(); i++) {
				int [] sites = getPattern(getPatternIndex(i));
    			for (int j = 0; j < getTaxonCount(); j++) {
    				counts.get(j).add(sites[j]);
    			}
			}
    	}
    	return counts;
    }
    
    @Override
    protected void calcPatterns() {
        int nrOfTaxa = counts.size();
        int nrOfSites = alignmentInput.get().getSiteCount();
                
        // convert data to transposed int array
        int[][] data = new int[nrOfSites][nrOfTaxa];
        String missingChar = Character.toString(DataType.MISSING_CHAR);
        String gapChar = Character.toString(DataType.GAP_CHAR);
        for (int i = 0; i < nrOfTaxa; i++) {
            List<Integer> sites = counts.get(i);
            for (int j = 0; j < nrOfSites; j++) {
                data[j][i] = sites.get(j);
            }
        }
        // make a new data where all constant sites are removed
        List<Integer> indexAfterRearrange = new ArrayList<>();
        for (int j = 0; j < nrOfSites; j++) {
			int[] noInState = new int[m_dataType.getStateCount()];
			for (int i = 0; i < nrOfTaxa; i++) {
				if (data[j][i] < m_dataType.getStateCount()) {
					noInState[data[j][i]]++;
				}
			}
			// if more than one noInState is larger than 1, then keep the state
			int largerZero = 0;
			for (int i = 0; i < m_dataType.getStateCount(); i++) {
				if (noInState[i] > 0) {
					largerZero++;
				}
			}
			if (largerZero > 1) {
				indexAfterRearrange.add(j);
			}
		}
        
        int[][] data2 = new int[indexAfterRearrange.size()][nrOfTaxa];
		for (int j = 0; j < indexAfterRearrange.size(); j++) {
			for (int i = 0; i < nrOfTaxa; i++) {
				data2[j][i] = data[indexAfterRearrange.get(j)][i];
			}
		}
        
        // sort data
        SiteComparator comparator = new SiteComparator();
        Arrays.sort(data2, comparator);
        

        // count patterns in sorted data
        int[] weights = new int[nrOfSites];
        int nrOfPatterns = 1;
        if (nrOfSites > 0) {
	        weights[0] = 1;
	        for (int i = 1; i < data2.length; i++) {
	        	// check if this patterns has variable sites that are not ambiguous
	        	
	            if (comparator.compare(data2[i - 1], data2[i]) != 0) {
	                nrOfPatterns++;
	                data2[nrOfPatterns - 1] = data2[i];
	            }
	            weights[nrOfPatterns - 1]++;
	        }
        } else {
            nrOfPatterns = 0;
        }
                        
        // reserve memory for patterns
        patternWeight = new int[nrOfPatterns];
        sitePatterns = new int[nrOfPatterns][nrOfTaxa];
        for (int i = 0; i < nrOfPatterns; i++) {
            patternWeight[i] = weights[i];
            sitePatterns[i] = data2[i];
        }

        // find patterns for the sites
        patternIndex = new int[nrOfSites];
        for (int i = 0; i < nrOfSites; i++) {
            int[] sites = new int[nrOfTaxa];
            for (int j = 0; j < nrOfTaxa; j++) {
                sites[j] = counts.get(j).get(i);
            }
            patternIndex[i] = Arrays.binarySearch(sitePatterns, sites, comparator);
        }

        if (siteWeights != null) {
        	// TODO: fill in weights with siteweights.
        	throw new RuntimeException("Cannot handle site weights in FilteredAlignment. Remove \"weights\" from data input.");
        }

        // determine maximum state count
        // Usually, the state count is equal for all sites,
        // though for SnAP analysis, this is typically not the case.
        maxStateCount = 0;
        for (int stateCount1 : stateCounts) {
            maxStateCount = Math.max(maxStateCount, stateCount1);
        }
        Log.info.println(getPatternCount() + " patterns");
        
        // counts are not valid any more -- better set to null in case
        // someone gets bitten by this.
        this.counts = null;
    }
    
}
