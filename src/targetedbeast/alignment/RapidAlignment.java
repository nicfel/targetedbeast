/*
* File Alignment.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package targetedbeast.alignment;

import java.util.*;


import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.pkgmgmt.PackageManager;

@Description("Adapted version of alignment where the first 4 bases are always 0, 1, 2, 3, etc. and the rest are the actual data. This is used for rapid calculations.")
public class RapidAlignment extends Alignment {

	public RapidAlignment() {
    	if (types.size() == 0) {
    		findDataTypes();    		
    	}
    }

    /**
     * Constructor for testing purposes.
     *
     * @param sequences
     * @param stateCount
     * @param dataType
     * @deprecated This is the deprecated legacy form and will be removed
     * at some point. Use {@link #Alignment(List, String)} instead.
     */
    @Deprecated
    public RapidAlignment(List<Sequence> sequences, Integer stateCount, String dataType) {
        this(sequences, dataType);
    }

    /**
     * Constructor for testing purposes.
     *
     * @param sequences
     * @param dataType
     */
    public RapidAlignment(List<Sequence> sequences, String dataType) {
        for (Sequence sequence : sequences) {
            sequenceInput.setValue(sequence, this);
        }
        dataTypeInput.setValue(dataType, this);
    	if (types.size() == 0) {
    		findDataTypes();    		
    	}
        initAndValidate();
    }

    @Override
    public void initAndValidate() {

        if (sequenceInput.get().size() == 0 && defaultInput.get().size() == 0) {
            throw new IllegalArgumentException("Either a sequence input must be specified, or a map of strings must be specified");
        }

        if (siteWeightsInput.get() != null) {
            String str = siteWeightsInput.get().trim();
            String[] strs = str.split(",");
            siteWeights = new int[strs.length];
            for (int i = 0; i < strs.length; i++) {
                siteWeights[i] = Integer.parseInt(strs[i].trim());
            }
        }

        // determine data type, either user defined or one of the standard ones
        if (userDataTypeInput.get() != null) {
            m_dataType = userDataTypeInput.get();
        } else {
            initDataType();
        }

        // initialize the sequence list
        if (sequenceInput.get().size() > 0) {
            sequences = sequenceInput.get();
        } else {
            // alignment defined by a map of id -> sequence
            List<String> taxa = new ArrayList<>();
            taxa.addAll(map.keySet());
            sequences.clear();
            for (String key : taxa) {
                String sequence = map.get(key);
                sequences.add(new Sequence(key, sequence));
            }
        }

        // initialize the alignment from the given list of sequences
        initializeWithSequenceList(sequences, true);

        if (taxonSetInput.get() != null && taxonSetInput.get().getTaxonCount() > 0) {
            sortByTaxonSet(taxonSetInput.get());
        }
        Log.info.println(toString(false));
    }

    /**
     * Initializes data types using
     * {@link PackageManager#find(Class, String[]) PackageManager.find}
     */
    @Override
    protected void initDataType() {
        if (types.get(dataTypeInput.get()) == null) {
            throw new IllegalArgumentException("data type + '" + dataTypeInput.get() + "' cannot be found. " +
                    "Choose one of " + Arrays.toString(types.keySet().toArray(new String[0])));
        }
        
        m_dataType = types.get(dataTypeInput.get());
    }

    /**
     * Initializes the alignment given the provided list of sequences and no other information.
     * It site weights and/or data type have been previously set up with initAndValidate then they
     * remain in place. This method is used mainly to re-order the sequences to a new taxon order
     * when an analysis of multiple alignments on the same taxa are undertaken.
     *
     * @param sequences
     */
    protected void initializeWithSequenceList(List<Sequence> sequences, boolean log) {
        this.sequences = sequences;
        taxaNames.clear();
        stateCounts.clear();
        counts.clear();
        try {
            for (Sequence seq : sequences) {

                counts.add(seq.getSequence(m_dataType));
                if (taxaNames.contains(seq.getTaxon())) {
                    throw new RuntimeException("Duplicate taxon found in alignment: " + seq.getTaxon());
                }
                taxaNames.add(seq.getTaxon());
                tipLikelihoods.add(seq.getLikelihoods());
                // if seq.isUncertain() == false then the above line adds 'null'
	            // to the list, indicating that this particular sequence has no tip likelihood information
                usingTipLikelihoods |= (seq.getLikelihoods() != null);	            

                if (seq.totalCountInput.get() != null) {
                    stateCounts.add(seq.totalCountInput.get());
                } else {
                    stateCounts.add(m_dataType.getStateCount());
                }
            }
            if (counts.size() == 0) {
                // no sequence data
                throw new RuntimeException("Sequence data expected, but none found");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
        sanityCheckCalcPatternsSetUpAscertainment(log);
    }

    /**
     * SiteComparator is used for ordering the sites,
     * which makes it easy to identify patterns.
     */
    public class SiteComparator implements Comparator<int[]> {
        @Override
		public int compare(int[] o1, int[] o2) {
        	// first check if one or the other has all the same elements
        	boolean allSame1 = true;
        	boolean allSame2 = true;
        	for (int i = 1; i < o1.length; i++) {
				if (o1[i] != o1[0]) {
					allSame1 = false;
					break;
				}
        	}
        	for (int i = 1; i < o2.length; i++) {
                if (o2[i] != o2[0]) {
					allSame2 = false;
					break;
                }
        	}
        	// if both are the same, check the first element
        	if (allSame1 && allSame2) {
				if (o1[0] > o2[0]) {
					return 1;
				}
				if (o1[0] < o2[0]) {
					return -1;
				}
				return 0;
        	}
        	// if one is all the same, it is always less than the other
        	if (allSame1) {
        		return -1;
        	}
        	if (allSame2) {
        		return 1;
        	}
        	
        	
            for (int i = 0; i < o1.length; i++) {
                if (o1[i] > o2[i]) {
                    return 1;
                }
                if (o1[i] < o2[i]) {
                    return -1;
                }
            }
            return 0;
        }
    } // class SiteComparator

    @Override
    protected void calcPatterns() {
        calcPatterns(true);
    }

    /**
     * calculate patterns from sequence data
     * *
     */
    @Override
    protected void calcPatterns(boolean log) {
        int taxonCount = counts.size();
        int siteCount = counts.get(0).size();
        
        // convert data to transposed int array
        int[][] data = new int[siteCount+stateCounts.get(0)][taxonCount];
        for (int i = 0; i < taxonCount; i++) {
            List<Integer> sites = counts.get(i);
            for (int j = 0; j < siteCount; j++) {
                data[j][i] = sites.get(j);
            }
            // add 0 to stateCounts for the first stateCounts
			for (int j = 0; j < stateCounts.get(0); j++) {
				data[j+siteCount][i] = j;
			}
        }
        
        // sort data
        SiteComparator comparator = new SiteComparator();
        Arrays.sort(data, comparator);

        // count patterns in sorted data
        // if (siteWeights != null) the weights are recalculated below
        int patterns = 1;
        int[] weights = new int[data.length];
        weights[0] = 1;
        for (int i = 1; i < data.length; i++) {
        	// check if all bases are ambiguous, in which case we don't want to count them as a pattern
        	boolean allAmbiguous = true;
        	for (int j = 0; j < taxonCount; j++) {
				if (data[i][j] < stateCounts.get(j)) {
					allAmbiguous = false;
					break;
				}
        	}
        	        	
        	if (!allAmbiguous) {
	            if (usingTipLikelihoods || comparator.compare(data[i - 1], data[i]) != 0) {
	            	// In the case where we're using tip probabilities, we need to treat each 
	            	// site as a unique pattern, because it could have a unique probability vector.
	                patterns++;
	                data[patterns - 1] = data[i];
	            }
	            weights[patterns - 1]++;	           
        	}
        }
        
        // subtract 1 for all the first patterns
        for (int i = 0; i < stateCounts.get(0); i++) {
        	weights[i]--;
        }

        // reserve memory for patterns
        patternWeight = new int[patterns];
        sitePatterns = new int[patterns][taxonCount];
        for (int i = 0; i < patterns; i++) {
            patternWeight[i] = weights[i];
            sitePatterns[i] = data[i];
        }
        

        // find patterns for the sites
        patternIndex = new int[siteCount];
        for (int i = 0; i < siteCount; i++) {
            int[] sites = new int[taxonCount];
            for (int j = 0; j < taxonCount; j++) {
                sites[j] = counts.get(j).get(i);
            }
            patternIndex[i] = Arrays.binarySearch(sitePatterns, sites, comparator);
        }

        if (siteWeights != null) {
            Arrays.fill(patternWeight, 0);
            for (int i = 0; i < siteCount; i++) {
                patternWeight[patternIndex[i]] += siteWeights[i];
            }
        }

        // determine maximum state count
        // Usually, the state count is equal for all sites,
        // though for SnAP analysis, this is typically not the case.
        maxStateCount = 0;
        for (int m_nStateCount1 : stateCounts) {
            maxStateCount = Math.max(maxStateCount, m_nStateCount1);
        }
        // report some statistics
        if (log && taxaNames.size() < 30) {
            for (int i = 0; i < taxaNames.size(); i++) {
                Log.info.println(taxaNames.get(i) + ": " + counts.get(i).size() + " " + stateCounts.get(i));
            }
        }

        if (stripInvariantSitesInput.get()) {
            // don't add patterns that are invariant, e.g. all gaps
            if (log) Log.info.println("Stripping invariant sites");

            int removedSites = 0;
            for (int i = 0; i < patterns; i++) {
                int[] pattern = sitePatterns[i];
                int value = pattern[0];
                boolean isInvariant = true;
                for (int k = 1; k < pattern.length; k++) {
                    if (pattern[k] != value) {
                        isInvariant = false;
                        break;
                    }
                }
                if (isInvariant) {
                    removedSites += patternWeight[i];
                    patternWeight[i] = 0;

                    if (log) Log.info.print(" <" + value + "> ");
                }
            }
            if (log) Log.info.println(" removed " + removedSites + " sites ");
        }
    } // calcPatterns

} // class Data
