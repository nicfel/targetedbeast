package targetedbeast.likelihood;

import java.util.Arrays;
import java.util.List;

import beast.base.evolution.likelihood.LikelihoodCore;

/**
 * standard likelihood core, uses no caching *
 */
public class RapidLikelihoodCore extends LikelihoodCore {
	protected int nrOfStates;
	protected int nrOfNodes;
	protected int nrOfPatterns;
	protected int partialsSize;
	protected int matrixSize;
	protected int nrOfMatrices;

	protected boolean integrateCategories;

	protected double[][][] partials;

	// keeps track of the intermediate partials
	protected double[] intermediatePartials1;
	protected double[] intermediatePartials2;

	
	// array to keep track of the partials after multiplications with the matrices
	protected double[][][] edgePartials;

	
	// repeated calculation, can be saved
	int correction;

	protected int[][] states;

	protected double[][][] matrices;

	protected int[] currentMatrixIndex;
	protected int[] storedMatrixIndex;
	protected int[] currentPartialsIndex;
	protected int[] storedPartialsIndex;
	protected int[] currentEdgeIndex;
	protected int[] storedEdgeIndex;
	

	protected boolean useScaling = false;

	protected double[][][] scalingFactors;
	
	double mutcount = 0.0;
	int count = 0;

	private double scalingThreshold = 1.0E-100;
	double SCALE = 2;

	public RapidLikelihoodCore(int nrOfStates) {
		this.nrOfStates = nrOfStates;
	} // c'tor

	protected void calculateStatesPruning(int[] stateIndex, double[] matrices, double[] intermediatePartials) {

		int v;
		
		for (int l = 0; l < nrOfMatrices; l++) {
			v = l * nrOfPatterns * nrOfStates;
			for (int k = 0; k < nrOfStates; k++) {			
				int state1 = stateIndex[k];
				int w = l * matrixSize;
				for (int i = 0; i < nrOfStates; i++) {
					intermediatePartials[v] = matrices[w + state1];	
					v++;
					w+= nrOfStates;
				}
			}
		}
	}

	protected void calculatePartialPruning(double[] partials, double[] matrices, double[] intermediatePartials,
			int[] calcForPatterns) {
		
		double sum;
		int u,v,w;		
		
		for (int l = 0; l < nrOfMatrices; l++) {
			u = l * nrOfPatterns * nrOfStates;
			v=u;
			for (int k = 0; k < nrOfStates; k++) {
				w = l * matrixSize; // pick the correct matrix index
				for (int i = 0; i < nrOfStates; i++) {
					sum = 0.0;
					for (int j = 0; j < nrOfStates; j++) {
						sum += matrices[w] * partials[v + j];
						w++;
					}

					intermediatePartials[u] = sum;
					u++;
				}
				v += nrOfStates;
			}
		}	
		for (int l = 0; l < nrOfMatrices; l++) {
			for (int k = 0; k < nrOfPatterns; k++) {
				if (calcForPatterns[k] == -1) {
					break;
				}
				
				u = l * nrOfPatterns * nrOfStates + (calcForPatterns[k]) * nrOfStates;
				v = u;
				w = 0;
			
				w = l * matrixSize;
				


				for (int i = 0; i < nrOfStates; i++) {
					sum = 0.0;
					for (int j = 0; j < nrOfStates; j++) {
						sum += matrices[w] * partials[v + j];
						w++;
					}
					intermediatePartials[u] = sum;
					u++;
				}
				v += nrOfStates;
			}
		}
		

	}

	protected void updatePartials(double[] partials, 
			double[] intermediatePartials1, double[] intermediatePartials2,
			int[] mutations1, 
			int[] mutations2,
			int[] calcForPatterns) {
		

		// calculate the product for the first nrOfStates
		for (int l = 0; l < nrOfMatrices; l++) {
			int u = l * nrOfPatterns * nrOfStates;
			for (int k = 0; k < nrOfStates; k++) {
				for (int i = 0; i < nrOfStates; i++) {
					partials[u] = intermediatePartials1[u] * intermediatePartials2[u];
					u++;
				}
			}
		}	
		
		// calculate the product of the two intermediate matrices to get the partials at the parent node
		for (int m = 0; m < calcForPatterns.length; m++) {
			if (calcForPatterns[m] == -1) {
				break;
			}
			int mutation1 = mutations1[calcForPatterns[m]];
			int mutation2 = mutations2[calcForPatterns[m]];
			int parentMut = calcForPatterns[m];
			
						
			for (int l = 0; l < nrOfMatrices; l++) {
				int v1 = mutation1 * nrOfStates + l * nrOfPatterns * nrOfStates;
				int v2 = mutation2 * nrOfStates + l * nrOfPatterns * nrOfStates;
				int u = parentMut * nrOfStates + l * nrOfPatterns * nrOfStates;

				for (int i = 0; i < nrOfStates; i++) {
					partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
					u++;
					v1++;
					v2++;
				}
			}
		}	
	}
	

	/**
	 * Integrates partials across categories.
	 *
	 * @param inPartials  the array of partials to be integrated
	 * @param proportions the proportions of sites in each category
	 * @param outPartials an array into which the partials will go
	 */
	@Override
	protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {
		
		int u = 0;
		int v = 0;
		for (int k = 0; k < nrOfPatterns; k++) {

			for (int i = 0; i < nrOfStates; i++) {

				outPartials[u] = inPartials[v] * proportions[0];

				u++;
				v++;
			}
		}


		for (int l = 1; l < nrOfMatrices; l++) {
			u = 0;

			for (int k = 0; k < nrOfPatterns; k++) {

				for (int i = 0; i < nrOfStates; i++) {

					outPartials[u] += inPartials[v] * proportions[l];
					u++;
					v++;
				}
			}
		}
	}
	

	/**
	 * Calculates pattern log likelihoods at a node.
	 *
	 * @param partials          the partials used to calculate the likelihoods
	 * @param frequencies       an array of state frequencies
	 * @param outLogLikelihoods an array into which the likelihoods will go
	 */
	@Override
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
		int v = 0;
		for (int k = 0; k < nrOfPatterns; k++) {

			double sum = 0.0;
			for (int i = 0; i < nrOfStates; i++) {

				sum += frequencies[i] * partials[v];
				v++;
			}
			outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
		}
	}

	/**
	 * initializes partial likelihood arrays.
	 *
	 * @param nodeCount           the number of nodes in the tree
	 * @param patternCount        the number of patterns
	 * @param matrixCount         the number of matrices (i.e., number of
	 *                            categories)
	 * @param integrateCategories whether sites are being integrated over all
	 *                            matrices
	 */
	@Override
	public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories,
			boolean useAmbiguities) {

		this.nrOfNodes = nodeCount;
		this.nrOfPatterns = patternCount;
		this.nrOfMatrices = matrixCount;

		this.integrateCategories = integrateCategories;

		if (integrateCategories) {
			partialsSize = patternCount * nrOfStates * matrixCount;
		} else {
			partialsSize = patternCount * nrOfStates;
		}

		partials = new double[2][nodeCount][];
		edgePartials = new double[2][nodeCount][];

		currentMatrixIndex = new int[nodeCount];
		storedMatrixIndex = new int[nodeCount];

		currentPartialsIndex = new int[nodeCount];
		storedPartialsIndex = new int[nodeCount];
		
		currentEdgeIndex = new int[nodeCount];
		storedEdgeIndex = new int[nodeCount];

		states = new int[nodeCount][];

		for (int i = 0; i < nodeCount; i++) {
			partials[0][i] = null;
			partials[1][i] = null;
			edgePartials[0][i] = null;
			edgePartials[1][i] = null;

			states[i] = null;
		}

		matrixSize = nrOfStates * nrOfStates;

		matrices = new double[2][nodeCount][matrixCount * matrixSize];

		intermediatePartials1 = new double[partialsSize];
		intermediatePartials2 = new double[partialsSize];

		correction = nrOfMatrices * nrOfStates;
	}

	/**
	 * cleans up and deallocates arrays.
	 */
	@Override
	public void finalize() throws java.lang.Throwable {
		nrOfNodes = 0;
		nrOfPatterns = 0;
		nrOfMatrices = 0;

		partials = null;
		currentPartialsIndex = null;
		storedPartialsIndex = null;
		states = null;
		matrices = null;
		currentMatrixIndex = null;
		storedMatrixIndex = null;

		scalingFactors = null;
	}

	@Override
	public void setUseScaling(double scale) {
		useScaling = (scale != 1.0);

		if (useScaling) {
			scalingFactors = new double[2][nrOfNodes][nrOfPatterns];
		}
	}

	/**
	 * Allocates partials for a node
	 */
	@Override
	public void createNodePartials(int nodeIndex) {

		this.partials[0][nodeIndex] = new double[partialsSize];
		this.partials[1][nodeIndex] = new double[partialsSize];
	}

	public void createEdgePartials(int nodeIndex){
		this.edgePartials[0][nodeIndex] = new double[partialsSize];
		this.edgePartials[1][nodeIndex] = new double[partialsSize];
	}
	
	
	/**
	 * Sets partials for a node
	 */
	@Override
	public void setNodePartials(int nodeIndex, double[] partials) {

		if (this.partials[0][nodeIndex] == null) {
			createNodePartials(nodeIndex);
		}
		if (partials.length < partialsSize) {
			int k = 0;
			for (int i = 0; i < nrOfMatrices; i++) {
				System.arraycopy(partials, 0, this.partials[0][nodeIndex], k, partials.length);
				k += partials.length;
			}
		} else {
			System.arraycopy(partials, 0, this.partials[0][nodeIndex], 0, partials.length);
		}
	}
	
	public void setEdgePartials(int nodeIndex, double edgePartials) {
		if (this.edgePartials[0][nodeIndex] == null) {
			createEdgePartials(nodeIndex);
		}
		if (edgePartials < partialsSize) {
			int k = 0;
			for (int i = 0; i < nrOfMatrices; i++) {
				System.arraycopy(partials, 0, this.edgePartials[0][nodeIndex], k, partials.length);
				k += partials.length;
			}
		} else {
			System.arraycopy(partials, 0, this.edgePartials[0][nodeIndex], 0, partials.length);
		}
	}

	@Override
	public void getNodePartials(int nodeIndex, double[] partialsOut) {
		System.arraycopy(partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0, partialsOut, 0, partialsOut.length);
	}

	/**
	 * Allocates states for a node
	 */
	public void createNodeStates(int nodeIndex) {

		this.states[nodeIndex] = new int[nrOfStates];
	}

	/**
	 * Sets states for a node
	 */
	@Override
	public void setNodeStates(int nodeIndex, int[] states) {

		if (this.states[nodeIndex] == null) {
			createNodeStates(nodeIndex);
		}
		System.arraycopy(states, 0, this.states[nodeIndex], 0, nrOfStates);
	}
	
	public void setEdgeStates(int nodeIndex, int[] states) {

        if (this.states[nodeIndex] == null) {
            createNodeStates(nodeIndex);
        }
        System.arraycopy(states, 0, this.states[nodeIndex], 0, nrOfStates);
    }

	/**
	 * Gets states for a node
	 */
	@Override
	public void getNodeStates(int nodeIndex, int[] states) {
		System.arraycopy(this.states[nodeIndex], 0, states, 0, nrOfStates);
	}

	@Override
	public void setNodeMatrixForUpdate(int nodeIndex) {
		currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
	}
	
	public void setEdgeForUpdate(int nodeIndex) {
		currentEdgeIndex[nodeIndex] = 1 - currentEdgeIndex[nodeIndex];
	}
	/**
	 * Sets probability matrix for a node
	 */
	@Override
	public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
		System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex], matrixIndex * matrixSize,
				matrixSize);
	}

	public void setPaddedNodeMatrices(int nodeIndex, double[] matrix) {
		System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex], 0, nrOfMatrices * matrixSize);
	}

	/**
	 * Gets probability matrix for a node
	 */
	@Override
	public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
		System.arraycopy(matrices[currentMatrixIndex[nodeIndex]][nodeIndex], matrixIndex * matrixSize, matrix, 0,
				matrixSize);
	}

	@Override
	public void setNodePartialsForUpdate(int nodeIndex) {
		currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
	}

	/**
	 * Sets the currently updating node partials for node nodeIndex. This may need
	 * to repeatedly copy the partials for the different category partitions
	 */
	public void setCurrentNodePartials(int nodeIndex, double[] partials) {
		if (partials.length < partialsSize) {
			int k = 0;
			for (int i = 0; i < nrOfMatrices; i++) {
				System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], k,
						partials.length);
				k += partials.length;
			}
		} else {
			System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0,
					partials.length);
		}
	}

	/**
	 * Calculates partial likelihoods at a node.
	 *
	 * @param nodeIndex1       the 'child 1' node
	 * @param nodeIndex2       the 'child 2' node
	 * @param nodeIndex3       the 'parent' node
	 * @param calcForPatterns
	 * @param mutations2
	 * @param mutations1
	 * @param calcForPatterns2
	 */
	public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] mutations1,
			int[] mutations2, int[] calcForPatterns1, int[] calcForPatterns2, int[] calcForPatternsParent,
			boolean recalculate1, boolean recalculate2) {
		
		
		if (recalculate1) {
			if (states[nodeIndex1] != null) {
				calculateStatesPruning(states[nodeIndex1], 
						matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
						edgePartials[currentEdgeIndex[nodeIndex1]][nodeIndex1]);
			} else {
				 calculatePartialPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],
						matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1], 
						edgePartials[currentEdgeIndex[nodeIndex1]][nodeIndex1], calcForPatterns1);
			}
		}
		
		if (recalculate2) {
			if (states[nodeIndex2] != null) {
				calculateStatesPruning(states[nodeIndex2], 
						matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
						edgePartials[currentEdgeIndex[nodeIndex2]][nodeIndex2]);
	        } else {
	        	calculatePartialPruning(partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],
						matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2], 
						edgePartials[currentEdgeIndex[nodeIndex2]][nodeIndex2], calcForPatterns2);
	        }
		}
		
		updatePartials(partials[currentPartialsIndex[nodeIndex3]][nodeIndex3], 
				edgePartials[currentEdgeIndex[nodeIndex1]][nodeIndex1], edgePartials[currentEdgeIndex[nodeIndex2]][nodeIndex2],
				mutations1, mutations2, calcForPatternsParent);		
		
		
		
//		if (states[nodeIndex1] != null) {
//			calculateStatesPruning(states[nodeIndex1], 
//					matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
//					intermediatePartials1);
//		} else {
//			 calculatePartialPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],
//					matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1], 
//					intermediatePartials1, calcForPatterns1);
//		}
//		
//		if (states[nodeIndex2] != null) {
//			calculateStatesPruning(states[nodeIndex2], 
//					matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
//					intermediatePartials2);
//        } else {
//        	calculatePartialPruning(partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],
//					matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2], 
//					intermediatePartials2, calcForPatterns2);
//        }
//		
//		updatePartials(partials[currentPartialsIndex[nodeIndex3]][nodeIndex3], 
//				intermediatePartials1, intermediatePartials2,
//				mutations1, mutations2, calcForPatternsParent);		
        
		// print out all the targets fo mutations1
		if (useScaling) {
			scalePartials(nodeIndex3);
		}
	}

	@Override
	public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
		calculateIntegratePartials(partials[currentPartialsIndex[nodeIndex]][nodeIndex], proportions, outPartials);
	}

	/**
	 * Scale the partials at a given node. This uses a scaling suggested by Ziheng
	 * Yang in Yang (2000) J. Mol. Evol. 51: 423-432
	 * <p/>
	 * This function looks over the partial likelihoods for each state at each
	 * pattern and finds the largest. If this is less than the scalingThreshold
	 * (currently set to 1E-40) then it rescales the partials for that pattern by
	 * dividing by this number (i.e., normalizing to between 0, 1). It then stores
	 * the log of this scaling. This is called for every internal node after the
	 * partials are calculated so provides most of the performance hit. Ziheng
	 * suggests only doing this on a proportion of nodes but this sounded like a
	 * headache to organize (and he doesn't use the threshold idea which improves
	 * the performance quite a bit).
	 *
	 * @param nodeIndex
	 */
	protected void scalePartials(int nodeIndex) {
//        int v = 0;
//    	double [] partials = m_fPartials[m_iCurrentPartialsIndices[nodeIndex]][nodeIndex];
//        for (int i = 0; i < m_nPatternCount; i++) {
//            for (int k = 0; k < m_nMatrixCount; k++) {
//                for (int j = 0; j < m_nStateCount; j++) {
//                	partials[v] *= SCALE;
//                	v++;
//                }
//            }
//        }
		int u = 0;

		for (int i = 0; i < nrOfPatterns; i++) {

			double scaleFactor = 0.0;
			int v = u;
			for (int k = 0; k < nrOfMatrices; k++) {
				for (int j = 0; j < nrOfStates; j++) {
					if (partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] > scaleFactor) {
						scaleFactor = partials[currentPartialsIndex[nodeIndex]][nodeIndex][v];
					}
					v++;
				}
				v += (nrOfPatterns - 1) * nrOfStates;
			}

			if (scaleFactor < scalingThreshold) {

				v = u;
				for (int k = 0; k < nrOfMatrices; k++) {
					for (int j = 0; j < nrOfStates; j++) {
						partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] /= scaleFactor;
						v++;
					}
					v += (nrOfPatterns - 1) * nrOfStates;
				}
				scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = Math.log(scaleFactor);

			} else {
				scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = 0.0;
			}
			u += nrOfStates;

		}
	}

	/**
	 * This function returns the scaling factor for that pattern by summing over the
	 * log scalings used at each node. If scaling is off then this just returns a 0.
	 *
	 * @return the log scaling factor
	 */
	@Override
	public double getLogScalingFactor(int patternIndex_) {
//    	if (m_bUseScaling) {
//    		return -(m_nNodeCount/2) * Math.log(SCALE);
//    	} else {
//    		return 0;
//    	}        
		double logScalingFactor = 0.0;
		if (useScaling) {
			for (int i = 0; i < nrOfNodes; i++) {
				logScalingFactor += scalingFactors[currentPartialsIndex[i]][i][patternIndex_];
			}
		}
		return logScalingFactor;
	}

	/**
	 * Gets the partials for a particular node.
	 *
	 * @param nodeIndex   the node
	 * @param outPartials an array into which the partials will go
	 */
	public void getPartials(int nodeIndex, double[] outPartials) {
		double[] partials1 = partials[currentPartialsIndex[nodeIndex]][nodeIndex];

		System.arraycopy(partials1, 0, outPartials, 0, partialsSize);
	}

	/**
	 * Store current state
	 */
	@Override
	public void restore() {
		// Rather than copying the stored stuff back, just swap the pointers...
		int[] tmp1 = currentMatrixIndex;
		currentMatrixIndex = storedMatrixIndex;
		storedMatrixIndex = tmp1;

		int[] tmp2 = currentPartialsIndex;
		currentPartialsIndex = storedPartialsIndex;
		storedPartialsIndex = tmp2;
		
		int[] tmp3 = currentEdgeIndex;
		currentEdgeIndex = storedEdgeIndex;
		storedEdgeIndex = tmp3;
	}

	@Override
	public void unstore() {
		System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
		System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);
		System.arraycopy(storedEdgeIndex, 0, currentEdgeIndex, 0, nrOfNodes);
	}

	/**
	 * Restore the stored state
	 */
	@Override
	public void store() {
		System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
		System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);
		System.arraycopy(currentEdgeIndex, 0, storedEdgeIndex, 0, nrOfNodes);
	}

//	@Override
//    public void calcRootPsuedoRootPartials(double[] frequencies, int nodeIndex, double [] pseudoPartials) {
//		int u = 0;
//		double [] inPartials = m_fPartials[m_iCurrentPartials[nodeIndex]][nodeIndex];
//		for (int k = 0; k < m_nPatterns; k++) {
//			for (int l = 0; l < m_nMatrices; l++) {
//				for (int i = 0; i < m_nStates; i++) {
//					pseudoPartials[u] = inPartials[u] * frequencies[i];
//					u++;
//				}
//			}
//		}
//    }
//	@Override
//    public void calcNodePsuedoRootPartials(double[] inPseudoPartials, int nodeIndex, double [] outPseudoPartials) {
//		double [] partials = m_fPartials[m_iCurrentPartials[nodeIndex]][nodeIndex];
//		double [] oldPartials = m_fPartials[m_iStoredPartials[nodeIndex]][nodeIndex];
//		int maxK = m_nPatterns * m_nMatrices * m_nStates; 
//		for (int k = 0; k < maxK; k++) {
//			outPseudoPartials[k] = inPseudoPartials[k] * partials[k] / oldPartials[k];
//		}
//	}
//    
//	@Override
//    public void calcPsuedoRootPartials(double [] parentPseudoPartials, int nodeIndex, double [] pseudoPartials) {
//		int v = 0;
//		int u = 0;
//		double [] matrices = m_fMatrices[m_iCurrentMatrices[nodeIndex]][nodeIndex];
//		for (int k = 0; k < m_nPatterns; k++) {
//			for (int l = 0; l < m_nMatrices; l++) {
//				for (int i = 0; i < m_nStates; i++) {
//					int w = 0;
//					double sum = 0;
//					for (int j = 0; j < m_nStates; j++) {
//					      sum += parentPseudoPartials[u+j] * matrices[w + i];
//					      w+=m_nStates;
//					}
//					pseudoPartials[v] = sum;
//					v++;
////					int w = l * m_nMatrixSize;
////					double sum = 0;
////					for (int j = 0; j < m_nStates; j++) {
////					      sum += parentPseudoPartials[u+j] * matrices[w+j];
////					}
////					pseudoPartials[v] = sum;
////					v++;
//				}
//				u += m_nStates;
//			}
//		}
//    }
//
//
//    @Override
//    void integratePartialsP(double [] inPartials, double [] proportions, double [] m_fRootPartials) {
//		int maxK = m_nPatterns * m_nStates;
//		for (int k = 0; k < maxK; k++) {
//			m_fRootPartials[k] = inPartials[k] * proportions[0];
//		}
//
//		for (int l = 1; l < m_nMatrices; l++) {
//			int n = maxK * l;
//			for (int k = 0; k < maxK; k++) {
//				m_fRootPartials[k] += inPartials[n+k] * proportions[l];
//			}
//		}
//    } // integratePartials
//
//	/**
//	 * Calculates pattern log likelihoods at a node.
//	 * @param partials the partials used to calculate the likelihoods
//	 * @param frequencies an array of state frequencies
//	 * @param outLogLikelihoods an array into which the likelihoods will go
//	 */
//    @Override
//	public void calculateLogLikelihoodsP(double[] partials,double[] outLogLikelihoods)
//	{
//        int v = 0;
//		for (int k = 0; k < m_nPatterns; k++) {
//            double sum = 0.0;
//			for (int i = 0; i < m_nStates; i++) {
//				sum += partials[v];
//				v++;
//			}
//            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
//		}
//	}
//	
//	
//	//    @Override
////    LikelihoodCore feelsGood() {return null;}

	@Override
	public boolean getUseScaling() {
		return useScaling;
	}

	@Override
	public void calculatePartials(int node1, int node2Index, int node3) {
		throw new RuntimeException("Not implemented");
	}

} // class BeerLikelihoodCore
