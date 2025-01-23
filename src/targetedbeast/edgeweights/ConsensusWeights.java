package targetedbeast.edgeweights;

import java.io.PrintStream;
import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import beast.base.inference.State;
import targetedbeast.alignment.RapidAlignment;
import targetedbeast.likelihood.RapidLikelihoodCore;
import targetedbeast.likelihood.RapidLikelihoodCore4;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using "
		+ "a variant of the 'peeling algorithm'. For details, see"
		+ "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class ConsensusWeights extends CalculationNode implements EdgeWeights {
	
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);
    
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	
	protected int hasDirt;


	/**
	 * memory allocation for likelihoods for each of the patterns *
	 */
	protected double[] patternLogLikelihoods;
	/**
	 * memory allocation for the root partials *
	 */
	protected double[] m_fRootPartials;
	/**
	 * memory allocation for probability tables obtained from the SiteModel *
	 */
	protected double[] probabilities;

	protected int matrixSize;

	private int[][][] calcForPatterns;

//	private List<Integer[]>[][] mutations;
	private int[][][] mutations;

	private int[][][] map;

	private double[][][] consensus;

	private boolean[] changed;
	private boolean[] changedChildren;

	private int[] activeIndex;
	private int[] storedActiveIndex;

	private int[] activeMutationsIndex;
	private int[] storedActiveMutationsIndex;

	public double[][] edgeMutations;

	private boolean operatorUpdated = false;

	int stateCount;
	public int patterns;

	/**
	 * flag to indicate ascertainment correction should be applied *
	 */
	protected boolean useAscertainedSitePatterns = false;

	/**
	 * dealing with proportion of site being invariant *
	 */
	protected double proportionInvariant = 0;

	public double getProportionInvariant() {
		return proportionInvariant;
	}

	public void setProportionInvariant(double proportionInvariant) {
		this.proportionInvariant = proportionInvariant;
	}

	List<Integer> constantPattern = null;

	public List<Integer> getConstantPattern() {
		return constantPattern;
	}

	public void setConstantPattern(List<Integer> constantPattern) {
		this.constantPattern = constantPattern;
	}
	
	long countwith=0;
	long countwithout=0;
	
	int count=0;

	@Override
	public void initAndValidate() {
		// sanity check: alignment should have same #taxa as tree

		if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
			String leaves = "?";
			if (treeInput.get() instanceof Tree) {
				leaves = String.join(", ", ((Tree) treeInput.get()).getTaxaNames());
			}
			throw new IllegalArgumentException(String.format(
					"The number of leaves in the tree (%d) does not match the number of sequences (%d). "
							+ "The tree has leaves [%s], while the data refers to taxa [%s].",
					treeInput.get().getLeafNodeCount(), dataInput.get().getTaxonCount(), leaves,
					String.join(", ", dataInput.get().getTaxaNames())));
		}

		stateCount = dataInput.get().getMaxStateCount();
		patterns = dataInput.get().getPatternCount();

		edgeMutations = new double[2][treeInput.get().getNodeCount()];

		consensus = new double[2][treeInput.get().getNodeCount()][patterns * stateCount];

		activeIndex = new int[treeInput.get().getNodeCount()];
		storedActiveIndex = new int[treeInput.get().getNodeCount()];

		activeMutationsIndex = new int[treeInput.get().getNodeCount()];
		storedActiveMutationsIndex = new int[treeInput.get().getNodeCount()];

		changed = new boolean[treeInput.get().getNodeCount()];
		changedChildren = new boolean[treeInput.get().getNodeCount()];

		getLeaveMutations(treeInput.get().getRoot());
		
		updateMutations();
	}

	protected RapidLikelihoodCore createLikelihoodCore(int stateCount) {
		if (stateCount == 4) {
			return new RapidLikelihoodCore4();
		} else {
			throw new IllegalArgumentException("Only stateCount=4 is supported at the moment");
//			return new RapidLikelihoodCore(stateCount);
		}
	}

	double ratio = 0;
	int nrratio = 0;
	int nrratio2 = 0;
	double totrati=0;
	
	public void updateMutations() {
		Arrays.fill(changed, true);
		Arrays.fill(changedChildren, true);
		getMutations(treeInput.get().getRoot());
	}

	public void updateByOperator() {
		operatorUpdated = true;
		Arrays.fill(changed, false);
		Arrays.fill(changedChildren, false);
		getFilthyNodes(treeInput.get().getRoot());
		getMutations(treeInput.get().getRoot());
	}

	public void updateByOperatorWithoutNode(int ignore, List<Integer> nodes) {
		Arrays.fill(changed, false);
		Arrays.fill(changedChildren, false);
		for (Integer nodeNo : nodes) {
			changed[nodeNo] = true;
		}
		getConsensusWithoutNode(treeInput.get().getRoot(), ignore);
	}

	public void fakeUpdateByOperator() {
		operatorUpdated = true;
		// used for operators that change the tree, but without affecting the order of
		// the patterns
	}


	private boolean getFilthyNodes(Node node) {
		// compute the number of patterns for each node.
		if (node.isLeaf()) {
			return false;
		} else {
			boolean left = getFilthyNodes(node.getLeft());
			boolean right = getFilthyNodes(node.getRight());
			if (left || right || node.isDirty() > 1) {
				changed[node.getNr()] = true;
				if (node.isDirty() == 3) {
					return false;
				} else {
				}
			}
		}

		return changed[node.getNr()];
	}

	/** CalculationNode methods **/

	private void getLeaveMutations(Node n) {
		if (n.isLeaf()) {
			int[] patterns = new int[dataInput.get().getPatternCount()];
			for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
				patterns[i] = dataInput.get().getPattern(n.getNr(), i);
				if (patterns[i] > 3) {
					for (int j = 0; j < dataInput.get().getMaxStateCount(); j++) {
						consensus[0][n.getNr()][i * stateCount + j] = 1.0/dataInput.get().getMaxStateCount();
					}
					patterns[i] = -2;
				}else {				
					consensus[0][n.getNr()][i * stateCount + patterns[i]] = 1;
				}
				
				
			}
			calcForPatterns[0][n.getNr()] = new int[dataInput.get().getPatternCount()];
			calcForPatterns[0][n.getNr()][0] = -1;

			System.arraycopy(patterns, 0, map[0][n.getNr()], 0, patterns.length);
		} else {
			getLeaveMutations(n.getLeft());
			getLeaveMutations(n.getRight());
		}
	}

	public int[] getMutations(Node n) {
		if (n.isLeaf()) {
			// the active index for leaves is always 0, could be changed to make the arrays
			// shorter
			return map[0][n.getNr()];
		} else {
			if (changed[n.getNr()]) {
				// compare the patterns of the two lineages
				int[] left = getMutations(n.getLeft());
				int[] right = getMutations(n.getRight());

				// set this node index to the active index
				activeIndex[n.getNr()] = 1 - activeIndex[n.getNr()];
				activeMutationsIndex[n.getLeft().getNr()] = 1 - activeMutationsIndex[n.getLeft().getNr()];
				activeMutationsIndex[n.getRight().getNr()] = 1 - activeMutationsIndex[n.getRight().getNr()];

				int activeInd = activeIndex[n.getNr()];
				int activeIndLeft = activeIndex[n.getLeft().getNr()];
				int activeIndRight = activeIndex[n.getRight().getNr()];
				for (int i = 0; i < consensus[0][0].length; i++) {
					double val = (consensus[activeIndLeft][n.getLeft().getNr()][i]
							+ consensus[activeIndRight][n.getRight().getNr()][i]) / 2;
					if (val > 0.5) {
						val = 1;
					} else if (val < 0.5) {
						val = 0;
					}
					consensus[activeInd][n.getNr()][i] = val;

				}

				// copy the left pattern to the new map
				System.arraycopy(left, 0, map[activeIndex[n.getNr()]][n.getNr()], 0, left.length);
				calcForPatterns[activeIndex[n.getNr()]][n.getNr()] = new int[dataInput.get().getPatternCount()];

				updatePatterns(left, right, mutations[activeMutationsIndex[n.getLeft().getNr()]][n.getLeft().getNr()],
						mutations[activeMutationsIndex[n.getRight().getNr()]][n.getRight().getNr()],
						calcForPatterns[activeIndex[n.getNr()]][n.getNr()], map[activeIndex[n.getNr()]][n.getNr()],
						n.getLeft().getNr(), n.getRight().getNr(), n.getNr());

			} else {
				getMutations(n.getLeft());
				getMutations(n.getRight());
			}
			// print the left and the right

			return map[activeIndex[n.getNr()]][n.getNr()];
		}
	}

	public void getConsensusWithoutNode(Node n, int ignore) {
		if (n.isLeaf()) {
			return;
		} else {
			if (changed[n.getNr()]) {

				// set this node index to the active index
				activeIndex[n.getNr()] = 1 - activeIndex[n.getNr()];

				int activeInd = activeIndex[n.getNr()];
				int activeIndLeft = activeIndex[n.getLeft().getNr()];
				int activeIndRight = activeIndex[n.getRight().getNr()];
				
				if (n.getLeft().getNr() != ignore && n.getRight().getNr() != ignore) {
					for (int i = 0; i < consensus[0][0].length; i++) {
						double val = (consensus[activeIndLeft][n.getLeft().getNr()][i]
								+ consensus[activeIndRight][n.getRight().getNr()][i]) / 2;
						if (val > 0.5) {
							val = 1;
						} else if (val < 0.5) {
							val = 0;
						}
						consensus[activeInd][n.getNr()][i] = val;
					}
				} else if (n.getLeft().getNr() == ignore) {
					System.arraycopy(consensus[activeIndRight][n.getRight().getNr()], 0,
							consensus[activeInd][n.getNr()], 0, consensus[0][0].length); // cop
				} else if (n.getRight().getNr() == ignore) {
					System.arraycopy(consensus[activeIndLeft][n.getLeft().getNr()], 0, consensus[activeInd][n.getNr()],
							0, consensus[0][0].length);
				}
			}
			getConsensusWithoutNode(n.getLeft(), ignore);
			getConsensusWithoutNode(n.getRight(), ignore);
			return;
		}
	}

	private void updatePatterns(int[] left, int[] right, int[] mutations_left, int[] mutations_right, int[] pattern,
			int[] map, int nodeLeft, int nodeRight, int node) {

		// get the number of leaves below left and right
		double sumleft = 0.01;
		double sumright = 0.01;

		int c = 0;
		for (int i = 0; i < left.length; i++) {
			// check if the consensus patterns differ
			if (left[i] == right[i]) {
				if (left[i] == -1) { // no new pattern has been added
					pattern[c] = i;
					c++;
					mutations_left[i] = i;
					mutations_right[i] = i;
					for (int j = 0; j < stateCount; j++) {
						double parent = consensus[activeIndex[node]][node][i * stateCount + j];
						if (parent > 0) {
							sumleft += Math
									.abs(parent - consensus[activeIndex[nodeLeft]][nodeLeft][i * stateCount + j]);
							sumright += Math
									.abs(parent - consensus[activeIndex[nodeRight]][nodeRight][i * stateCount + j]);
						}
					}
				}
			} else {
				for (int j = 0; j < stateCount; j++) {
					double parent = consensus[activeIndex[node]][node][i * stateCount + j];
					if (parent > 0) {
						sumleft += Math.abs(parent - consensus[activeIndex[nodeLeft]][nodeLeft][i * stateCount + j]);
						sumright += Math.abs(parent - consensus[activeIndex[nodeRight]][nodeRight][i * stateCount + j]);
					}
				}

				pattern[c] = i;
				c++;
				if (left[i] != -1) {
					// get the current pattern for this nodes
					mutations_left[i] = left[i];
				} else { // the left edge already has this pattern
					mutations_left[i] = i;
				}
				if (right[i] != -1) {
					mutations_right[i] = right[i];

				} else { // the right edge already has this pattern
					mutations_right[i] = i;
				}
				map[i] = -1;
			}
		}
		pattern[c] = -1;
		edgeMutations[activeMutationsIndex[nodeLeft]][nodeLeft] = sumleft / 2;
		edgeMutations[activeMutationsIndex[nodeRight]][nodeRight] = sumright / 2;
	}

	
	/**
	 * check state for changed variables and update temp results if necessary *
	 */
	@Override
	protected boolean requiresRecalculation() {
		hasDirt = Tree.IS_CLEAN;

		if (dataInput.get().isDirtyCalculation()) {
			hasDirt = Tree.IS_FILTHY;
			return true;
		}
		return treeInput.get().somethingIsDirty();
	}

	@Override
	public void store() {
//    	System.err.println("store ");
		super.store();
		if (!operatorUpdated) { // avoid storing again if the operator has already done it
			if (operatorUpdated) {
				operatorUpdated = false;
				return;
			}
			System.arraycopy(activeIndex, 0, storedActiveIndex, 0, activeIndex.length);
			System.arraycopy(activeMutationsIndex, 0, storedActiveMutationsIndex, 0, activeMutationsIndex.length);
		}
	}

	public void prestore() {
//		System.out.println("prestore ");
		System.arraycopy(activeIndex, 0, storedActiveIndex, 0, activeIndex.length);
		System.arraycopy(activeMutationsIndex, 0, storedActiveMutationsIndex, 0, activeMutationsIndex.length);
	}

	public void reset() {
//		System.out.println("reset ");
		// undoes any previous calculation
		System.arraycopy(storedActiveIndex, 0, activeIndex, 0, activeIndex.length);
	}

	public void unstore() {
	}

	@Override
	public void restore() {
//		System.out.println("restore ");


		super.restore();

		System.arraycopy(storedActiveIndex, 0, activeIndex, 0, activeIndex.length);
		System.arraycopy(storedActiveMutationsIndex, 0, activeMutationsIndex, 0, activeMutationsIndex.length);
	}


	public double getEdgeMutations(int i) {
		return edgeMutations[activeMutationsIndex[i]][i];
	}
	
	public int[] getCalcPatterns(int i) {
		return calcForPatterns[activeIndex[i]][i];
	}

	public boolean getChanged(int i) {
		return changed[i];
	}

	public double[] getConsensus(int nr) {
		return consensus[activeIndex[nr]][nr];
	}

	@Override
	public double getEdgeWeights(int nodeNr) {
		return getEdgeMutations(nodeNr);
	}

	@Override
	public double[] getTargetWeights(int fromNodeNr, List<Node> toNodeNrs) {
		double[] distances = new double[toNodeNrs.size()];
		double[] currConsensus = getConsensus(fromNodeNr);
		
		for (int k = 0; k < toNodeNrs.size(); k++) {
			int nodeNo = toNodeNrs.get(k).getNr();
			double[] consensus = getConsensus(nodeNo);
			// calculate the distance between the two consensus
			double sum = 0.1;
			for (int l = 0; l < consensus.length; l++) {
				sum += Math.abs(currConsensus[l] - consensus[l]);
			}
			distances[k] = 1 / (sum);
		}		
		return distances;
	}

} // class TreeLikelihood
