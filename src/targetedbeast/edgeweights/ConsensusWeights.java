package targetedbeast.edgeweights;

import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;

@Description("Keeps track of the consensus sequences and the number of mutations between consensus sequences along edges")
public class ConsensusWeights extends CalculationNode implements EdgeWeights {
	
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);
    
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    
    final public Input<Double> maxWeightInput = new Input<>("maxWeight", "maximum weight for an edge", 10.0);
    
    final public Input<Double> minWeightInput = new Input<>("minWeight", "maximum weight for an edge", 0.01);

	protected int hasDirt;

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
	int patternCount;	
	int maxStateCount;
	
	double maxWeight;
	double minWeight;

	@Override
	public void initAndValidate() {

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
		patternCount = dataInput.get().getPatternCount();
		maxStateCount = dataInput.get().getMaxStateCount();

		edgeMutations = new double[2][treeInput.get().getNodeCount()];
		// should probably be changes to a non double
		consensus = new double[2][treeInput.get().getNodeCount()][patternCount * stateCount];

		activeIndex = new int[treeInput.get().getNodeCount()];
		storedActiveIndex = new int[treeInput.get().getNodeCount()];

		activeMutationsIndex = new int[treeInput.get().getNodeCount()];
		storedActiveMutationsIndex = new int[treeInput.get().getNodeCount()];

		changed = new boolean[treeInput.get().getNodeCount()];
		changedChildren = new boolean[treeInput.get().getNodeCount()];
			
		maxWeight = maxWeightInput.get();
		minWeight = minWeightInput.get();

		initLeaveConsensus(treeInput.get().getRoot());
		
		updateWeights();
	}

	
	public void updateWeights() {
		Arrays.fill(changed, true);
		Arrays.fill(changedChildren, true);
		getNodeConsensusSequences(treeInput.get().getRoot());
	}

	public void updateByOperator() {
		operatorUpdated = true;
		Arrays.fill(changed, false);
		Arrays.fill(changedChildren, false);
		getFilthyNodes(treeInput.get().getRoot());
		getNodeConsensusSequences(treeInput.get().getRoot());
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

	private void initLeaveConsensus(Node n) {
		if (n.isLeaf()) {
			int[] patterns = new int[patternCount];
			for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
				patterns[i] = dataInput.get().getPattern(n.getNr(), i);
				if (patterns[i] >= maxStateCount) {
					for (int j = 0; j < dataInput.get().getMaxStateCount(); j++) {
						consensus[0][n.getNr()][i * stateCount + j] = 1.0/maxStateCount;
					}
					patterns[i] = -2;
				}else {				
					consensus[0][n.getNr()][i * stateCount + patterns[i]] = 1;
				}			
			}
		} else {
			initLeaveConsensus(n.getLeft());
			initLeaveConsensus(n.getRight());
		}
	}

	public void getNodeConsensusSequences(Node n) {
		if (n.isLeaf()) {
			// the active index for leaves is always 0, could be changed to make the arrays
			// shorter
			return;
		} else {
			if (changed[n.getNr()]) {
				// compare the patterns of the two lineages
				getNodeConsensusSequences(n.getLeft());
				getNodeConsensusSequences(n.getRight());

				// set this node index to the active index
				activeIndex[n.getNr()] = 1 - activeIndex[n.getNr()];
				activeMutationsIndex[n.getLeft().getNr()] = 1 - activeMutationsIndex[n.getLeft().getNr()];
				activeMutationsIndex[n.getRight().getNr()] = 1 - activeMutationsIndex[n.getRight().getNr()];

				int activeInd = activeIndex[n.getNr()];
				int activeIndLeft = activeIndex[n.getLeft().getNr()];
				int activeIndRight = activeIndex[n.getRight().getNr()];
				double sumMuts = minWeight;
				for (int i = 0; i < consensus[0][0].length; i++) {
					double val = (consensus[activeIndLeft][n.getLeft().getNr()][i]
							+ consensus[activeIndRight][n.getRight().getNr()][i]) / 2;
					if (val > 0.5) {
						val = 1;
					} else if (val < 0.5) {
						val = 0;
					} else {
						sumMuts += 1;
					}
					consensus[activeInd][n.getNr()][i] = val;
				}
				edgeMutations[activeMutationsIndex[n.getLeft().getNr()]][n.getLeft().getNr()] = Math.min(maxWeight, sumMuts);
				edgeMutations[activeMutationsIndex[n.getRight().getNr()]][n.getRight().getNr()] = Math.min(maxWeight, sumMuts);
				
			} else {
				getNodeConsensusSequences(n.getLeft());
				getNodeConsensusSequences(n.getRight());
			}
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
		
		if (!operatorUpdated) {
			Arrays.fill(changed, false);
			Arrays.fill(changedChildren, false);
			getFilthyNodes(treeInput.get().getRoot());
			getNodeConsensusSequences(treeInput.get().getRoot());
		}
		
		operatorUpdated = false;
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
