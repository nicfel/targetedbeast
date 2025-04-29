package targetedbeast.edgeweights;

import java.io.PrintStream;
import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;

@Description("Keeps track of the consensus sequences and the number of mutations between consensus sequences along edges"
		+ "Consensus weights is a distribution to ensure that it is updated correctly")
public class ParsimonyWeights extends Distribution implements EdgeWeights, Loggable {
	
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);
    
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    
    final public Input<Double> maxWeightInput = new Input<>("maxWeight", "maximum weight for an edge", 10.0);
    
    final public Input<Double> minWeightInput = new Input<>("minWeight", "maximum weight for an edge", 0.01);

	protected int hasDirt;

	private byte[][][] consensus;

	private boolean[] changed;
	private boolean[] changedChildren;

	private int[] activeIndex;
	private int[] storedActiveIndex;

	private int[] activeMutationsIndex;
	private int[] storedActiveMutationsIndex;

	public double[][] edgeMutations;
	
	public List<Mutation>[][] mutations;

	private boolean operatorUpdated = false;

	int stateCount;
	int patternCount;	
	int maxStateCount;
	
	double maxWeight;
	double minWeight;
	
	double totalMuts[];

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

		
		if (stateCount != 4) {
			throw new IllegalArgumentException("Only 4 state sequences are supported at this point");
		}
		
		edgeMutations = new double[2][treeInput.get().getNodeCount()];
		// should probably be changes to a non double
		consensus = new byte[2][treeInput.get().getNodeCount()][patternCount];

		activeIndex = new int[treeInput.get().getNodeCount()];
		storedActiveIndex = new int[treeInput.get().getNodeCount()];

		activeMutationsIndex = new int[treeInput.get().getNodeCount()];
		storedActiveMutationsIndex = new int[treeInput.get().getNodeCount()];
		
		mutations = new ArrayList[2][treeInput.get().getNodeCount()];

		changed = new boolean[treeInput.get().getNodeCount()];
		changedChildren = new boolean[treeInput.get().getNodeCount()];
			
		maxWeight = maxWeightInput.get();
		minWeight = minWeightInput.get();
		
		totalMuts = new double[patternCount];

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
			int patterns;
			for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
				patterns = dataInput.get().getPattern(n.getNr(), i);			
				if (patterns >= maxStateCount) {
					consensus[0][n.getNr()][i] = 10;
				}else {				
					consensus[0][n.getNr()][i] = (byte) patterns;
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
				byte left, right, val;
				
				mutations[activeIndLeft][n.getLeft().getNr()] = new ArrayList<>();
				mutations[activeIndRight][n.getRight().getNr()] = new ArrayList<>();
				
				double sumLeft = 0;
				double sumRight = 0;
				
				for (int i = 0; i < consensus[0][0].length; i++) {
					left = consensus[activeIndLeft][n.getLeft().getNr()][i];
					right = consensus[activeIndRight][n.getRight().getNr()][i];
										
					if (left!=right) {
						val = getCombination(left, right);	
//						System.out.println("Left: " + left + " Right: " + right + " Val: " + val);
//						mutations[activeIndLeft][n.getLeft().getNr()].add(new Mutation(left, val, i));
//						mutations[activeIndRight][n.getRight().getNr()].add(new Mutation(right, val, i));
						sumLeft += getDiff(left, val);
						sumRight += getDiff(right, val);
					} else {
						val = left;
					}
					consensus[activeInd][n.getNr()][i] = val;
				}
				edgeMutations[activeMutationsIndex[n.getLeft().getNr()]][n.getLeft().getNr()] = Math.min(maxWeight, sumLeft);
				edgeMutations[activeMutationsIndex[n.getRight().getNr()]][n.getRight().getNr()] = Math.min(maxWeight, sumRight);				
			} else {
				getNodeConsensusSequences(n.getLeft());
				getNodeConsensusSequences(n.getRight());
			}
		}
	}

	private double getDiff(byte child, byte parent) {
		if (parent == child) {
			return 0;
		}else if (child == 10) {
			return 0;
		}else if (child < 4 && parent <4) {
			return 1;
        }else {
        	return 0.5;
        }		
	}


	private byte getCombination(byte left, byte right) {
		if (left > 3 && right > 3) {
			if (left < right) {
				if (left == 4) {
					if (right == 5) {
                        return 0;
                    } else if (right == 6) {
                        return 0;
                    } else if (right == 10) {
                    	return 4;                    	
                    } else {
                        return 10;
                    }
				}else if (left == 5) {
                    if (right == 6) {
                        return 0;
                    } else if(right == 7) {
                        return 2;
                    } else if(right == 10) {
                    	return 5;
                    }else {
                    	return 10;
                    }
                    
                } else if(left == 6) {
					if (right == 8) {
						return 3;
					} else if (right == 9) {
						return 3;
					} else if (right == 10) {
						return 6;
					} else {
						return 10;
					}
                } else if (left ==7) {
                    if (right == 8) {
                        return 1;
                    } else if (right == 10) {
                    	return 7;
                    } else {
                        return 10;
                    }
                } else if (left == 8) {
					if (right == 9) {
						return 3;
					} else if (right == 10) {
						return 8;
					} else {
						return 10;
					}
                } else if (left == 9) {
					return 9; // can only happen if right==10
                } else {
                    throw new IllegalArgumentException("Error in getCombination");
                }
			} else { // left > right
                if (right == 4) {
                    if (left == 5) {
                        return 0;
                    } else if (left == 6) {
                        return 0;
                    } else if (left == 10) {
                        return 4;
                    } else {
                        return 10;
                    }
                } else if (right == 5) {
                    if (left == 6) {
                        return 0;
                    } else if (left == 7) {
                        return 2;
                    } else if (left == 10) {
                        return 5;
                    } else {
                        return 10;
                    }
                } else if (right == 6) {
                    if (left == 8) {
                        return 3;
                    } else if (left == 9) {
                        return 3;
                    } else if (left == 10) {
                        return 6;
                    } else {
                        return 10;
                    }
                } else if (right == 7) {
                    if (left == 8) {
                        return 1;
                    } else if (left == 10) {
                        return 7;
                    } else {
                        return 10;
                    }
                } else if (right == 8) {
                    if (left == 9) {
                        return 3;
                    } else if (left == 10) {
                        return 8;
                    } else {
                        return 10;
                    }
                } else if (right == 9) {
                    return 9; // can only happen if left==10
                } else {
                    throw new IllegalArgumentException("Error in getCombination");
                }				
			}		
		}
		
		if (left < right) {
			if (left==0) {
				if (right==1) {
                    return 4;
				}else if (right==2) {
                    return 5;
				} else if (right == 3) {
					return 6;
				} else {
					return 0;
				}
			} else if (left == 1) {
				if (right == 2) {
					return 7;
				} else if (right == 3) {
					return 8;
				} else {
					return 1;
				}
			} else if (left == 2) {
				if (right == 3) {
					return 9;
				} else {
					return 2;
				}
			} else {
				return 3;
			}
		}else {
			// do the opposite of above
			if (right == 0) {
				if (left == 1) {
					return 4;
				} else if (left == 2) {
					return 5;
				} else if (left == 3) {
					return 6;
				} else {
					return 0;
				}
			} else if (right == 1) {
				if (left == 2) {
					return 7;
				} else if (left == 3) {
					return 8;
				} else {
					return 1;
				}
			} else if (right == 2) {
				if (left == 3) {
					return 9;
				} else {
					return 2;
				}
			} else {
				return 3;
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
				byte left, right, val;

				
				if (n.getLeft().getNr() != ignore && n.getRight().getNr() != ignore) {
					for (int i = 0; i < consensus[0][0].length; i++) {
						left = consensus[activeIndLeft][n.getLeft().getNr()][i];
						right = consensus[activeIndRight][n.getRight().getNr()][i];
											
						if (left!=right) {
							val = getCombination(left, right);	
						} else {
							val = left;
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
		System.arraycopy(activeIndex, 0, storedActiveIndex, 0, activeIndex.length);
		System.arraycopy(activeMutationsIndex, 0, storedActiveMutationsIndex, 0, activeMutationsIndex.length);
	}

	public void reset() {
		// undoes any previous calculation
		System.arraycopy(storedActiveIndex, 0, activeIndex, 0, activeIndex.length);
		System.arraycopy(storedActiveMutationsIndex, 0, activeMutationsIndex, 0, activeMutationsIndex.length);
		operatorUpdated = false;
	}

	public void unstore() {
	}

	@Override
	public void restore() {
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

	public byte[] getConsensus(int nr) {
		return consensus[activeIndex[nr]][nr];
	}

	@Override
	public double getEdgeWeights(int nodeNr) {
		return getEdgeMutations(nodeNr);
	}

	
	@Override	
	public double[] getTargetWeights(int fromNodeNr, List<Node> toNodeNrs) {
		double[] distances = new double[toNodeNrs.size()];
		byte[] currConsensus = getConsensus(fromNodeNr);
		
		for (int k = 0; k < toNodeNrs.size(); k++) {
			int nodeNo = toNodeNrs.get(k).getNr();
			byte[] consensus = getConsensus(nodeNo);
			// calculate the distance between the two consensus
			double sum = 0.1;
			for (int l = 0; l < consensus.length; l++) {
//				if (consensus[l] == 1 || currConsensus[l] == 1)
//					continue;
				sum += getDiff(currConsensus[l], consensus[l]);
			}
			distances[k] = 1 / (sum);
		}		
		return distances;
	}
	
	
	
	public double[] getTargetWeightsInteger(int fromNodeNr, List<Integer> toNodeNrs) {
		double[] distances = new double[toNodeNrs.size()];
		byte[] currConsensus = getConsensus(fromNodeNr);
		
		for (int k = 0; k < toNodeNrs.size(); k++) {
			int nodeNo = toNodeNrs.get(k);
			byte[] consensus = getConsensus(nodeNo);
			// calculate the distance between the two consensus
			double sum = 0.1;
			for (int l = 0; l < consensus.length; l++) {
				sum += getDiff(currConsensus[l], consensus[l]);
			}
			distances[k] = 1 / (sum);
		}		
		return distances;
	}


	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
	
	@Override
	public void init(PrintStream out) {
//    	out.print("mutations\t");
//        Node node = treeInput.get().getRoot();
		out.println("#NEXUS\n");
		out.println("Begin trees;");
	}

	@Override
	public void log(long sample, PrintStream out) {
		Tree tree = (Tree) treeInput.get();
		out.print("tree STATE_" + sample + " = ");
		// Don't sort, this can confuse CalculationNodes relying on the tree
		// tree.getRoot().sort();
//        final int[] dummy = new int[1];
		final String newick = toNewick(tree.getRoot());
		out.print(newick);
		out.print(";");
		
		// calculate the total number of mutations
		double totalMutations = 0;
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isRoot())
				continue;
			totalMutations += edgeMutations[activeMutationsIndex[i]][i];
		}
//		System.out.println("Total mutations: " + totalMutations + " number of patters " + patternCount);
//		System.exit(0);

	}
	
	public String getTree() {
		Tree tree = (Tree) treeInput.get();
		return toNewick(tree.getRoot());
	}

	public String toNewick(Node n) {
		final StringBuilder buf = new StringBuilder();
		if (!n.isLeaf()) {
			buf.append("(");
			boolean isFirst = true;
			for (Node child : n.getChildren()) {
				if (isFirst)
					isFirst = false;
				else
					buf.append(",");
				buf.append(toNewick(child));
			}
			buf.append(")");

			if (n.getID() != null)
				buf.append(n.getID());
		} else {
			if (n.getID() != null)
				buf.append(n.getID());

		}

		double diff = 0.0;

		if (!n.isRoot() && !n.getParent().isRoot()) {

			byte[] grandParentConsensus = getConsensus(n.getParent().getParent().getNr());
			byte[] parentConsensus = getConsensus(n.getParent().getNr());
			Node otherChild = n.getParent().getLeft() == n ? n.getParent().getRight() : n.getParent().getLeft();
			byte[] otherChildConsensus = getConsensus(otherChild.getNr());
			// get any mutations from otherChild to parentConsensus, then check if that site
			// mutated to grandParentConsensus
			int doublemuts = 0;
			for (int i = 0; i < grandParentConsensus.length; i++) {
				if (parentConsensus[i] != otherChildConsensus[i]
						&& (otherChildConsensus[i] == grandParentConsensus[i])) {
					doublemuts++;
					break;
				}
			}
			diff = doublemuts;
		}

		// format diff to only use 4 decimals
		String diffStr = String.format("%.1f", diff);
		
//		System.out.println("consensus: " + Arrays.toString(getConsensus(n.getNr())));
//		System.out.println("mutations: " + edgeMutations[activeMutationsIndex[n.getNr()]][n.getNr()]);
		
		if (!n.isRoot()) {
			String str = "null";
			if (mutations[activeMutationsIndex[n.getNr()]][n.getNr()] != null) {
				str = mutations[activeMutationsIndex[n.getNr()]][n.getNr()].toString();
				str = str.replace(" ", "");
				str = str.replace("[", "{");
				str = str.replace("]", "}");
			}			
			buf.append(
					"[&sum=" +  edgeMutations[activeMutationsIndex[n.getNr()]][n.getNr()] + "]");
		} else {
			buf.append("[&muts=" + edgeMutations[activeMutationsIndex[n.getNr()]][n.getNr()] + "]");
		}
		buf.append(":").append(n.getLength());

		return buf.toString();
	}

	/**
	 * @see beast.base.core.Loggable *
	 */
	@Override
	public void close(PrintStream out) {
		out.print("End;");
	}


	@Override
	public double minEdgeWeight() {
		return minWeight;
	}


	@Override
	public byte[] getNodeConsensus(int NodeNo) {		
		return getConsensus(NodeNo);
	}
	
	
	private class Mutation{
	    byte from;
		byte to;
		int postion;
		
		Mutation(byte from, byte to, int postion) {
			this.from = from;
			this.to = to;
			this.postion = postion;
		}
		
		@Override
		public String toString() {
            return from +":" + postion + ":" + to;
        }
		
	}
		
	


} // class TreeLikelihood
