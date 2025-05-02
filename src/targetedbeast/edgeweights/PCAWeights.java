package targetedbeast.edgeweights;

import java.io.PrintStream;
import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import targetedbeast.util.Alignment2PCA;

@Description("Keeps track of the distances in PCA space. "
		+ "PCAWeights is a distribution to ensure that it is updated correctly")
public class PCAWeights extends Distribution implements EdgeWeights, Loggable {
	
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);
    
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    
    final public Input<Double> maxWeightInput = new Input<>("maxWeight", "maximum weight for an edge", 10.0);
    
    final public Input<Double> minWeightInput = new Input<>("minWeight", "maximum weight for an edge", 0.01);

	final public Input<Integer> dimensionInput = new Input<>("dimension", "dimension of PCA points", 2);

    final public Input<String> valueInput = new Input<>("value", "comma separated list of taxon=x1 x2 x3 pairs, "
    		+ "where <dimension> number of dimensions are specified, e.g. taxon1=0.2 0.4, taxon2=0.4 0.3, taxon3=0.9 0.1. "
    		+ "If value is specified, data is ignored");

	final public Input<Boolean> distanceBasedInput = new Input<>("distanceBased", "flag to indicate PCA is done based on distance matrix (if true) or normalised alignment (if false)", true);

	protected int hasDirt;

//	private byte[][][] consensus;
	// 2 x nr of taxa x dimensions
	private double[][][] points;

	private boolean[] changed;
	private boolean[] changedChildren;

	private int[] activeIndex;
	private int[] storedActiveIndex;

	private int[] activeMutationsIndex;
	private int[] storedActiveMutationsIndex;
	
	
	// for every leaf node, generate a sequence number
	// leaf nodes with the same sequence have the same number
	private int [] sequenceID;

	public double[][] edgeMutations;
	
	public List<Mutation>[][] mutations;

	private boolean operatorUpdated = false;

	int stateCount;
	int patternCount;	
	int maxStateCount;
	
	double maxWeight;
	double minWeight;
	
	double totalMuts[];

	private double [] diff;
	private int range;
	final static byte UNKNOWN = 0xf;
	private int dim;

	
	@Override
	public void initAndValidate() {
		dim = dimensionInput.get();
		
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
		points = new double[2][treeInput.get().getNodeCount()][dimensionInput.get()];
		if (valueInput.get() != null) {
			parseValue();
		} else {
			calcValue();
		}

		activeIndex = new int[treeInput.get().getNodeCount()];
		storedActiveIndex = new int[treeInput.get().getNodeCount()];

		activeMutationsIndex = new int[treeInput.get().getNodeCount()];
		storedActiveMutationsIndex = new int[treeInput.get().getNodeCount()];
		
		mutations = new ArrayList[2][treeInput.get().getNodeCount()];
		for (int i = 0; i < treeInput.get().getNodeCount(); i++) {
			mutations[0][i] = new ArrayList<>();
			mutations[1][i] = new ArrayList<>();
		}

		changed = new boolean[treeInput.get().getNodeCount()];
		changedChildren = new boolean[treeInput.get().getNodeCount()];
			
		maxWeight = maxWeightInput.get();
		minWeight = minWeightInput.get();
		
		totalMuts = new double[patternCount];

		
		initDiff();
		updateWeights();
		initSequenceID();
	}

	
	private void calcValue() {
		Map<String, double[]> map = Alignment2PCA.getPoints(dataInput.get(), dimensionInput.get(), distanceBasedInput.get());

		int taxonCount = treeInput.get().getLeafNodeCount();
		int dim = dimensionInput.get();
		List<String> taxa = treeInput.get().getTaxonset().asStringList();
		for (int i = 0; i < taxonCount; i++) {
			String taxon = taxa.get(i);
			double [] p = map.get(taxon);
			for (int j = 0; j < dim; j++) {
				points[0][i][j] = p[j];
			}
		}
	}


	private void parseValue() {
		String [] strs0 = valueInput.get().split(",");
        List<String> labels = treeInput.get().getTaxonset().asStringList();
        int dimension = dimensionInput.get();
		if (strs0.length != labels.size()) {
			Log.warning("Number of points specified (" + strs0.length +" should equal number of taxa (" + labels.size()+")");
		}
		
		for  (int i = 0; i < strs0.length; i++) {
			String trait = strs0[i];
            trait = trait.replaceAll("\\s+", " ");
            String[] strs = trait.split("=");
            if (strs.length != 2) {
                throw new IllegalArgumentException("could not parse trait: " + trait);
            }
            String taxonID = normalize(strs[0]);
            int taxonNr = labels.indexOf(taxonID);
            
            String [] point = strs[1].split("\\s+");
            if (point.length != dimension) {
                throw new IllegalArgumentException("could not parse trait: " + trait + " since dimension is not " + dimension);
            }
            for (int j = 0; j < dimension; j++) {
            	points[0][taxonNr][j] = Double.parseDouble(point[j]);
            }
		}
	}

	/**
     * remove start and end spaces
     */
    protected String normalize(String str) {
        if (str.charAt(0) == ' ') {
            str = str.substring(1);
        }
        if (str.endsWith(" ")) {
            str = str.substring(0, str.length() - 1);
        }
        return str;
    }


	private void initSequenceID() {
		TreeInterface tree = treeInput.get();
		Alignment data = dataInput.get();
		
		sequenceID = new int[tree.getNodeCount()];
		Map<String,Integer> sequenceMap = new HashMap<>();
		int k = 0;
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			String taxon = tree.getNode(i).getID();
			String seq = data.getSequenceAsString(taxon); 
			if (!sequenceMap.containsKey(seq)) {
				sequenceMap.put(seq, k++);
			}
			sequenceID[i] = sequenceMap.get(seq);
		}
		for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
			sequenceID[i] = i;			
		}
	}


	// diff is 0 if states are the same
	// diff is 1 of there is no overlap between states
	// diff is 0.5 if there is overlap in states, but states differ
	// diff is 0 if one of the states is fully ambiguous 
	private void initDiff() {
		int n = (int) Math.pow(2,stateCount);
		range = n;
		diff = new double[n*n];
		// only go to n-1, since state n-1 represents UNKNOWN states
		for (int i = 0; i < n - 1; i++) {
			for (int j = 0; j < i; j++) {
				int intersection = i & j;
				if (intersection == 0) {
					// no overlap between i and j
					diff[i*n+j] = 1;
				} if (intersection == i) {
					// diff from i to j may not need a mutation
					diff[i*n+j] = 0.25;
				} else {
					// some overlap between i and j
					diff[i*n+j] = 0.5;
				}
			}
			for (int j = i+1; j < n - 1; j++) {
				int intersection = i & j;
				if (intersection == 0) {
					// no overlap between i and j
					diff[i*n+j] = 1;
				} if (intersection == i) {
					// diff from i to j may not need a mutation
					diff[i*n+j] = 0.25;
				} else {
					// some overlap between i and j
					diff[i*n+j] = 0.5;
				}
			}
		}
	}


	private void updateWeights() {
		Arrays.fill(changed, true);
		Arrays.fill(changedChildren, true);
		getNodeConsensusSequences(treeInput.get().getRoot());
	}

	@Override
	public void updateByOperator() {
		operatorUpdated = true;
		Arrays.fill(changed, false);
		Arrays.fill(changedChildren, false);
		getFilthyNodes(treeInput.get().getRoot());
		getNodeConsensusSequences(treeInput.get().getRoot());
	}

	@Override
	public void updateByOperatorWithoutNode(int ignore, List<Integer> nodes) {
		Arrays.fill(changed, false);
		Arrays.fill(changedChildren, false);
		for (Integer nodeNo : nodes) {
			changed[nodeNo] = true;
		}
		getConsensusWithoutNode(treeInput.get().getRoot(), ignore);
	}

	@Override
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

	private void getNodeConsensusSequences(Node n) {
		if (n.isLeaf()) {
			// the active index for leaves is always 0, could be changed to make the arrays
			// shorter
			return;
		} else {
			final int nodeNr = n.getNr();
			if (changed[nodeNr]) {
				// compare the patterns of the two lineages
				getNodeConsensusSequences(n.getLeft());
				getNodeConsensusSequences(n.getRight());
				
				final int leftNr = n.getLeft().getNr();
				final int rightNr = n.getRight().getNr();

				// set this node index to the active index
				activeIndex[nodeNr] = 1 - activeIndex[nodeNr];
				activeMutationsIndex[leftNr] = 1 - activeMutationsIndex[leftNr];
				activeMutationsIndex[rightNr] = 1 - activeMutationsIndex[rightNr];

				int activeInd = activeIndex[nodeNr];
				int activeIndLeft = activeIndex[leftNr];
				int activeIndRight = activeIndex[rightNr];
				double left, right, val;
				
				mutations[activeIndLeft][leftNr].clear();// = new ArrayList<>();
				mutations[activeIndRight][rightNr].clear();// = new ArrayList<>();
				
				double sumLeft = 0;
				double sumRight = 0;
				
				final double[] leftconsensus = points[activeIndLeft][leftNr]; 
				final double[] rightconsensus = points[activeIndRight][rightNr]; 
				final double [] currentconsensus = points[activeInd][nodeNr];

				for (int i = 0; i < dim; i++) {
					
					left = leftconsensus[i];
					right = rightconsensus[i];
					val = (left + right) / 2.0; 
										
					sumLeft += Math.abs(left - val);
					sumRight += Math.abs(right - val);
					currentconsensus[i] = val;
				}
				edgeMutations[activeMutationsIndex[leftNr]][leftNr] = Math.min(maxWeight, sumLeft);
				edgeMutations[activeMutationsIndex[rightNr]][rightNr] = Math.min(maxWeight, sumRight);				
			} else {
				getNodeConsensusSequences(n.getLeft());
				getNodeConsensusSequences(n.getRight());
			}
		}
	}

//	private double getDiffX(byte child, byte parent) {
//		if (parent == child) {
//			return 0;
//		} else if (child == 0xf) {
//			return 0;
//		} else if (child < 4 && parent <4) {
//			return 1;
//        } else {
//        	return 0.5;
//        }		
//	}
//
//
//	private byte getCombination(byte left, byte right) {
//		return (byte) (left | right);
//	}

	private void getConsensusWithoutNode(Node n, int ignore) {
		if (n.isLeaf()) {
			return;
		} else {
			final int nodeNr = n.getNr();
			if (changed[nodeNr]) {

				// set this node index to the active index
				activeIndex[nodeNr] = 1 - activeIndex[nodeNr];

				final int leftNr = n.getLeft().getNr();
				final int rightNr = n.getRight().getNr();

				int activeInd = activeIndex[nodeNr];
				int activeIndLeft = activeIndex[leftNr];
				int activeIndRight = activeIndex[rightNr];
				double left, right, val;

				
				final double[] leftconsensus = points[activeIndLeft][leftNr]; 
				final double[] rightconsensus = points[activeIndRight][rightNr];
				final double [] currentconsensus = points[activeInd][nodeNr];
				
				if (leftNr != ignore && rightNr != ignore) {
					for (int i = 0; i < dim; i++) {
						left = leftconsensus[i];
						right = rightconsensus[i];
											
						val = (left + right) / 2.0;	
						currentconsensus[i] = val;
					}
				} else if (leftNr == ignore) {
					System.arraycopy(rightconsensus, 0,
							currentconsensus, 0, dim); 
				} else if (rightNr == ignore) {
					System.arraycopy(leftconsensus, 0,
							currentconsensus, 0, dim); 
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

	@Override
	public void prestore() {
		System.arraycopy(activeIndex, 0, storedActiveIndex, 0, activeIndex.length);
		System.arraycopy(activeMutationsIndex, 0, storedActiveMutationsIndex, 0, activeMutationsIndex.length);
	}

	@Override
	public void reset() {
		// undoes any previous calculation
		System.arraycopy(storedActiveIndex, 0, activeIndex, 0, activeIndex.length);
		System.arraycopy(storedActiveMutationsIndex, 0, activeMutationsIndex, 0, activeMutationsIndex.length);
		operatorUpdated = false;
	}

//	public void unstore() {
//	}

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

//	public byte[] getConsensus(int nr) {
//		return consensus[activeIndex[nr]][nr];
//	}


	public double[] getPoints(int nr) {
		return points[activeIndex[nr]][nr];
	}
	
	@Override
	public double getEdgeWeights(int nodeNr) {
		return getEdgeMutations(nodeNr);
	}

	
	@Override	
	public double[] getTargetWeights(int fromNodeNr, List<Node> toNodeNrs) {
		double[] distances = new double[toNodeNrs.size()];
		double[] currConsensus = getPoints(fromNodeNr);
		
		for (int k = 0; k < toNodeNrs.size(); k++) {
			int nodeNo = toNodeNrs.get(k).getNr();
			if (sequenceID[nodeNo] == sequenceID[fromNodeNr]) {
				distances[k] = 0;				
			} else {
				double[] consensus = getPoints(nodeNo);
				// calculate the distance between the two consensus
				double sum = 0.1;
				for (int l = 0; l < dim; l++) {
					sum += Math.abs(currConsensus[l] - consensus[l]);
				}
				distances[k] = 1 / (sum);
			}
		}
		return distances;
	}
	
	@Override
	public double[] getTargetWeightsInteger(int fromNodeNr, List<Integer> toNodeNrs) {
		double[] distances = new double[toNodeNrs.size()];
		double[] currConsensus = getPoints(fromNodeNr);
		
		for (int k = 0; k < toNodeNrs.size(); k++) {
			int nodeNo = toNodeNrs.get(k);
			double sum = 0.1;
			if (sequenceID[nodeNo] != sequenceID[fromNodeNr]) {
				double[] consensus = getPoints(nodeNo);
				// calculate the distance between the two consensus
				for (int l = 0; l < dim; l++) {
					sum += Math.abs(currConsensus[l] - consensus[l]);
				}
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

//			byte[] grandParentConsensus = getConsensus(n.getParent().getParent().getNr());
//			byte[] parentConsensus = getConsensus(n.getParent().getNr());
//			Node otherChild = n.getParent().getLeft() == n ? n.getParent().getRight() : n.getParent().getLeft();
//			byte[] otherChildConsensus = getConsensus(otherChild.getNr());
//			// get any mutations from otherChild to parentConsensus, then check if that site
//			// mutated to grandParentConsensus
//			int doublemuts = 0;
//			for (int i = 0; i < grandParentConsensus.length; i++) {
//				if (parentConsensus[i] != otherChildConsensus[i]
//						&& (otherChildConsensus[i] == grandParentConsensus[i])) {
//					doublemuts++;
//					break;
//				}
//			}
//			diff = doublemuts;
		}

		// format diff to only use 4 decimals
		String diffStr = String.format("%.1f", diff);
		
//		System.out.println("consensus: " + Arrays.toString(getConsensus(nodeNr)));
//		System.out.println("mutations: " + edgeMutations[activeMutationsIndex[nodeNr]][nodeNr]);
		
		final int nodeNr = n.getNr();
		if (!n.isRoot()) {
			String str = null;
			if (mutations[activeMutationsIndex[nodeNr]][nodeNr] != null) {
				str = mutations[activeMutationsIndex[nodeNr]][nodeNr].toString();
				str = str.replace(" ", "");
				str = str.replace("[", "{");
				str = str.replace("]", "}");
			}			
			buf.append(
					"[&sum=" +  edgeMutations[activeMutationsIndex[nodeNr]][nodeNr]);
			if (str!=null) {
				buf.append(",muts=" + str);
			}
			buf.append("]");
		} else {
			buf.append("[&sum=" + edgeMutations[activeMutationsIndex[nodeNr]][nodeNr] + "]");
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


//	@Override
//	public byte[] getNodeConsensus(int NodeNo) {		
//		return getConsensus(NodeNo);
//	}
	
	
	private class Mutation {
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
