package targetedbeast.likelihood;

import java.io.PrintStream;
import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
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
import beast.base.inference.State;
import targetedbeast.alignment.RapidAlignment;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using "
		+ "a variant of the 'peeling algorithm'. For details, see"
		+ "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class RapidTreeLikelihood extends RapidGenericTreeLikelihood {

	final public Input<Boolean> m_useAmbiguities = new Input<>("useAmbiguities",
			"flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)",
			false);
	final public Input<Boolean> m_useTipLikelihoods = new Input<>("useTipLikelihoods",
			"flag to indicate that partial likelihoods are provided at the tips", false);
	final public Input<String> implementationInput = new Input<>("implementation",
			"name of class that implements this treelikelihood potentially more efficiently. "
					+ "This class will be tried first, with the TreeLikelihood as fallback implementation. "
					+ "When multi-threading, multiple objects can be created.",
			"beast.evolution.likelihood.BeagleTreeLikelihood");

	public static enum Scaling {
		none, always, _default
	};

	final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of "
			+ Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling.none,
			Scaling.values());

	final public Input<Frequencies> rootFrequenciesInput = new Input<>("rootFrequencies",
			"prior state frequencies at root, optional", Input.Validate.OPTIONAL);

	/**
	 * calculation engine *
	 */
	protected RapidLikelihoodCore likelihoodCore;

	public RapidLikelihoodCore getLikelihoodCore() {
		return likelihoodCore;
	}

	/**
	 * BEASTObject associated with inputs. Since none of the inputs are StateNodes,
	 * it is safe to link to them only once, during initAndValidate.
	 */
	protected SubstitutionModel substitutionModel;
	protected SiteModel.Base m_siteModel;
	protected BranchRateModel.Base branchRateModel;

	public SubstitutionModel getSubstitutionModel() {
		return substitutionModel;
	}

	/**
	 * flag to indicate the // when CLEAN=0, nothing needs to be recalculated for
	 * the node // when DIRTY=1 indicates a node partial needs to be recalculated //
	 * when FILTHY=2 indicates the indices for the node need to be recalculated //
	 * (often not necessary while node partial recalculation is required)
	 */
	protected int hasDirt;

	/**
	 * Lengths of the branches in the tree associated with each of the nodes in the
	 * tree through their node numbers. By comparing whether the current branch
	 * length differs from stored branch lengths, it is tested whether a node is
	 * dirty and needs to be recomputed (there may be other reasons as well...).
	 * These lengths take branch rate models in account.
	 */
	protected double[] m_branchLengths;
	protected double[] storedBranchLengths;

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

		int nodeCount = treeInput.get().getNodeCount();
		if (!(siteModelInput.get() instanceof SiteModel.Base)) {
			throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
		}
		m_siteModel = (SiteModel.Base) siteModelInput.get();
		m_siteModel.setDataType(dataInput.get().getDataType());
		substitutionModel = m_siteModel.substModelInput.get();

		if (branchRateModelInput.get() != null) {
			branchRateModel = branchRateModelInput.get();
		} else {
			branchRateModel = new StrictClockModel();
		}
		m_branchLengths = new double[nodeCount];
		storedBranchLengths = new double[nodeCount];

		stateCount = dataInput.get().getMaxStateCount();
		patterns = dataInput.get().getPatternCount();
		likelihoodCore = createLikelihoodCore(stateCount);

		mutations = new int[2][treeInput.get().getNodeCount()][patterns];
		edgeMutations = new double[2][treeInput.get().getNodeCount()];

		calcForPatterns = new int[2][treeInput.get().getNodeCount()][patterns];
//        storedCalcForPatterns = new int[treeInput.get().getNodeCount()][patterns];

		map = new int[2][treeInput.get().getNodeCount()][patterns];
		consensus = new double[2][treeInput.get().getNodeCount()][patterns * stateCount];
//        storeMap = new int[treeInput.get().getNodeCount()][patterns];

		activeIndex = new int[treeInput.get().getNodeCount()];
		storedActiveIndex = new int[treeInput.get().getNodeCount()];

		activeMutationsIndex = new int[treeInput.get().getNodeCount()];
		storedActiveMutationsIndex = new int[treeInput.get().getNodeCount()];

		changed = new boolean[treeInput.get().getNodeCount()];
		changedChildren = new boolean[treeInput.get().getNodeCount()];

		getLeaveMutations(treeInput.get().getRoot());

		String className = getClass().getSimpleName();

		RapidAlignment alignment = dataInput.get();

		Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
		Log.info.println("  " + alignment.toString(true));
		// print startup messages via Log.print*
		
		updateMutations();

		proportionInvariant = m_siteModel.getProportionInvariant();
		m_siteModel.setPropInvariantIsCategory(false);
		if (proportionInvariant > 0) {
			calcConstantPatternIndices(patterns, stateCount);
		}

		initCore();

		patternLogLikelihoods = new double[patterns];
		m_fRootPartials = new double[patterns * stateCount];
		matrixSize = (stateCount + 1) * (stateCount + 1);
		probabilities = new double[(stateCount + 1) * (stateCount + 1)];
		Arrays.fill(probabilities, 1.0);

		if (dataInput.get().isAscertained) {
			useAscertainedSitePatterns = true;
		}
	}

	protected RapidLikelihoodCore createLikelihoodCore(int stateCount) {
		if (stateCount == 4) {
			return new RapidLikelihoodCore4();
		} else {
			throw new IllegalArgumentException("Only stateCount=4 is supported at the moment");
//			return new RapidLikelihoodCore(stateCount);
		}
	}

	/**
	 * Determine indices of m_fRootProbabilities that need to be≠ updates≠ // due to
	 * sites being invariant. If none of the sites are invariant, // the 'site
	 * invariant' category does not contribute anything to the // root probability.
	 * If the site IS invariant for a certain character, // taking ambiguities in
	 * account, there is a contribution of 1 from // the 'site invariant' category.
	 */
	protected void calcConstantPatternIndices(final int patterns, final int stateCount) {
		constantPattern = new ArrayList<>();
		for (int i = 0; i < patterns; i++) {
			final int[] pattern = dataInput.get().getPattern(i);
			final boolean[] isInvariant = new boolean[stateCount];
			Arrays.fill(isInvariant, true);
			for (final int state : pattern) {
				final boolean[] isStateSet = dataInput.get().getStateSet(state);
				if (m_useAmbiguities.get() || !dataInput.get().getDataType().isAmbiguousCode(state)) {
					for (int k = 0; k < stateCount; k++) {
						isInvariant[k] &= isStateSet[k];
					}
				}
			}
			for (int k = 0; k < stateCount; k++) {
				if (isInvariant[k]) {
					constantPattern.add(i * stateCount + k);
				}
			}
		}
	}

	protected void initCore() {
		final int nodeCount = treeInput.get().getNodeCount();
		likelihoodCore.initialize(nodeCount, dataInput.get().getPatternCount(), m_siteModel.getCategoryCount(), true,
				m_useAmbiguities.get());

		final int extNodeCount = nodeCount / 2 + 1;
		final int intNodeCount = nodeCount / 2;

		if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
			setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
			throw new UnsupportedOperationException("Ambiguities and tip likelihoods not supported");
		} else {
			setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
		}
		hasDirt = Tree.IS_FILTHY;
		for (int i = 0; i < intNodeCount; i++) {
			likelihoodCore.createNodePartials(extNodeCount + i);
		}
	}

	/**
	 * This method samples the sequences based on the tree and site model.
	 */
	@Override
	public void sample(State state, Random random) {
		throw new UnsupportedOperationException("Can't sample a fixed alignment!");
	}

	/**
	 * set leaf states in likelihood core *
	 */
	protected void setStates(Node node, int patternCount) {
		if (node.isLeaf()) {
			int[] states = new int[dataInput.get().getMaxStateCount()];
			for (int i = 0; i < dataInput.get().getMaxStateCount(); i++) {
				states[i] = i;
			}
			likelihoodCore.setNodeStates(node.getNr(), states);
		} else {
			setStates(node.getLeft(), patternCount);
			setStates(node.getRight(), patternCount);
		}
	}

	/**
	 *
	 * @param taxon the taxon name as a string
	 * @param data  the alignment
	 * @return the taxon index of the given taxon name for accessing its sequence
	 *         data in the given alignment, or -1 if the taxon is not in the
	 *         alignment.
	 */
	private int getTaxonIndex(String taxon, RapidAlignment data) {
		int taxonIndex = data.getTaxonIndex(taxon);
		if (taxonIndex == -1) {
			if (taxon.startsWith("'") || taxon.startsWith("\"")) {
				taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
			}
			if (taxonIndex == -1) {
				throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
			}
		}
		return taxonIndex;
	}

	/**
	 * set leaf partials in likelihood core *
	 */
	protected void setPartials(Node node, int patternCount) {
		if (node.isLeaf()) {
			RapidAlignment data = dataInput.get();
			int states = data.getDataType().getStateCount();
			double[] partials = new double[patternCount * states];
			int k = 0;
			int taxonIndex = getTaxonIndex(node.getID(), data);
			for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {
				double[] tipLikelihoods = data.getTipLikelihoods(taxonIndex, patternIndex_);
				if (tipLikelihoods != null) {
					for (int state = 0; state < states; state++) {
						partials[k++] = tipLikelihoods[state];
					}
				} else {
					int stateCount = data.getPattern(taxonIndex, patternIndex_);
					boolean[] stateSet = data.getStateSet(stateCount);
					for (int state = 0; state < states; state++) {
						partials[k++] = (stateSet[state] ? 1.0 : 0.0);
					}
				}
			}
			likelihoodCore.setNodePartials(node.getNr(), partials);

		} else {
			setPartials(node.getLeft(), patternCount);
			setPartials(node.getRight(), patternCount);
		}
	}

	// for testing
	public double[] getRootPartials() {
		return m_fRootPartials.clone();
	}

	/**
	 * Calculate the log likelihood of the current state.
	 *
	 * @return the log likelihood.
	 */
	double m_fScale = 1.01;
	int m_nScale = 0;
	int X = 100;

	double ratio = 0;
	int nrratio = 0;
	int nrratio2 = 0;
	double totrati=0;
	
	public void updateMutations() {
		Arrays.fill(changed, true);
		Arrays.fill(changedChildren, true);
		getMutations(treeInput.get().getRoot());
	}

	@Override
	public double calculateLogP() {
		final TreeInterface tree = treeInput.get();

		if (!operatorUpdated) {
			Arrays.fill(changed, false);
			Arrays.fill(changedChildren, false);
			getFilthyNodes(tree.getRoot());
			getMutations(tree.getRoot());
		}

		operatorUpdated = false;
		ratio=0;
		nrratio=0;
		
		

		try {
			if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
				calcLogP();
		} catch (ArithmeticException e) {
			return Double.NEGATIVE_INFINITY;
		}
		
		if (count%1000==0) {
			System.out.println("Ratio: " + countwith/(0.0+countwithout));
			countwith=0;
			countwithout=0;
		}
		
		count++;
		
		
//		totrati+=(ratio/nrratio)/calcForPatterns[0][0].length;
//		nrratio2++;
//		System.out.println("Ratio: " + totrati/nrratio2 + " " + (ratio/nrratio)/calcForPatterns[0][0].length);

		m_nScale++;
		if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
			throw new RuntimeException("LogP is positive, scaling is needed but not implementes");
		} else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // &&
																												// !m_likelihoodCore.getUseScaling())
																												// {
			m_nScale = 0;
			m_fScale *= 1.01;
			Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
			likelihoodCore.setUseScaling(m_fScale);
			likelihoodCore.unstore();
			hasDirt = Tree.IS_FILTHY;
			traverse(tree.getRoot());
			calcLogP();
			return logP;
		}
		return logP;
	}

	protected void calcLogP() {
		logP = 0.0;
		if (useAscertainedSitePatterns) {
			final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
			for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
				logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
			}
		} else {
			for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
				logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
			}
		}
	}

	/* Assumes there IS a branch rate model as opposed to traverse() */
	protected int traverse(final Node node) {

		int update = (node.isDirty() | hasDirt);

		final int nodeIndex = node.getNr();

		final double branchRate = branchRateModel.getRateForBranch(node);
		final double branchTime = node.getLength() * branchRate;

		// First update the transition probability matrix(ices) for this branch
		// if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime !=
		// m_StoredBranchLengths[nodeIndex])) {
		if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
			m_branchLengths[nodeIndex] = branchTime;
			final Node parent = node.getParent();
			likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
			for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
				final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
				substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(),
						jointBranchRate, probabilities);
				// System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
				likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
			}
			update |= Tree.IS_DIRTY;
		}

		// If the node is internal, update the partial likelihoods.
		if (!node.isLeaf()) {

			// Traverse down the two child nodes
			final Node child1 = node.getLeft(); // Two children
			final int update1 = traverse(child1);

			final Node child2 = node.getRight();
			final int update2 = traverse(child2);

			// If either child node was updated then update this node too
			if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

				final int childNum1 = child1.getNr();
				final int childNum2 = child2.getNr();

				likelihoodCore.setNodePartialsForUpdate(nodeIndex);
				update |= (update1 | update2);
				if (update >= Tree.IS_FILTHY) {
					likelihoodCore.setNodeStatesForUpdate(nodeIndex);
				}

				if (m_siteModel.integrateAcrossCategories()) {
					try {
						likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex,
								mutations[activeMutationsIndex[childNum1]][childNum1],
								mutations[activeMutationsIndex[childNum2]][childNum2],
								calcForPatterns[activeIndex[childNum1]][childNum1],
								calcForPatterns[activeIndex[childNum2]][childNum2],
								calcForPatterns[activeIndex[nodeIndex]][nodeIndex]);
						
						countwithout+=calcForPatterns[0][0].length;
						countwith +=4;
						for (int i = 0; i < calcForPatterns[activeIndex[nodeIndex]][nodeIndex].length; i++) {
							if (calcForPatterns[activeIndex[nodeIndex]][nodeIndex][i] == -1)
								break;
							countwith++;
						}					
						
					} catch (Exception e) {
						System.out.println(childNum1 + " " + childNum2 + " " + nodeIndex);
						System.exit(0);
					}
				} else {
					throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
					// m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum,
					// siteCategories);
				}

				if (node.isRoot()) {
					// No parent this is the root of the beast.tree -
					// calculate the pattern likelihoods

					final double[] proportions = m_siteModel.getCategoryProportions(node);
					likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

					if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
						proportionInvariant = m_siteModel.getProportionInvariant();
						// some portion of sites is invariant, so adjust root partials for this
						for (final int i : constantPattern) {
							m_fRootPartials[i] += proportionInvariant;
						}
					}

					double[] rootFrequencies = substitutionModel.getFrequencies();
					if (rootFrequenciesInput.get() != null) {
						rootFrequencies = rootFrequenciesInput.get().getFreqs();
					}
					likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
				}

			}
		}
		return update;
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

	/*
	 * return copy of pattern log likelihoods for each of the patterns in the
	 * alignment
	 */
	public double[] getPatternLogLikelihoods() {
		return patternLogLikelihoods.clone();
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
		if (m_siteModel.isDirtyCalculation()) {
			hasDirt = Tree.IS_DIRTY;
			return true;
		}
		if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
			// m_nHasDirt = Tree.IS_DIRTY;
			return true;
		}
		if (rootFrequenciesInput.get() != null && rootFrequenciesInput.get().isDirtyCalculation()) {
			hasDirt = Tree.IS_DIRTY;
			return true;
		}
		return treeInput.get().somethingIsDirty();
	}

	@Override
	public void store() {
//    	System.err.println("store ");
		super.store();
		if (!operatorUpdated) { // avoid storing again if the operator has already done it
			if (likelihoodCore != null) {
				likelihoodCore.store();
			}
			System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);

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
		if (likelihoodCore != null) {
			likelihoodCore.store();
		}
		System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
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

		if (likelihoodCore != null) {
			likelihoodCore.restore();
		}

		super.restore();

		double[] tmp = m_branchLengths;
		m_branchLengths = storedBranchLengths;
		storedBranchLengths = tmp;
		System.arraycopy(storedActiveIndex, 0, activeIndex, 0, activeIndex.length);
		System.arraycopy(storedActiveMutationsIndex, 0, activeMutationsIndex, 0, activeMutationsIndex.length);
	}

	/**
	 * @return a list of unique ids for the state nodes that form the argument
	 */
	@Override
	public List<String> getArguments() {
		return Collections.singletonList(dataInput.get().getID());
	}

	/**
	 * @return a list of unique ids for the state nodes that make up the conditions
	 */
	@Override
	public List<String> getConditions() {
		return m_siteModel.getConditions();
	}

	public double getEdgeMutations(int i) {
//		System.err.println("getEdgeMutations " +i + " " + activeMutationsIndex[i] + " " + edgeMutations[activeMutationsIndex[i]][i]);
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

//	@Override
//	public void init(PrintStream out) {
////    	out.print("mutations\t");
////        Node node = treeInput.get().getRoot();
//		out.println("#NEXUS\n");
//		out.println("Begin trees;");
//	}
//
//	@Override
//	public void log(long sample, PrintStream out) {
//		double totalMut = 0;
//		for (int i = 0; i < edgeMutations[0].length; i++) {
//			totalMut += edgeMutations[activeMutationsIndex[i]][i] - 0.01;
//		}
//		avg_muts = totalMut;
////		out.print(totalMut + "\t");	
//
//		Tree tree = (Tree) treeInput.get();
//		out.print("tree STATE_" + sample + " = ");
//		// Don't sort, this can confuse CalculationNodes relying on the tree
//		// tree.getRoot().sort();
////        final int[] dummy = new int[1];
//		final String newick = getTree();
//		out.print(newick);
//		out.print(";");
//	}

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

			double[] grandParentConsensus = getConsensus(n.getParent().getParent().getNr());
			double[] parentConsensus = getConsensus(n.getParent().getNr());
			Node otherChild = n.getParent().getLeft() == n ? n.getParent().getRight() : n.getParent().getLeft();
			double[] otherChildConsensus = getConsensus(otherChild.getNr());
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

		double sum = 0;

		// format diff to only use 4 decimals
		String diffStr = String.format("%.1f", diff);
		double total_muts = 0;
		if (n.isLeaf() || n.isRoot()) {

		} else {
			total_muts = edgeMutations[activeMutationsIndex[n.getNr()]][n.getNr()]
					+ edgeMutations[activeMutationsIndex[n.getLeft().getNr()]][n.getLeft().getNr()]
					+ edgeMutations[activeMutationsIndex[n.getRight().getNr()]][n.getRight().getNr()];
		}
		if (!n.isRoot() && !n.getParent().isRoot()) {
			buf.append(
					"[&diffs=" + diffStr + ", muts=" + edgeMutations[activeMutationsIndex[n.getNr()]][n.getNr()] + "]");
		} else {
			buf.append("[&muts=" + edgeMutations[activeMutationsIndex[n.getNr()]][n.getNr()] + "]");
		}
		buf.append(":").append(n.getLength());

		return buf.toString();
	}

//	/**
//	 * @see beast.base.core.Loggable *
//	 */
//	@Override
//	public void close(PrintStream out) {
//		out.print("End;");
//	}

	double avg_muts;
	List<Double> prev_heights;

	public String getTree() {
		double totalMut = 0;
		double totalLength = 0;

		for (int i = 0; i < edgeMutations[0].length; i++) {
			totalMut += edgeMutations[activeMutationsIndex[i]][i] - 0.01;
			totalLength += treeInput.get().getNode(i).getLength();
		}
		avg_muts = totalMut / totalLength;

		// compute the standard deviation of the avg muts
		double deviation = 0;
		for (int i = 0; i < edgeMutations[0].length; i++) {
			if (treeInput.get().getNode(i).isLeaf()) {
				continue;
			}
			int left = treeInput.get().getNode(i).getLeft().getNr();
			int right = treeInput.get().getNode(i).getRight().getNr();
			double thisRate = edgeMutations[activeMutationsIndex[left]][left]
					+ edgeMutations[activeMutationsIndex[right]][right];
			double thisLength = treeInput.get().getNode(left).getLength() + treeInput.get().getNode(right).getLength();
			deviation += Math.abs(thisRate / thisLength - avg_muts);
		}
		double sd = deviation / edgeMutations[0].length;
//		System.out.println(totalMut);
//		System.out.println("avg muts = " + avg_muts + " sd " + sd + " clockRate = " + (branchRateModel.getRateForBranch(treeInput.get().getNode(0))*dataInput.get().getSiteCount()));

		return toNewick(treeInput.get().getRoot());
	}


} // class TreeLikelihood
