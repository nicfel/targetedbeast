
package targetedbeast.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import targetedbeast.edgeweights.ConsensusWeights;
import targetedbeast.edgeweights.EdgeWeights;
import targetedbeast.likelihood.RapidTreeLikelihood;

@Description("Performs narrow move on nodes that are close together in terms of the number of mutations between their consensus sequences")
public class WeightBasedNodeRandomizer extends TreeOperator {
	public final Input<Double> percentageInput = new Input<>("percentage",
			"percentage of nodes below limit to pick", 0.25);
    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights", "input of weights to be used for targetedn tree operations", Input.Validate.REQUIRED);
    final public Input<Boolean> optimiseInput = new Input<>("optimise",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			true);
    public Input<RealParameter> ratesInput = new Input<>("rates",
			"if provided, branch rates are co-scaled to keep genetic distance (rate x time) constant across the swap (ORC constant-distance NNI)");

    EdgeWeights edgeWeights;
    RealParameter branchRates;
    double limit;
    double percentage;
    double lower;

	@Override
	public void initAndValidate() {
		edgeWeights = edgeWeightsInput.get();
		branchRates = ratesInput.get();
		limit = edgeWeights.minEdgeWeight();
		percentage = percentageInput.get();
		lower = 1.0/treeInput.get().getInternalNodeCount();
	}

	/**
	 * override this for proposals,
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted *
	 */
	@Override
	public double proposal() {
		final Tree tree = treeInput.get();
//		System.out.println(tree + ";");

		double logHastingsRatio = 0;
		
		int numberOfAttempts = (int) Math.ceil(tree.getNodeCount() * percentage);
		
//		numberOfAttempts = Math.max(numberOfAttempts, 2);
		
		
		int totalCandidates = 0;
		boolean[] isCandidate = new boolean[tree.getNodeCount()];
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (isEligible(tree.getNode(i))) {
				isCandidate[i] = true;
				totalCandidates++;
			}
		}


		if (totalCandidates == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		// look for groups of nodes that have less than 0.1 mutations between them
		for (int k = 0; k < numberOfAttempts; k++) {
			
			double scaler = Randomizer.nextDouble() * totalCandidates;
			int nodeNr = -1;
			int currDev = 0;
			for (int j = 0; j < isCandidate.length; j++) {
				if (!isCandidate[j])
					continue;
				currDev++;
				
				if (currDev > scaler) {
					nodeNr = j;
					break;
				}
			}
			
			Node j = tree.getNode(nodeNr);
			 
			

			Node parent = j.getParent();
			Node sibling = parent.getLeft() == j ? parent.getRight() : parent.getLeft();
			Node grandParent = parent.getParent();
			Node uncle = grandParent.getLeft() == parent ? grandParent.getRight() : grandParent.getLeft();
			

			final int internalNodes = tree.getInternalNodeCount();
			if (internalNodes <= 1) {
				return Double.NEGATIVE_INFINITY;
			}

			
//			System.out.println(sibling.getNr() + " " + uncle.getNr() + " " + parent.getNr() + " " + grandParent.getNr());
//			System.out.println(tree +";");
			exchangeNodes(sibling, uncle, parent, grandParent);
//			System.out.println(tree +";");
			grandParent.makeDirty(3 - grandParent.isDirty());

			// ORC constant-distance rate co-scaling: heights are unchanged by the NNI, so only the
			// sibling (moves parent->grandParent) and uncle (moves grandParent->parent) branches change
			// duration. Rescale their rates to keep rate*time constant, adding the Jacobian to the HR.
			if (branchRates != null) {
				double tOldSib = parent.getHeight() - sibling.getHeight();      // sibling was below parent
				double tNewSib = grandParent.getHeight() - sibling.getHeight(); // now below grandParent
				double tOldUnc = grandParent.getHeight() - uncle.getHeight();   // uncle was below grandParent
				double tNewUnc = parent.getHeight() - uncle.getHeight();        // now below parent
				if (tNewSib <= 0 || tNewUnc <= 0) {
					return Double.NEGATIVE_INFINITY;
				}
				branchRates.setValue(sibling.getNr(), branchRates.getValue(sibling.getNr()) * tOldSib / tNewSib);
				branchRates.setValue(uncle.getNr(),   branchRates.getValue(uncle.getNr())   * tOldUnc / tNewUnc);
				logHastingsRatio += Math.log(tOldSib / tNewSib) + Math.log(tOldUnc / tNewUnc);
			}
			// calculate reverse contribution
			// check if the parentIndex.getNr() changes the isCandidate status
			int reverseTotalCandidates = totalCandidates;
			boolean siblingNowIsCandidate = false;
			boolean uncleNowIsCandidate = false;
//			System.out.println(sibling.getNr() + " " + uncle.getNr() + " " + parent.getNr() + " " + grandParent.getNr());
//			
//			System.exit(0);
			
			List<Node> nodesToRecheckEligibility = new ArrayList<>();
			nodesToRecheckEligibility.add(sibling);
			nodesToRecheckEligibility.add(uncle);
			if (!sibling.isLeaf()) {
				nodesToRecheckEligibility.add(sibling.getLeft());
				nodesToRecheckEligibility.add(sibling.getRight());
			}
			if (!uncle.isLeaf()) {
				nodesToRecheckEligibility.add(uncle.getLeft());
				nodesToRecheckEligibility.add(uncle.getRight());
			}
			if (!j.isLeaf()) {
				nodesToRecheckEligibility.add(j.getLeft());
				nodesToRecheckEligibility.add(j.getRight());
			}
			
			for (Node n : nodesToRecheckEligibility) {
				boolean isCandidateNow = isEligible(n);

				if (!isCandidate[n.getNr()] && isCandidateNow) {
					reverseTotalCandidates++;
					isCandidate[n.getNr()] = true;
				} else if (isCandidate[n.getNr()] && !isCandidateNow) {
					reverseTotalCandidates--;			
					isCandidate[n.getNr()] = false;
				}

			}
			
			
			logHastingsRatio += Math.log((float) totalCandidates / reverseTotalCandidates);
			totalCandidates = reverseTotalCandidates;

		}
//		System.out.println(tree + ";");
		return logHastingsRatio;
	}

	// Eligibility for the sibling<->uncle swap. MUST be identical in the forward scan and the
	// reverse recheck, otherwise reverseTotalCandidates diverges from |C(y)| and detailed balance
	// breaks (e.g. an all-N leaf has edge weight 0 and would be counted as a candidate the forward
	// scan never picks).
	private boolean isEligible(final Node n) {
		if (n.isLeaf() || n.isRoot() || n.getParent().isRoot()) {
			return false;
		}
		final Node parent = n.getParent();
		final Node grandParent = parent.getParent();
		final Node uncle = grandParent.getLeft() == parent ? grandParent.getRight() : grandParent.getLeft();
		if (uncle.getHeight() > parent.getHeight()) {
			return false;
		}
		return edgeWeights.getEdgeWeights(n.getNr()) <= limit
				&& edgeWeights.getEdgeWeights(parent.getNr()) <= limit;
	}

	private int isg(final Node n) {
		return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
	}

	private int sisg(final Node n) {
		return n.isLeaf() ? 0 : isg(n);
	}

	public double newNarrow(final Tree tree) {

		double limit = edgeWeights.minEdgeWeight();

//		System.out.println(((ConsensusWeights) edgeWeights).getTree() + ";");

		double logq = 0;
		// choose a random node avoiding root
		double totalDeviation = 0;
		double[] deviation = new double[tree.getNodeCount()];
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isLeaf())
				continue;
			
			Node left = tree.getNode(i).getLeft();
			Node right = tree.getNode(i).getRight();
			
			if (left.getHeight() < right.getHeight()) {
				Node tmp = left;
				left = right;
				right = tmp;
			}
			if (left.isLeaf())
				continue;

			// only pick nodes that have less than 0.01 mutations on left and at least one
			// child has less than 0.01 mutations
			if (edgeWeights.getEdgeWeights(left.getNr()) <= limit
					&& (edgeWeights.getEdgeWeights(right.getNr()) <= limit)) {
				// check the child node with the the higher
				deviation[i] = 1;
				totalDeviation += deviation[i];
			}
		}
		if (totalDeviation == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		double scaler = Randomizer.nextDouble() * totalDeviation;
		int nodeNr = -1;
		double currDev = 0;
		for (int j = 0; j < deviation.length; j++) {
			currDev += deviation[j];
			if (currDev > scaler) {
				nodeNr = j;
				break;
			}
		}

		Node grandParent = tree.getNode(nodeNr);

		final int internalNodes = tree.getInternalNodeCount();
		if (internalNodes <= 1) {
			return Double.NEGATIVE_INFINITY;
		}

		Node parentIndex = grandParent.getLeft();
		Node uncle = grandParent.getRight();
		if (parentIndex.getHeight() < uncle.getHeight()) {
			parentIndex = grandParent.getRight();
			uncle = grandParent.getLeft();
		}

		if (parentIndex.isLeaf()) {
			// tree with dated tips
			return Double.NEGATIVE_INFINITY;
		}

		int validGP = 0;
		{
			for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; ++i) {
				validGP += isg(tree.getNode(i));
			}
		}

		final int c2 = sisg(parentIndex) + sisg(uncle);

		// check if both children have less than 0.01 mutations
		Node i;
		if (edgeWeights.getEdgeWeights(parentIndex.getLeft().getNr()) > limit) {
			i = parentIndex.getRight();
		} else if (edgeWeights.getEdgeWeights(parentIndex.getRight().getNr()) > limit) {
			i = parentIndex.getLeft();
		} else {
			i = (Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight());
		}

		exchangeNodes(i, uncle, parentIndex, grandParent);

		grandParent.makeDirty(3 - grandParent.isDirty());

		final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);

		logq += Math.log((float) validGP / validGPafter);		
		return logq;
	}

	protected void exchangeNodes(Node i, Node j, Node p, Node jP) {
		// precondition p -> i & jP -> j
		replace(p, i, j);
		replace(jP, j, i);
		// postcondition p -> j & p -> i
	}
	

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(percentage);
            final double f = Math.exp(delta);
            double lower = 1/treeInput.get().getNodeCount();
        	percentage = Math.min(Math.max(f, lower), 0.1);
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return percentage;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
    	percentage = value;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = percentage * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }
}
