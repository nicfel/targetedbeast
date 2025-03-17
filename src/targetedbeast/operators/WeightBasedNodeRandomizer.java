
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

    EdgeWeights edgeWeights;
    double limit;
    double percentage;
    double lower;
    
	@Override
	public void initAndValidate() {
		edgeWeights = edgeWeightsInput.get();
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

		double logHastingsRatio = 0;
		
		int numberOfAttempts = (int) (tree.getNodeCount() * percentage);
		
		int totalCandidates = 0;
		boolean[] isCandidate = new boolean[tree.getNodeCount()];
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
					&& (edgeWeights.getEdgeWeights(right.getNr()) <= limit) &&
					(edgeWeights.getEdgeWeights(left.getLeft().getNr()) <= limit ||
					 edgeWeights.getEdgeWeights(left.getRight().getNr()) <= limit)) {
				
				// check the child node with the the higher
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
			// calculate reverse contribution
			// check if the parentIndex.getNr() changes the isCandidate status
			int reverseTotalCandidates = totalCandidates;
			int pIN = parentIndex.getNr();
			boolean nowIsCandidate = false;
			
			Node left = tree.getNode(pIN).getLeft();
			Node right = tree.getNode(pIN).getRight();

			if (left.getHeight() < right.getHeight()) {
				Node tmp = left;
				left = right;
				right = tmp;
			}
			if (!left.isLeaf()) {
				if (edgeWeights.getEdgeWeights(left.getNr()) <= limit
						&& (edgeWeights.getEdgeWeights(right.getNr()) <= limit) &&
						(edgeWeights.getEdgeWeights(left.getLeft().getNr()) <= limit ||
						 edgeWeights.getEdgeWeights(left.getRight().getNr()) <= limit)) {
	                nowIsCandidate = true;
	            }
			}

			///////////
			
			if (nowIsCandidate) {
				reverseTotalCandidates++;
			}
			if (isCandidate[pIN]) {
				reverseTotalCandidates--;
			}
			
			isCandidate[pIN] = nowIsCandidate;		
			

			
			logHastingsRatio += Math.log((float) reverseTotalCandidates / totalCandidates);		
			totalCandidates = reverseTotalCandidates;

		}
//		System.out.println(tree + ";");
		return logHastingsRatio;
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
