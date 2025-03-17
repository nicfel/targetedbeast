
package targetedbeast.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
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
public class HeightBasedNodeRandomizer extends TreeOperator {
	public final Input<Double> percentageInput = new Input<>("percentage",
			"percentage of nodes below limit to pick", 0.25);
    final public Input<Boolean> optimiseInput = new Input<>("optimise",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			true);

    double limit;
    double percentage;
    double lower;
    
	@Override
	public void initAndValidate() {
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
		
		List<Double> heights = new ArrayList<>();
		for (int i = 0; i < tree.getNodeCount(); i++) {			
			double heightDiff = tree.getNode(i).getLength();
			heights.add(heightDiff);
		}
		Collections.sort(heights);
		int numberOfAttempts = (int) (heights.size() * percentage)+1;
		double limit = heights.get(numberOfAttempts);
		while (limit == 0) {
			numberOfAttempts++;
			limit = heights.get(numberOfAttempts);
		}
		
//		System.out.println(numberOfAttempts + " " + limit);
//		System.out.println(tree + ";");

		// check which ones are candidate edges
		boolean[] isCandidate = new boolean[tree.getNodeCount()];
		int totalCandidates = 0;
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

			double heightDiff = tree.getNode(i).getHeight() - tree.getNode(i).getLeft().getHeight();

			// only pick nodes that have less than 0.01 mutations on left and at least one
			// child has less than 0.01 mutations
			if (heightDiff < limit) {
				// check the child node with the the higher
				isCandidate[i] = true;
				totalCandidates++;
			}
		}
		
		
		// look for groups of nodes that have less than 0.1 mutations between them
		for (int k = 0; k < numberOfAttempts; k++) {
			
			
			
			// choose a random node avoiding root
			if (totalCandidates == 0) {
				return Double.NEGATIVE_INFINITY;
			}

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
			// parent has to be older than uncle

			if (parentIndex.isLeaf()) {
				// tree with dated tips
				return Double.NEGATIVE_INFINITY;
			}

			// pick one of the children at random
			Node i = (Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight());
			// exchange nodes
			exchangeNodes(i, uncle, parentIndex, grandParent);
			// set dirty flags
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
				double heightDiff = tree.getNode(pIN).getHeight() - tree.getNode(pIN).getLeft().getHeight();
				if (heightDiff < limit) {
	                nowIsCandidate = true;
	            }
			}
			
			if (nowIsCandidate) {
				reverseTotalCandidates++;
			}
			if (isCandidate[pIN]) {
				reverseTotalCandidates--;
			}
			
			isCandidate[pIN] = nowIsCandidate;		
			

			
			logHastingsRatio += Math.log((float) reverseTotalCandidates / totalCandidates);		
			totalCandidates = reverseTotalCandidates;
			
			
			
//			logHastingsRatio += newNarrow(tree, limit, isCandidate, totalCandidates);
//			if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
//				return Double.NEGATIVE_INFINITY;
//			}
		}
//		System.out.println(tree + ";\n");
//		System.exit(0);
		return logHastingsRatio;
	}


	public double newNarrow(final Tree tree, double limit, boolean[] isCandidate, int totalCandidates) {

		System.out.println(tree + ";");
		double logq = 0;
		// choose a random node avoiding root
		if (totalCandidates == 0) {
			return Double.NEGATIVE_INFINITY;
		}

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
		// parent has to be older than uncle

		if (parentIndex.isLeaf()) {
			// tree with dated tips
			return Double.NEGATIVE_INFINITY;
		}

		// pick one of the children at random
		Node i = (Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight());
		// exchange nodes
		exchangeNodes(i, uncle, parentIndex, grandParent);
		// set dirty flags
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
			double heightDiff = tree.getNode(pIN).getHeight() - tree.getNode(pIN).getLeft().getHeight();
			if (heightDiff < limit) {
                nowIsCandidate = true;
            }
		}
		
		if (nowIsCandidate) {
			reverseTotalCandidates++;
		}
		if (isCandidate[pIN]) {
			reverseTotalCandidates--;
		}
		
		isCandidate[pIN] = nowIsCandidate;		
		

		
		logq += Math.log((float) reverseTotalCandidates / totalCandidates);		
		totalCandidates = reverseTotalCandidates;
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
        	percentage = Math.min(Math.max(f, lower), 0.25);
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
