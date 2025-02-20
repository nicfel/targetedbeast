
package targetedbeast.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import targetedbeast.edgeweights.EdgeWeights;
import targetedbeast.likelihood.RapidTreeLikelihood;

@Description("Performs narrow move on nodes that are close together in terms of the number of mutations between their consensus sequences")
public class NodeRandomizer extends TreeOperator {
	public final Input<Double> percentageInput = new Input<>("percentage",
			"percentage of nodes below limit to pick", 0.25);
    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights", "input of weights to be used for targetedn tree operations", Input.Validate.REQUIRED);

    EdgeWeights edgeWeights;
    double limit;
    double percentage;
    
	@Override
	public void initAndValidate() {
		edgeWeights = edgeWeightsInput.get();
		limit = edgeWeights.minEdgeWeight();
		percentage = percentageInput.get();
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
//		System.out.println();
//		System.out.println(tree + ";");

		double logHastingsRatio = 0;
		
		int numberOfAttempts = (int) (tree.getNodeCount() * percentage);

		// look for groups of nodes that have less than 0.1 mutations between them
		for (int i = 0; i < numberOfAttempts; i++) {
			logHastingsRatio += newNarrow(tree);
			if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
				return Double.NEGATIVE_INFINITY;
			}

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
}
