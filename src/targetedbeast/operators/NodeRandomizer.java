
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
import targetedbeast.likelihood.RapidTreeLikelihood;

@Description("Performs narrow move on nodes that are close together in terms of the number of mutations between their consensus sequences")
public class NodeRandomizer extends TreeOperator {
	final public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
			"if true (default) a narrow exchange is performed, otherwise a wide exchange", true);
	public final Input<Integer> numberOfAttemptsInput = new Input<>("attempts",
			"number of attempts to make before giving up", 50);

	public Input<RapidTreeLikelihood> rapidTreeLikelihoodInput = new Input<>("rapidTreeLikelihood",
			"The likelihood to be used for the tree proposal. If not specified, the tree likelihood is calculated from the tree.");
	final public Input<Boolean> isNewNarrowInput = new Input<>("isNewNarrow",
			"if true (default) a new narrow exchange is performed, otherwise the old one", false);

	@Override
	public void initAndValidate() {
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

		// look for groups of nodes that have less than 0.1 mutations between them
		for (int i = 0; i < numberOfAttemptsInput.get(); i++) {
			logHastingsRatio += newNarrow(tree);
			if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
				return Double.NEGATIVE_INFINITY;
			}

		}
		return logHastingsRatio;
	}

	private int isg(final Node n) {
		return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
	}

	private int sisg(final Node n) {
		return n.isLeaf() ? 0 : isg(n);
	}

	public double newNarrow(final Tree tree) {

		double limit = 0.1;

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
			if (rapidTreeLikelihoodInput.get().getEdgeMutations(left.getNr()) < limit
					&& (rapidTreeLikelihoodInput.get().getEdgeMutations(left.getLeft().getNr()) < limit
							|| rapidTreeLikelihoodInput.get().getEdgeMutations(left.getRight().getNr()) < limit)

			) {

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
		if (rapidTreeLikelihoodInput.get().getEdgeMutations(parentIndex.getLeft().getNr()) >= limit) {
			i = parentIndex.getRight();
		} else if (rapidTreeLikelihoodInput.get().getEdgeMutations(parentIndex.getRight().getNr()) >= limit) {
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

	public double newWide(final Tree tree) {
		double logHR = 0;

		// get the total amount of mutations
		double totalMutations = 0;
		for (int i = 0; i < tree.getNodeCount(); i++) {
			totalMutations += rapidTreeLikelihoodInput.get().getEdgeMutations(i);
		}

		// add plus one for each node tree.getNodeCount(), except the root node
		totalMutations += tree.getNodeCount() - 1;
		// pick a random numner
		double scaler = Randomizer.nextDouble() * totalMutations;
		int randomNode = -1;
		double currMuts = 0;
		for (int i = 0; i < tree.getNodeCount(); i++) {
			currMuts += rapidTreeLikelihoodInput.get().getEdgeMutations(i) + 1;
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
		}

//		if (tree.getNode(randomNode).getParent().isRoot()) {
//			return Double.NEGATIVE_INFINITY;
//		}

//    	System.out.println("before = " + (rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode)+1)/ totalMutations);
//		System.out.println("\n"+rapidTreeLikelihoodInput.get().toNewick(treeInput.get().getRoot()) + " ;");

		logHR -= Math.log((rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode) + 1) / totalMutations);

		// pick a random number between 0 and root height

		// pick the other node that has to co-exist with the parent node height
		double parentHeight = tree.getNode(randomNode).getParent().getHeight();
		List<Integer> nodes = new ArrayList<>();
		for (int i = 0; i < tree.getNodeCount(); i++) {
			Node node = tree.getNode(i);
			if (node.getHeight() < parentHeight && node.getParent().getHeight() > parentHeight) {
				nodes.add(i);
			}
		}

		// pick the re-attachment edge
		if (nodes.size() == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		Node i = tree.getNode(randomNode);
		Node j = tree.getNode(nodes.get(Randomizer.nextInt(nodes.size())));

//        System.out.println(i.getNr() + " " + j.getNr() + " " + i.getParent().getNr() + " " + j.getParent().getNr());

		final Node p = i.getParent();
		final Node jP = j.getParent();

		if ((p != jP) && (i != jP) && (j != p) && (j.getHeight() < p.getHeight()) && (i.getHeight() < jP.getHeight())) {
			exchangeNodes(i, j, p, jP);

			// All the nodes on the path from i/j to the common ancestor of i/j parents had
			// a topology change,
			// so they need to be marked FILTHY.
			if (markCladesInput.get()) {
				Node iup = p;
				Node jup = jP;
				while (iup != jup) {
					if (iup.getHeight() < jup.getHeight()) {
						assert !iup.isRoot();
						iup = iup.getParent();
						iup.makeDirty(Tree.IS_FILTHY);
					} else {
						assert !jup.isRoot();
						jup = jup.getParent();
						jup.makeDirty(Tree.IS_FILTHY);
					}
				}
			}

			rapidTreeLikelihoodInput.get().prestore();
			rapidTreeLikelihoodInput.get().updateByOperator();
			// recalculate hasting ratio
			// calculate tree length
			totalMutations = 0;
			for (int k = 0; k < tree.getNodeCount(); k++) {
				totalMutations += rapidTreeLikelihoodInput.get().getEdgeMutations(k);
			}
			totalMutations += tree.getNodeCount() - 1;

//    		System.out.println(rapidTreeLikelihoodInput.get().toNewick(treeInput.get().getRoot()) + " ;");
//        	System.out.println("after = " + (rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode)+1)/ totalMutations);

//        	System.out.println(i.getNr() + " " + randomNode + " " + j.getNr());

			logHR += Math.log((rapidTreeLikelihoodInput.get().getEdgeMutations(j.getNr()) + 1) / totalMutations);

			// get the next shared node
			while (i != j) {
				if (i.getHeight() < j.getHeight()) {
					i = i.getParent();
				} else {
					j = j.getParent();
				}
			}

			i.makeDirty(3 - i.isDirty());
//        	System.out.println("logHR " + logHR);
			rapidTreeLikelihoodInput.get().unstore();
			return logHR;
		}

		// Randomly selected nodes i and j are not valid candidates for a wide exchange.
		// reject instead of counting (like we do for narrow).
		return Double.NEGATIVE_INFINITY;
	}

	/* exchange sub-trees whose root are i and j */

	protected void exchangeNodes(Node i, Node j, Node p, Node jP) {
		// precondition p -> i & jP -> j
		replace(p, i, j);
		replace(jP, j, i);
		// postcondition p -> j & p -> i
	}
}
