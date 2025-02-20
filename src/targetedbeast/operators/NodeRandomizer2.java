
package targetedbeast.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import targetedbeast.edgeweights.EdgeWeights;

@Description("Performs narrow move on nodes that are close together in terms of the number of mutations between their consensus sequences")
public class NodeRandomizer2 extends TreeOperator {
	public final Input<Double> percentageInput = new Input<>("percentage",
			"percentage of nodes below limit to pick", 0.05);
    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights", "input of weights to be used for targetedn tree operations", Input.Validate.REQUIRED);
    public Input<Double> mutationLimit = new Input<>("mutLimit", "limit on the edge weight used for the narrow move");
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
		if (mutationLimit.get() != null) {
			limit = mutationLimit.get();
		}
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
			if (edgeWeights.getEdgeWeights(left.getNr()) <= limit) {
				// check the child node with the the higher
				deviation[i] = 1;
				totalDeviation += deviation[i];
			}else {
				deviation[i] = 0.01; // here to enable the reverse move if needed
				totalDeviation += deviation[i];
			}
		}
				
		List<Integer> nodesToOperateOn = new ArrayList<>();
		
		int noNodesToSelect = (int) (tree.getInternalNodeCount()*percentage);
		
		for (int i = 0; i < noNodesToSelect; i++) {
			// randomly pick nodes to operate on using deviation
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
            nodesToOperateOn.add(nodeNr);
		}
		
		// calculate the HR contribution for the forward move
		for (int nodeNr : nodesToOperateOn) {
			logHastingsRatio -= Math.log(deviation[nodeNr] / totalDeviation);
		}
		
		for (int nodeNr : nodesToOperateOn) {
			logHastingsRatio += newNarrow(tree, nodeNr);
			if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
				return Double.NEGATIVE_INFINITY;
			}
		}
		
		// calculate the HR contribution for the reverse move
        edgeWeights.prestore();
        edgeWeights.updateByOperator();
        
		// choose a random node avoiding root
		totalDeviation = 0;
		deviation = new double[tree.getNodeCount()];
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

			// only pick nodes that have less don't have any mutations on them
			if (edgeWeights.getEdgeWeights(left.getNr()) <= limit) {
				// check the child node with the the higher
				deviation[i] = 1;
				totalDeviation += deviation[i];
			}else {
				deviation[i] = 0.01; // here to enable the reverse move if needed
				totalDeviation += deviation[i];
			}
		}
		
		for (int nodeNr : nodesToOperateOn) {
			logHastingsRatio += Math.log(deviation[nodeNr] / totalDeviation);
		}
		return logHastingsRatio;
	}

	private int isg(final Node n) {
		return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
	}

	private int sisg(final Node n) {
		return n.isLeaf() ? 0 : isg(n);
	}

	public double newNarrow(final Tree tree, int nodeNr) {
		Node grandParent = tree.getNode(nodeNr);

		final int internalNodes = tree.getInternalNodeCount();

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

//		logq += Math.log((float) validGP / validGPafter);
//		System.out.println(Math.log((float) validGP / validGPafter));	
//		if (Math.log((float) validGP / validGPafter) != 0) {
//			System.out.println(tree + ";");
//			System.exit(0);
//		}
		return 0;//Math.log((float) validGP / validGPafter);
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
            delta += Math.log(1.0 / percentage - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return percentage;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
    	percentage = Math.max(Math.min(value, 1.0), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(percentage, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

	
}
