
package targetedbeast.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import targetedbeast.edgeweights.EdgeWeights;
import targetedbeast.likelihood.RapidTreeLikelihood;


/**
 * Implements the subtree slide move.
 */
@Description("Performs a scale move on the intervals between nodes.")
public class IntervalScaleOperator extends TreeOperator {

    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);

    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);

    final public Input<Boolean> optimiseInput = new Input<>("optimise",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			true);

	public Input<List<RealParameter>> downInput = new Input<>("down", "down parameter to scale", new ArrayList<>());

	public Input<List<RealParameter>> upInput = new Input<>("up", "up parameter to scale", new ArrayList<>());

	public Input<Boolean> scaleAllNodesIndependentlyInput = new Input<>("scaleAllNodesIndependently",
			"if true, all nodes are scaled with a different factor, otherwise a single factor is used", false);
	
    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights", "input of weights to be used for targetedn tree operations");

    private double scaleFactor;

    private double upper, lower;
    
    EdgeWeights edgeWeights = null;

	@Override
	public void initAndValidate() {
        scaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();
        
		if (edgeWeightsInput.get() != null) {
			edgeWeights = edgeWeightsInput.get();
		}
	}

	@Override
	public double proposal() {
		
		final Tree tree = (Tree) InputUtil.get(treeInput, this);
		if (scaleAllNodesIndependentlyInput.get()) {
			double logHR = resampleNodeHeight(tree.getRoot());			
			return logHR;
		}else {

			double scaler = getScaler();
			double lengthBefore = getTreeLength(tree.getRoot());
			int numbers = resampleNodeHeight(tree.getRoot(), scaler);
			double lengthAfter = getTreeLength(tree.getRoot());
			double actualScaler = lengthAfter / lengthBefore;
			
			double logHR = Math.log(scaler) * (numbers - 2);

			for (RealParameter down : downInput.get()) {
				down.setValue(down.getValue() / actualScaler);
				logHR -= Math.log(actualScaler);
			}
			for (RealParameter up : upInput.get()) {
				up.setValue(up.getValue() * actualScaler);
				logHR += Math.log(actualScaler);
			}
			return logHR;
		}
	}

	private int resampleNodeHeight(Node node, double scaler) {
		if (node.isLeaf()) {
			return 0;
		}
		double oldHeights = node.getHeight() - Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		int logHR = 1;
		logHR += resampleNodeHeight(node.getLeft(), scaler);
		logHR += resampleNodeHeight(node.getRight(), scaler);

		// resample the height
		double minHeight = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		double newHeight = oldHeights * scaler;
		node.setHeight(newHeight + minHeight);
		return logHR;
	}
	
	private double resampleNodeHeight(Node node) {
		if (node.isLeaf()) {
			return 0.0;
		}
		double oldHeights = node.getHeight() - Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		double logHR = 0.0;
		logHR += resampleNodeHeight(node.getLeft());
		logHR += resampleNodeHeight(node.getRight());
			

		// resample the height
		double scaler = -1;
		if (edgeWeights!=null) {
			// calculate the mutations above and below the node
			double total_muts = 0.0;
			total_muts += edgeWeights.getEdgeWeights(node.getLeft().getNr());
			total_muts += edgeWeights.getEdgeWeights(node.getRight().getNr());
			scaler = getScalerExp(1/total_muts);
		} else {
            scaler = getScalerExp();
		}
		double minHeight = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		double newHeight = oldHeights * scaler;
		node.setHeight(newHeight + minHeight);
		logHR += Math.log(scaler);
		return logHR;
	}

	private double getTreeLength(Node node) {
		double length = 0;
		if (!node.isRoot()) {
			length = node.getLength();
		}
		if (!node.isLeaf()) {
			length += getTreeLength(node.getLeft());
			length += getTreeLength(node.getRight());
		}

		return length;
	}
	
	
	protected double getScalerExp(double mutl_factor) {
		return Math.exp(Randomizer.nextGaussian()*(1-scaleFactor)*mutl_factor);
	}

	
	protected double getScalerExp() {
		return Math.exp(Randomizer.nextGaussian()*(1-scaleFactor));
	}

    protected double getScaler() {
        return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / scaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }
}
