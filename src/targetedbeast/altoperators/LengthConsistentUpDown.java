package targetedbeast.altoperators;


import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.inference.operator.kernel.KernelOperator;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;


@Description("Like the UpDownOperator, this element represents an operator that scales "
		+ "two (or more) parameters in different directions, but uses a Bactrian proposal distribution for the scale value. "
        + "The up parameter is multiplied by this scale and the down parameter is divided by this scale.")
public class LengthConsistentUpDown extends KernelOperator {
    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "magnitude factor used for scaling", Validate.REQUIRED);
    final public Input<List<StateNode>> upInput = new Input<>("up",
            "zero or more items to scale upwards", new ArrayList<>());
    final public Input<List<StateNode>> downInput = new Input<>("down",
            "zero or more items to scale downwards", new ArrayList<>());
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Boolean> elementWiseInput = new Input<>("elementWise", "flag to indicate that the scaling is applied to a random index in multivariate parameters (default false)", false);

    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 10.0);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 0.0);

    double scaleFactor;
    private double upper,lower;



    @Override
    public void initAndValidate() {
    	super.initAndValidate();
        scaleFactor = scaleFactorInput.get();
        // sanity checks
        if (upInput.get().size() + downInput.get().size() == 0) {
        	Log.warning.println("WARNING: At least one up or down item must be specified");
        }
        if (upInput.get().size() == 0 || downInput.get().size() == 0) {
        	Log.warning.println("WARNING: no " + (upInput.get().size() == 0 ? "up" : "down") + " item specified in UpDownOperator");
        }
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();
    }
    
	protected double getScaler(int i) {
		return kernelDistribution.getScaler(i, Double.NaN, getCoercableParameterValue());
	}
    

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal
     *         should not be accepted
     */
    @Override
    public final double proposal() {

        final double scale = getScaler(0);
        int goingUp = 0, goingDown = 0;


        if (elementWiseInput.get()) {
            int size = 0;
            for (StateNode up : upInput.get()) {
                if (size == 0) size = up.getDimension();
                if (size > 0 && up.getDimension() != size) {
                    throw new RuntimeException("elementWise=true but parameters of differing lengths!");
                }
                goingUp += 1;
            }

            for (StateNode down : downInput.get()) {
                if (size == 0) size = down.getDimension();
                if (size > 0 && down.getDimension() != size) {
                    throw new RuntimeException("elementWise=true but parameters of differing lengths!");
                }
                goingDown += 1;
            }

            int index = Randomizer.nextInt(size);

            for (StateNode up : upInput.get()) {
                if (up instanceof RealParameter) {
                    RealParameter p = (RealParameter) up;
                    p.setValue(index, p.getValue(index) * scale);
                }
                if (outsideBounds(up)) {
                    return Double.NEGATIVE_INFINITY;
                }
            }

            for (StateNode down : downInput.get()) {
                if (down instanceof RealParameter) {
                    RealParameter p = (RealParameter) down;
                    p.setValue(index, p.getValue(index) / scale);
                }
                if (outsideBounds(down)) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        } else {

            try {
                for (StateNode up : upInput.get()) {
                    up = up.getCurrentEditable(this);
					if (up instanceof Tree) {
						int val = scaleTree((Tree) up, scale);
						if (val == -1) {
							return Double.NEGATIVE_INFINITY;
						}
						goingUp += val;
					} else {
						goingUp += up.scale(scale);
					}
                }
                // separated this into second loop because the outsideBounds might return true transiently with
                // related variables which would be BAD. Note current implementation of outsideBounds isn't dynamic,
                // so not currently a problem, but this became a problem in BEAST1 so this is preemptive strike.
                // Same below for down
                for (StateNode up : upInput.get()) {
                    if (outsideBounds(up)) {
                        return Double.NEGATIVE_INFINITY;
                    }
                }

                for (StateNode down : downInput.get()) {
                    down = down.getCurrentEditable(this);
					if (down instanceof Tree) {
						int val = scaleTree((Tree) down, 1.0 / scale);
						if (val == -1) {
							return Double.NEGATIVE_INFINITY;
						}
						goingDown += val;
						System.out.println(down);
					} else {
						goingUp += down.scale(1.0 / scale);
					}

                }
                for (StateNode down : downInput.get()) {
                    if (outsideBounds(down)) {
                        return Double.NEGATIVE_INFINITY;
                    }
                }
            } catch (Exception e) {
                // scale resulted in invalid StateNode, abort proposal
                return Double.NEGATIVE_INFINITY;
            }
        }
        return (goingUp - goingDown) * Math.log(scale);
    }

    private int scaleTree(Tree tree, double scale) {
    	System.out.println("scaling tree");
    	System.out.println(scale);
    	System.out.println("=============");
    	System.out.println(tree);
    	
    	double[] tmp = scaleSubtreeLength(tree.getRoot(), scale);
    	System.exit(0);
		if (tmp[0] == Double.NEGATIVE_INFINITY)
			return -1;
    	return tree.getInternalNodeCount();
	}
    
	private double[] scaleSubtreeLength(Node node, double scale) {
		double[] oldNewLengths = new double[] {0.0, 0.0};
		// keep track of the edge length before any scaling
		if (node.isLeaf()) {
			oldNewLengths[0] = node.getParent().getHeight() - node.getHeight();
			return oldNewLengths;		
		}
		
		if (!node.isRoot()) 
			oldNewLengths[0] = node.getParent().getHeight() - node.getHeight();
		
		// get the old and the 		
		double subtreeLength = 0.0;
		double subsubtreeLength = 0.0;
		System.out.println(":::::::::::");
		for (Node child : node.getChildren()) {
			double[] tmp = scaleSubtreeLength(child, scale);
			subtreeLength += tmp[0];
			subsubtreeLength += tmp[1];
			if (node.getNr()==26)
				System.out.println("tmp: " + tmp[0] + " " + tmp[1]);
		}
		System.out.println("=============");
		System.out.println(node.getHeight() + " " + node.getNr());
		System.out.println(subtreeLength + " " + subsubtreeLength);
		System.out.println(subsubtreeLength*scale);
		System.out.println((subtreeLength+subsubtreeLength)*scale);
		System.out.println(node.getLeft().getHeight() + " " + node.getRight().getHeight());
		
		
		// get the new length 
		double remainingNewLength = subtreeLength*scale;
		
		System.out.println(remainingNewLength);
		// compute the new node height to match that
		remainingNewLength += node.getLeft().getHeight() + node.getRight().getHeight();
		System.out.println(node.getHeight() + " " + remainingNewLength/2);
		
		// set the new height
		node.setHeight(remainingNewLength/2);
		
		if (node.getHeight() < node.getLeft().getHeight() || node.getHeight() < node.getRight().getHeight()) {
			System.out.println("reject");
			System.out.println(node.getHeight() + " " + node.getLeft().getHeight() + " " + node.getRight().getHeight());
			System.exit(0);
			oldNewLengths[0] = Double.NEGATIVE_INFINITY;
		}
		
//		oldNewLengths[0] = subtreeLength;	
		oldNewLengths[1] += subtreeLength;
		
//		System.exit(0);
		
		return oldNewLengths;
		
	}

	private boolean outsideBounds(final StateNode node) {
        if (node instanceof Parameter<?>) {
            final Parameter<?> p = (Parameter<?>) node;
            final Double lower = (Double) p.getLower();
            final Double upper = (Double) p.getUpper();
            final Double value = (Double) p.getValue();
            if (value < lower || value > upper) {
                return true;
            }
        }
        return false;
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
    	if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
    	}
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value,upper),lower);
    }

    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }

    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = getCoercableParameterValue() * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scale factor to about " + formatter.format(newWindowSize);
        } else return "";
    }
    
} // class UpDownOperator
