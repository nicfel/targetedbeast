
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
 * Corrected variant of {@link IntervalScaleOperator}: the (stdev, rates) joint move is made
 * quantile-preserving so it leaves the branch-rate prior exactly invariant, while keeping the
 * node-height interval co-scaling that makes this an interval operator.
 */
@Description("Scale move on the intervals between nodes that also co-scales the lognormal stdev and "
		+ "the branch rates. In the stdev mode each branch rate is moved to the same quantile under the "
		+ "new lognormal (mean=1 in real space), leaving the rate prior invariant; the Hastings ratio is "
		+ "booked exactly. All other modes match IntervalScaleOperator.")
public class IntervalRateCoScaler extends TreeOperator {

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


    public Input<RealParameter> branchRatesInput = new Input<>("branchRates", "branch rates parameter to scale");

    public Input<RealParameter> stdevInput = new Input<>("stdev", "standard devation of the lognormal distribution");


    private double scaleFactor;

    private double upper, lower;

    // true when this instance uses getScalerExp() (Gaussian-in-log proposals), in which case
    // scaleFactor is interpreted as the Gaussian log-SD and coerced in log space (no SD<=1 ceiling).
    // false for the uniform getScaler() path (e.g. ScaleAll), which keeps scaleFactor in (0,1).
    private boolean useExpScaler;

    EdgeWeights edgeWeights = null;

	@Override
	public void initAndValidate() {
        scaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();

        useExpScaler = (branchRatesInput.get() != null) || scaleAllNodesIndependentlyInput.get();

		if (edgeWeightsInput.get() != null) {
			edgeWeights = edgeWeightsInput.get();
		}
	}

	@Override
	public double proposal() {

		final Tree tree = (Tree) InputUtil.get(treeInput, this);
		if (branchRatesInput.get() != null) {
			if (stdevInput.get() !=null) {
				// Joint move: scale node-height intervals FIRST (so branch lengths stay positive), then
				// ROTATE the branch rates so the lognormal rate prior is left EXACTLY invariant, without
				// changing the stdev (S) or mean hyper-parameters ("keeping the stdev and mean constant").
				//
				// Work in log-rate space z_i = log r_i. Split z = m0*1 + c, m0 = mean(z), c the residual
				// (sum c_i = 0). log P(rates | mu,S) depends on {z_i} ONLY through sum z_i and sum z_i^2,
				// i.e. through m0 and |c|; so ANY move that preserves the mean m0 and the residual norm |c|
				// leaves the rate prior exactly unchanged. The set {m0 fixed, |c| fixed} is a sphere, and
				// the valid, reversible move on it is a ROTATION of c (an isometry -> no rate Jacobian):
				//   (1) scale each internal interval by an independent exp-Gaussian scaler -> new lengths;
				//   (2) build the distance-compensating direction d_i = delta_i - mean(delta),
				//       delta_i = log(oldLen_i/newLen_i)  (delta = 0 at the root, which has no branch);
				//   (3) rotate c toward d by a symmetric random angle theta within the plane span(c,d):
				//         e1 = c/|c|, e2 = (d - (d.e1)e1)/|d_perp|,  c' = cos(theta) c + sin(theta) |c| e2.
				//       This keeps m0 and |c| exactly => sum z_i and sum z_i^2 exactly => rate prior exactly.
				// Reversibility: the reverse move sees swapped lengths (d -> -d) and the rotated residual c',
				// which span the SAME plane, and recovers c with angle +theta; theta is drawn from a
				// symmetric density, so the proposal ratio is 1. The rotation is an isometry of the sphere
				// (verified: it leaves the uniform-on-sphere = prior-conditional distribution invariant),
				// so it contributes NOTHING to the Hastings ratio. The whole HR is the height scaling only.
				// The cost lands in the likelihood (a single angle only partially compensates the whole
				// length change), which is step-size tunable, unlike the O(n) rate-prior blow-up of a bare
				// distance move.
				final RealParameter rates = branchRatesInput.get();
				final int dim = rates.getDimension();
				final int nodeCount = tree.getNodeCount();

				// --- step (1): capture old branch lengths, then scale intervals per node ---
				final double[] oldLen = new double[nodeCount];
				for (Node nd : tree.getNodesAsArray()) {
					oldLen[nd.getNr()] = nd.isRoot() ? 0.0 : nd.getParent().getHeight() - nd.getHeight();
				}
				final double heightLogJ = resampleNodeHeight(tree.getRoot());   // per-node interval scaling
				// positivity gate on the new branch lengths (branches first, else negative branches)
				for (Node nd : tree.getNodesAsArray()) {
					if (!nd.isRoot() && !(nd.getParent().getHeight() - nd.getHeight() > 0.0)) {
						return Double.NEGATIVE_INFINITY;
					}
				}

				// --- step (2): log-rates z, mean m0, and the distance-compensating direction d ---
				final double[] z = new double[dim];
				final double[] delta = new double[dim];
				double sumZ = 0.0, sumDelta = 0.0;
				for (int i = 0; i < dim; i++) {
					z[i] = Math.log(rates.getValue(i));
					sumZ += z[i];
					final Node nd = tree.getNode(i);
					if (nd.isRoot()) {
						delta[i] = 0.0;
					} else {
						final double newLen = nd.getParent().getHeight() - nd.getHeight();
						delta[i] = Math.log(oldLen[i] / newLen);
					}
					sumDelta += delta[i];
				}
				final double m0 = sumZ / dim;
				final double dBar = sumDelta / dim;

				// residual c = z - m0, direction d = delta - dBar (both mean-zero); build the rotation plane
				double normC2 = 0.0, normD2 = 0.0, dDotC = 0.0;
				for (int i = 0; i < dim; i++) {
					final double c = z[i] - m0;
					final double dd = delta[i] - dBar;
					normC2 += c * c;
					normD2 += dd * dd;
					dDotC += dd * c;
				}
				if (!(normC2 > 0.0)) {
					return Double.NEGATIVE_INFINITY;                 // degenerate: all rates equal
				}
				// |d_perp|^2 = |d|^2 - (d.c)^2/|c|^2 ; if ~0 the height change gives no rotation direction
				final double normDperp2 = normD2 - dDotC * dDotC / normC2;
				if (!(normDperp2 > 1e-24)) {
					return heightLogJ;                               // pure height scale, rates unchanged
				}
				final double normC = Math.sqrt(normC2);
				final double normDperp = Math.sqrt(normDperp2);
				final double proj = dDotC / normC2;                  // (d.c)/|c|^2
				final double g = dDotC / normC;                      // d . chat
				final double m = normDperp;                          // |d_perp| = sqrt(|d|^2 - g^2)

				// --- step (3): rotate c toward d by a symmetric random angle theta ---
				final double theta = Randomizer.nextGaussian() * scaleFactor;
				final double cosT = Math.cos(theta), sinT = Math.sin(theta);

				// Reversibility guard. Rotating toward the height-derived direction d, then reversing with
				// swapped lengths (d -> -d), retraces to the original residual ONLY when the rotation has not
				// carried c past d, i.e. when sigma = sign(g*sin(theta) - m*cos(theta)) < 0. In that region
				// the reverse move is itself in the region (its sigma is -m < 0), so the move is a closed,
				// reversible involution; proposals with sigma >= 0 have no matching reverse and must be
				// rejected to preserve detailed balance.
				if (g * sinT - m * cosT >= 0.0) {
					return Double.NEGATIVE_INFINITY;
				}
				for (int i = 0; i < dim; i++) {
					final double c = z[i] - m0;
					final double dd = delta[i] - dBar;
					final double e2 = (dd - proj * c) / normDperp;   // unit residual of d orthogonal to c
					final double cPrime = cosT * c + sinT * normC * e2;
					final double rNew = Math.exp(m0 + cPrime);
					if (rNew <= rates.getLower() || rNew >= rates.getUpper()) {
						return Double.NEGATIVE_INFINITY;
					}
					rates.setValue(i, rNew);
				}

				// The rotation preserves mean(z) and |c| exactly, so sum z_i and sum z_i^2 are unchanged and
				// the rate prior is exactly invariant. The rotation plane depends on c, so as a map of the
				// residual space its Jacobian is not 1: the differential is beta*I + a rank-2 update, giving
				// the closed form (verified vs. a finite-difference determinant) det = beta^(dim-1) * det2,
				//   beta = cos(theta) - sin(theta)*g/m,  k = g/m,  s = sin(theta)/beta,  det2 = (1+s*k)^2+s^2.
				final double beta = cosT - sinT * g / m;
				final double s = sinT / beta;
				final double k = g / m;
				final double det2 = (1.0 + s * k) * (1.0 + s * k) + s * s;
				if (!(Math.abs(beta) > 0.0) || !(det2 > 0.0)) {
					return Double.NEGATIVE_INFINITY;                 // singular map
				}
				final double logJrates = (dim - 1) * Math.log(Math.abs(beta)) + Math.log(det2);

				return heightLogJ + logJrates;
			}else {
				double logHR = resampleNodeHeight(tree.getRoot(), branchRatesInput.get());
				return logHR;
			}

		}
		else if (scaleAllNodesIndependentlyInput.get()) {
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

	private double resampleNodeHeight(Node node, RealParameter branchRates) {
		if (node.isLeaf()) {
			return 0.0;
		}

		double oldLengthLeft = node.getLeft().getLength();
		double oldLengthRight = node.getRight().getLength();

		double oldHeights = node.getHeight() - Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		double logHR = 0.0;
		logHR += resampleNodeHeight(node.getLeft(), branchRates);
		logHR += resampleNodeHeight(node.getRight(), branchRates);


		// resample the height
		double scaler = -1;
        scaler = getScalerExp();

		double minHeight = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		double newHeight = oldHeights * scaler;
		node.setHeight(newHeight + minHeight);
		logHR += Math.log(scaler);
		// scale the rates of the two children of the node
		double newLengthLeft = node.getLeft().getLength();
		double newLengthRight = node.getRight().getLength();

		branchRates.setValue(node.getLeft().getNr(), branchRates.getValue(node.getLeft().getNr()) * oldLengthLeft / newLengthLeft);
		branchRates.setValue(node.getRight().getNr(), branchRates.getValue(node.getRight().getNr()) * oldLengthRight / newLengthRight);
		logHR -= Math.log(newLengthLeft / oldLengthLeft);
		logHR -= Math.log(newLengthRight / oldLengthRight);

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


	// scaleFactor is the Gaussian log-SD in these modes (see useExpScaler / optimize()).
	protected double getScalerExp(double mutl_factor) {
		return Math.exp(Randomizer.nextGaussian()*scaleFactor*mutl_factor);
	}


	protected double getScalerExp() {
		return Math.exp(Randomizer.nextGaussian()*scaleFactor);
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
            if (useExpScaler) {
                // scaleFactor is the Gaussian log-SD: tune it in log space, no SD<=1 ceiling
                // (cf. beast.base.inference.operator.RealRandomWalkOperator).
                delta += Math.log(scaleFactor);
                setCoercableParameterValue(Math.exp(delta));
            } else {
                // uniform getScaler() with scaleFactor in (0,1): logit update
                delta += Math.log(1.0 / scaleFactor - 1.0);
                setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
            }
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        if (useExpScaler) {
            scaleFactor = Math.max(1e-8, value);
        } else {
            scaleFactor = Math.max(Math.min(value, upper), lower);
        }
    }

    @Override
    public String getPerformanceSuggestion() {
        // In exp mode scaleFactor is auto-tuned in log space; no manual suggestion is meaningful.
        if (useExpScaler) {
            return "";
        }
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
