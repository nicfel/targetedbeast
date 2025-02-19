/*
* File SubtreeSlide.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
/*
 * SubtreeSlideOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package targetedbeast.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import targetedbeast.edgeweights.EdgeWeights;

/**
 * Implements the subtree slide move.
 */
@Description("")
public class RangeSlide extends TreeOperator {

	final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());
	final public Input<Boolean> optimiseInput = new Input<>("optimise",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			true);
	final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, "
			+ "i.e. -1. (when positive, gets multiplied by tree-height/log2(n-taxa).", -1.0);

	public Input<Boolean> useWeightedStepInput = new Input<>("useWeightedStep", "Use weighted step", false);

    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights", "input of weights to be used for targetedn tree operations", Input.Validate.REQUIRED);

	// shadows size
	protected double size;
	private double limit;
	EdgeWeights edgeWeights;
    KernelDistribution kernelDistribution;

	@Override
	public void initAndValidate() {
		size = sizeInput.get();
		limit = limitInput.get();
		edgeWeights = edgeWeightsInput.get();
        kernelDistribution = kernelDistributionInput.get();
	}

	/**
	 * Do a probabilistic subtree slide move.
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted *
	 */
	@Override
	public double proposal() {
		final Tree tree = (Tree) InputUtil.get(treeInput, this);

		double val = getDelta();

		if (useWeightedStepInput.get())
			return doWeightedStep(tree, val);
		else
			return doStep(tree, val);
	}

	private double getDelta() {
    	return kernelDistribution.getRandomDelta(0, Double.NaN, size);
	}

	private double doStep(Tree tree, double val) {

		double logHastingsRatio = 0.0;
		// choose a random node avoiding root
		double totalWeight = 0;
		double[] weight = new double[tree.getNodeCount()];
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isRoot())
				continue;
			weight[i] = edgeWeights.getEdgeWeights(i);
			totalWeight += weight[i];
		}

		double scaler = Randomizer.nextDouble() * totalWeight;
		int nodeNr = -1;
		double currDev = 0;
		for (int i = 0; i < weight.length; i++) {
			currDev += weight[i];
			if (currDev > scaler) {
				nodeNr = i;
				break;
			}
		}

		logHastingsRatio -= Math.log(weight[nodeNr] / totalWeight);

		Node i = tree.getNode(nodeNr);

		if (i.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}
		final Node p = i.getParent();

		Node CiP = getOtherChild(p, i);
		Node PiP = p.getParent();

		// 2. choose a delta to move
		final double delta = Math.abs(val);
		// keeps track of all target nodes
		Map<Node, Double> targets = new HashMap<>();
		targets.putAll(collectTargets(p, delta, true));
		targets.putAll(collectTargets(CiP, delta, false));

		// calculate the minimal possible heights
		double minHeight = i.getHeight();
		// remove all targets that are
		List<Node> keys = new ArrayList<>(targets.keySet());

		for (Node node : keys) {
			if (targets.get(node) < minHeight) {
				targets.remove(node);
			}
		}

		if (targets.size() == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		// randomly choose a target edge to attach p to
		int targetIndex = Randomizer.nextInt(targets.size());
		Node target = targets.keySet().toArray(new Node[0])[targetIndex];
		double newHeight = targets.get(target);

		logHastingsRatio -= Math.log(1.0 / targets.size());

		// check if target is p or CiP
		if (target != p && target != CiP) {
			if (target.getParent() == null) { // new root
				// Attach the other child of p directly to the parent of p, replacing p
				replace(PiP, p, CiP);
				// Attach p to the targets parent, replacing the target
				replace(p, CiP, target);
				p.setParent(null);
				tree.setRoot(p);
				p.setHeight(newHeight);
			} else if (p.isRoot()) { // was root
				// new root is CiP
				replace(target.getParent(), target, p);
				replace(p, CiP, target);
				CiP.setParent(null);
				tree.setRoot(CiP);
			} else {
				// remove p from its current position
				// Attach the other child of p to the parent of p
				replace(PiP, p, CiP);
				// Attach p to the targets parent, replacing the target
				replace(target.getParent(), target, p);
				// Attach the target to p
				replace(p, CiP, target);

				// set the highest node to

				// find the common ancestor of PiP and target.getParent() and set it to make
				// dirty 3
				Node tP = target.getParent();
				while (tP != PiP) {
					if (tP.getHeight() < PiP.getHeight()) {
						tP = tP.getParent();
					} else {
						PiP = PiP.getParent();
					}
				}
				PiP.makeDirty(3 - i.isDirty());
			}
		}

		p.setHeight(newHeight);

		Map<Node, Double> newTargets = new HashMap<>();
		newTargets.putAll(collectTargets(p, delta, true));
		newTargets.putAll(collectTargets(getOtherChild(p, i), delta, false));

		List<Node> newKeys = new ArrayList<>(newTargets.keySet());

		for (Node node : newKeys) {
			if (newTargets.get(node) < minHeight) {
				newTargets.remove(node);
			}
		}
		logHastingsRatio += Math.log(1.0 / newTargets.size());

		// check if any of the branch lengths are negative
		for (int j = 0; j < tree.getNodeCount(); j++) {
			if (tree.getNode(j).getLength() < 0) {
				throw new RuntimeException("negative branch length in Range Slide proposal, should not occur");
			}
		}

		edgeWeights.prestore();
		edgeWeights.updateByOperator();
		
		// choose a random node avoiding root
		totalWeight = 0;
		weight = new double[tree.getNodeCount()];
		for (int j = 0; j < tree.getNodeCount(); j++) {
			if (tree.getNode(j).isRoot())
				continue;
			weight[j] = edgeWeights.getEdgeWeights(j);
			totalWeight += weight[j];
		}
		logHastingsRatio += Math.log(weight[i.getNr()] / totalWeight);

		return logHastingsRatio;
	}

	private double doWeightedStep(Tree tree, double val) {
		
		double logHastingsRatio = 0.0;
		// choose a random node avoiding root
		double totalWeight = 0;
		double[] weight = new double[tree.getNodeCount()];
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isRoot())
				continue;
			weight[i] = edgeWeights.getEdgeWeights(i);
			totalWeight += weight[i];
		}

		double scaler = Randomizer.nextDouble() * totalWeight;
		int nodeNr = -1;
		double currDev = 0;
		for (int i = 0; i < weight.length; i++) {
			currDev += weight[i];
			if (currDev > scaler) {
				nodeNr = i;
				break;
			}
		}
		logHastingsRatio -= Math.log(weight[nodeNr] / totalWeight);


		Node i = tree.getNode(nodeNr);

		if (i.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}

		final Node p = i.getParent();
		Node CiP = getOtherChild(p, i);
		Node PiP = p.getParent();

		// 2. choose a delta to move
		double delta = Math.abs(val);
		// keeps track of all target nodes
		Map<Node, Double> targets = new HashMap<>();
		targets.putAll(collectTargets(p, delta, true));
		targets.putAll(collectTargets(CiP, delta, false));

		// calculate the minimal possible heights
		double minHeight = i.getHeight();
		// remove all targets that are
		List<Node> keys = new ArrayList<>(targets.keySet());

		for (Node node : keys) {
			if (targets.get(node) < minHeight) {
				targets.remove(node);
			}
		}

		if (targets.size() == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		// get all nodes in target that are ancestors to p
		List<Integer> ancestors = new ArrayList<>();
		// make a copy of p
		Node ancestor = i.getParent();
		while (ancestor != null) {
			// check if p is in the target
			if (targets.containsKey(ancestor)) {
				ancestors.add(ancestor.getNr());
			}
			ancestor = ancestor.getParent();
		}

		// calculate the consensus sequences without the node i
		edgeWeights.prestore();
		edgeWeights.updateByOperatorWithoutNode(i.getNr(), ancestors);

		double[] distance = edgeWeights.getTargetWeights(i.getNr(), new ArrayList<>(targets.keySet()));
		
		double totalDistance = 0;		
		for (int k=0; k < distance.length; k++) {
            totalDistance += distance[k];
        }

		double scaler2 = Randomizer.nextDouble() * totalDistance;
		double currDist = 0;
		int targetIndex = -1;
		int k = 0;
		for (k = 0; k < distance.length; k++) {
			currDist += distance[k];
			if (currDist > scaler2) {
				targetIndex = k;
				break;
			}
		}
		
		// pick a random node
		logHastingsRatio -= Math.log(distance[targetIndex] / totalDistance);

		Node target = targets.keySet().toArray(new Node[0])[targetIndex];
		
		double newHeight = targets.get(target);

		// check if target is p or CiP
		if (target != p && target != CiP) {
			if (target.getParent() == null) { // new root
				// Attach the other child of p directly to the parent of p, replacing p
				replace(PiP, p, CiP);
				// Attach p to the targets parent, replacing the target
				replace(p, CiP, target);
				p.setParent(null);
				tree.setRoot(p);
				p.setHeight(newHeight);
			} else if (p.isRoot()) { // was root
				// new root is CiP
				replace(target.getParent(), target, p);
				replace(p, CiP, target);
				CiP.setParent(null);
				tree.setRoot(CiP);
			} else {
				// remove p from its current position
				// Attach the other child of p to the parent of p
				replace(PiP, p, CiP);
				// Attach p to the targets parent, replacing the target
				replace(target.getParent(), target, p);
				// Attach the target to p
				replace(p, CiP, target);
				// find the common ancestor of PiP and target.getParent() and set it to make
				// dirty 3
				Node tP = target.getParent();
				while (tP != PiP) {
					if (tP.getHeight() < PiP.getHeight()) {
						tP = tP.getParent();
					} else {
						PiP = PiP.getParent();
					}
				}
				PiP.makeDirty(3 - i.isDirty());
			}
		}

		p.setHeight(newHeight);

		Map<Node, Double> newTargets = new HashMap<>();
		// 0 branch lengths can cause issues here in that the edge where the node is from may not be found again for the 
		// reverse calculation of the HR;
		newTargets.putAll(collectTargets(p, delta, true));
		newTargets.putAll(collectTargets(getOtherChild(p, i), delta, false));

		List<Node> newKeys = new ArrayList<>(newTargets.keySet());

		for (Node node : newKeys) {
			if (newTargets.get(node) < minHeight) {
				newTargets.remove(node);
			}
		}

		if (logHastingsRatio == Double.NEGATIVE_INFINITY) {
			return logHastingsRatio;
		}
		
		// calculate the distance of the previous edge without node i, prior to updating
		// everything
		k = 0;
		distance = edgeWeights.getTargetWeights(i.getNr(), new ArrayList<>(newTargets.keySet()));
		edgeWeights.reset();

		totalDistance = 0;
		int oldNodeIndex = -1;

		for (Node newTarget : newTargets.keySet()) {
			int nodeNo = newTarget.getNr();
			if (newTarget.getNr() == p.getNr()) { // if the new target is p, use the old node number
				nodeNo = target.getNr();
			}
			if (newTarget == CiP) { // if the new target is CiP, use the old node number
				oldNodeIndex = k;
			}
			
			totalDistance += distance[k];
			k++;
		}
	

		// check if CiP is in the new targets
		if (oldNodeIndex == -1) {
			// get the index of p in the new targets
			k = 0;
			for (Node newTarget : newTargets.keySet()) {
				if (newTarget == p) {
					oldNodeIndex = k;
					break;
				}
				k++;
			}
		}
		
		if (oldNodeIndex == -1) {
			System.err.println("couldn't find original node for backwards move,\n"
					+ "return negative infinity instead, can occur when the other child edge of p has length 0");
			return Double.NEGATIVE_INFINITY;
		}
		
		logHastingsRatio += Math.log(distance[oldNodeIndex] / totalDistance);

		// reset the consensus calculated without the node i		
		edgeWeights.prestore();
		edgeWeights.updateByOperator();
		
		// choose a random node avoiding root
		totalWeight = 0;
		weight = new double[tree.getNodeCount()];
		for (int j = 0; j < tree.getNodeCount(); j++) {
			if (tree.getNode(j).isRoot())
				continue;
			weight[j] = edgeWeights.getEdgeWeights(j);
			totalWeight += weight[j];
		}
		logHastingsRatio += Math.log(weight[i.getNr()] / totalWeight);
		return logHastingsRatio;
	}

	private Map<Node, Double> collectTargets(Node start, double distanceLeft, boolean isUp) {
		Map<Node, Double> targets = new HashMap<>();

		if (start.isLeaf()) {
			if (start.getLength() > distanceLeft) {
				Double newHeight = start.getParent().getHeight() - distanceLeft;
				// add start and new height to the targets
				targets.put(start, newHeight);
				return targets;

			} else {
				return targets;
			}
		}

		if (isUp) {
			if (start.isRoot()) {
				targets.put(start, start.getHeight() + distanceLeft);
				return targets;
			}
			// check if we need to move to the next node
			if (start.getLength() > distanceLeft) {
				targets.put(start, start.getHeight() + distanceLeft);
				return targets;
			} else {
				distanceLeft -= start.getLength();
				targets.putAll(collectTargets(start.getParent(), distanceLeft, true));
				// move to the other child
				targets.putAll(collectTargets(getOtherChild(start.getParent(), start), distanceLeft, false));
			}
		} else {
			// check if we need to move to the next node
			if (start.getLength() > distanceLeft) {
				targets.put(start, start.getParent().getHeight() - distanceLeft);
				return targets;
			} else {
				distanceLeft -= start.getLength();
				targets.putAll(collectTargets(start.getLeft(), distanceLeft, false));
				targets.putAll(collectTargets(start.getRight(), distanceLeft, false));
			}
		}
		return targets;

	}

	@Override
	public void replace(final Node node, final Node child, final Node replacement) {
		node.removeChild(child);
		node.addChild(replacement);
		node.makeDirty(2);
		replacement.makeDirty(2);
	}

	/**
	 * automatic parameter tuning *
	 */
	@Override
	public void optimize(final double logAlpha) {
		if (optimiseInput.get()) {
			double delta = calcDelta(logAlpha);
			delta += Math.log(size);
			final double f = Math.exp(delta);
//            double f = Math.exp(delta);
			if (limit > 0) {
				final Tree tree = treeInput.get();
				final double h = tree.getRoot().getHeight();
				final double k = Math.log(tree.getLeafNodeCount()) / Math.log(2);
				final double lim = (h / k) * limit;
				if (f <= lim) {
					size = f;
				}
			} else {
				size = f;
			}
		}
	}

	@Override
	public double getCoercableParameterValue() {
		return size;
	}

	@Override
	public void setCoercableParameterValue(final double value) {
		size = value;
	}

	@Override
	public String getPerformanceSuggestion() {
		final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
		final double targetProb = getTargetAcceptanceProbability();

		double ratio = prob / targetProb;

		if (ratio > 2.0)
			ratio = 2.0;
		if (ratio < 0.5)
			ratio = 0.5;

		final double newDelta = size * ratio;

		final DecimalFormat formatter = new DecimalFormat("#.###");
		if (prob < 0.10) {
			return "Try decreasing size to about " + formatter.format(newDelta);
		} else if (prob > 0.40) {
			return "Try increasing size to about " + formatter.format(newDelta);
		} else
			return "";
	}

}
