/*
* File Uniform.java
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
 * UniformOperator.java
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

package targetedbeast.altoperators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import targetedbeast.likelihood.RapidTreeLikelihood;

@Description("Randomly selects true internal tree node (i.e. not the root) and move node height uniformly in interval "
		+ "restricted by the nodes parent and children.")
public class Uniform extends TreeOperator {
	public final Input<Integer> numberOfAttemptsInput = new Input<>("attempts",
			"number of attempts to make before giving up", 50);
	public Input<RapidTreeLikelihood> rapidTreeLikelihoodInput = new Input<>("rapidTreeLikelihood",
			"The likelihood to be used for the tree proposal. If not specified, the tree likelihood is calculated from the tree.");

	// empty constructor to facilitate construction by XML + initAndValidate
	public Uniform() {
	}

	public Uniform(Tree tree) {
		try {
			initByName(treeInput.getName(), tree);
		} catch (Exception e) {
			e.printStackTrace(); // To change body of catch statement use File | Settings | File Templates.
			throw new RuntimeException("Failed to construct Uniform Tree Operator.");
		}
	}

	@Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		Tree tree = treeInput.get();
		tree.startEditing(this);
		
//		System.out.println(rapidTreeLikelihoodInput.get().getTree() + ";");

		// Abort if no non-root internal nodes
		if (tree.getInternalNodeCount() == 1)
			return Double.NEGATIVE_INFINITY;
		
		// check potential nodes based on below and above mutations being <0.01
		List<Node> potentialNodes = new ArrayList<>();
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isLeaf() || tree.getNode(i).isRoot())
				continue;
			Node node = tree.getNode(i);
			if (rapidTreeLikelihoodInput.get().getEdgeMutations(node.getNr()) < 0.01 && 
					(rapidTreeLikelihoodInput.get().getEdgeMutations(node.getLeft().getNr()) < 0.01 ||
					rapidTreeLikelihoodInput.get().getEdgeMutations(node.getRight().getNr()) < 0.01)) {
				potentialNodes.add(node);
			}
		}
		
        // Abort if no internal nodes
		double logHastingsRatio = 0.0;

		for (int a = 0; a < potentialNodes.size(); a++) {

			Node node = potentialNodes.get(Randomizer.nextInt(potentialNodes.size()));
			

			if (node.isRoot()) {
				return Double.NEGATIVE_INFINITY;
			}

			final double upper = node.getParent().getHeight();
			final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
			final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;

			node.setHeight(newValue);
		}

//		System.out.println(rapidTreeLikelihoodInput.get().getTree() + ";\n");
		
		return logHastingsRatio;

	}

}
