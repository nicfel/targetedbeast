/*
 * WilsonBalding.java
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
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import targetedbeast.likelihood.RapidTreeLikelihood;

/**
 * WILSON, I. J. and D. J. BALDING, 1998  Genealogical inference from microsatellite data.
 * Genetics 150:499-51
 * http://www.genetics.org/cgi/ijlink?linkType=ABST&journalCode=genetics&resid=150/1/499
 */
@Description("Implements the unweighted Wilson-Balding branch swapping move. " +
        "This move is similar to one proposed by WILSON and BALDING 1998  " +
        "and involves removing a subtree and re-attaching it on a new parent branch. " +
        "See <a href='http://www.genetics.org/cgi/content/full/161/3/1307/F1'>picture</a>.")
public class Uniform2 extends TreeOperator {
	
    public Input<RapidTreeLikelihood> rapidTreeLikelihoodInput = new Input<>("rapidTreeLikelihood", "The likelihood to be used for the tree proposal. If not specified, the tree likelihood is calculated from the tree.");


    @Override
    public void initAndValidate() {
    }

    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     */
    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        Tree tree = (Tree) InputUtil.get(treeInput, this);

        double logHastingsRatio = 0.0;

        // choose a random node avoiding root
        double totalMutations = 0;
    	for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isRoot())
				continue;			
			totalMutations += rapidTreeLikelihoodInput.get().getEdgeMutations(i);		
    	}
    	
        double scaler = Randomizer.nextDouble() * totalMutations;
        int randomNode = -1;
        double currMuts = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
        	currMuts += rapidTreeLikelihoodInput.get().getEdgeMutations(i);
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
        }
        
        logHastingsRatio -= Math.log(rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode) / totalMutations);

        
        Node i = tree.getNode(randomNode);
        Node p = i.getParent();
        // choose another random node to insert i above
        Node j = tree.getRoot();
        Node jP;
        
		if (p.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}		
        
        double currHeight = p.getHeight();
        double minHeight = i.getHeight();
        
        List<Integer> coExistingNodes = new ArrayList<>();
		for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;
			
//			if (tree.getNode(k).getHeight() < currHeight && tree.getNode(k).getParent().getHeight() > currHeight)
//				coExistingNodes.add(k);
			if (tree.getNode(k).getHeight()  > minHeight)
				coExistingNodes.add(k);

		}
		if (coExistingNodes.size() == 0)
			return Double.NEGATIVE_INFINITY;
		
		double[] distance = new double[coExistingNodes.size()];
		double totalDistance = 0;
		double[] currConsensus = rapidTreeLikelihoodInput.get().getConsensus(i.getNr());
		for (int k = 0; k < coExistingNodes.size(); k++) {
			double[] consensus = rapidTreeLikelihoodInput.get().getConsensus(coExistingNodes.get(k));
			// calculate the distance between the two consensus
			distance[k] = 0;
			for (int l = 0; l < consensus.length; l++) {
				distance[k] += Math.abs(currConsensus[l] - consensus[l]);
			}			
			distance[k] = 1/(distance[k]+0.1);
			totalDistance += distance[k];				
		}
		
		
		double scaler2 = Randomizer.nextDouble() * totalDistance;
		double currDist = 0;
		int nodeNr = -1;
		for (int k = 0; k < coExistingNodes.size(); k++) {
			currDist += distance[k];
			if (currDist > scaler2) {
				nodeNr = k;
				j = treeInput.get().getNode(coExistingNodes.get(nodeNr));
				break;
			}
		}
		
		try {
			logHastingsRatio -= Math.log(distance[nodeNr] / totalDistance);
		} catch (Exception e) {
			System.out.println("nodeNr = " + nodeNr);
			System.out.println("totalDistance = " + totalDistance);
			System.out.println("distance = " + Arrays.toString(distance));
			System.out.println("scaler2 = " + distance[nodeNr]);
		}
        final Node CiP = getOtherChild(p, i);
        
        totalDistance -= distance[nodeNr];        
		double[] consensus = rapidTreeLikelihoodInput.get().getConsensus(CiP.getNr());
		distance[nodeNr]=0;
		for (int l = 0; l < consensus.length; l++) {
			distance[nodeNr] += Math.abs(currConsensus[l] - consensus[l]);
		}			
		distance[nodeNr] = 1/(distance[nodeNr]+0.1);		
		totalDistance += distance[nodeNr];
		logHastingsRatio += Math.log(distance[nodeNr] / totalDistance);

		
		jP = j.getParent();
    
        // disallow moves that change the root.
        if (j.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        assert jP != null;  // j != root tested above
        final int pnr = p.getNr();
        final int jPnr = jP.getNr();
        
        if ( jPnr == pnr || j.getNr() == pnr || jPnr == i.getNr())
            return Double.NEGATIVE_INFINITY;


        final Node PiP = p.getParent();
        
//    	System.out.println(rapidTreeLikelihoodInput.get().toNewick(treeInput.get().getRoot()) + " ;");

        double newMinAge = Math.max(i.getHeight(), j.getHeight());
        double newRange = jP.getHeight() - newMinAge;
        double newAge = newMinAge + (Randomizer.nextDouble() * newRange);
        double oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
        double oldRange = PiP.getHeight() - oldMinAge;
        logHastingsRatio += Math.log(newRange / Math.abs(oldRange));

//        if (oldRange == 0 || newRange == 0) {
            // This happens when some branch lengths are zero.
            // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
            // node i can be catapulted anywhere in the tree, resulting in
            // very bad trees that are always accepted.
            // For symmetry, newRange = 0 should therefore be ruled out as well
//            return Double.NEGATIVE_INFINITY;
//        }

        // root changing moves disallowed earlier. Comment out just like in BEAST 1

        //update
//        if (j.isRoot()) {
//            // 1. remove edges <p, CiP>
//            // 2. add edges <k, p>, <p, j>, <PiP, CiP>
//
//            replace(p, CiP, j);
//            replace(PiP, p, CiP);
//
//            // p is the new root
//            p.setParent(null);
//            tree.setRoot(p);
//
//        } else if (p.isRoot()) {
//            // 1. remove edges <k, j>, <p, CiP>, <PiP, p>
//            // 2. add edges <k, p>, <p, j>, <PiP, CiP>
//
//            replace(jP, j, p);
//            //replace(p, CiP, p);
//            replace(p, CiP, j);
//
//            // CiP is the new root
//            CiP.setParent(null);
//            tree.setRoot(CiP);
//
//        } else {
//            // 1. remove edges <k, j>, <p, CiP>, <PiP, p>
            // 2. add edges <k, p>, <p, j>, <PiP, CiP>

            // disconnect p
//            final Node pP = p.getParent();
            // remove p from PiP, add the other child of p to PiP, frees p
            replace(PiP, p, CiP);
            // add j as child to p, removing CiP as child of p
            replace(p, CiP, j);
            // then parent node of j to p
            replace(jP, j, p);

        // mark paths to common ancestor as changed
            if( markCladesInput.get() ) {
                Node iup = PiP;
                Node jup = p;
                while (iup != jup) {
                    if( iup.getHeight() < jup.getHeight() ) {
                        assert !iup.isRoot();
                        iup = iup.getParent();
                        iup.makeDirty(Tree.IS_FILTHY);
                    } else {
                        assert !jup.isRoot();
                        jup = jup.getParent();
                        jup.makeDirty(Tree.IS_FILTHY);
                    }
                }
                jup.makeDirty(3-i.isDirty());
            }

//        }

        p.setHeight(newAge);
            
//        while (CiP!=j) {
//			if (CiP.getHeight() < j.getHeight()) {
//				CiP = i.getParent();
//			} else {
//				j = j.getParent();
//			}
//        }
//        
//        i.makeDirty(3-i.isDirty());
        
        
        rapidTreeLikelihoodInput.get().prestore();
        rapidTreeLikelihoodInput.get().updateByOperator();
        // recalculate hasting ratio
        // calculate tree length
        totalMutations = 0;
    	for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;			
			totalMutations += rapidTreeLikelihoodInput.get().getEdgeMutations(k);	
    	}
    	
    	logHastingsRatio += Math.log(rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode) / totalMutations);

    	
    	rapidTreeLikelihoodInput.get().unstore();
    	    	
    	
        
        return logHastingsRatio;
    }


} // class WilsonBalding
