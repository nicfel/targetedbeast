
package targetedbeast.operators;


import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import targetedbeast.likelihood.RapidTreeLikelihood;


@Description("This operator changes a beast.tree.")
abstract public class TargetedTreeOperator extends TreeOperator {
    public Input<RapidTreeLikelihood> rapidTreeLikelihoodInput = new Input<>("rapidTreeLikelihood", "The likelihood to be used for the tree proposal. If not specified, the tree likelihood is calculated from the tree.");
    public Input<Double> mutationLimitInput = new Input<>("mutationLimit", "Input of the number of mutations to be used as a limit", 5.0);

    double mutationLimit;
    double logHastingsRatio;
    
    @Override
    public void initAndValidate() {
    	mutationLimit = mutationLimitInput.get();    
    }
    
    protected Node getNodeToOperate() {
		double totalMutations = 0;
		for (int i = 0; i < treeInput.get().getNodeCount(); i++) {
			if (treeInput.get().getNode(i).isRoot())
				continue;
			totalMutations += Math.min(mutationLimit, rapidTreeLikelihoodInput.get().getEdgeMutations(i) + 0.25);
		}
		
		
		double scaler = Randomizer.nextDouble() * totalMutations;
		int randomNode = -1;
		double currMuts = 0;
		for (int i = 0; i < treeInput.get().getNodeCount(); i++) {
			currMuts += Math.min(mutationLimit, rapidTreeLikelihoodInput.get().getEdgeMutations(i) + 0.25);
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
		}
		logHastingsRatio -= Math.log(Math.min(mutationLimit, rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode) + 0.25) / totalMutations);
		return treeInput.get().getNode(randomNode);
    }

    protected void getReverseNodeToOperateContribution(int nodeNo) {
		rapidTreeLikelihoodInput.get().updateByOperator();

		// recalculate hasting ratio
		// calculate tree length
		double totalMutations = 0;
		for (int k = 0; k < treeInput.get().getNodeCount(); k++) {
			if (treeInput.get().getNode(k).isRoot())
				continue;
			totalMutations += Math.min(mutationLimit, rapidTreeLikelihoodInput.get().getEdgeMutations(k) + 0.25);
		}

		logHastingsRatio += Math
				.log(Math.min(mutationLimit, rapidTreeLikelihoodInput.get().getEdgeMutations(nodeNo) + 0.25) / totalMutations);

    }
    
    protected void doReplacementMove(Node PiP, Node p, Node CiP, Node j, Node jP) {
    
		// remove p from PiP, add the other child of p to PiP, frees p
		replace(PiP, p, CiP);
		// add j as child to p, removing CiP as child of p
		replace(p, CiP, j);
		// then parent node of j to p
		replace(jP, j, p);

    	
	    Node iup = PiP;
	    Node jup = p;
	    while (iup != jup) {
	        if( iup.getHeight() < jup.getHeight() ) {
	            assert !iup.isRoot();
	            iup = iup.getParent();
	        } else {
	            assert !jup.isRoot();
	            jup = jup.getParent();
	        }
	    }
	    jup.makeDirty(3-jup.isDirty());
    }

	protected double computeDistances(List<Integer> potentialTargets, double[] distance, Node i) {
		double[] currConsensus = rapidTreeLikelihoodInput.get().getConsensus(i.getNr());
		double totalDistance = 0;
		for (int k = 0; k < potentialTargets.size(); k++) {
			double[] consensus = rapidTreeLikelihoodInput.get().getConsensus(potentialTargets.get(k));
			// calculate the distance between the two consensus
			double sum = 0;
			for (int l = 0; l < consensus.length; l++) {
				sum += Math.abs(currConsensus[l] - consensus[l]);
			}
			distance[k] = 1 / (sum + 0.5);
			totalDistance += distance[k];
		}

		return totalDistance;
	}



}
