
package targetedbeast.operators;

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
@Description("Picks a node, looks for all coexisting lineages and then performs a weighted move")
public class WeightedWideOperator extends TreeOperator {
	
    public Input<RapidTreeLikelihood> rapidTreeLikelihoodInput = new Input<>("rapidTreeLikelihood", "The likelihood to be used for the tree proposal. If not specified, the tree likelihood is calculated from the tree.");
    public Input<Double> mutationLimitInput = new Input<>("mutationLimit", "Input of the number of mutations to be used as a limit", 5.0);


    double limit;
    
    @Override
    public void initAndValidate() {
		limit = mutationLimitInput.get();
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
			totalMutations += Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(i)+0.25) ;		
    	}
        double scaler = Randomizer.nextDouble() * totalMutations;
        int randomNode = -1;
        double currMuts = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
        	currMuts +=  Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(i)+0.25);
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
        }
        
        logHastingsRatio -= Math.log(Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode)+0.25) / totalMutations);

        
        Node i = tree.getNode(randomNode);
        Node p = i.getParent();
        // choose another random node to insert i above
        Node j = tree.getRoot();
        Node jP;
        final Node CiP = getOtherChild(p, i);

        
		if (p.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}		
        
        double currHeight = p.getHeight();
        
        List<Integer> coExistingNodes = new ArrayList<>();
		for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;
			
			if (tree.getNode(k).getHeight() < currHeight && tree.getNode(k).getParent().getHeight() > currHeight)
				coExistingNodes.add(k);
		}
		if (coExistingNodes.size() == 0)
			return Double.NEGATIVE_INFINITY;
		
		// add the current node to the list
		coExistingNodes.add(CiP.getNr());
		
		double[] distance = new double[coExistingNodes.size()];
		double totalDistance = 0;
		double[] currConsensus = rapidTreeLikelihoodInput.get().getConsensus(i.getNr());
		for (int k = 0; k < coExistingNodes.size(); k++) {
			double[] consensus = rapidTreeLikelihoodInput.get().getConsensus(coExistingNodes.get(k));
			// calculate the distance between the two consensus
			double sum = 0;
			for (int l = 0; l < consensus.length; l++) {
				sum += Math.abs(currConsensus[l] - consensus[l]);
			}			
			distance[k] = 1/(sum+1);
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
		
//		if (j.getNr() == CiP.getNr())
//			return Double.NEGATIVE_INFINITY;
		
		// calculate HR contribution
		logHastingsRatio -= Math.log(distance[nodeNr] / totalDistance);
		logHastingsRatio += Math.log(distance[coExistingNodes.indexOf(CiP.getNr())] / totalDistance);

		
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

        // disconnect p
//            final Node pP = p.getParent();
        // remove p from PiP, add the other child of p to PiP, frees p
        replace(PiP, p, CiP);
        // add j as child to p, removing CiP as child of p
        replace(p, CiP, j);
        // then parent node of j to p
        replace(jP, j, p);

        // mark paths to common ancestor as changed
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
        
        
        rapidTreeLikelihoodInput.get().prestore();
        rapidTreeLikelihoodInput.get().updateByOperator();
        // recalculate hasting ratio
        // calculate tree length
        totalMutations = 0;
    	for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;			
			totalMutations +=  Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(k)+0.25);	
    	}
    	
    	logHastingsRatio += Math.log(Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode)+0.25)/ totalMutations);

    	
    	rapidTreeLikelihoodInput.get().unstore();
    	    	        
        return logHastingsRatio;
    }


} // class WilsonBalding
