
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
import targetedbeast.edgeweights.ConsensusWeights;
import targetedbeast.edgeweights.EdgeWeights;
import targetedbeast.likelihood.RapidTreeLikelihood;
@Description("Picks a node, looks for all coexisting lineages and then performs a weighted move")
public class WeightedWideOperator extends TreeOperator {
	
    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights", "input of weights to be used for targetedn tree operations");

    public Input<Double> mutationLimitInput = new Input<>("mutationLimit", "Input of the number of mutations to be used as a limit", 15.0);

    public 

    double limit;
    
    EdgeWeights edgeWeights;
    
    @Override
    public void initAndValidate() {
		limit = mutationLimitInput.get();
		edgeWeights = edgeWeightsInput.get();
    }
    
    @Override
    public double proposal() {
		if (edgeWeights != null) {
			return weightedProposal();
		}else {
			return unweightedProposal();
		}
    }


    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     */
    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    private double weightedProposal() {
        Tree tree = (Tree) InputUtil.get(treeInput, this);

        double logHastingsRatio = 0.0;
//        System.out.println(tree + ";");
//        System.out.println(" ");
//        System.out.println(((ConsensusWeights) edgeWeights).getTree() + ";");

        // choose a random node avoiding root
        double totalMutations = 0;
    	for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isRoot())
				continue;			
			totalMutations += edgeWeights.getEdgeWeights(i);		
    	}
        double scaler = Randomizer.nextDouble() * totalMutations;
        int randomNode = -1;
        double currMuts = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
        	currMuts +=  edgeWeights.getEdgeWeights(i);
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
        }
        
        logHastingsRatio -= Math.log(edgeWeights.getEdgeWeights(randomNode) / totalMutations);       
        
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
		
		double[] distance = edgeWeights.getTargetWeightsInteger(i.getNr(), coExistingNodes);
		double totalDistance = 0;
		for (int k = 0; k < coExistingNodes.size(); k++) {
			distance[k] = distance[k]* distance[k];
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
                
        edgeWeights.prestore();
        edgeWeights.updateByOperator();
        // recalculate hasting ratio
        // calculate tree length
        totalMutations = 0;
    	for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;			
			totalMutations += edgeWeights.getEdgeWeights(k);	
    	}
    	
//        System.out.println(((ConsensusWeights) edgeWeights).getTree() + ";");

//		System.out.println("B " + mutationsBefore + " A " + edgeWeights.getEdgeWeights(randomNode) + " "
//				+" B "+ 1/distance[coExistingNodes.indexOf(CiP.getNr())] + " A " + 1/distance[nodeNr]);
    	
    	logHastingsRatio += Math.log(edgeWeights.getEdgeWeights(randomNode)/ totalMutations);
//    	System.out.println(logHastingsRatio);
    	return logHastingsRatio;
    }

    private double unweightedProposal() {
        Tree tree = (Tree) InputUtil.get(treeInput, this);

        double logHastingsRatio = 0.0;

        // choose a random node avoiding root
        double scaler = Randomizer.nextDouble() * tree.getNodeCount();
        int randomNode = -1;
        double currMuts = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
        	currMuts ++;
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
        }
        
        Node i = tree.getNode(randomNode);
        
        if (i.isRoot())
            return Double.NEGATIVE_INFINITY;
        
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
		
		double scaler2 = Randomizer.nextDouble() * coExistingNodes.size();
		double currDist = 0;
		int nodeNr = -1;
		for (int k = 0; k < coExistingNodes.size(); k++) {
			currDist++;
			if (currDist > scaler2) {				
				nodeNr = k;
				j = treeInput.get().getNode(coExistingNodes.get(nodeNr));
				break;
			}
		}
				
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
                
    	return logHastingsRatio;
    }


} // class WilsonBalding
