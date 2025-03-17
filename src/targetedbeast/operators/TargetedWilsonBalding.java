
package targetedbeast.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import targetedbeast.edgeweights.ConsensusWeights;
import targetedbeast.edgeweights.EdgeWeights;
import targetedbeast.likelihood.RapidTreeLikelihood;

/**
 * WILSON, I. J. and D. J. BALDING, 1998 Genealogical inference from
 * microsatellite data. Genetics 150:499-51
 * http://www.genetics.org/cgi/ijlink?linkType=ABST&journalCode=genetics&resid=150/1/499
 */
@Description("Implements the unweighted Wilson-Balding branch swapping move. "
		+ "This move is similar to one proposed by WILSON and BALDING 1998  "
		+ "and involves removing a subtree and re-attaching it on a new parent branch. "
		+ "See <a href='http://www.genetics.org/cgi/content/full/161/3/1307/F1'>picture</a>.")
public class TargetedWilsonBalding extends TreeOperator {

    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights", "input of weights to be used for targetedn tree operations", Input.Validate.REQUIRED);

    public Input<Double> mutationLimitInput = new Input<>("mutationLimit", "Input of the number of mutations to be used as a limit", 15.0);
	
    public Input<Boolean> useRepeatedMutationsInput = new Input<>("useRepeatedMutations", "flag to indicate if repeated mutations should be used as weights for operations", false);
    
    double limit;
    
    EdgeWeights edgeWeights;
    
    @Override
    public void initAndValidate() {
		limit = mutationLimitInput.get();
		edgeWeights = edgeWeightsInput.get();
    }

	/**
	 * WARNING: Assumes strictly bifurcating beast.tree.
	 */
	/**
	 * override this for proposals,
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted *
	 */
//	@Override
	public double proposal() {
        Tree tree = (Tree) InputUtil.get(treeInput, this);
        
        
//        System.out.println(((ConsensusWeights) edgeWeights).getTree() + ";");

        double logHastingsRatio = 0.0;
//        System.out.println(tree + ";");

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
//        System.out.println(totalMutations + " " + logHastingsRatio);

		
		Node i = tree.getNode(randomNode);
		
		Node p = i.getParent();
		if (p.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}

		double minHeight = i.getHeight();

		List<Integer> coExistingNodes = new ArrayList<>();
		for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;

			if (tree.getNode(k).getParent().getHeight() > minHeight)
				coExistingNodes.add(k);
		}
		if (coExistingNodes.size() == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		final Node CiP = getOtherChild(p, i);

		// get all nodes in target that are ancestors to p
		List<Integer> ancestors = new ArrayList<>();
		// make a copy of p
		Node ancestor = i.getParent();
		while (ancestor != null) {
			// check if p is in the target
			if (coExistingNodes.contains(ancestor.getNr())) {
				ancestors.add(ancestor.getNr());
			}
			ancestor = ancestor.getParent();
		}

		// calculate the consensus sequences without the node i
		edgeWeights.prestore();
		edgeWeights.updateByOperatorWithoutNode(i.getNr(), ancestors);

		// remove p as potential targets
		coExistingNodes.remove(coExistingNodes.indexOf(p.getNr()));
		try {
			coExistingNodes.remove(coExistingNodes.indexOf(i.getNr()));
		} catch (Exception e) {
			System.err.println("couldn't find node " + i.getNr() + " among coexisting nodes ");
			System.err.println(tree +";");
			return Double.NEGATIVE_INFINITY;
		}
//		coExistingNodes.remove(coExistingNodes.indexOf(CiP.getNr()));

		double[] distance = edgeWeights.getTargetWeightsInteger(i.getNr(), coExistingNodes);
		
		double siblingDistance = distance[coExistingNodes.indexOf(CiP.getNr())];
		siblingDistance = Math.pow(siblingDistance, 2);
		distance[coExistingNodes.indexOf(CiP.getNr())] = 0;
		
		
		double totalDistance = 0;
		for (int k = 0; k < coExistingNodes.size(); k++) {
			distance[k] = Math.pow(distance[k], 2);
			totalDistance += distance[k];				
		}
		
		Node j = tree.getRoot();
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
		
		logHastingsRatio -= Math.log(distance[nodeNr] / totalDistance);
		logHastingsRatio += Math.log(siblingDistance / (totalDistance + siblingDistance - distance[nodeNr]));	      
      
		edgeWeights.reset();

		Node jP = j.getParent();

		// disallow moves that change the root.
		if (j.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}

		final int pnr = p.getNr();
		final int jPnr = jP.getNr();

		if (jPnr == pnr || j.getNr() == pnr || jPnr == i.getNr()) {
			return Double.NEGATIVE_INFINITY;
		}

		final Node PiP = p.getParent();

		double newMinAge = Math.max(i.getHeight(), j.getHeight());
		double newRange = jP.getHeight() - newMinAge;
		double newAge = newMinAge + (Randomizer.nextDouble() * newRange);
		double oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
		double oldRange = PiP.getHeight() - oldMinAge;
		logHastingsRatio += Math.log(newRange / Math.abs(oldRange));
		
		p.setHeight(newAge);

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
			if (iup.getHeight() < jup.getHeight()) {
				assert !iup.isRoot();
				iup = iup.getParent();
			} else {
				assert !jup.isRoot();
				jup = jup.getParent();
			}
		}
		jup.makeDirty(3 - jup.isDirty());
		
		
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
    	
    	logHastingsRatio += Math.log(edgeWeights.getEdgeWeights(randomNode)/ totalMutations);
    	
        return logHastingsRatio;
	}

//	@Override
	public double proposal2() {
        Tree tree = (Tree) InputUtil.get(treeInput, this);
        
        double[] weights = getWeights();
        
//        System.out.println(((ConsensusWeights) edgeWeights).getTree() + ";");

        double logHastingsRatio = 0.0;
//        System.out.println(tree + ";");

        // choose a random node avoiding root
        double totalMutations = 0;
    	for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isRoot())
				continue;			
			totalMutations += weights[i];		
    	}
    	
        double scaler = Randomizer.nextDouble() * totalMutations;
        int randomNode = -1;
        double currMuts = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
        	currMuts +=  weights[i];
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
        }
                
        logHastingsRatio -= Math.log(weights[randomNode] / totalMutations);
//        System.out.println(totalMutations + " " + logHastingsRatio);

		
		Node i = tree.getNode(randomNode);
		
		Node p = i.getParent();
		if (p.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}

		double minHeight = i.getHeight();

		List<Integer> coExistingNodes = new ArrayList<>();
		for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;

			if (tree.getNode(k).getParent().getHeight() > minHeight)
				coExistingNodes.add(k);
		}
		if (coExistingNodes.size() == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		final Node CiP = getOtherChild(p, i);

		// get all nodes in target that are ancestors to p
		List<Integer> ancestors = new ArrayList<>();
		// make a copy of p
		Node ancestor = i.getParent();
		while (ancestor != null) {
			// check if p is in the target
			if (coExistingNodes.contains(ancestor.getNr())) {
				ancestors.add(ancestor.getNr());
			}
			ancestor = ancestor.getParent();
		}

		// calculate the consensus sequences without the node i
		edgeWeights.prestore();
		edgeWeights.updateByOperatorWithoutNode(i.getNr(), ancestors);

		// remove p as potential targets
		coExistingNodes.remove(coExistingNodes.indexOf(p.getNr()));
		coExistingNodes.remove(coExistingNodes.indexOf(i.getNr()));
//		coExistingNodes.remove(coExistingNodes.indexOf(CiP.getNr()));

		double[] distance = edgeWeights.getTargetWeightsInteger(i.getNr(), coExistingNodes);
		
		double siblingDistance = distance[coExistingNodes.indexOf(CiP.getNr())];
		siblingDistance *= siblingDistance;
		distance[coExistingNodes.indexOf(CiP.getNr())] = 0;
		
		
		double totalDistance = 0;
		for (int k = 0; k < coExistingNodes.size(); k++) {
			distance[k] *= distance[k]; // square the dist
			totalDistance += distance[k];				
		}
		
		
		
		Node j = tree.getRoot();
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
		
		

		logHastingsRatio -= Math.log(distance[nodeNr] / totalDistance);
//		System.out.println("hr2: " + logHastingsRatio);
		logHastingsRatio += Math.log(siblingDistance / (totalDistance + siblingDistance - distance[nodeNr]));
//		System.out.println("hr3: " + logHastingsRatio);
		
		
		double diff = 1/distance[nodeNr] - 1/siblingDistance;

		
//      System.out.println("hr3: " + logHastingsRatio);
      
//      System.out.println(Arrays.toString(distance));
//      System.out.println(nodeNr + " " + coExistingNodes.indexOf(CiP.getNr()));

		      
      
		edgeWeights.reset();

		Node jP = j.getParent();

		// disallow moves that change the root.
		if (j.isRoot()) {
			return Double.NEGATIVE_INFINITY;
		}

		final int pnr = p.getNr();
		final int jPnr = jP.getNr();

		if (jPnr == pnr || j.getNr() == pnr || jPnr == i.getNr()) {
			return Double.NEGATIVE_INFINITY;
		}

		final Node PiP = p.getParent();

		double newMinAge = Math.max(i.getHeight(), j.getHeight());
		double newRange = jP.getHeight() - newMinAge;
		double newAge = newMinAge + (Randomizer.nextDouble() * newRange);
		double oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
		double oldRange = PiP.getHeight() - oldMinAge;
		logHastingsRatio += Math.log(newRange / Math.abs(oldRange));
		
//		System.out.println("hr4: " + logHastingsRatio);

		p.setHeight(newAge);
//		System.out.println(i.getLength());
//
//    	System.out.println(tree +";");

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
			if (iup.getHeight() < jup.getHeight()) {
				assert !iup.isRoot();
				iup = iup.getParent();
			} else {
				assert !jup.isRoot();
				jup = jup.getParent();
			}
		}
		jup.makeDirty(3 - jup.isDirty());
		
		
        edgeWeights.prestore();
        edgeWeights.updateByOperator();
        
        weights = getWeights();
        // recalculate hasting ratio
        // calculate tree length
        totalMutations = 0;
    	for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;			
			totalMutations += weights[k];	
    	}
    	
    	logHastingsRatio += Math.log(weights[randomNode]/ totalMutations);
    	    	
        return logHastingsRatio;
	}


	
//	@Override
//	public void accept() {
//		accepted.add(lastNode);
//		acceptedKids.add(lastNoKids);
////		System.out.println("A  " + (posteriorInput.get().calculateLogP()-postBefore) + " " + (accepted.size()/(0.0 + accepted.size() + rejected.size())) + " " + lastNode );
////		System.out.println();
////		System.out.println(treeBefore +";");
////		System.out.println(treeAfter +";");
//	}
//	
//	@Override
//	public void reject(int reason) {
//		rejected.add(lastNode);
//		rejectedKids.add(lastNoKids);
////		System.out.println(reason + " " + (posteriorInput.get().calculateLogP()-postBefore) + " " + (accepted.size()/(0.0 + accepted.size() + rejected.size())) + " " + lastNode );
//
//	}
	

	double[] getWeights() {
		double[] sumMuts = new double[edgeWeights.getNodeConsensus(0).length / 4]; 
		for (int i = 0; i < treeInput.get().getNodeCount(); i++) {
			Node n = treeInput.get().getNode(i);
			if (treeInput.get().getNode(i).isRoot())
				continue;
			
			byte[] parentSeq = edgeWeights.getNodeConsensus(n.getParent().getNr());
			byte[] seq = edgeWeights.getNodeConsensus(n.getNr());
			
			for (int j = 0; j < seq.length; j++) {
				if (seq[j]!=1) {
					sumMuts[j/4] += Math.abs(parentSeq[j] - seq[j]);					
				}			
			}		
		}
		
		double sums = 0;
		for (int j = 0; j < sumMuts.length; j++) {
			sumMuts[j]= Math.max(sumMuts[j]/12-1, 0.1);
			sums += sumMuts[j];
		}		
		
		double[] weights = new double[treeInput.get().getNodeCount()];
		// calculate the weights for each edge based on the number of dublicate mutations
		
		for (int i = 0; i < treeInput.get().getNodeCount(); i++) {
			Node n = treeInput.get().getNode(i);
			if (treeInput.get().getNode(i).isRoot())
				continue;
			
			byte[] parentSeq = edgeWeights.getNodeConsensus(n.getParent().getNr());
			byte[] seq = edgeWeights.getNodeConsensus(n.getNr());
//			byte[] otherSeq = edgeWeights.getNodeConsensus(getOtherChild(n.getParent(), n).getNr());
			weights[i] = 0.1;
			for (int j = 0; j < seq.length; j++) {
				if (sumMuts[j/4]>=1) {
					weights[i] += Math.abs(parentSeq[j] - seq[j]);
				}			
				
			}		
				
		}
//		System.out.println(treeInput.get() + ";");
		
//		System.exit(0);
		return weights;
	}
	
	@Override
	public void replace(final Node node, final Node child, final Node replacement) {
		node.removeChild(child);
		node.addChild(replacement);
		node.makeDirty(2);
		replacement.makeDirty(2);
	}

} // class WilsonBalding
