
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

	public Input<RapidTreeLikelihood> rapidTreeLikelihoodInput = new Input<>("rapidTreeLikelihood",
			"The likelihood to be used for the tree proposal. If not specified, the tree likelihood is calculated from the tree.");
	
//	public Input<CompoundDistribution> posteriorInput = new Input<>("posterior", "The posterior distribution to be used for the tree proposal. If not specified, the tree likelihood is calculated from the tree.");

    double limit = 5.0;
	
	@Override
	public void initAndValidate() {

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
	@Override
	public double proposal() {
		Tree tree = (Tree) InputUtil.get(treeInput, this);
		double logHastingsRatio = 0.0;
		
        // choose a random node avoiding root
        double totalMutations = 0;
    	for (int i = 0; i < tree.getNodeCount(); i++) {
			if (tree.getNode(i).isRoot())
				continue;			
			totalMutations += Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(i)+0.1) ;		
    	}
        double scaler = Randomizer.nextDouble() * totalMutations;
        int randomNode = -1;
        double currMuts = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
        	currMuts +=  Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(i)+0.1);
			if (currMuts > scaler) {
				randomNode = i;
				break;
			}
        }
        
        logHastingsRatio -= Math.log(Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode)+0.1) / totalMutations);
//		
//		double limit = 1;
//		
//		double totsum = 0;
//		// choose a random node avoiding root
//		double[] deviation = new double[tree.getNodeCount()];
//		for (int k = 0; k < tree.getNodeCount(); k++) {
//			Node n = tree.getNode(k);
//			if (n.isRoot() || n.isLeaf())
//				continue;
//			sum += rapidTreeLikelihoodInput.get().getEdgeMutations(n.getNr()) > limit ? 1 : 0;
//			sum += rapidTreeLikelihoodInput.get().getEdgeMutations(n.getLeft().getNr()) > limit ? 1 : 0;
//			sum += rapidTreeLikelihoodInput.get().getEdgeMutations(n.getRight().getNr()) > limit ? 1 : 0;
//            deviation[k] = sum;
//			totsum += ;
//		}
//				
//		double scaler = Randomizer.nextDouble() * totsum;
//		int randomNode = -1;
//		double currMuts = 0;
//		for (int k = 0; k < tree.getNodeCount(); k++) {
//			currMuts += deviation[k];
//			if (currMuts > scaler) {
//				randomNode = k;
//				break;
//			}
//		}
//		logHastingsRatio -= Math
//				.log(deviation[randomNode] / totsum);

//		double ssum = -1;
//		ssum += rapidTreeLikelihoodInput.get().getEdgeMutations(tree.getNode(randomNode).getNr()) > limit ? 1 : 0;
//		ssum += rapidTreeLikelihoodInput.get().getEdgeMutations(tree.getNode(randomNode).getLeft().getNr()) > limit ? 1 : 0;
//		ssum += rapidTreeLikelihoodInput.get().getEdgeMutations(tree.getNode(randomNode).getRight().getNr()) > limit ? 1 : 0;
//
//		lastNode += Math.min(ssum, 0.1);
//
//		lastNoKids = getChildNodeCount(tree.getNode(randomNode));
		
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
		rapidTreeLikelihoodInput.get().prestore();
		rapidTreeLikelihoodInput.get().updateByOperatorWithoutNode(i.getNr(), ancestors);

		// remove p as potential targets
		coExistingNodes.remove(coExistingNodes.indexOf(p.getNr()));
//		coExistingNodes.remove(coExistingNodes.indexOf(CiP.getNr()));

		// calculate the distance between the consensus of i and the other nodes after
		// removing i from the consensus
		double[] distance = new double[coExistingNodes.size()];
		double totalDistance = 0;
		double[] currConsensus = rapidTreeLikelihoodInput.get().getConsensus(i.getNr());
		for (int k = 0; k < coExistingNodes.size(); k++) {
			double[] consensus = rapidTreeLikelihoodInput.get().getConsensus(coExistingNodes.get(k));
			// get the common ancestor node between i and coExistingNodes.get(k)
//			Node commonAncestor = getCommonAncestor(tree, i, tree.getNode(coExistingNodes.get(k)));
//			int[] calcForPatterns = rapidTreeLikelihoodInput.get().getCalcPatterns(commonAncestor.getNr());			
			
			// calculate the distance between the two consensus
			double sum = 1;
			for (int l = 0; l < consensus.length; l++) {				
				sum += Math.abs(currConsensus[l] - consensus[l]);
			}
			
			distance[k] = 1 / (sum);
			totalDistance += distance[k];
		}

		// choose another random node to insert i above
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

		// pick a random node
		logHastingsRatio -= Math.log(distance[nodeNr] / totalDistance);
		logHastingsRatio += Math.log(distance[coExistingNodes.indexOf(CiP.getNr())] / totalDistance);
		rapidTreeLikelihoodInput.get().reset();

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

		rapidTreeLikelihoodInput.get().updateByOperator();
		
        // recalculate hasting ratio
        // calculate tree length
        totalMutations = 0;
    	for (int k = 0; k < tree.getNodeCount(); k++) {
			if (tree.getNode(k).isRoot())
				continue;			
			totalMutations += Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(k)+0.1);	
    	}
    	
    	logHastingsRatio += Math.log(Math.min(limit, rapidTreeLikelihoodInput.get().getEdgeMutations(randomNode)+0.1)/ totalMutations);

    	return logHastingsRatio;
	}

	private int getChildNodeCount(Node node) {
		int children = 0;
		if (node.isLeaf()) {
			return 1;
		}else {
			children += getChildNodeCount(node.getLeft());
			children += getChildNodeCount(node.getRight());
		}
		return children;
	}

	private Node getCommonAncestor(Tree tree, Node i, Node j) {
		Node ii = i;
		Node jj = j;
		while (ii != jj) {
            if (ii.getHeight() < jj.getHeight()) {
                ii = ii.getParent();
            } else {
                jj = jj.getParent();
            }
        }
		return ii;
	}
	
	
	// calculate double mutations
	private double[] calculateDoubleMuts(Tree tree) {
		double[] totalMuts;
		
		double[] c = rapidTreeLikelihoodInput.get().getConsensus(0);
		totalMuts = new double[c.length];

		for (int i = 0; i < tree.getNodeCount(); i++) {
			Node n = tree.getNode(i);
			if (n.isRoot())
                continue;
			
			int[] calcForPattern = rapidTreeLikelihoodInput.get().getCalcPatterns(n.getParent().getNr());
			
			double[] consensus1 = rapidTreeLikelihoodInput.get().getConsensus(n.getNr());
			double[] consensus2 = rapidTreeLikelihoodInput.get().getConsensus(n.getParent().getNr());
			
			for (int j = 0; j < consensus1.length; j++) {
				totalMuts[j] += Math.abs(consensus1[j] - consensus2[j]);
			}
		}
		// get all positions with total Muts > 2
		List<Integer> doubleMuts = new ArrayList<>();
		for (int i = 0; i < totalMuts.length; i++) {
            if (totalMuts[i] > 2) {
                doubleMuts.add(i);
            }
        }
		
		double[] offset = new double[tree.getNodeCount()];
		for (Integer i = 0; i < tree.getNodeCount(); i++) {
			Node n = tree.getNode(i);
			if (n.isRoot())
				continue;
			double sum = 0.1;
			
			double[] consensus1 = rapidTreeLikelihoodInput.get().getConsensus(n.getNr());
			double[] consensus2 = rapidTreeLikelihoodInput.get().getConsensus(n.getParent().getNr());

			for (Integer j = 0; j < doubleMuts.size(); j++) {
				sum += Math.abs(consensus1[doubleMuts.get(j)] - consensus2[doubleMuts.get(j)]);
			}
			offset[i] = sum;
		}
		return offset;
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
	
	
	@Override
	public void replace(final Node node, final Node child, final Node replacement) {
		node.removeChild(child);
		node.addChild(replacement);
		node.makeDirty(2);
		replacement.makeDirty(2);
	}

} // class WilsonBalding
