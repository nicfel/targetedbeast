package targetedbeast.coalescent;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.CalculationNode;
import beast.base.inference.State;
import beast.base.util.Binomial;


/**
 * @author Alexei Drummond
 */

@Description("Calculates the probability of a beast.tree conditional on a population size function. " +
        "Note that this does not take the number of possible tree interval/tree topology combinations " +
        "in account, in other words, the constant required for making this a proper distribution that integrates " +
        "to unity is not calculated (partly, because we don't know how for sequentially sampled data).")
public class RapidCoalescent extends TreeDistribution {
    final public Input<PopulationFunction> popSizeInput = new Input<>("populationModel", "A population size model", Validate.REQUIRED);


    @Override
    public void initAndValidate() {
        calculateLogP();
    }


    /**
     * do the actual calculation *
     */
    @Override
    public double calculateLogP() {

        logP = calculateLogLikelihood(popSizeInput.get());

        if (Double.isInfinite(logP)) {
        	logP = Double.NEGATIVE_INFINITY;
        }

        return logP;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO this should eventually sample a coalescent tree conditional on population size function
        throw new UnsupportedOperationException("This should eventually sample a coalescent tree conditional on population size function.");
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return Collections.singletonList(treeIntervalsInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
    public List<String> getConditions() {
        return popSizeInput.get().getParameterIds();
    }


    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     *
     * @param intervals       the intervals whose likelihood is computed
     * @param popSizeFunction the population size function
     * @return the log likelihood of the intervals given the population size function
     */
    public double calculateLogLikelihood(PopulationFunction popSizeFunction) {
        return calculateLogLikelihood(popSizeFunction, 0.0);
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a population size function.
     *
     * @param intervals       the intervals whose likelihood is computed
     * @param popSizeFunction the population size function
     * @param threshold       the minimum allowable coalescent interval size; negative infinity will be returned if
     *                        any non-zero intervals are smaller than this
     * @return the log likelihood of the intervals given the population size function
     */
    public double calculateLogLikelihood(PopulationFunction popSizeFunction, double threshold) { 	
    	
    	if (treeIntervalsInput.get().treeInput.get().somethingIsDirty()) {
    	}
    	
        double logL = 0.0;
        
        // get the tree
    	Tree tree = treeIntervalsInput.get().treeInput.get();
    	
        List<Node> nodes = new ArrayList<>();
        nodes.add(tree.getRoot().getLeft());
        nodes.add(tree.getRoot().getRight());
//        List<Boolean> isDirty = new ArrayList<>();
//        isDirty.add(treeInput.get().getRoot().isDirty()!=Tree.IS_CLEAN);     
        
        final int nodeCount = tree.getNodeCount();
        int c = nodeCount - 2; // skip the first interval
        double prev_time = tree.getRoot().getHeight();
        double curr_time;
        int old_lineageCount = 2;
        int new_lineageCount = 2;
        
        logL -= Math.log(popSizeFunction.getPopSize(prev_time));

        while (c > 0) {       
        	int index = 0;
			curr_time = nodes.get(index).getHeight();


			if (nodes.get(index).isLeaf()) {
				new_lineageCount--;
				nodes.remove(index);
			}else {								
				Node left = nodes.get(index).getLeft();
				Node right = nodes.get(index).getRight();
				nodes.remove(index);
				insertOrderedNode(nodes, left);
				insertOrderedNode(nodes, right);
				new_lineageCount++;
				
				logL -= Math.log(popSizeFunction.getPopSize(curr_time));				
			}			
			c--;
			
			final double intervalArea = popSizeFunction.getIntegral(curr_time, prev_time);
			
            final double kChoose2 = Binomial.choose2(old_lineageCount);
            // common part
            logL += -kChoose2 * intervalArea;
            
			if (old_lineageCount == 2) {
				logL -= Math.log(popSizeFunction.getPopSize(curr_time));
			}
			old_lineageCount = new_lineageCount;
			prev_time = curr_time;
			
		}
        return logL;
    }
    
    private void insertOrderedNode(List<Node> nodes, Node node) {
        for (int i = 0; i < nodes.size(); i++) {
			if (node.getHeight() > nodes.get(i).getHeight()) {
				nodes.add(i, node);
				return;
			}
        }    	
        nodes.add(node);
    }

    private int getMaxNextTime(List<Node> nodes) {
		double maxTime = -1.0;
		int minIndex = -1;
		for (int i = 0; i < nodes.size(); i++) {
			Node node = nodes.get(i);
			if (node.getHeight() > maxTime) {
				maxTime = node.getHeight();
				minIndex = i;
			}
		}
		if (minIndex == -1) {
			System.out.println(treeIntervalsInput.get().treeInput.get());
			System.out.println(nodes.size());
			for (Node n : nodes) {
				System.out.println(n.getHeight());
			}
            throw new RuntimeException("No max time found ");
		}
		return minIndex;
	}

    
    @Override
    protected boolean requiresRecalculation() {
        return ((CalculationNode) popSizeInput.get()).isDirtyCalculation() || super.requiresRecalculation();
    }
}
