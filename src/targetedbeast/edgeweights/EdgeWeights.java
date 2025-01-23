package targetedbeast.edgeweights;

import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;

@Description("Interface for weights to be used to use weighted tree operations")
public interface EdgeWeights {
	
	
	void updateMutations ();

    void updateByOperator();

	void updateByOperatorWithoutNode(int ignore, List<Integer> nodes);

	
	void fakeUpdateByOperator();

	void prestore();

	void reset();

	/**
	 * Get the edge weights for a node for operations
	 * 
	 * @param nodeNr
	 * @return
	 */
	double getEdgeWeights(int nodeNr);
	
	/**
	 * Get the distance from one node to a list of nodes
	 * @param fromNodeNr
	 * @param toNodeNrs
	 * @return
	 */
	public double[] getTargetWeights(int fromNodeNr , List<Node> toNodeNrs);
}
