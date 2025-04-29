package targetedbeast.edgeweights;

import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;

@Description("Interface for weights to be used to use weighted tree operations")
public interface EdgeWeights {
	
	
//	void updateWeights ();

    void updateByOperator();

	void updateByOperatorWithoutNode(int ignore, List<Integer> nodes);
	
	// byte[] getNodeConsensus(int NodeNo);

	
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
	double[] getTargetWeights(int fromNodeNr, List<Node> toNodeNrs);
	
	/**
	 * Get the distance from one node to a list of nodes
	 * 
	 * @param fromNodeNr
	 * @param toNodeNrs
	 * @return
	 */
	double[] getTargetWeightsInteger(int fromNodeNr, List<Integer> toNodeNrs);

	double minEdgeWeight();

}
