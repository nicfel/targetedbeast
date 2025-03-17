package targetedbeast.init;

import java.util.List;

import beast.base.evolution.tree.ClusterTree;
import beast.base.evolution.tree.Node;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;

public class ensureNoZeroEdgeLength extends ClusterTree implements StateNodeInitialiser {

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		double rootHeight = getRoot().getHeight();
		
		double minEdgeLength = 1/(rootHeight*(getInternalNodeCount()/10.0));
		
//		System.out.println("minEdgeLength: " + minEdgeLength);
		
		resampleNodeHeight(getRoot(), minEdgeLength);
		
		// check if any edge length is below 0 or NaN
		for (int i = 0; i < getInternalNodeCount(); i++) {
			Node node = getNode(i);
			if (node.getLength() == Double.NaN || node.getLength()<0.0) {
				System.out.println("Edge length is below 0 or NaN");
				System.exit(0);
			}
		}
//		
//		
//		System.out.println(this);
//		System.exit(0);
		

	}
	
	private void resampleNodeHeight(Node node, double addEdgeLength) {
		if (node.isLeaf()) {
			return;
		}
		double oldHeights = node.getHeight() - Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		resampleNodeHeight(node.getLeft(), addEdgeLength);
		resampleNodeHeight(node.getRight(), addEdgeLength);
			
		double minHeight = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		double newHeight = oldHeights + addEdgeLength;
		node.setHeight(newHeight + minHeight);
		return;
	}

	@Override
    public void initStateNodes() {
        super.initStateNodes();
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
       super.getInitialisedStateNodes(stateNodes);
    }

}
