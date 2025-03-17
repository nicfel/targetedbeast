package targetedbeast.logger;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

public class CladeLogger extends CalculationNode implements Loggable {
	
	public Input<Tree> treeInput = new Input<>("tree", "The tree to log clades for", Input.Validate.REQUIRED);

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
			out.print("# " + i + "\t" + treeInput.get().getNode(i).getID() + "\n");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		logClades(treeInput.get().getRoot(), out);
		
	}

	private List<Integer> logClades(Node node, PrintStream out) {
		List<Integer> clade = new ArrayList<>();
		if (node.isLeaf()) {
			clade.add(node.getNr());
			return clade;
		} else {
			for (Node child : node.getChildren()) {
				clade.addAll(logClades(child, out));
			}
		}
		// 
		Collections.sort(clade);
		String cls = "";
		for (int i =0; i < clade.size()-1; i++) {
            cls = cls + clade.get(i) + ",";
            }
		cls = cls + clade.get(clade.size()-1);
		out.print(cls + "\t");
		return clade;
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

}
