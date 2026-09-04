package targetedbeast.logger;

import java.io.PrintStream;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import targetedbeast.edgeweights.EdgeWeights;
import targetedbeast.edgeweights.ParsimonyWeights2;

/**
 * Tree logger that runs the full Fitch parsimony algorithm over the actual sequence bases and annotates
 * every branch with the number of substitutions (mutations) reconstructed on it: {@code [&mutations=N]}.
 *
 * Up-pass: each node's state set is the intersection of its children's sets, or their union if that
 * intersection is empty. Leaves are their observed base (A/C/G/T = a single state; N, gap and other
 * ambiguities = fully ambiguous {A,C,G,T}). Down-pass: the root is resolved to a base in its set, and
 * each node keeps its parent's base if that base is in the node's set, else takes a base from its set;
 * a mutation sits on a branch wherever the resolved base differs from the parent's. A fully-ambiguous
 * leaf (N) therefore resolves to its parent's base and carries no spurious mutation. Counts are weighted
 * by pattern weight (so a column shared by k sites contributes k). Colour/size branches by
 * {@code mutations} in FigTree/baltic to plot where substitutions fall.
 */
@Description("Full-Fitch tree logger: annotates each branch with [&mutations=0/1] (any parsimony substitution on it) and, if given, [&edgeweight=...] from the parsimony edge weights.")
public class FitchMutationTreeLogger extends CalculationNode implements Loggable {

    public Input<Tree> treeInput = new Input<>("tree", "tree to log", Input.Validate.REQUIRED);
    public Input<Alignment> dataInput = new Input<>("data",
            "alignment providing the actual bases for the Fitch reconstruction", Input.Validate.REQUIRED);
    public Input<EdgeWeights> edgeWeightsInput = new Input<>("edgeWeights",
            "optional parsimony edge weights, also logged per branch as [&edgeweight=...]");

    private Tree tree;
    private Alignment data;
    private EdgeWeights edgeWeights;
    private int patternCount;
    private int[] patternWeight;
    private int[][] leafBits;    // [leaf nr][pattern] = state set bitmask (single base, or 0xf for N/ambiguous)
    private int[][] up;          // [node nr][pattern] Fitch up-pass state set (internal nodes; leaves reuse leafBits)
    private int[][] fin;         // [node nr][pattern] resolved single-base bitmask
    private long[] mut;          // [node nr] substitutions on the branch above the node

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        data = dataInput.get();
        edgeWeights = edgeWeightsInput.get();
        patternCount = data.getPatternCount();
        patternWeight = new int[patternCount];
        int leafCount = tree.getLeafNodeCount();
        int[][] columns = new int[patternCount][];
        for (int p = 0; p < patternCount; p++) {
            columns[p] = data.getPattern(p);
            patternWeight[p] = data.getPatternWeight(p);
        }
        leafBits = new int[leafCount][patternCount];
        for (int leaf = 0; leaf < leafCount; leaf++) {
            int ti = data.getTaxonIndex(tree.getNode(leaf).getID());
            for (int p = 0; p < patternCount; p++) {
                int st = columns[p][ti];
                leafBits[leaf][p] = (st >= 0 && st < 4) ? (1 << st) : 0xf;   // definite base, else fully ambiguous
            }
        }
        int nodeCount = tree.getNodeCount();
        up = new int[nodeCount][];
        fin = new int[nodeCount][];
        for (int i = 0; i < nodeCount; i++) {
            fin[i] = new int[patternCount];
            if (i >= leafCount) up[i] = new int[patternCount];   // internal buffers; leaves point at leafBits
        }
        mut = new long[nodeCount];
    }

    @Override
    public void init(PrintStream out) {
        tree.init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (edgeWeights instanceof ParsimonyWeights2) {
            ((ParsimonyWeights2) edgeWeights).forceRecompute();   // make weights current for this tree
        }
        Node root = ((Tree) tree.getCurrent()).getRoot();
        upPass(root);
        downPass(root);
        out.print("tree STATE_" + sample + " = [&R] ");
        StringBuilder buf = new StringBuilder();
        toNewick(buf, root);
        buf.append(";");
        out.print(buf);
    }

    private void upPass(Node node) {
        if (node.isLeaf()) {
            up[node.getNr()] = leafBits[node.getNr()];
            return;
        }
        Node l = node.getLeft(), r = node.getRight();
        upPass(l);
        upPass(r);
        int[] sl = up[l.getNr()], sr = up[r.getNr()], out = up[node.getNr()];
        for (int p = 0; p < patternCount; p++) {
            int inter = sl[p] & sr[p];
            out[p] = inter != 0 ? inter : (sl[p] | sr[p]);
        }
    }

    private void downPass(Node node) {
        int nr = node.getNr();
        int[] s = up[nr], f = fin[nr];
        if (node.isRoot()) {
            for (int p = 0; p < patternCount; p++) {
                f[p] = Integer.lowestOneBit(s[p]);
            }
        } else {
            int[] pf = fin[node.getParent().getNr()];
            long m = 0;
            for (int p = 0; p < patternCount; p++) {
                int pb = pf[p];
                int fb = ((s[p] & pb) != 0) ? pb : Integer.lowestOneBit(s[p]);
                f[p] = fb;
                if (fb != pb) {
                    m += patternWeight[p];
                }
            }
            mut[nr] = m;
        }
        for (Node c : node.getChildren()) {
            downPass(c);
        }
    }

    private void toNewick(StringBuilder buf, Node node) {
        if (node.isLeaf()) {
            buf.append(node.getNr() + 1);
        } else {
            buf.append("(");
            List<Node> children = node.getChildren();
            for (int i = 0; i < children.size(); i++) {
                if (i > 0) buf.append(",");
                toNewick(buf, children.get(i));
            }
            buf.append(")");
        }
        if (!node.isRoot()) {
            buf.append("[&mutations=").append(mut[node.getNr()] > 0 ? 1 : 0);   // 0 = none, 1 = at least one
            if (edgeWeights != null) {
                buf.append(",edgeweight=").append(edgeWeights.getEdgeWeights(node.getNr()));
            }
            buf.append("]");
            buf.append(":").append(node.getLength());
        }
    }

    @Override
    public void close(PrintStream out) {
        tree.close(out);
    }
}
