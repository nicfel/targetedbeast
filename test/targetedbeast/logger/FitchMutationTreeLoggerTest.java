package targetedbeast.logger;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.TreeParser;

public class FitchMutationTreeLoggerTest {

    private String logTree(Alignment aln, String newick) {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", aln, "newick", newick, "IsLabelledNewick", true, "adjustTipHeights", false);
        FitchMutationTreeLogger logger = new FitchMutationTreeLogger();
        logger.initByName("tree", tree, "data", aln);
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(bos);
        logger.init(ps);
        logger.log(0L, ps);
        logger.close(ps);
        ps.flush();
        return bos.toString();
    }

    private int totalMutations(String out) {
        int sum = 0;
        Matcher m = Pattern.compile("\\[&mutations=(\\d+)\\]").matcher(out);
        while (m.find()) sum += Integer.parseInt(m.group(1));
        return sum;
    }

    @Test
    public void testSingleSplitIsOneMutation() {
        // ((0,1),(2,3)); site A A T T -> parsimony score 1 (one A<->T change on the internal split)
        Alignment aln = new Alignment();
        aln.initByName(
                "sequence", new Sequence("0", "A"),
                "sequence", new Sequence("1", "A"),
                "sequence", new Sequence("2", "T"),
                "sequence", new Sequence("3", "T"),
                "dataType", "nucleotide");
        String out = logTree(aln, "((0,1),(2,3))");
        assertEquals("total substitutions = parsimony score", 1, totalMutations(out));
        assertTrue("exactly one branch carries the mutation", out.contains("[&mutations=1]"));
    }

    @Test
    public void testNCreatesNoMutation() {
        // ((0,1),(2,3)); site A A T N -> still score 1; the N leaf (taxon 3 -> label 4) resolves to its
        // parent and carries no mutation.
        Alignment aln = new Alignment();
        aln.initByName(
                "sequence", new Sequence("0", "A"),
                "sequence", new Sequence("1", "A"),
                "sequence", new Sequence("2", "T"),
                "sequence", new Sequence("3", "N"),
                "dataType", "nucleotide");
        String out = logTree(aln, "((0,1),(2,3))");
        assertEquals("N is ignored: still one real substitution", 1, totalMutations(out));
        assertTrue("N leaf carries zero mutations\n" + out, out.contains("4[&mutations=0]"));
    }

    @Test
    public void testMultipleMutationsOnOneBranchIsStillOne() {
        // two sites both A A T T -> both changes fall on the SAME branch; binary must read 1, not 2
        Alignment aln = new Alignment();
        aln.initByName(
                "sequence", new Sequence("0", "AA"),
                "sequence", new Sequence("1", "AA"),
                "sequence", new Sequence("2", "TT"),
                "sequence", new Sequence("3", "TT"),
                "dataType", "nucleotide");
        String out = logTree(aln, "((0,1),(2,3))");
        assertEquals("binary: one branch carries mutations regardless of how many", 1, totalMutations(out));
    }
}
