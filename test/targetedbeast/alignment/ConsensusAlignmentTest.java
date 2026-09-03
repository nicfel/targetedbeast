package targetedbeast.alignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;

/**
 * Tests for {@link ConsensusAlignment}, which filters an alignment down to the parsimony-informative
 * sites. The rule (nucleotide data): count only A/C/G/T — anything else (N, gap, ambiguity codes) is
 * treated as missing — and keep a site iff at least two of those four states each occur in at least
 * two taxa. So a site needs "more than one mutation" that agree, not a lone singleton.
 */
public class ConsensusAlignmentTest {

    private static Alignment nucAlignment(String[][] taxonSeq) {
        Alignment aln = new Alignment();
        List<Object> args = new ArrayList<>();
        for (String[] ts : taxonSeq) {
            args.add("sequence");
            args.add(new Sequence(ts[0], ts[1]));
        }
        args.add("dataType");
        args.add("nucleotide");
        aln.initByName(args.toArray());
        return aln;
    }

    /**
     * 8 sites, one per rule case (columns read down the taxa t0..t3):
     *   site 0  A A G G  informative (A x2, G x2)                     -> KEPT
     *   site 1  A A A A  constant                                      -> dropped
     *   site 2  A A A G  singleton G                                   -> dropped
     *   site 3  A A N G  G singleton once N is ignored                 -> dropped
     *   site 4  A A - G  G singleton once gap is ignored               -> dropped
     *   site 5  C C T T  informative (C x2, T x2)                      -> KEPT
     *   site 6  A C G T  four singletons                              -> dropped
     *   site 7  A A G G  duplicate of site 0 (same pattern)            -> KEPT
     */
    @Test
    public void testParsimonyInformativeFiltering() {
        Alignment aln = nucAlignment(new String[][]{
                {"t0", "AAAAACAA"},
                {"t1", "AAAAACCA"},
                {"t2", "GAAN-TGG"},
                {"t3", "GAGGGTTG"},
        });
        ConsensusAlignment ca = new ConsensusAlignment();
        ca.initByName("data", aln);

        // sites 0, 5, 7 survive; 0 and 7 share a pattern, so 2 distinct patterns
        assertEquals("only two distinct informative patterns should survive", 2, ca.getPatternCount());

        int[] kept = {0, 5, 7};
        int[] dropped = {1, 2, 3, 4, 6};
        for (int s : kept)    assertTrue("site " + s + " should be kept",    ca.getPatternIndex(s) >= 0);
        for (int s : dropped) assertFalse("site " + s + " should be dropped", ca.getPatternIndex(s) >= 0);

        // pattern weights: AAGG twice (sites 0 & 7), CCTT once (site 5)
        int[] weights = new int[ca.getPatternCount()];
        for (int p = 0; p < ca.getPatternCount(); p++) weights[p] = ca.getPatternWeight(p);
        Arrays.sort(weights);
        assertArrayEqualsMsg(new int[]{1, 2}, weights);

        // reconstructed sequence spans only the retained sites, in original order (0,5,7)
        assertEquals("ACA", ca.getSequenceAsString("t0"));   // A(site0) C(site5) A(site7)
        assertEquals("GTG", ca.getSequenceAsString("t2"));   // G       T        G
    }

    /**
     * Non-ACGT symbols (N, gap) are missing, not states:
     *   site 0  A A N N G G  A x2 and G x2 (N ignored) -> KEPT (informative despite the Ns)
     *   site 1  A A A A N G  G is a singleton once N is ignored -> dropped
     */
    @Test
    public void testMissingBasesIgnored() {
        Alignment aln = nucAlignment(new String[][]{
                {"t0", "AA"},
                {"t1", "AA"},
                {"t2", "NA"},
                {"t3", "NA"},
                {"t4", "GN"},
                {"t5", "GG"},
        });
        ConsensusAlignment ca = new ConsensusAlignment();
        ca.initByName("data", aln);

        assertEquals("only the A/G 2-2 site survives; the N-masked singleton does not", 1, ca.getPatternCount());
        assertTrue("site 0 (AANNGG) informative", ca.getPatternIndex(0) >= 0);
        assertFalse("site 1 (AAAANG) singleton", ca.getPatternIndex(1) >= 0);
    }

    private static void assertArrayEqualsMsg(int[] expected, int[] actual) {
        assertEquals("pattern-weight multiset size", expected.length, actual.length);
        for (int i = 0; i < expected.length; i++)
            assertEquals("pattern weight " + i, expected[i], actual[i]);
    }
}
