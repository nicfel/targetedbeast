package targetedbeast;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.operator.Exchange;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.operator.SubtreeSlide;
import beast.base.evolution.speciation.YuleModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeStatLogger;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.DirectSimulator;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.distribution.Exponential;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;

import targetedbeast.alignment.ConsensusAlignment;
import targetedbeast.edgeweights.ParsimonyWeights2;
import targetedbeast.operators.BactrianIntervalScaleOperator;
import targetedbeast.operators.HeightBasedNodeRandomizer;
import targetedbeast.operators.IntervalRateStdevScaler;
import targetedbeast.operators.IntervalScaleOperator;
import targetedbeast.operators.RangeSlide;
import targetedbeast.operators.TargetedWilsonBalding;
import targetedbeast.operators.TargetedWilsonBaldingRates;
import targetedbeast.operators.WeightBasedNodeRandomizer;
import targetedbeast.operators.WeightedWideOperator;

/**
 * Full MCMC sampling UNDER THE PRIOR for EVERY targeted operator.
 *
 * Target = prior only (no likelihood):
 *   tree      ~ Yule(birthDiffRate)
 *   rate_i    ~ LogNormal(mean=1, S=ucldStdev)
 *   ucldStdev ~ Exponential(mean=1)
 *
 * A correct operator must leave that prior invariant. The reference is a standard-operator MCMC
 * (the control), which shares the harness's finite-sample behaviour; each targeted operator is run
 * as the dominant operator on top of the same supporting set and must recover the SAME ucldStdev and
 * Yule tree stats (height, length) as the control (a broken Hastings ratio deviates). A DirectSimulator
 * baseline (independent draws) is printed as the ground truth for reference.
 *
 * The alignment is realistic and actually exercised: taxa {0,1,2} form a zero-mutation cluster (taxon 2
 * carries N's, so it stays flat) which gives WeightBasedNodeRandomizer real candidates, while taxa 3,4
 * create several parsimony-informative splits (with N's AT informative sites) so the edge-weight
 * targeting of the other operators is non-trivial. After each run the operator's own accept/reject
 * counters are asserted non-zero, proving the MCMC actually invoked it and it made accepted moves.
 */
public class RatePriorMCMCTest {

    static final int NUM_TAXA = 8;
    static final double BIRTH_DIFF_RATE = 2.0;
    static final long CHAIN_LENGTH = 5_000_000L;
    static final int LOG_EVERY = 500;
    static final int BASELINE_DRAWS = 20_000;
    static final String[] STATS = {"ucldStdev", "ucldMean", "tree.height", "tree.treeLength", "tree.imbalance"};

    // Where the per-operator sample logs are written so the distributions can be plotted afterwards.
    static final String LOG_DIR = System.getProperty("ratePriorLogDir", "/tmp/rateprior_logs");
    // Every targeted operator kind, used by the combined "ALL" run.
    static final String[] ALL_KINDS = {
        "WBNR", "WBNR_rates", "HeightBased", "RangeSlide", "RangeSlideSqrt",
        "WideTargeted", "WideUntargeted", "WideSqrt", "TWB", "TWB_rates", "TWB_rates_sqrt", "TWB_rates_edgeLength",
        "IntervalRandom", "IntervalRates", "IntervalRatesStdevFixed",
        "BactrianInterval", "BactrianUpDown", "IntervalScaleAll", "IntervalUpDown"
    };

    static double[][] cachedControl = null;
    static double[] cachedBaseline = null;

    static class Model {
        Tree tree; RealParameter rates, stdev, mean; CompoundDistribution posterior; ParsimonyWeights2 ew;
    }

    TaxonSet taxonSet() {
        TaxonSet ts = new TaxonSet();
        for (int i = 0; i < NUM_TAXA; i++) ts.initByName("taxon", new Taxon(String.valueOf(i)));
        return ts;
    }

    /**
     * A realistic 8-taxon alignment, built like real sequence data so the targeted operators behave
     * as they do in production (no degenerate zero-weight states).
     *
     * The parsimony-informative structure is a set of splits that PAIRWISE-SEPARATE every taxon
     * (a 3-bit code: split A=bit2, B=bit1, C=bit0), so no two taxa have identical consensus
     * sequences and hence no pair has zero parsimony distance -- the collapse that made the tiny
     * 5-taxon toy produce totalWeight/totalDistance == 0 and crash the weighted node picks.
     * Split A is deep (wide), C is shallow (narrow) so cherries {0,1},{2,3},{4,5},{6,7} are tight,
     * low-weight edges for WeightBasedNodeRandomizer to pin. Splits D,E add nesting/homoplasy.
     * N's are placed AT informative columns of splits A and C (so the edge-weight calculation must
     * treat them as fully ambiguous), and a gap/singleton column exercises the rest of the range.
     */
    static Alignment richAlignment() {
        StringBuilder[] s = new StringBuilder[NUM_TAXA];
        for (int i = 0; i < NUM_TAXA; i++) s[i] = new StringBuilder();
        addBlock(s, new int[]{0, 1, 2, 3}, 8, 'A', 'C');   // A (bit2): deep split, wide
        addColumn(s, "AANACCNC");                          // A-split PI column, N in t2 (left) & t6 (right)
        addColumn(s, "AAAACNCC");                          // A-split PI column, N in t5 (right)
        addBlock(s, new int[]{0, 1, 4, 5}, 4, 'G', 'T');   // B (bit1)
        addBlock(s, new int[]{0, 2, 4, 6}, 3, 'A', 'G');   // C (bit0): shallow split, narrow
        addColumn(s, "ANAGAGAG");                          // C-split PI column, N in t1 (right)
        addBlock(s, new int[]{0, 1}, 3, 'C', 'A');         // D: tight cherry {0,1}
        addBlock(s, new int[]{6, 7}, 3, 'T', 'C');         // E: clade {6,7}
        addBlock(s, new int[]{}, 3, 'A', 'A');             // constant (filtered by ConsensusAlignment)
        addColumn(s, "ACGTAC-T");                          // gap + singletons (filtered)

        Alignment d = new Alignment();
        List<Object> init = new ArrayList<>();
        for (int i = 0; i < NUM_TAXA; i++) { init.add("sequence"); init.add(new Sequence(String.valueOf(i), s[i].toString())); }
        init.add("dataType"); init.add("nucleotide");
        d.initByName(init.toArray());
        return d;
    }

    /** Append {@code width} identical columns: taxa in {@code leftSet} get {@code left}, the rest {@code right}. */
    static void addBlock(StringBuilder[] s, int[] leftSet, int width, char left, char right) {
        boolean[] isLeft = new boolean[s.length];
        for (int t : leftSet) isLeft[t] = true;
        for (int w = 0; w < width; w++)
            for (int i = 0; i < s.length; i++) s[i].append(isLeft[i] ? left : right);
    }

    /** Append one explicit column given as an 8-char string (one base/N/gap per taxon, in taxon order). */
    static void addColumn(StringBuilder[] s, String column) {
        for (int i = 0; i < s.length; i++) s[i].append(column.charAt(i));
    }

    ParsimonyWeights2 parsimonyWeights(Tree tree, Alignment aln) {
        ConsensusAlignment ca = new ConsensusAlignment();
        ca.initByName("data", aln);
        ParsimonyWeights2 ew = new ParsimonyWeights2();
        ew.initByName("tree", tree, "data", ca);
        return ew;
    }

    Model buildModel() {
        Model m = new Model();
        m.tree = new Tree();
        m.tree.initByName("taxonset", taxonSet());
        m.tree.setID("tree");
        m.rates = new RealParameter();
        m.rates.initByName("value", "1.0", "dimension", 2 * NUM_TAXA - 1, "lower", "0.0");
        m.rates.setID("rate");
        m.stdev = new RealParameter();
        m.stdev.initByName("value", "0.5", "lower", "0.0");
        m.stdev.setID("ucldStdev");
        m.mean = new RealParameter();
        m.mean.initByName("value", "1.0", "lower", "0.0");
        m.mean.setID("ucldMean");
        m.ew = parsimonyWeights(m.tree, richAlignment());

        // rate prior: LogNormal with real-space mean = ucldMean (a sampled hyper-parameter, so the
        // scale-all up/down move that co-scales the mean rate has something to move and be tested on).
        LogNormalDistributionModel rateDistr = new LogNormalDistributionModel();
        rateDistr.initByName("S", m.stdev, "M", m.mean, "meanInRealSpace", true);
        Prior ratesPrior = new Prior();
        ratesPrior.initByName("x", m.rates, "distr", rateDistr);
        Exponential stdevDistr = new Exponential();
        stdevDistr.initByName("mean", "1.0");
        Prior stdevPrior = new Prior();
        stdevPrior.initByName("x", m.stdev, "distr", stdevDistr);
        Exponential meanDistr = new Exponential();
        meanDistr.initByName("mean", "1.0");
        Prior meanPrior = new Prior();
        meanPrior.initByName("x", m.mean, "distr", meanDistr);
        YuleModel treePrior = new YuleModel();
        treePrior.initByName("tree", m.tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);

        m.posterior = new CompoundDistribution();
        m.posterior.initByName("distribution", treePrior, "distribution", ratesPrior,
                "distribution", stdevPrior, "distribution", meanPrior, "distribution", m.ew);
        return m;
    }

    Operator makeOperator(String kind, Model m) {
        switch (kind) {
            case "STD": return null;
            case "WBNR": {
                WeightBasedNodeRandomizer o = new WeightBasedNodeRandomizer();
                o.initByName("tree", m.tree, "weight", 10.0, "percentage", 0.25, "edgeWeights", m.ew);
                return o;
            }
            case "WBNR_rates": {
                WeightBasedNodeRandomizer o = new WeightBasedNodeRandomizer();
                o.initByName("tree", m.tree, "weight", 10.0, "percentage", 0.25, "edgeWeights", m.ew, "rates", m.rates);
                return o;
            }
            case "HeightBased": {
                HeightBasedNodeRandomizer o = new HeightBasedNodeRandomizer();
                o.initByName("tree", m.tree, "weight", 10.0, "percentage", 0.25);
                return o;
            }
            case "RangeSlide": {
                RangeSlide o = new RangeSlide();
                o.initByName("tree", m.tree, "weight", 10.0, "size", 0.1, "edgeWeights", m.ew);
                return o;
            }
            case "WideTargeted": {
                WeightedWideOperator o = new WeightedWideOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "edgeWeights", m.ew);
                return o;
            }
            case "WideUntargeted": {
                WeightedWideOperator o = new WeightedWideOperator();
                o.initByName("tree", m.tree, "weight", 10.0);
                return o;
            }
            case "WideSqrt": {
                WeightedWideOperator o = new WeightedWideOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "edgeWeights", m.ew, "sqrtWeights", true);
                return o;
            }
            case "TWB": {
                TargetedWilsonBalding o = new TargetedWilsonBalding();
                o.initByName("tree", m.tree, "weight", 10.0, "edgeWeights", m.ew);
                return o;
            }
            case "TWB_rates": {
                TargetedWilsonBaldingRates o = new TargetedWilsonBaldingRates();
                o.initByName("tree", m.tree, "weight", 10.0, "edgeWeights", m.ew, "rates", m.rates);
                return o;
            }
            case "IntervalRandom": {
                IntervalScaleOperator o = new IntervalScaleOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "scaleFactor", 0.1, "scaleAllNodesIndependently", true);
                return o;
            }
            case "IntervalRates": {
                IntervalScaleOperator o = new IntervalScaleOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "scaleFactor", 0.25, "branchRates", m.rates);
                return o;
            }
            case "IntervalRatesStdevFixed": {
                IntervalRateStdevScaler o = new IntervalRateStdevScaler();
                o.initByName("tree", m.tree, "weight", 10.0, "scaleFactor", 0.25, "branchRates", m.rates, "stdev", m.stdev);
                return o;
            }
            case "BactrianInterval": {   // pure single scaler, no co-scaling (isolates the tree Jacobian)
                BactrianIntervalScaleOperator o = new BactrianIntervalScaleOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "scaleFactor", 0.1);
                return o;
            }
            case "BactrianUpDown": {     // scale-all co-scaling the mean rate (down)
                BactrianIntervalScaleOperator o = new BactrianIntervalScaleOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "scaleFactor", 0.1, "down", m.mean);
                return o;
            }
            case "IntervalScaleAll": {   // pure single scaler, no co-scaling (isolates the tree Jacobian)
                IntervalScaleOperator o = new IntervalScaleOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "scaleFactor", 0.9);
                return o;
            }
            case "RangeSlideSqrt": {
                RangeSlide o = new RangeSlide();
                o.initByName("tree", m.tree, "weight", 10.0, "size", 0.1, "edgeWeights", m.ew, "sqrtWeights", true);
                return o;
            }
            case "TWB_rates_sqrt": {
                TargetedWilsonBaldingRates o = new TargetedWilsonBaldingRates();
                o.initByName("tree", m.tree, "weight", 10.0, "edgeWeights", m.ew, "rates", m.rates, "sqrtWeights", true);
                return o;
            }
            case "TWB_rates_edgeLength": {
                TargetedWilsonBaldingRates o = new TargetedWilsonBaldingRates();
                o.initByName("tree", m.tree, "weight", 10.0, "edgeWeights", m.ew, "rates", m.rates, "useEdgeLength", true);
                return o;
            }
            case "IntervalUpDown": {     // scale-all co-scaling the mean rate (down)
                IntervalScaleOperator o = new IntervalScaleOperator();
                o.initByName("tree", m.tree, "weight", 10.0, "scaleFactor", 0.9, "down", m.mean);
                return o;
            }
        }
        throw new IllegalArgumentException("unknown operator kind " + kind);
    }

    double stat(Model m, String name) {
        switch (name) {
            case "ucldStdev": return m.stdev.getValue();
            case "ucldMean": return m.mean.getValue();
            case "tree.height": return m.tree.getRoot().getHeight();
            case "tree.treeLength": return TreeUtils.getTreeLength(m.tree, m.tree.getRoot());
            case "tree.imbalance": {
                beast.base.evolution.tree.Node r = m.tree.getRoot();
                return Math.abs(r.getLeft().getLeafNodeCount() - r.getRight().getLeafNodeCount());
            }
        }
        return Double.NaN;
    }

    /** Logs the root imbalance |#left leaves - #right leaves| so the MCMC prior test also covers tree shape. */
    public static class ImbalanceLogger extends beast.base.core.BEASTObject implements beast.base.core.Loggable {
        final Tree tree;
        public ImbalanceLogger(Tree tree) { this.tree = tree; }
        @Override public void initAndValidate() {}
        @Override public void init(java.io.PrintStream out) { out.print("tree.imbalance\t"); }
        @Override public void log(long sample, java.io.PrintStream out) {
            beast.base.evolution.tree.Node r = tree.getRoot();
            out.print(Math.abs(r.getLeft().getLeafNodeCount() - r.getRight().getLeafNodeCount()) + "\t");
        }
        @Override public void close(java.io.PrintStream out) {}
    }

    double[] baseline() throws Exception {
        if (cachedBaseline != null) return cachedBaseline;
        Model m = buildModel();
        DirectSimulator sim = new DirectSimulator();
        sim.initByName("distribution", m.posterior, "nSamples", 1);
        double[] sum = new double[STATS.length];
        new File(LOG_DIR).mkdirs();
        System.out.println("[RatePriorMCMCTest] per-operator sample logs (+ ratePrior_BASELINE.log) written to: " + LOG_DIR);
        java.io.PrintWriter bw = new java.io.PrintWriter(new File(LOG_DIR, "ratePrior_BASELINE.log"));
        bw.print("Sample"); for (String st : STATS) bw.print("\t" + st); bw.println();
        for (int k = 0; k < BASELINE_DRAWS; k++) {
            sim.run();
            bw.print(k);
            for (int s = 0; s < STATS.length; s++) { double v = stat(m, STATS[s]); sum[s] += v; bw.print("\t" + v); }
            bw.println();
        }
        bw.close();
        cachedBaseline = new double[STATS.length];
        for (int s = 0; s < STATS.length; s++) cachedBaseline[s] = sum[s] / BASELINE_DRAWS;
        return cachedBaseline;
    }

    /** Returns per-stat {mean, blockSE} recovered by an MCMC using the given operator kind. */
    double[][] runMCMC(String kind) throws Exception {
        Model m = buildModel();

        List<Operator> ops = new ArrayList<>();
        ScaleOperator stdevScaler = new ScaleOperator();
        stdevScaler.initByName("parameter", m.stdev, "scaleFactor", 0.6, "weight", 3.0);
        ScaleOperator ratesScaler = new ScaleOperator();
        ratesScaler.initByName("parameter", m.rates, "scaleFactor", 0.75, "weight", 3.0);
        ScaleOperator meanScaler = new ScaleOperator();
        meanScaler.initByName("parameter", m.mean, "scaleFactor", 0.75, "weight", 3.0);
        SubtreeSlide subtreeSlide = new SubtreeSlide();
        subtreeSlide.initByName("tree", m.tree, "weight", 5.0);
        Exchange narrow = new Exchange();
        narrow.initByName("tree", m.tree, "isNarrow", true, "weight", 3.0);
        ScaleOperator treeScaler = new ScaleOperator();
        treeScaler.initByName("tree", m.tree, "scaleFactor", 0.8, "weight", 1.0);
        ops.add(stdevScaler); ops.add(ratesScaler); ops.add(meanScaler); ops.add(subtreeSlide); ops.add(narrow); ops.add(treeScaler);

        Operator underTest = null;
        if (kind.equals("ALL")) {
            for (String k : ALL_KINDS) { Operator op = makeOperator(k, m); if (op != null) ops.add(op); }
        } else {
            underTest = makeOperator(kind, m);
            if (underTest != null) ops.add(underTest);
        }

        DirectSimulator sim = new DirectSimulator();
        sim.initByName("distribution", m.posterior, "nSamples", 1);
        sim.run();

        State state = new State();
        state.initByName("stateNode", m.tree, "stateNode", m.rates, "stateNode", m.stdev, "stateNode", m.mean);

        new File(LOG_DIR).mkdirs();
        File logFile = new File(LOG_DIR, "ratePrior_" + kind + ".log");
        TreeStatLogger treeStats = new TreeStatLogger();
        treeStats.initByName("tree", m.tree);
        Logger tracelog = new Logger();
        tracelog.initByName("fileName", logFile.getAbsolutePath(), "logEvery", LOG_EVERY,
                "log", m.stdev, "log", m.mean, "log", treeStats, "log", new ImbalanceLogger(m.tree));

        Logger.FILE_MODE = Logger.LogFileMode.overwrite;
        MCMC mcmc = new MCMC();
        List<Object> init = new ArrayList<>(List.of(
                "chainLength", CHAIN_LENGTH, "state", state, "distribution", m.posterior, "logger", tracelog));
        for (Operator op : ops) { init.add("operator"); init.add(op); }
        mcmc.initByName(init.toArray());
        mcmc.run();

        // Fire-check from the REAL MCMC: the operator must actually have been invoked and made
        // accepted moves. For the edge-weight operators this proves the alignment was exercised
        // (their proposals are driven entirely by the parsimony edge weights derived from it).
        if (underTest != null) {
            int acc = underTest.get_m_nNrAccepted(), rej = underTest.get_m_nNrRejected();
            System.out.printf("[%s] MCMC usage: %d accepted / %d rejected (%d proposals)%n",
                    kind, acc, rej, acc + rej);
            Assert.assertTrue(kind + " was never invoked by the MCMC", acc + rej > 0);
            Assert.assertTrue(kind + " never made an accepted move (alignment not productively exercised?)", acc > 0);
        }

        List<double[]> rows = new ArrayList<>();
        int[] cidx = new int[STATS.length];
        java.util.Arrays.fill(cidx, -1);
        boolean header = false;
        for (String line : java.nio.file.Files.readAllLines(logFile.toPath())) {
            if (line.startsWith("#") || line.trim().isEmpty()) continue;
            String[] f = line.split("\t");
            if (!header) {
                for (int i = 0; i < f.length; i++)
                    for (int s = 0; s < STATS.length; s++)
                        if (f[i].equals(STATS[s])) cidx[s] = i;
                header = true;
                continue;
            }
            double[] r = new double[STATS.length];
            boolean ok = true;
            for (int s = 0; s < STATS.length; s++) {
                try { r[s] = Double.parseDouble(f[cidx[s]]); } catch (Exception e) { ok = false; }
            }
            if (ok) rows.add(r);
        }
        int burn = rows.size() / 5;
        List<double[]> post = rows.subList(burn, rows.size());
        double[][] out = new double[STATS.length][2];
        int B = 20, blk = post.size() / B;
        for (int s = 0; s < STATS.length; s++) {
            double mean = 0; for (double[] r : post) mean += r[s]; mean /= post.size();
            double[] bm = new double[B];
            for (int b = 0; b < B; b++) { double t = 0; for (int k = 0; k < blk; k++) t += post.get(b * blk + k)[s]; bm[b] = t / blk; }
            double bMean = 0; for (double x : bm) bMean += x; bMean /= B;
            double bVar = 0; for (double x : bm) bVar += (x - bMean) * (x - bMean); bVar /= (B - 1);
            out[s][0] = mean; out[s][1] = Math.sqrt(bVar / B);
        }
        return out;
    }

    double[][] control() throws Exception {
        if (cachedControl == null) cachedControl = runMCMC("STD");
        return cachedControl;
    }

    void checkOperator(String kind) throws Exception {
        double[] base = baseline();
        double[][] ctrl = control();
        double[][] op = runMCMC(kind);
        System.out.printf("%n=== %s  (control-relative; baseline = DirectSimulator truth) ===%n", kind);
        boolean[] fail = new boolean[STATS.length];
        for (int s = 0; s < STATS.length; s++) {
            double tol = Math.max(4 * Math.sqrt(ctrl[s][1] * ctrl[s][1] + op[s][1] * op[s][1]), 0.03 * Math.abs(ctrl[s][0]));
            double diff = op[s][0] - ctrl[s][0];
            fail[s] = Math.abs(diff) > tol;
            System.out.printf("  %-14s baseline=%9.4f  control=%9.4f  %s=%9.4f +-%7.4f  diff=%+8.4f tol=%7.4f %s%n",
                    STATS[s], base[s], ctrl[s][0], kind, op[s][0], op[s][1], diff, tol, fail[s] ? "<== FAIL" : "ok");
        }
        for (int s = 0; s < STATS.length; s++) {
            double tol = Math.max(4 * Math.sqrt(ctrl[s][1] * ctrl[s][1] + op[s][1] * op[s][1]), 0.03 * Math.abs(ctrl[s][0]));
            Assert.assertEquals(kind + " must recover the control's " + STATS[s], ctrl[s][0], op[s][0], tol);
        }
    }

    @Test public void testWeightBasedNodeRandomizer()        throws Exception { checkOperator("WBNR"); }
    @Test public void testWeightBasedNodeRandomizerRates()   throws Exception { checkOperator("WBNR_rates"); }
    @Test public void testHeightBasedNodeRandomizer()        throws Exception { checkOperator("HeightBased"); }
    @Test public void testRangeSlide()                       throws Exception { checkOperator("RangeSlide"); }
    @Test public void testWeightedWideTargeted()             throws Exception { checkOperator("WideTargeted"); }
    @Test public void testWeightedWideUntargeted()           throws Exception { checkOperator("WideUntargeted"); }
    @Test public void testWeightedWideSqrt()                 throws Exception { checkOperator("WideSqrt"); }
    @Test public void testTargetedWilsonBalding()            throws Exception { checkOperator("TWB"); }
    @Test public void testTargetedWilsonBaldingRates()       throws Exception { checkOperator("TWB_rates"); }
    @Test public void testIntervalScaleRandom()              throws Exception { checkOperator("IntervalRandom"); }
    @Test public void testIntervalScaleWithRates()           throws Exception { checkOperator("IntervalRates"); }
    @Test public void testIntervalRateStdevScalerFixed()     throws Exception { checkOperator("IntervalRatesStdevFixed"); }
    @Test public void testBactrianIntervalScaleOperator()        throws Exception { checkOperator("BactrianInterval"); }
    @Test public void testBactrianUpDown()                       throws Exception { checkOperator("BactrianUpDown"); }
    @Test public void testIntervalScaleAll()                     throws Exception { checkOperator("IntervalScaleAll"); }
    @Test public void testRangeSlideSqrt()                       throws Exception { checkOperator("RangeSlideSqrt"); }
    @Test public void testRangeSlideWeighted()                   throws Exception { checkOperator("RangeSlideWeighted"); }
    @Test public void testTargetedWilsonBaldingRatesSqrt()       throws Exception { checkOperator("TWB_rates_sqrt"); }
    @Test public void testTargetedWilsonBaldingRatesEdgeLength() throws Exception { checkOperator("TWB_rates_edgeLength"); }
    @Test public void testIntervalScaleUpDown()                  throws Exception { checkOperator("IntervalUpDown"); }
    @Test public void testAllOperators()                         throws Exception { checkOperator("ALL"); }
}
