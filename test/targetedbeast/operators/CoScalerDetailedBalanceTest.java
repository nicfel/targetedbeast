package targetedbeast.operators;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestName;

import beast.base.evolution.speciation.YuleModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeUtils;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.DirectSimulator;
import beast.base.inference.Distribution;
import beast.base.inference.distribution.Gamma;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;

import targetedbeast.util.Counter;
import targetedbeast.util.DefaultHashMap;
import targetedbeast.DetailedBalanceTest;
import targetedbeast.alignment.ConsensusAlignment;
import targetedbeast.edgeweights.ParsimonyWeights2;

/**
 * Detailed-balance tests for the co-scaling behaviour of {@link IntervalScaleOperator}.
 *
 * The methodology is the same as {@link targetedbeast.DetailedBalanceTest}: states are drawn iid
 * from a known target with a {@link DirectSimulator}, one operator step is applied, and the
 * accept-weighted flow between coarse-grained state buckets is checked for symmetry
 * (detailed balance: p_i * q_ij = p_j * q_ji). A wrong Hastings ratio breaks the symmetry.
 *
 * The only difference is the target. The co-scaling moves change real-valued parameters (up/down,
 * branchRates, stdev) in addition to the tree, and the co-scaling Jacobian is only balanced if those
 * parameters are part of the sampled target. Each test therefore builds a JOINT prior and buckets on
 * the co-scaled quantities, not just on the tree. Because those parameters cannot be reached through
 * {@code DetailedBalanceTest}'s tree-only {@code OperatorFactory}, the harness is inlined here rather
 * than reused.
 */
public class CoScalerDetailedBalanceTest {

    protected static final int NUM_SAMPLES = 200_000;
    protected static final int NUM_TAXA = 5;
    protected static final double BIRTH_DIFF_RATE = 2.0;

    /** Maps the current (mutated) state to a coarse-grained bucket label. */
    public interface StateGroupMapper {
        public String op();
    }

    @Rule public TestName testName = new TestName();

    // Append one Markdown row per test to the shared validation table (same file as DetailedBalanceTest;
    // header written only when the file is empty/absent). Path: -DdetBalReport=<path>.
    static synchronized void writeReportRow(String suite, String test, String worst,
                                            double diff, double tol, int pairs, boolean pass) {
        String path = System.getProperty("detBalReport", "validation/detailed_balance.md");
        try {
            java.io.File f = new java.io.File(path);
            if (f.getParentFile() != null) f.getParentFile().mkdirs();
            boolean header = !f.exists() || f.length() == 0;
            try (java.io.PrintWriter w = new java.io.PrintWriter(new java.io.OutputStreamWriter(
                    new java.io.FileOutputStream(f, true), java.nio.charset.StandardCharsets.UTF_8))) {
                if (header) {
                    w.println("| Suite | Test / operator | Worst bucket pair | worst Δflow | tol (3σ) | pairs | result |");
                    w.println("|---|---|---|---:|---:|---:|:---:|");
                }
                String result = pairs == 0 ? "⚠️ vacuous"
                        : (!Double.isFinite(diff) || !Double.isFinite(tol)) ? "⚠️ non-finite"
                        : (pass ? "✅ pass" : "❌ FAIL");
                w.printf("| %s | `%s` | %s | %.3f | %.3f | %d | %s |%n", suite, test, worst, diff, tol, pairs, result);
            }
        } catch (java.io.IOException e) {
            throw new RuntimeException("could not write detailed-balance report to " + path, e);
        }
    }

    // Append one CSV row per bucket pair (compared or not) for plotting all buckets. Path: -DdetBalPairs=<path>.
    static synchronized void writePairRow(String suite, String test, String mapper, String pair,
                                          double fwd, double bwd, double tol, int nFwd, int nBwd, boolean compared) {
        String path = System.getProperty("detBalPairs", "validation/detailed_balance_pairs.csv");
        try {
            java.io.File f = new java.io.File(path);
            if (f.getParentFile() != null) f.getParentFile().mkdirs();
            boolean header = !f.exists() || f.length() == 0;
            try (java.io.PrintWriter w = new java.io.PrintWriter(new java.io.OutputStreamWriter(
                    new java.io.FileOutputStream(f, true), java.nio.charset.StandardCharsets.UTF_8))) {
                if (header) w.println("suite,test,mapper,pair,forward,backward,tol,nFwd,nBwd,compared");
                w.printf("%s,%s,%s,%s,%.6f,%.6f,%.6f,%d,%d,%b%n",
                        suite, test, mapper, pair.replace(",", ";"), fwd, bwd, tol, nFwd, nBwd, compared);
            }
        } catch (java.io.IOException e) {
            throw new RuntimeException("could not write detailed-balance pairs to " + path, e);
        }
    }

    @Test
    public void testUpDownCoScaler() throws Exception {
        // Co-scaled parameter with its own prior: log(down) ~ N(0, 1).
        RealParameter down = new RealParameter();
        down.initByName("value", "1.0", "lower", 0.0);

        LogNormalDistributionModel downDist = new LogNormalDistributionModel();
        downDist.initByName("M", "0.0", "S", "1.0", "meanInRealSpace", false);
        Prior downPrior = new Prior();
        downPrior.initByName("x", down, "distr", downDist);

        // Joint target: Yule(tree) x LogNormal(down). DirectSimulator draws both each iteration.
        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));
        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);

        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", yule, "distribution", downPrior);

        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        IntervalScaleOperator operator = new IntervalScaleOperator();
        operator.initByName("tree", tree, "weight", 1.0, "down", down);

        // The tree is scaled up and down is scaled inversely; down is lognormal and scaled
        // multiplicatively, so bucket it on log(down) (additive shifts).
        Map<String, StateGroupMapper> stateGroupers = new LinkedHashMap<>();
        addTreeBuckets(stateGroupers, tree);
        stateGroupers.put("Down",       () -> String.format("%.1f", Math.log(down.getValue())));

        testDetailedBalance(prior, simulator, operator, stateGroupers);
    }

    @Test
    public void testBactrianUpDownCoScaler() throws Exception {
        // Same up/down co-scaling as testUpDownCoScaler but for the Bactrian-kernel interval scaler.
        RealParameter down = new RealParameter();
        down.initByName("value", "1.0", "lower", 0.0);
        LogNormalDistributionModel downDist = new LogNormalDistributionModel();
        downDist.initByName("M", "0.0", "S", "1.0", "meanInRealSpace", false);
        Prior downPrior = new Prior();
        downPrior.initByName("x", down, "distr", downDist);

        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));
        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", yule, "distribution", downPrior);
        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        BactrianIntervalScaleOperator operator = new BactrianIntervalScaleOperator();
        operator.initByName("tree", tree, "weight", 1.0, "scaleFactor", 0.1, "down", down);

        Map<String, StateGroupMapper> stateGroupers = new LinkedHashMap<>();
        addTreeBuckets(stateGroupers, tree);
        stateGroupers.put("Down",       () -> String.format("%.1f", Math.log(down.getValue())));

        testDetailedBalance(prior, simulator, operator, stateGroupers);
    }

    @Test
    public void testWeightBasedNodeRandomizerRatesCoScaler() throws Exception {
        // WeightBasedNodeRandomizer with rate co-scaling, against Yule x LogNormal(rates | fixed S).
        int dim = 2 * NUM_TAXA - 2;
        RealParameter rates = new RealParameter();
        rates.initByName("value", "1.0", "dimension", dim, "lower", 0.0);
        LogNormalDistributionModel ratesDist = new LogNormalDistributionModel();
        ratesDist.initByName("M", "0.0", "S", "0.5", "meanInRealSpace", false);
        Prior ratesPrior = new Prior();
        ratesPrior.initByName("x", rates, "distr", ratesDist);

        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));
        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", yule, "distribution", ratesPrior);
        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        // flat-cluster alignment so the randomizer has real candidates to fire on
        ConsensusAlignment ca = new ConsensusAlignment();
        ca.initByName("data", DetailedBalanceTest.createFlatClusterAlignment());
        ParsimonyWeights2 ew = new ParsimonyWeights2();
        ew.initByName("tree", tree, "data", ca);

        WeightBasedNodeRandomizer operator = new WeightBasedNodeRandomizer();
        operator.initByName("tree", tree, "weight", 1.0, "percentage", 0.5, "edgeWeights", ew, "rates", rates);

        // WBNR reshuffles topology (root imbalance) and co-scales the moved branches' rates (log-rate
        // spread); it barely changes total tree length, so bucket on what it actually moves.
        Map<String, StateGroupMapper> stateGroupers = new LinkedHashMap<>();
        addTreeBuckets(stateGroupers, tree);
        stateGroupers.put("RateLogSD", () -> String.format("%.1f", Math.log(sdLogRates(rates))));

        testDetailedBalance(prior, simulator, operator, stateGroupers, ew);
    }

    @Test
    public void testBranchRatesLengthCoScaler() throws Exception {
        // No-stdev branchRates path: each interval is scaled and the two child rates are co-scaled by
        // oldLength/newLength so that r_e * L_e stays constant (likelihood-invariant). Here it is tested
        // against the joint prior Yule(tree) x LogNormal(rates | fixed S).
        int dim = 2 * NUM_TAXA - 2;   // one rate per non-root node; assumes root is the highest node
        RealParameter rates = new RealParameter();
        rates.initByName("value", "1.0", "dimension", dim, "lower", 0.0);

        LogNormalDistributionModel ratesDist = new LogNormalDistributionModel();
        ratesDist.initByName("M", "0.0", "S", "0.5", "meanInRealSpace", false);   // fixed S
        Prior ratesPrior = new Prior();
        ratesPrior.initByName("x", rates, "distr", ratesDist);

        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));
        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);

        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", yule, "distribution", ratesPrior);

        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        IntervalRateCoScaler operator = new IntervalRateCoScaler();
        operator.initByName("tree", tree, "weight", 1.0, "branchRates", rates);   // no stdev

        // The per-branch co-scaling changes the rates heterogeneously and the interval scaling changes
        // the tree, so bucket on both the spread of the log-rates and the tree length.
        Map<String, StateGroupMapper> stateGroupers = new LinkedHashMap<>();
        stateGroupers.put("RateLogSD",  () -> String.format("%.1f", Math.log(sdLogRates(rates))));
        addTreeBuckets(stateGroupers, tree);

        testDetailedBalance(prior, simulator, operator, stateGroupers);
    }

    @Test
    public void testBranchRatesStdevCoScaler() throws Exception {
        // Variance move: the hyperparameter stdev has its own prior and log(rate) ~ N(0, stdev).
        RealParameter stdev = new RealParameter();
        stdev.initByName("value", "0.5", "lower", 0.0);
        Gamma stdevDist = new Gamma();
        stdevDist.initByName("alpha", "2.0", "beta", "0.5");   // shape/scale, mean = 1
        Prior stdevPrior = new Prior();
        stdevPrior.initByName("x", stdev, "distr", stdevDist);

        int dim = 2 * NUM_TAXA - 2;   // one rate per non-root node; assumes root is the highest node
        RealParameter rates = new RealParameter();
        rates.initByName("value", "1.0", "dimension", dim, "lower", 0.0);
        LogNormalDistributionModel ratesDist = new LogNormalDistributionModel();
        ratesDist.initByName("M", "0.0", "S", stdev, "meanInRealSpace", false);   // S references stdev
        Prior ratesPrior = new Prior();
        ratesPrior.initByName("x", rates, "distr", ratesDist);

        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));
        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);

        // Joint target: DirectSimulator draws stdev ~ Gamma, then rates ~ LogNormal(stdev), then tree.
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", yule, "distribution", ratesPrior, "distribution", stdevPrior);

        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        // IntervalRateCoScaler in stdev mode is the branchRates+stdev operator (the rotation move): it
        // co-scales node-height intervals and ROTATES the branch rates, keeping the rate prior (and hence
        // the realized log-rate SD and the stdev hyper-parameter) exactly fixed. So it moves the tree and
        // the individual rates, but not sd(log rates)/stdev - bucket on what it actually changes.
        IntervalRateCoScaler operator = new IntervalRateCoScaler();
        operator.initByName("tree", tree, "weight", 1.0, "scaleFactor", 0.25, "branchRates", rates, "stdev", stdev);

        Map<String, StateGroupMapper> stateGroupers = new LinkedHashMap<>();
        addTreeBuckets(stateGroupers, tree);
        stateGroupers.put("Rate0",      () -> String.format("%.1f", Math.log(rates.getValue(0))));

        testDetailedBalance(prior, simulator, operator, stateGroupers);
    }

    @Test
    public void testTargetedWilsonBaldingRatesLengthPath() throws Exception {
        wbRatesTest(false, false);   // useEdgeLength=false -> EdgeWeightsTreeProposal, raw weights
    }

    @Test
    public void testTargetedWilsonBaldingRatesEdgePath() throws Exception {
        wbRatesTest(true, false);    // useEdgeLength=true  -> LengthWeightedTreeProposal
    }

    @Test
    public void testTargetedWilsonBaldingRatesSqrt() throws Exception {
        wbRatesTest(false, true);    // sqrtWeights node selection (the production setting)
    }

    private void wbRatesTest(boolean useEdgeLength, boolean sqrtWeights) throws Exception {
        // joint target: Yule(tree) x LogNormal(rates | fixed S). The moved-stem rate co-scaling
        // Jacobian is only balanced if rates are part of the sampled target.
        int dim = 2 * NUM_TAXA - 2;   // one rate per non-root node (assumes root = highest-numbered node)
        RealParameter rates = new RealParameter();
        rates.initByName("value", "1.0", "dimension", dim, "lower", 0.0);
        LogNormalDistributionModel ratesDist = new LogNormalDistributionModel();
        ratesDist.initByName("M", "0.0", "S", "0.5", "meanInRealSpace", false);
        Prior ratesPrior = new Prior();
        ratesPrior.initByName("x", rates, "distr", ratesDist);

        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));
        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);

        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", yule, "distribution", ratesPrior);

        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        // edge weights the operator targets on (parsimony weights over a small fixed alignment)
        ConsensusAlignment ca = new ConsensusAlignment();
        ca.initByName("data", DetailedBalanceTest.createDummyAlignment());
        ParsimonyWeights2 ew = new ParsimonyWeights2();
        ew.initByName("tree", tree, "data", ca);

        TargetedWilsonBaldingRates operator = new TargetedWilsonBaldingRates();
        operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", ew, "rates", rates,
                "useEdgeLength", useEdgeLength, "sqrtWeights", sqrtWeights);

        // The move regrafts a subtree (tree changes) and co-scales the moved stem's rate (rate spread
        // changes), so bucket on both the tree length and the log-rate spread.
        Map<String, StateGroupMapper> stateGroupers = new LinkedHashMap<>();
        addTreeBuckets(stateGroupers, tree);
        stateGroupers.put("RateLogSD",  () -> String.format("%.1f", Math.log(sdLogRates(rates))));

        testDetailedBalance(prior, simulator, operator, stateGroupers, ew);
    }

    @Test
    public void testRangeSlideRatesEdge() throws Exception {
        rangeSlideRatesTest(false, false);   // edge-weight node selection + rate co-scaling
    }

    @Test
    public void testRangeSlideRatesLength() throws Exception {
        rangeSlideRatesTest(true, false);    // branch-length node selection + rate co-scaling
    }

    @Test
    public void testRangeSlideRatesSqrt() throws Exception {
        rangeSlideRatesTest(false, true);    // sqrt edge-weight node selection + rate co-scaling
    }

    private void rangeSlideRatesTest(boolean weightByBranchLength, boolean sqrtWeights) throws Exception {
        // joint target: Yule(tree) x LogNormal(rates | fixed S). RangeSlide co-scales every edge whose
        // length changes, so the co-scaling Jacobian only balances if rates are part of the target.
        int dim = 2 * NUM_TAXA - 2;
        RealParameter rates = new RealParameter();
        rates.initByName("value", "1.0", "dimension", dim, "lower", 0.0);
        LogNormalDistributionModel ratesDist = new LogNormalDistributionModel();
        ratesDist.initByName("M", "0.0", "S", "0.5", "meanInRealSpace", false);
        Prior ratesPrior = new Prior();
        ratesPrior.initByName("x", rates, "distr", ratesDist);

        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));
        YuleModel yule = new YuleModel();
        yule.initByName("tree", tree, "birthDiffRate", "" + BIRTH_DIFF_RATE);

        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", yule, "distribution", ratesPrior);

        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        ConsensusAlignment ca = new ConsensusAlignment();
        ca.initByName("data", DetailedBalanceTest.createDummyAlignment());
        ParsimonyWeights2 ew = new ParsimonyWeights2();
        ew.initByName("tree", tree, "data", ca);

        RangeSlide operator = new RangeSlide();
        operator.initByName("tree", tree, "weight", 1.0, "size", 0.5, "edgeWeights", ew, "rates", rates,
                "weightByBranchLength", weightByBranchLength, "sqrtWeights", sqrtWeights);

        Map<String, StateGroupMapper> stateGroupers = new LinkedHashMap<>();
        addTreeBuckets(stateGroupers, tree);
        stateGroupers.put("RateLogSD",  () -> String.format("%.1f", Math.log(sdLogRates(rates))));

        testDetailedBalance(prior, simulator, operator, stateGroupers, ew);
    }

    /**
     * Draw states iid from the target, apply one operator step, and verify the accept-weighted flow
     * between buckets is symmetric for each mapper. Also asserts the test is non-vacuous: at least one
     * bucket-pair with enough counts must actually be compared, otherwise the operator never moved a
     * state across a boundary and the check is meaningless.
     */
    public void testDetailedBalance(Distribution prior, DirectSimulator simulator, TreeOperator operator,
                                    Map<String, StateGroupMapper> stateGroupers) throws Exception {
        testDetailedBalance(prior, simulator, operator, stateGroupers, null);
    }

    public void testDetailedBalance(Distribution prior, DirectSimulator simulator, TreeOperator operator,
                                    Map<String, StateGroupMapper> stateGroupers, ParsimonyWeights2 ew) throws Exception {
        // Create a map from mapper name to data structures
        Map<String, Counter<String>> groupCounters = new HashMap<>();
        Map<String, Counter<String>> proposalCounters = new HashMap<>();
        Map<String, DefaultHashMap<String, Double>> flows = new HashMap<>();

        // Initialize data structures for each mapper
        for (String mapperName : stateGroupers.keySet()) {
            groupCounters.put(mapperName, new Counter<>());
            proposalCounters.put(mapperName, new Counter<>());
            flows.put(mapperName, new DefaultHashMap<>(0.0));
        }

        for (int i = 0; i < NUM_SAMPLES; i++) {

            simulator.run();
            // DirectSimulator redraws the state but leaves the cached edge weights stale; refresh them
            // so the operator's forward selection reads weights that match the current tree.
            if (ew != null) ew.forceRecompute();

            // Compute all beforeKeys for all mappers
            Map<String, String> beforeKeys = new HashMap<>();
            for (String mapperName : stateGroupers.keySet()) {
                String beforeKey = stateGroupers.get(mapperName).op();
                beforeKeys.put(mapperName, beforeKey);
                groupCounters.get(mapperName).increment(beforeKey);
            }

            double logPBefore = freshLogP(prior);

            double logHR = operator.proposal();

            double logPAfter = freshLogP(prior);

            double pAccept = Math.min(1, Math.exp(logPAfter - logPBefore + logHR));

            // Compute all afterKeys for all mappers and update flow counters
            for (String mapperName : stateGroupers.keySet()) {
                String afterKey = stateGroupers.get(mapperName).op();
                String beforeKey = beforeKeys.get(mapperName);

                String transitionKey = beforeKey + "-" + afterKey;

                proposalCounters.get(mapperName).increment(transitionKey);
                flows.get(mapperName).put(transitionKey, flows.get(mapperName).get(transitionKey) + pAccept);
            }
        }

        // Verify detailed balance; accumulate the worst pair + overall pass for the validation table.
        double worstRatio = -1, worstDiff = 0, worstTol = 0; String worstLabel = "-"; int totalPairs = 0;
        boolean pass = true;
        for (String mapperName : stateGroupers.keySet()) {
            System.out.println("\n=== Detailed balance test for mapper: " + mapperName + " ===");
            Counter<String> groupCounter = groupCounters.get(mapperName);
            Counter<String> proposalCounter = proposalCounters.get(mapperName);
            DefaultHashMap<String, Double> flow = flows.get(mapperName);

            System.out.println(groupCounter);

            int comparedPairs = 0;
            for (String group : groupCounter.keySet()) {
                for (String toGroup : groupCounter.keySet()) {
                    if (group.compareTo(toGroup) < 0) continue;   // check each pair once
                    if (group.equals(toGroup)) continue;          // skip self-loops

                    String keyForward = group + "-" + toGroup;
                    String keyBackward = toGroup + "-" + group;
                    double fwd = flow.get(keyForward);
                    double bwd = flow.get(keyBackward);
                    int nFwd = proposalCounter.getCount(keyForward);
                    int nBwd = proposalCounter.getCount(keyBackward);

                    boolean compared = Math.min(groupCounter.getCount(group), groupCounter.getCount(toGroup)) >= 20
                                    && Math.min(nFwd, nBwd) >= 20;
                    double tolerance = Double.NaN;
                    if (compared) {
                        double forwVariance = estimateFlowVariance(groupCounter.getCount(group), nFwd, fwd);
                        double backVariance = estimateFlowVariance(groupCounter.getCount(toGroup), nBwd, bwd);
                        tolerance = 3 * Math.sqrt(forwVariance + backVariance);
                        double diff = Math.abs(fwd - bwd);
                        comparedPairs++; totalPairs++;
                        if (diff > tolerance) pass = false;
                        double ratio = tolerance > 0 ? diff / tolerance : (diff > 0 ? Double.POSITIVE_INFINITY : 0);
                        if (ratio > worstRatio) { worstRatio = ratio; worstDiff = diff; worstTol = tolerance; worstLabel = mapperName + ":" + keyForward; }
                        if (Math.max(fwd, bwd) > 0)
                            System.out.println(String.format("%-8s:  %8.4f <-> %8.4f    (tol=%8.4f     diff=%8.4f)", keyForward, fwd, bwd, tolerance, diff));
                    }
                    writePairRow("Co-scale", testName.getMethodName(), mapperName, keyForward, fwd, bwd, tolerance, nFwd, nBwd, compared);
                }
            }
        }

        // Record the row first (so failures still appear in the table), then assert. A single vacuous
        // bucket (e.g. a scaler never changing topology) is expected and shows as "vacuous" in the matrix;
        // only fail if the operator crossed no bucket boundary at all.
        boolean vacuous = totalPairs == 0;
        writeReportRow("Co-scale", testName.getMethodName(), worstLabel, worstDiff, worstTol, totalPairs, pass && !vacuous);
        Assert.assertFalse("Vacuous detailed-balance test (no comparable cross-bucket transitions at all): " + testName.getMethodName(), vacuous);
        Assert.assertTrue("Detailed balance FAILED: worst " + worstLabel + " |Δflow|=" + worstDiff + " > tol=" + worstTol, pass);
    }

    /**
     * CompoundDistribution.calculateLogP() caches per child and only recomputes children flagged dirty;
     * without a State managing dirtiness the cache is stale (or NaN on the first call), which silently
     * broke the co-scale detailed-balance flows. Sum the children's freshly recomputed logP instead, so
     * the acceptance ratio sees the actual current state.
     */
    static double freshLogP(Distribution prior) {
        if (prior instanceof CompoundDistribution) {
            double s = 0;
            for (Distribution d : ((CompoundDistribution) prior).pDistributions.get()) s += d.calculateLogP();
            return s;
        }
        return prior.calculateLogP();
    }

    /** The three tree-shape buckets (length, imbalance, root-child height), matching DetailedBalanceTest,
     *  so every co-scale operator is also checked on the tree dimensions it moves. */
    static void addTreeBuckets(Map<String, StateGroupMapper> m, Tree tree) {
        m.put("TreeLength", () -> String.format("%.0f", 2 * TreeUtils.getTreeLength(tree, tree.getRoot())));
        m.put("TreeImbalance", () -> String.valueOf(DetailedBalanceTest.rootImbalance(tree.getRoot())));
        m.put("RootFirstChildHeight", () -> {
            java.util.List<beast.base.evolution.tree.Node> rc = tree.getRoot().getChildren();
            return String.format("%.0f", 2 * Math.max(rc.get(0).getHeight(), rc.get(1).getHeight()));
        });
    }

    double estimateFlowVariance(int stateCounts, int proposalCounts, double totalFlow) {
        double proposalFrequency = (double) proposalCounts / stateCounts;
        double meanAcceptRate = (totalFlow / proposalCounts);
        return NUM_SAMPLES * meanAcceptRate * proposalFrequency * (1 - proposalFrequency);
    }

    /** Population standard deviation of the log branch rates (the quantity the variance move scales). */
    static double sdLogRates(RealParameter rates) {
        int d = rates.getDimension();
        double m = 0;
        for (int i = 0; i < d; i++) m += Math.log(rates.getValue(i));
        m /= d;
        double v = 0;
        for (int i = 0; i < d; i++) {
            double x = Math.log(rates.getValue(i)) - m;
            v += x * x;
        }
        return Math.sqrt(v / d);
    }

    protected TaxonSet getTaxonSet(int numTaxa) {
        TaxonSet taxonSet = new TaxonSet();
        for (int i = 0; i < numTaxa; i++) {
            taxonSet.initByName("taxon", new Taxon(String.valueOf(i)));
        }
        return taxonSet;
    }
}
