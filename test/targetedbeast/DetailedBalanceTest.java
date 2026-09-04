package targetedbeast;

import targetedbeast.alignment.ConsensusAlignment;
import targetedbeast.edgeweights.ConstantWeights;
import targetedbeast.edgeweights.FelsensteinWeights;
import targetedbeast.edgeweights.PCAWeights;
import targetedbeast.edgeweights.ParsimonyWeights2;
import targetedbeast.likelihood.SlowTreeLikelihood;
import targetedbeast.operators.BactrianIntervalScaleOperator;
import targetedbeast.operators.HeightBasedNodeRandomizer;
import targetedbeast.operators.IntervalScaleOperator;
import targetedbeast.operators.RangeSlide;
import targetedbeast.operators.TargetedWilsonBalding;
import targetedbeast.operators.WeightBasedNodeRandomizer;
import targetedbeast.operators.WeightedWideOperator;
//import targetedbeast.operators.TargetedWilsonBaldingFixedHeight;
import targetedbeast.util.Counter;
import targetedbeast.util.DefaultHashMap;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestName;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.speciation.YuleModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.DirectSimulator;
import beast.base.inference.Distribution;
import java.util.*;


public class DetailedBalanceTest {

    protected static final int NUM_SAMPLES = 200_000;
    protected static final int NUM_TAXA = 5;
    protected static final double BIRTH_DIFF_RATE = 2.0;
    protected static final double TOLERANCE = 0.001;
    protected static Alignment alignment;

    public interface TreeGroupMapper {
        public String op(Tree tree);
    }

    public static Map<String, TreeGroupMapper> treeGroupers;

    @Rule public TestName testName = new TestName();

    // Append one Markdown row per detailed-balance test to a shared validation table. Both this class
    // and CoScalerDetailedBalanceTest write the same file; the header is written only when the file is
    // empty/absent, so the runner just deletes the file before the run. Path: -DdetBalReport=<path>.
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

    // Append one CSV row per bucket pair (compared or not) so the render script can plot forward vs
    // backward flow for ALL buckets - showing which pairs the operator actually moves across (hit) and
    // which are vacuous (low count). Path: -DdetBalPairs=<path>.
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

    @BeforeClass
    public static void setUpClass() {
        alignment = createDummyAlignment();

        treeGroupers = new HashMap<>();

        // treeGroupers.put("RootHeight1Decimal", tree -> String.format("%.1f", (1 * tree.getRoot().getHeight())));
        treeGroupers.put("RootFirstChildHeight", tree -> {
            List<Node> rootChildren = tree.getRoot().getChildren();
            double firstChildHeight = Math.max(rootChildren.get(0).getHeight(),
                                               rootChildren.get(1).getHeight());
            return String.format("%.0f", 2 * firstChildHeight);
        });

        treeGroupers.put("TreeLength", tree -> {
            double treeLength = TreeUtils.getTreeLength(tree, tree.getRoot());
            return String.format("%.0f", 2 * treeLength);
        });

        treeGroupers.put("TreeImbalance", tree -> String.valueOf(rootImbalance(tree.getRoot())));
    }

    // reset before every test so a WBNR test's alignment override can't leak into other tests
    @Before
    public void resetAlignment() {
        alignment = createDummyAlignment();
    }


    @Test
    public void testWilsonBalding() throws Exception {
        TargetedWilsonBaldingConstFactory opFactory = new TargetedWilsonBaldingConstFactory();
        testDetailedBalance(opFactory);
    }

    @Test
    public void testTargWilsonBalding() throws Exception {
        TargetedWilsonBaldingFactory opFactory = new TargetedWilsonBaldingFactory();
        testDetailedBalance(opFactory);
    }

    @Test
    public void testTargWilsonBaldingFelsenstein() throws Exception {
        TargetedWilsonBaldingFelsensteinFactory opFactory = new TargetedWilsonBaldingFelsensteinFactory();
        testDetailedBalance(opFactory);
    }

    @Test
    public void testTargWilsonBaldingFixedHeight() throws Exception {
        TargetedWilsonBaldingFactory opFactory = new TargetedWilsonBaldingFactory();
        testDetailedBalance(opFactory);
    }

    @Test
    public void testTargWilsonBaldingEdgeWeightsPath() throws Exception {   // useEdgeLength=true -> EdgeWeightsTreeProposal
        testDetailedBalance(new TargetedWilsonBaldingEdgeWeightsFactory());
    }

    @Test
    public void testTargWilsonBaldingEdgeConst() throws Exception {   // EdgeWeightsTreeProposal with CONSTANT weights
        testDetailedBalance(new TargetedWilsonBaldingEdgeConstFactory());
    }

    @Test
    public void testTargWilsonBaldingEdgeUniformNode() throws Exception {   // uniform NODE pick + parsimony TARGET
        testDetailedBalance(new TargetedWilsonBaldingEdgeUniformNodeFactory());
    }

    @Test
    public void testTargWilsonBaldingEdgeFullRecompute() throws Exception {   // parsimony NODE, but full post-move recompute
        testDetailedBalance(new TargetedWilsonBaldingEdgeFullRecomputeFactory());
    }

    @Test
    public void testTargWilsonBaldingEdgeFreshWeights() throws Exception {   // brand-new edge-weights every iteration (no cache)
        TargetedWilsonBaldingEdgeWeightsFactory f = new TargetedWilsonBaldingEdgeWeightsFactory();
        f.freshEachIteration = true;
        testDetailedBalance(f);
    }

    @Test
    public void testTargWilsonBaldingEdgeUniformTarget() throws Exception {   // parsimony NODE + FLAT target distances
        testDetailedBalance(new TargetedWilsonBaldingEdgeUniformTargetFactory());
    }

    /** Parsimony weights but with a flat (constant) edge weight, so node selection is uniform while the
     *  regraft-target distances stay parsimony-based. Isolates the node-selection term from the target term. */
    static class UniformNodeParsimony extends ParsimonyWeights2 {
        @Override public double getEdgeWeights(int nodeNr) { return 1.0; }
    }

    /** Parsimony weights whose post-move recompute is a full (non-incremental) recompute. If the
     *  edge-weight node selection now balances, the bug is the incremental updateByOperator recompute. */
    static class FullRecomputeParsimony extends ParsimonyWeights2 {
        @Override public void updateByOperator() { forceRecompute(); }
    }

    /** Parsimony edge weights (so NODE selection stays parsimony) but FLAT regraft-target distances,
     *  to test whether the bug needs the distance calculation at all. */
    static class UniformTargetParsimony extends ParsimonyWeights2 {
        @Override public double[] getTargetWeights(int from, java.util.List<Node> to) {
            double[] d = new double[to.size()]; java.util.Arrays.fill(d, 1.0); return d;
        }
        @Override public double[] getTargetWeightsInteger(int from, java.util.List<Integer> to) {
            double[] d = new double[to.size()]; java.util.Arrays.fill(d, 1.0); return d;
        }
    }
    
    @Test
	public void testWeightedWideOperator() throws Exception {
		WeightedWideOperatorFactory opFactory = new WeightedWideOperatorFactory();
		testDetailedBalance(opFactory);
	}
    
    @Test
	public void testUntargetedWideOperator() throws Exception {
		UntargetedWideOperatorFactory opFactory = new UntargetedWideOperatorFactory();
		testDetailedBalance(opFactory);
	}

    
    @Test
	public void testRangeSlideOperator() throws Exception {
    	testRangeSlideOperatorFactory opFactory = new testRangeSlideOperatorFactory();
		testDetailedBalance(opFactory);
	}
    
    @Test
	public void testIntervalScaleOperator() throws Exception {
    	testIntervalScaleOperatorFactory opFactory = new testIntervalScaleOperatorFactory();
		testDetailedBalance(opFactory);
	}
    
    @Test
	public void testIntervalRandomScaleOperator() throws Exception {
    	testIntervalRandomScaleOperatorFactory opFactory = new testIntervalRandomScaleOperatorFactory();
		testDetailedBalance(opFactory);
	}

    @Test
	public void testBactrianIntervalScaleOperator() throws Exception {
		testDetailedBalance(new testBactrianIntervalScaleOperatorFactory());
	}

    @Test
	public void testBactrianIntervalRandomScaleOperator() throws Exception {
		testDetailedBalance(new testBactrianIntervalRandomScaleOperatorFactory());
	}

    @Test
	public void testRangeSlideSqrtOperator() throws Exception {
		testDetailedBalance(new testRangeSlideSqrtOperatorFactory());
	}

    @Test
		public void testRangeSlideBranchLengthOperator() throws Exception {
			testDetailedBalance(new testRangeSlideBranchLengthOperatorFactory());
		}

    @Test
		public void testRangeSlideUniformOperator() throws Exception {
			testDetailedBalance(new testRangeSlideUniformOperatorFactory());
		}

    @Test
	public void testWeightedWideSqrtOperator() throws Exception {
		testDetailedBalance(new testWeightedWideSqrtOperatorFactory());
	}

    
    @Test
	public void testHeightBasedRandomizer() throws Exception {
    	testHeightBasedRandomizerFactory opFactory = new testHeightBasedRandomizerFactory();
        testDetailedBalance(opFactory);    	        
    }
    
    @Test
    public void testWeightBasedNodeRandomizer() throws Exception {
    	alignment = createFlatClusterAlignment();
    	testWeightBasedNodeRandomizerFactory opFactory = new testWeightBasedNodeRandomizerFactory();
    	testDetailedBalance(opFactory);
    }

    @Test
    public void testWeightBasedNodeRandomizerWithNs() throws Exception {
    	alignment = createAlignmentWithNs();
    	testWeightBasedNodeRandomizerFactory opFactory = new testWeightBasedNodeRandomizerFactory();
    	testDetailedBalance(opFactory);
    }

    

    public void testDetailedBalance(OperatorFactory operatorFactory) throws Exception {
        // Create a map from mapper name to data structures
        Map<String, Counter<String>> groupCounters = new HashMap<>();
        Map<String, Counter<String>> proposalCounters = new HashMap<>();
        Map<String, DefaultHashMap<String, Double>> flows = new HashMap<>();

        // Initialize data structures for each mapper
        for (String mapperName : treeGroupers.keySet()) {
            groupCounters.put(mapperName, new Counter<>());
            proposalCounters.put(mapperName, new Counter<>());
            flows.put(mapperName, new DefaultHashMap<>(0.0));
        }

        Tree tree = new Tree();
        tree.initByName("taxonset", getTaxonSet(NUM_TAXA));

        Distribution prior = getPrior(tree);
        DirectSimulator simulator = new DirectSimulator();
        simulator.initByName("distribution", prior, "nSamples", 1);

        // Create operator once outside the loop to avoid memory issues
        TreeOperator operator = operatorFactory.getOperator(tree);

        for (int i = 0; i < NUM_SAMPLES; i++) {

            simulator.run();
            // DirectSimulator redraws the tree but does not refresh the (cached) edge weights, so make
            // them current before the operator reads them for its forward node selection.
            if (operatorFactory.ew != null) operatorFactory.ew.forceRecompute();
            // Decisive caching control: rebuild the operator with a brand-new (zero-history) edge-weights
            // object each iteration, so no cross-proposal cache state can survive.
            if (operatorFactory.freshEachIteration) operator = operatorFactory.getOperator(tree);

            // Compute all beforeKeys for all mappers
            Map<String, String> beforeKeys = new HashMap<>();
            for (String mapperName : treeGroupers.keySet()) {
                TreeGroupMapper mapper = treeGroupers.get(mapperName);
                String beforeKey = mapper.op(tree);
                beforeKeys.put(mapperName, beforeKey);
                groupCounters.get(mapperName).increment(beforeKey);
            }

            double logPBefore = prior.calculateLogP();

            double logHR = operator.proposal();

            double logPAfter = prior.calculateLogP();

            double pAccept = Math.min(1, Math.exp(logPAfter - logPBefore + logHR));

            // Compute all afterKeys for all mappers and update flow counters
            for (String mapperName : treeGroupers.keySet()) {
                TreeGroupMapper mapper = treeGroupers.get(mapperName);
                String afterKey = mapper.op(tree);
                String beforeKey = beforeKeys.get(mapperName);

                String transitionKey = beforeKey + "-" + afterKey;
                // String keyNoTransition = beforeKey + "-" + beforeKey;

                // Count number of proposals
                proposalCounters.get(mapperName).increment(transitionKey);

                // Update the expected flow
                flows.get(mapperName).put(transitionKey, flows.get(mapperName).get(transitionKey) + pAccept);
                // flows.get(mapperName).put(keyNoTransition, flows.get(mapperName).get(keyNoTransition) + (1 - pAccept));

            }
        }

        // Verify detailed balance; accumulate the worst pair + overall pass for the validation table.
        double worstRatio = -1, worstDiff = 0, worstTol = 0; String worstLabel = "-"; int totalPairs = 0; boolean pass = true;
        for (String mapperName : treeGroupers.keySet()) {
            System.out.println("\n=== Detailed balance test for mapper: " + mapperName + " ===");
            Counter<String> groupCounter = groupCounters.get(mapperName);
            Counter<String> proposalCounter = proposalCounters.get(mapperName);
            DefaultHashMap<String, Double> flow = flows.get(mapperName);

            System.out.println(groupCounter);

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

                    // A pair is only *compared* if both directions have enough counts; the rest are
                    // recorded too (compared=false) so the plot shows the vacuous buckets.
                    boolean compared = Math.min(groupCounter.getCount(group), groupCounter.getCount(toGroup)) >= 20
                                    && Math.min(nFwd, nBwd) >= 20;
                    double tolerance = Double.NaN;
                    if (compared) {
                        double forwVariance = estimateFlowVariance(groupCounter.getCount(group), nFwd, fwd);
                        double backVariance = estimateFlowVariance(groupCounter.getCount(toGroup), nBwd, bwd);
                        tolerance = 3 * Math.sqrt(forwVariance + backVariance);
                        double diff = Math.abs(fwd - bwd);
                        totalPairs++;
                        if (diff > tolerance) pass = false;
                        double ratio = tolerance > 0 ? diff / tolerance : (diff > 0 ? Double.POSITIVE_INFINITY : 0);
                        if (ratio > worstRatio) { worstRatio = ratio; worstDiff = diff; worstTol = tolerance; worstLabel = mapperName + ":" + keyForward; }
                        if (Math.max(fwd, bwd) > 0)
                            System.out.println(String.format("%-8s:  %8.4f <-> %8.4f    (tol=%8.4f     diff=%8.4f)", keyForward, fwd, bwd, tolerance, diff));
                    }
                    writePairRow("Tree", testName.getMethodName(), mapperName, keyForward, fwd, bwd, tolerance, nFwd, nBwd, compared);
                }
            }
        }

        // Record the row first (so failures still appear in the table), then assert.
        boolean vacuous = totalPairs == 0;
        writeReportRow("Tree", testName.getMethodName(), worstLabel, worstDiff, worstTol, totalPairs, pass && !vacuous);
        Assert.assertFalse("Vacuous detailed-balance test (no comparable cross-bucket transitions): " + testName.getMethodName(), vacuous);
        Assert.assertTrue("Detailed balance FAILED: worst " + worstLabel + " |Δflow|=" + worstDiff + " > tol=" + worstTol, pass);
    }

    double estimateFlowVariance(int stateCounts, int proposalCounts, double totalFlow) {
        double proposalFrequency = (double) proposalCounts / stateCounts;
        double meanAcceptRate = (totalFlow / proposalCounts);
        System.err.println(meanAcceptRate);
        return NUM_SAMPLES * meanAcceptRate * proposalFrequency * (1 - proposalFrequency);
    }

    protected TaxonSet getTaxonSet(int numTaxa) {
        TaxonSet taxonSet = new TaxonSet();
        for (int i = 0; i < numTaxa; i++) {
            taxonSet.initByName("taxon", new Taxon(String.valueOf(i)));
        }
        return taxonSet;
    }

    protected TreeDistribution getPrior(Tree tree) {
        YuleModel treePrior = new YuleModel();
        treePrior.initByName(
                "tree", tree,
                "birthDiffRate", "" + BIRTH_DIFF_RATE
        );
        return treePrior;
    }

    static public Alignment createDummyAlignment() {
        Sequence seq0 = new Sequence("0", "TGATAAAGAGTTACTAGAGTAAATAATAGGAGCTCCCCCTAGACTATG");
        Sequence seq1 = new Sequence("1", "CGATACAGAATTACTAGAGTAAATAATAGGAGTATCCCCCTGACTATA");
        Sequence seq2 = new Sequence("2", "TGATAAAGAAATACTAGAGTAAATAATAGGAGTTTCCCCTTGACTAAG");
        Sequence seq3 = new Sequence("3", "AGATATAGAGTTACTAGAGTAAATAATAGAGGTACCCGCTTGACAATG");
        Sequence seq4 = new Sequence("4", "TGACA-AGAGTTACTAGAGTAAAAAATAGAGGTCTCCCCTTCAGTATG");
        // Sequence seq5 = new Sequence("5", "CGACGAAGAGTTACTAGAGTAAATAACAGGGGTTTCCCCTTAACCATAGGAGTCGAACCCATCCTTGAGAATCCCTGCCACCCGTCGCACCCTGTTCTAAGTAAGGGGTTATACCCTTCCCATACTAAGAAATTTAGGTTAAACACAGACCAAGAGCC");

        Alignment data = new Alignment();
        data.initByName(
            "sequence", seq0,
            "sequence", seq1,
            "sequence", seq2,
            "sequence", seq3,
            "sequence", seq4,
            // "sequence", seq5,
            "dataType", "nucleotide"
        );
        return data;
    }

    /**
     * Flat-cluster alignment where WeightBasedNodeRandomizer CAN fire.
     * Taxa 0,1,2 are identical -> the edges inside that block are zero-mutation, so a candidate
     * chain grandparent->parent->j exists whenever the sampled tree groups them; the swaps stay
     * inside the block (uncle is another A-consensus node) and are parsimony-neutral.
     * Taxa 3,4 are identical-but-different, so every column has two states each in >=2 taxa
     * (parsimony-informative) -> ConsensusAlignment is non-empty and edge weights are defined.
     */
    static public Alignment createFlatClusterAlignment() {
        Sequence s0 = new Sequence("0", "AAAAAAAAAAAA");
        Sequence s1 = new Sequence("1", "AAAAAAAAAAAA");
        Sequence s2 = new Sequence("2", "AAAAAAAAAAAA");
        Sequence s3 = new Sequence("3", "CCCCCCCCCCCC");
        Sequence s4 = new Sequence("4", "CCCCCCCCCCCC");
        Alignment data = new Alignment();
        data.initByName("sequence", s0, "sequence", s1, "sequence", s2,
                        "sequence", s3, "sequence", s4, "dataType", "nucleotide");
        return data;
    }

    /**
     * Same idea but with ambiguous characters. Taxon 2 is all-N (fully ambiguous), so its edge
     * weight is 0 against any parent (diff = 0 when one state is fully ambiguous) -> it is always
     * "flat" and extends the candidate region wherever it attaches. The A/C split of {0,1} vs {3,4}
     * still supplies parsimony-informative sites. Exercises the operator on ambiguous consensus.
     */
    static public Alignment createAlignmentWithNs() {
        Sequence s0 = new Sequence("0", "AAAAAAAAAAAA");
        Sequence s1 = new Sequence("1", "AAAAAAAAAAAA");
        Sequence s2 = new Sequence("2", "NNNNNNNNNNNN");
        Sequence s3 = new Sequence("3", "CCCCCCCCCCCC");
        Sequence s4 = new Sequence("4", "CCCCCCCCCCCC");
        Alignment data = new Alignment();
        data.initByName("sequence", s0, "sequence", s1, "sequence", s2,
                        "sequence", s3, "sequence", s4, "dataType", "nucleotide");
        return data;
    }

    /*
     * Factories to pass to the test function to create various operators
     */

    abstract class OperatorFactory {

        ParsimonyWeights2 ew;   // set by factories that use parsimony weights, so the harness can refresh them
        boolean freshEachIteration = false;   // if true, rebuild the operator (+ a fresh, cache-less edge-weights) every iteration

        abstract public TreeOperator getOperator(Tree tree);

        SiteModel getSiteModel() {
            SiteModel siteModel = new SiteModel();
            siteModel.initByName("substModel", new JukesCantor());
            return siteModel;
        }

        ConstantWeights getConstantWeights(Tree tree) {
            ConstantWeights edgeWeights = new ConstantWeights();
            edgeWeights.initByName("tree", tree);
            return edgeWeights;
        }

        ParsimonyWeights2 getParsimonyWeights(Tree tree) {
            ConsensusAlignment consAlignment = new ConsensusAlignment();
            consAlignment.initByName("data", alignment);
            ParsimonyWeights2 edgeWeights = new ParsimonyWeights2();
            edgeWeights.initByName("tree", tree, "data", consAlignment);
            return edgeWeights;
        }

        PCAWeights getPCAWeights(Tree tree) {
            PCAWeights edgeWeights = new PCAWeights();
            edgeWeights.initByName("tree", tree, "data", alignment, "dimension", 3);
            return edgeWeights;
        }
        
        FelsensteinWeights getLikelihoodWeights(Tree tree) {
            SlowTreeLikelihood treeLikelihood = new SlowTreeLikelihood();
            treeLikelihood.initByName("tree", tree, "siteModel", getSiteModel(), "data", alignment, "implementation", "SlowBeerLikelihoodCore4");
            FelsensteinWeights edgeWeights = new FelsensteinWeights();
            edgeWeights.initByName("tree", tree, "likelihood", treeLikelihood, "data", alignment);
            return edgeWeights;
        }
        
    }

    class TargetedWilsonBaldingConstFactory extends OperatorFactory {

        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            // alignment = new SimulatedAlignment();
            // alignment.initByName("tree", tree, "siteModel", getSiteModel(), "sequenceLength", 5, "dataType", "nucleotide");
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getConstantWeights(tree));
            return operator;
        }
    }

    class TargetedWilsonBaldingFactory extends OperatorFactory {

        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            ParsimonyWeights2 ew = getParsimonyWeights(tree);
            this.ew = ew;   // length path selects its target by parsimony weight too — keep it current
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", ew);
            return operator;
        }
    }

    class TargetedWilsonBaldingEdgeWeightsFactory extends OperatorFactory {
        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            ParsimonyWeights2 ew = getParsimonyWeights(tree);
            this.ew = ew;   // let the harness keep the weights current after each DirectSimulator redraw
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", ew, "useEdgeLength", true);
            return operator;
        }
    }

    class TargetedWilsonBaldingEdgeConstFactory extends OperatorFactory {
        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getConstantWeights(tree), "useEdgeLength", true);
            return operator;
        }
    }

    class TargetedWilsonBaldingEdgeUniformNodeFactory extends OperatorFactory {
        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            ConsensusAlignment ca = new ConsensusAlignment();
            ca.initByName("data", alignment);
            UniformNodeParsimony ew = new UniformNodeParsimony();
            ew.initByName("tree", tree, "data", ca);
            this.ew = ew;   // keep the (parsimony) target weights current after each DirectSimulator redraw
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", ew, "useEdgeLength", true);
            return operator;
        }
    }

    class TargetedWilsonBaldingEdgeFullRecomputeFactory extends OperatorFactory {
        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            ConsensusAlignment ca = new ConsensusAlignment();
            ca.initByName("data", alignment);
            FullRecomputeParsimony ew = new FullRecomputeParsimony();
            ew.initByName("tree", tree, "data", ca);
            this.ew = ew;
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", ew, "useEdgeLength", true);
            return operator;
        }
    }

    class TargetedWilsonBaldingEdgeUniformTargetFactory extends OperatorFactory {
        { freshEachIteration = true; }
        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            ConsensusAlignment ca = new ConsensusAlignment();
            ca.initByName("data", alignment);
            UniformTargetParsimony ew = new UniformTargetParsimony();
            ew.initByName("tree", tree, "data", ca);
            this.ew = ew;
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", ew, "useEdgeLength", true);
            return operator;
        }
    }

    class TargetedWilsonBaldingFelsensteinFactory extends OperatorFactory {

        @Override
        public TargetedWilsonBalding getOperator(Tree tree) {
            TargetedWilsonBalding operator = new TargetedWilsonBalding();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getLikelihoodWeights(tree));
            return operator;
        }
    }


    class WeightedWideOperatorFactory extends OperatorFactory {

        @Override
        public WeightedWideOperator getOperator(Tree tree) {
        	WeightedWideOperator operator = new WeightedWideOperator();
            operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getPCAWeights(tree));
            return operator;
        }
    }
    
    class UntargetedWideOperatorFactory extends OperatorFactory {

        @Override
        public WeightedWideOperator getOperator(Tree tree) {
        	WeightedWideOperator operator = new WeightedWideOperator();
            operator.initByName("tree", tree, "weight", 1.0);
            return operator;
        }
    }

    
	class testRangeSlideOperatorFactory extends OperatorFactory {

		@Override
		public RangeSlide getOperator(Tree tree) {
			RangeSlide operator = new RangeSlide();
			operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getPCAWeights(tree));
			return operator;
		}
	}
	
	class testIntervalScaleOperatorFactory extends OperatorFactory {

		@Override
		public IntervalScaleOperator getOperator(Tree tree) {
			IntervalScaleOperator operator = new IntervalScaleOperator();
			operator.initByName("tree", tree, "weight", 1.0);
			return operator;
		}
	}
	
	class testIntervalRandomScaleOperatorFactory extends OperatorFactory {

		@Override
		public IntervalScaleOperator getOperator(Tree tree) {
			IntervalScaleOperator operator = new IntervalScaleOperator();
			operator.initByName("tree", tree, "weight", 1.0, "scaleAllNodesIndependently", true);
			return operator;
		}
	}

	
	class testHeightBasedRandomizerFactory extends OperatorFactory {

		@Override
		public HeightBasedNodeRandomizer getOperator(Tree tree) {
			HeightBasedNodeRandomizer operator = new HeightBasedNodeRandomizer();
			operator.initByName("tree", tree, "weight", 1.0);
			return operator;
		}
	}
	
	class testWeightBasedNodeRandomizerFactory extends OperatorFactory {

		@Override
		public WeightBasedNodeRandomizer getOperator(Tree tree) {
			WeightBasedNodeRandomizer operator = new WeightBasedNodeRandomizer();
			operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getParsimonyWeights(tree));
			return operator;
		}
	}

	class testBactrianIntervalScaleOperatorFactory extends OperatorFactory {
		@Override
		public BactrianIntervalScaleOperator getOperator(Tree tree) {
			BactrianIntervalScaleOperator operator = new BactrianIntervalScaleOperator();
			operator.initByName("tree", tree, "weight", 1.0, "scaleFactor", 0.1);
			return operator;
		}
	}

	class testBactrianIntervalRandomScaleOperatorFactory extends OperatorFactory {
		@Override
		public BactrianIntervalScaleOperator getOperator(Tree tree) {
			BactrianIntervalScaleOperator operator = new BactrianIntervalScaleOperator();
			operator.initByName("tree", tree, "weight", 1.0, "scaleFactor", 0.1, "scaleAllNodesIndependently", true);
			return operator;
		}
	}

	class testRangeSlideSqrtOperatorFactory extends OperatorFactory {
		@Override
		public RangeSlide getOperator(Tree tree) {
			RangeSlide operator = new RangeSlide();
			operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getPCAWeights(tree), "sqrtWeights", true);
			return operator;
		}
	}

	class testRangeSlideBranchLengthOperatorFactory extends OperatorFactory {
		@Override
		public RangeSlide getOperator(Tree tree) {
			RangeSlide operator = new RangeSlide();
			// node selection proportional to branch length (time), not edge weight
			operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getPCAWeights(tree), "weightByBranchLength", true);
			return operator;
		}
	}

	class testRangeSlideUniformOperatorFactory extends OperatorFactory {
		@Override
		public RangeSlide getOperator(Tree tree) {
			RangeSlide operator = new RangeSlide();
			// uniform node selection, but edge weights still supplied so the shared cache stays current
			operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getPCAWeights(tree), "uniform", true);
			return operator;
		}
	}

	class testWeightedWideSqrtOperatorFactory extends OperatorFactory {
		@Override
		public WeightedWideOperator getOperator(Tree tree) {
			WeightedWideOperator operator = new WeightedWideOperator();
			operator.initByName("tree", tree, "weight", 1.0, "edgeWeights", getPCAWeights(tree), "sqrtWeights", true);
			return operator;
		}
	}



    /**
     * Compute tree-level imbalance as weighted sum of node imbalances
     * @param tree The tree to compute imbalance for
     * @return Weighted sum of node imbalances
     */
    static public double treeImbalance(Tree tree) {
        double weightedSum = 0.0;
        double totalWeight = 0.0;
        
        // Traverse all internal nodes and compute weighted imbalance
        for (Node node : tree.getNodesAsArray()) {
            if (!node.isLeaf()) {
                ImbalanceResult result = nodeImbalance(node);
                if (!Double.isNaN(result.imbalance)) {
                    weightedSum += result.imbalance * result.weight;
                    totalWeight += result.weight;
                }
            }
        }
        
        // Return weighted average if there are any internal nodes
        return totalWeight > 0 ? weightedSum / totalWeight : 0.0;
    }

    static public int rootImbalance(Node root) {
        return Math.abs(root.getLeft().getLeafNodeCount() - root.getRight().getLeafNodeCount());
    }

    /**
     * Helper class to hold imbalance and weight values
     */
    static class ImbalanceResult {
        public final double imbalance;
        public final double weight;
        
        public ImbalanceResult(double imbalance, double weight) {
            this.imbalance = imbalance;
            this.weight = weight;
        }
    }
    
    /**
     * Compute node imbalance metric
     * @param node The node to compute imbalance for
     * @return ImbalanceResult containing imbalance and weight
     */
    static public ImbalanceResult nodeImbalance(Node node) {
        int size = node.getLeafNodeCount();
        double I;
        double w = 1.0;

        if (node.getChildCount() < 2 || size < 4) {
            I = Double.NaN;
        } else {
            Node c1 = node.getChild(0);
            Node c2 = node.getChild(1);

            int bigger = Math.max(c1.getLeafNodeCount(), c2.getLeafNodeCount());
            double m = Math.ceil(size / 2.0);

            I = (bigger - m) / (size - m - 1);
        }

        if (size % 2 == 1) {  // odd
            w = 1.0;
        } else {              // even
            w = 1.0 - 1.0 / size;
            if (I == 0) {
                w *= 2;
            } else {
                assert Double.isNaN(I) || I > 0 : "I should be NaN or positive, but was: " + I;
            }
        }
        if (Double.isNaN(I)) {
            I = 0.0;
            w = 0.0;
        }
        return new ImbalanceResult(I, w);
    }

}
