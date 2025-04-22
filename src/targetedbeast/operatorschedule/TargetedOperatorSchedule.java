package targetedbeast.operatorschedule;


import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.operator.AdaptableOperatorSampler;
import beast.base.evolution.operator.EpochFlexOperator;
import beast.base.evolution.operator.Exchange;
import beast.base.evolution.operator.SubtreeSlide;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.operator.TreeStretchOperator;
import beast.base.evolution.operator.WilsonBalding;
import beast.base.evolution.operator.kernel.BactrianSubtreeSlide;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.operator.UpDownOperator;
import beast.base.inference.operator.kernel.BactrianUpDownOperator;
import targetedbeast.edgeweights.ConsensusWeights;
import targetedbeast.operators.HeightBasedNodeRandomizer;
import targetedbeast.operators.IntervalScaleOperator;
import targetedbeast.operators.RangeSlide;
import targetedbeast.operators.TargetedWilsonBalding;
import targetedbeast.operators.WeightBasedNodeRandomizer;
import targetedbeast.operators.WeightedWideOperator;

@Description("Operator schedule that replaces operators with Targeted operators")
@Citation("Bouckaert RR and Müller NF, Accelerated Bayesian phylodynamic inference through "
		+ "novel operators. In preparation, 2025")
public class TargetedOperatorSchedule extends OperatorSchedule {
	final public Input<Boolean> reweightInput = new Input<>("reweight", "if true, total weight of operators on trees remain the same", true);
	
	
	
	public TargetedOperatorSchedule() {
		super();
	}

	
	private ConsensusWeights weights = null;
	private Tree tree = null;
	private boolean intervalScaleOperatorAdded = false;
	private double weightFactor = 1.0;
	
	@Override
	public void initAndValidate() {
		weightFactor = 1.0;
		if (reweightInput.get()) {
			// get list of operators
			List<Operator> operators = null;
			for (BEASTInterface o : getOutputs()) {
				if (o instanceof Runnable) {
					operators = (List) o.getInput("operator");
				}
			}
			if (operators != null) {
				double weight = 0;
				for (Operator op : operators) {
					// tree operator?
					if (op instanceof TreeOperator) {
						weight += op.getWeight();
					}
				}
				weightFactor = weight / 75.0;
			}
		}

		super.initAndValidate();
	}
	
/*
NarrowExchange =>
	<operator id="WeightBasedNodeRandomizer" spec="targetedbeast.operators.WeightBasedNodeRandomizer" percentage="0.01" optimise="true" tree="@Tree.t:dna" weight="0.5" edgeWeights="@consensusWeights"/>
    <operator id="HeightBasedNodeRandomizer" spec="targetedbeast.operators.HeightBasedNodeRandomizer" percentage="0.01" optimise="true" tree="@Tree.t:dna" weight="0.5"/>

WideExchange =>
    <operator id="UnTargetedWide" spec="targetedbeast.operators.WeightedWideOperator" tree="@Tree.t:dna" weight="5.0"/>
    <operator id="TargetedWide" spec="targetedbeast.operators.WeightedWideOperator" tree="@Tree.t:dna" weight="15.0" edgeWeights="@consensusWeights"/>

SubTreeSlide =>    
    <operator id="CoalescentConstantSubtreeSlide.t:dna" spec="kernel.BactrianSubtreeSlide" tree="@Tree.t:dna" weight="20.0"/>
    <operator id="Range" spec="targetedbeast.operators.RangeSlide" tree="@Tree.t:dna" weight="30.0" size="0.1" edgeWeights="@consensusWeights"/>

    
WilsonBalding =>
    <operator id="TargetedWilsonBalding" spec="targetedbeast.operators.TargetedWilsonBalding" tree="@Tree.t:dna" weight="5" edgeWeights="@consensusWeights"/>
    <operator id="TargetedWilsonBalding2" spec="targetedbeast.operators.TargetedWilsonBalding" tree="@Tree.t:dna" weight="1" edgeWeights="@consensusWeights" useEdgeLength="true"/>

UpDownOperator
    <operator id="ScaleAll" spec="targetedbeast.operators.IntervalScaleOperator" tree="@Tree.t:dna" weight="0.25" scaleFactor="0.9" up="@popSize.t:dna">
            <down idref="clockRate.c:dna"/>
    </operator>
Flexer + Stretcher =>
    <operator id="ScaleIntervalsRandomly" spec="targetedbeast.operators.IntervalScaleOperator" tree="@Tree.t:dna" weight="0.25" scaleFactor="0.9" scaleAllNodesIndependently="true"/>

TreeRootScaler + UniformOperator =>
    ???    
*/
	
	
	@Override
	public void addOperator(Operator p) {
		if (p.getClass() == Exchange.class) {
			Exchange op = (Exchange) p;
			if (op.isNarrowInput.get()) {
				// Replace Narrow Exchange
				WeightBasedNodeRandomizer op1 = new WeightBasedNodeRandomizer();
				op1.initByName(op1.edgeWeightsInput.getName(), getConsensusWeights(),
						op1.percentageInput.getName(), 0.01,
						op1.treeInput.getName(), getTree(),
						op1.m_pWeight.getName(), 0.5 * weightFactor
						);
				op1.setID(op.getID() + "WeightBased");
				super.addOperator(op1);

				HeightBasedNodeRandomizer op2 = new HeightBasedNodeRandomizer();
				op2.initByName(
						op2.percentageInput.getName(), 0.01,
						op1.treeInput.getName(), getTree(),
						op2.m_pWeight.getName(), 0.5 * weightFactor
						);
				op1.setID(op.getID() + "HeightBased");
				super.addOperator(op2);
				
				Log.warning("replacing " + p.getID() + " with " + op1.getClass().getSimpleName() + " and " + op2.getClass().getSimpleName());
				
			} else {
				// Replace Wide Exchange
				WeightedWideOperator op1 = new WeightedWideOperator();
				op1.initByName(op1.edgeWeightsInput.getName(), getConsensusWeights(),
					op1.treeInput.getName(), getTree(),
					op1.m_pWeight.getName(), 0.25 * p.getWeight() * weightFactor
					);
				op1.setID(op.getID() + "UntargetedWide");
				super.addOperator(op1);
			
				WeightedWideOperator op2 = new WeightedWideOperator();
				op2.initByName(
					op2.edgeWeightsInput.getName(), getConsensusWeights(),
					op2.treeInput.getName(), getTree(),
					op2.m_pWeight.getName(), 0.75 * p.getWeight() * weightFactor
					);
				op2.setID(op.getID() + "TargetedWide");
				super.addOperator(op2);

				Log.warning("replacing " + p.getID() + " with " + op1.getClass().getSimpleName() + " and " + op2.getClass().getSimpleName());
			}
			
		} else if (p.getClass() == SubtreeSlide.class || p.getClass() == BactrianSubtreeSlide.class) {
			// SubTreeSlide => RangeSlide
			RangeSlide op1 = new RangeSlide();
			op1.initByName(
					op1.edgeWeightsInput.getName(), getConsensusWeights(),
					op1.treeInput.getName(), getTree(),
					op1.m_pWeight.getName(), 30.0 * weightFactor);
			op1.setID(p.getID() + "RangeSlide");
			super.addOperator(p);
			super.addOperator(op1);
			
			Log.warning("replacing " + p.getID() + " with " + p.getClass().getSimpleName() + " and " + op1.getClass().getSimpleName());

		} else if (p.getClass() == WilsonBalding.class) {
			// WilsonBalding => TargetedWilsonBalding
			TargetedWilsonBalding op1 = new TargetedWilsonBalding();
			op1.initByName(op1.edgeWeightsInput.getName(), getConsensusWeights(),
				op1.treeInput.getName(), getTree(),
				op1.m_pWeight.getName(), 5.0 * weightFactor
				);
			op1.setID(p.getID() + "Targeted");
			super.addOperator(op1);
		
			TargetedWilsonBalding op2 = new TargetedWilsonBalding();
			op2.initByName(
				op2.edgeWeightsInput.getName(), getConsensusWeights(),
				op2.treeInput.getName(), getTree(),
				op2.m_pWeight.getName(), 1.0 * weightFactor,
				op2.useEdgeLengthInput.getName(), true
				);
			op2.setID(p.getID() + "Targeted2");
			super.addOperator(op2);

			Log.warning("replacing " + p.getID() + " with " + op1.getClass().getSimpleName() + " and " + op2.getClass().getSimpleName() + "using edge lengths");

		} else if (p.getClass() == EpochFlexOperator.class || p.getClass() == TreeStretchOperator.class) {
			// EpochFlexOperator + TreeStretchOperator => IntervalScaleOperator
			if (!intervalScaleOperatorAdded) {
				IntervalScaleOperator op = new IntervalScaleOperator();
				op.initByName(
						op.edgeWeightsInput.getName(), getConsensusWeights(),
						op.treeInput.getName(), getTree(),
						op.m_pWeight.getName(), 0.25 * weightFactor
						);
				op.setID(p.getID() + "IntervalsScaler");
				intervalScaleOperatorAdded = true;
				super.addOperator(op);
				Log.warning("replacing " + p.getID() + " with " + op.getClass().getSimpleName());
			} else {
				Log.warning("removing " + p.getID() + " (replaced by IntervalScaleOperator)");
			}
		} else if (p.getClass() == UpDownOperator.class || p.getClass() == BactrianUpDownOperator.class) {
			// check if there is a tree in the UpDownOperator
			if (((List)p.getInput("down").get()).contains(getTree()) || 
					((List)p.getInput("up").get()).contains(getTree())) {
				IntervalScaleOperator op1 = new IntervalScaleOperator();
				
				List<StateNode> up = ((List<StateNode>)p.getInput("up").get());
				boolean treeGoesUp = up.contains(getTree()); 
				up.remove(getTree());
				List<StateNode> down = ((List<StateNode>)p.getInput("down").get());
				down.remove(getTree());
				op1.initByName(
						op1.edgeWeightsInput.getName(), getConsensusWeights(),
						op1.treeInput.getName(), getTree(),
						op1.downInput.getName(), treeGoesUp ? down : up,
						op1.upInput.getName(),   treeGoesUp ? up : down,
						op1.m_pWeight.getName(), 0.25 * weightFactor
						);
				op1.setID(p.getID());
				super.addOperator(op1);
				Log.warning("replacing " + p.getID() + " with " + op1.getClass().getSimpleName());
			} else {
				super.addOperator(p);
			}				
		} else if (p.getClass() == AdaptableOperatorSampler.class) {
			// may have an UpDownOperator hidden in there that needs changing
			AdaptableOperatorSampler op = (AdaptableOperatorSampler) p;
			IntervalScaleOperator op1 = null;
			Operator udop = null;
			for (Operator op2 : op.operatorsInput.get()) {
				if (op2 instanceof UpDownOperator || op2 instanceof BactrianUpDownOperator) {
					udop = op2;
					if (((List<StateNode>)udop.getInput("down").get()).contains(getTree()) || 
							((List<StateNode>)udop.getInput("up").get()).contains(getTree())) {
						// found one!
						op1 = new IntervalScaleOperator();
						List<StateNode> up = ((List<StateNode>)udop.getInput("up").get());
						boolean treeGoesUp = up.contains(getTree()); 
						up.remove(getTree());
						List<StateNode> down = ((List<StateNode>)udop.getInput("down").get());
						down.remove(getTree());
						op1.initByName(
								op1.edgeWeightsInput.getName(), getConsensusWeights(),
								op1.treeInput.getName(), getTree(),
								op1.downInput.getName(), treeGoesUp ? down : up,
								op1.upInput.getName(),   treeGoesUp ? up : down,
								op1.m_pWeight.getName(), 0.25 * weightFactor
								);
						op1.setID(op2
								.getID());
						
						Log.warning("replacing " + udop.getID() + " with " + op1.getClass().getSimpleName());
					}
				}
			}
			if (op1 != null) {
				op.operatorsInput.get().remove(udop);
				op.operatorsInput.get().add(op1);
				op.initAndValidate();
			}
			super.addOperator(p);
		} else {
			super.addOperator(p);
		}
	}

	/* Get ConsensusWeight object, if it already exists
	 * If it does not exist, create such object and add it to the posterior
	 * so it keeps up to date.
	 */
	private Object getConsensusWeights() {
		if (weights == null) {
			weights = new ConsensusWeights();
			
			tree = null;
			Alignment data = null;
			for (BEASTInterface o : getOutputs()) {
				if (o instanceof MCMC mcmc) {
					State state = mcmc.startStateInput.get();
					// find the tree
					for (StateNode sn : state.stateNodeInput.get()) {
						if (sn instanceof Tree) {
							if (tree != null) {
								throw new IllegalArgumentException("Expected only a single tree in the state");
							}
							tree = (Tree) sn;
							// find alignment
							for (BEASTInterface o2 : tree.getOutputs()) {
								if (o2 instanceof TreeLikelihood tl) {
									data = tl.dataInput.get();
								}
							}
						}
					}
				}
			}
			if (tree == null || data == null) {
				throw new IllegalArgumentException("Could not find tree or alignment");
			}
			weights.initByName("tree", tree, "data", data);

			for (BEASTInterface o : getOutputs()) {
				if (o instanceof MCMC mcmc) {
					Distribution d = mcmc.posteriorInput.get();
					if (d instanceof CompoundDistribution cp) {
						cp.pDistributions.get().add(weights);
					}
				}
			}
		}
		return weights;
	}
	
	/* Get tree, if it exists */
	private Tree getTree() {
		if (tree == null) {
			getConsensusWeights();
		}
		return tree;
	}

}
