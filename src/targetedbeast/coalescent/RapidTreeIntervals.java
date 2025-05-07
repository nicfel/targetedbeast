package targetedbeast.coalescent;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.IntervalList;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.inference.CalculationNode;
import beast.base.util.HeapSort;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



/*
 * TreeIntervals.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

/**
 * Extracts the intervals from a beast.tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeIntervals.java,v 1.9 2005/05/24 20:25:56 rambaut Exp $
 */
@Description("Extracts the intervals from a tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class RapidTreeIntervals extends TreeIntervals {

    protected double[] samplingTimes;
    protected int[] samplingNodeNo;

    
    public RapidTreeIntervals() {
        super();
    }

    public RapidTreeIntervals(Tree tree) {
        init(tree);
    }

    @Override
    public void initAndValidate() {
        // this initialises data structures that store/restore might need
        calculateSamplingTimes();

        calculateIntervals();
        intervalsKnown = false;
        
    }

    private void calculateSamplingTimes() {
    	Tree tree = treeInput.get();
    	
        final int nodeCount = tree.getLeafNodeCount();

        double[] samplingTimes = new double[nodeCount];
        int[] samplingNodeNo = new int[nodeCount];

        collectSamplingTimes(tree, samplingTimes, samplingNodeNo);

        indices = new int[nodeCount];

        HeapSort.sort(samplingTimes, indices);
        
        this.samplingTimes = new double[nodeCount];
        this.samplingNodeNo = new int[nodeCount];
        
        
        this.samplingTimes = new double[nodeCount];
		for (int i = 0; i < nodeCount; i++) {
			this.samplingTimes[i] = samplingTimes[indices[i]];
			this.samplingNodeNo[i] = samplingNodeNo[indices[i]];
		}

	}

	/**
     * CalculationNode methods *
     */
    @Override
    protected boolean requiresRecalculation() {
        // we only get here if the tree is dirty, which is a StateNode
        // since the StateNode can only become dirty through an operation,
        // we need to recalculate tree intervals
        intervalsKnown = false;
        return true;
    }

    @Override
    protected void restore() {
        super.restore();
    }

    @Override
    protected void store() {
        super.store();
    }


    @Override
	public int getSampleCount() {
        // Assumes a binary tree!
        return treeInput.get().getInternalNodeCount();
    }

    /**
     * get number of intervals
     */
    @Override
	public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalCount;
    }

    /**
     * Gets an interval.
     */
    @Override
	public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < 0 || i >= intervalCount) throw new IllegalArgumentException();
        return intervals[i];
    }

    /**
     * Defensive implementation creates copy
     *
     * @return
     */
    @Override
    public double[] getIntervals(double[] inters) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (inters == null) inters = new double[intervals.length];
        System.arraycopy(intervals, 0, inters, 0, intervals.length);
        return inters;
    }

    @Override
    public double[] getCoalescentTimes(double[] coalescentTimes) {

        if (!intervalsKnown) {
            calculateIntervals();
        }

        if (coalescentTimes == null) coalescentTimes = new double[getSampleCount()];

        double time = 0;
        int coalescentIndex = 0;
        for (int i = 0; i < intervals.length; i++) {
            time += intervals[i];
            for (int j = 0; j < getCoalescentEvents(i); j++) {
                coalescentTimes[coalescentIndex] = time;
                coalescentIndex += 1;
            }
        }
        return coalescentTimes;
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     */
    @Override
	public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        return lineageCounts[i];
    }


    /**
     * Returns the number of coalescent events in an interval
     */
    @Override
	public int getCoalescentEvents(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        if (i < intervalCount - 1) {
            return lineageCounts[i] - lineageCounts[i + 1];
        } else {
            return lineageCounts[i] - 1;
        }
    }
    
    @Override
	public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        int numEvents = getCoalescentEvents(i);

        if (numEvents > 0) return IntervalType.COALESCENT;
        else if (numEvents < 0) return IntervalType.SAMPLE;
        else return IntervalType.NOTHING;
    }



    /**
     * get the total height of the genealogy represented by these
     * intervals.
     */
    @Override
	public double getTotalDuration() {

        if (!intervalsKnown) {
            calculateIntervals();
        }
        double height = 0.0;
        for (int j = 0; j < intervalCount; j++) {
            height += intervals[j];
        }
        return height;
    }

    /**
     * Checks whether this set of coalescent intervals is fully resolved
     * (i.e. whether is has exactly one coalescent event in each
     * subsequent interval)
     */
    @Override
	public boolean isBinaryCoalescent() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) > 0) {
                if (getCoalescentEvents(i) != 1) return false;
            }
        }

        return true;
    }

    /**
     * Checks whether this set of coalescent intervals coalescent only
     * (i.e. whether is has exactly one or more coalescent event in each
     * subsequent interval)
     */
    @Override
	public boolean isCoalescentOnly() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) < 1) return false;
        }

        return true;
    }

    /**
     * Recalculates all the intervals for the given beast.tree.
     */
    @Override
    protected void calculateIntervals() {
    	try {
    		calculateIntervalsRapid();
    	} catch (IndexOutOfBoundsException e) {
    		// caught unsafe operation in calculateIntervalsRapid()
    		// so try it the slow way
    		super.calculateIntervals();
    	}
    }
    
    private void calculateIntervalsRapid() {
        Tree tree = treeInput.get();

        int nodeCount = tree.getNodeCount();
        int internalNodeCount = tree.getInternalNodeCount();

        times = new double[internalNodeCount];
        int[] nodeNo = new int[internalNodeCount];

        collectCoalescentTimes(tree, times, nodeNo);

        indices = new int[internalNodeCount];

        HeapSort.sort(times, indices);

        if (intervals == null || intervals.length != nodeCount) {
            intervals = new double[nodeCount];
            lineageCounts = new int[nodeCount];
            lineagesAdded = new List[nodeCount];
            lineagesRemoved = new List[nodeCount];
//            lineages = new List[nodeCount];

            storedIntervals = new double[nodeCount];
            storedLineageCounts = new int[nodeCount];

        } else {
            for (List<Node> l : lineagesAdded) {
                if (l != null) {
                    l.clear();
                }
            }
            for (List<Node> l : lineagesRemoved) {
                if (l != null) {
                    l.clear();
                }
            }
        }

        // start is the time of the first tip
        double nextCoal = times[indices[0]];
        double nextSampling = samplingTimes[0];
        int numLines = 0;
        int coalNo = 0;
        int sampNo = 0;
        double currTime = 0.0;
        
        intervalCount = 0;
        while (intervalCount < nodeCount) {
        	// check if the next event is a coalescent event
        	if (nextCoal < nextSampling) {
        		Node node = tree.getNode(nodeNo[indices[coalNo]]);
                removeLineage(intervalCount, node.getLeft());
                removeLineage(intervalCount, node.getRight());
                addLineage(intervalCount, node);
                intervals[intervalCount] = nextCoal - currTime;
                lineageCounts[intervalCount] = numLines;
                currTime = nextCoal;
                numLines -= 1;
                coalNo++;
                if (coalNo < times.length)
                	nextCoal = times[indices[coalNo]];
        	}else{
                addLineage(intervalCount, tree.getNode(samplingNodeNo[sampNo]));
                intervals[intervalCount] = nextSampling - currTime;
                lineageCounts[intervalCount] = numLines;
                currTime = nextSampling;
                numLines += 1;
                sampNo++;
                nextSampling = sampNo < samplingTimes.length ? samplingTimes[sampNo] : Double.POSITIVE_INFINITY;
        	}
        	intervalCount++;
        }
        
        intervalsKnown = true;
    }

    /**
     * Returns the time of the start of an interval
     *
     * @param i which interval
     * @return start time
     */
    public double getIntervalTime(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return times[indices[i]];
    }
    
    protected void addLineage(int interval, Node node) {
        if (lineagesAdded[interval] == null) lineagesAdded[interval] = new ArrayList<>();
        lineagesAdded[interval].add(node);
    }

    protected void removeLineage(int interval, Node node) {
        if (lineagesRemoved[interval] == null) lineagesRemoved[interval] = new ArrayList<>();
        lineagesRemoved[interval].add(node);
    }

    /**
     * @return the delta parameter of Pybus et al (Node spread statistic)
     */
    public double getDelta() {

        return IntervalList.Utils.getDelta(this);
    }

    /**
     * extract coalescent times and tip information into array times from beast.tree.
     *
     * @param tree        the beast.tree
     * @param times       the times of the nodes in the beast.tree
     * @param childCounts the number of children of each node
     */
    protected static void collectCoalescentTimes(Tree tree, double[] times, int[] nodeNo) {
        Node[] nodes = tree.getNodesAsArray();
        int c = 0;
        for (int i = 0; i < nodes.length; i++) {
            Node node = nodes[i];
            if (!node.isLeaf()) {
            	times[c] = node.getHeight();
            	nodeNo[c] = node.getNr();
            	c++;
            }
        }
    }
    
    protected static void collectSamplingTimes(Tree tree, double[] times, int[] nodeNo) {
        Node[] nodes = tree.getNodesAsArray();
        int c = 0;
        for (int i = 0; i < nodes.length; i++) {
            Node node = nodes[i];
            if (node.isLeaf()) {
            	times[c] = node.getHeight();
            	nodeNo[c] = node.getNr();
            	c++;
            }
        }
    }

    
}