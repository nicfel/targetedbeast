package targetedbeast.util;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.distance.JukesCantorDistance;
import beast.base.inference.Runnable;
import beast.pkgmgmt.BEASTClassLoader;
import beastfx.app.inputeditor.AlignmentImporter;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.Utils;
import cern.colt.Arrays;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import targetedbeast.alignment.ConsensusAlignment;

@Description("Convert alignment to n-dimensional points where distances between points resemble pairwise sequence distances")
public class Alignment2PCA extends Runnable {
	final public Input<File> dataFileInput = new Input<>("datafile", "fasta or nexus file with alignment");
	final public Input<Alignment> dataInput = new Input<>("data", "alignment object that PCA will be applied to");
	final public Input<OutFile> matrixFileInput = new Input<>("matrixfile", "if datafile is specified, "
			+ "file where distance matrix is written (if matrixfile is specified). "
			+ "if datafile is not specified, file where distance matrix is read from");
	final public Input<Integer> dimensionInput = new Input<>("dimension", "dimension of output points", 2);
	final public Input<File> outputInput = new Input<>("out", "output file with results in comma delimited format -- if not specified results are on standard output");
	final public Input<Boolean> distanceBasedInput = new Input<>("distanceBased", "flag to indicate PCA is done based on distance matrix (if true) or normalised alignment (if false)", true);
	
	final static boolean debug = false;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		// get data
		Alignment data = null;
		int n = 0;
		List<String> taxa = null;
		double [][] distance;
		int [] sequenceID = null;
		if (dataFileInput.get() != null) {
			data = getAlignment(dataFileInput.get());
			n = data.getTaxonCount();
			taxa = data.getTaxaNames();

			// remove noninformative sites 
			ConsensusAlignment c = new ConsensusAlignment();
			c.initByName("data", data);
			data = c;
			
			List<Integer> usedSequenceIDs = new ArrayList<>();
			sequenceID = calcSequenceIDs(data, usedSequenceIDs);
			
			// create pairwise distance matrix
			if (distanceBasedInput.get()) {
				distance = calcDistanceMatrix(sequenceID, usedSequenceIDs, data, n);
			} else {
				distance = calcMatrix(sequenceID, usedSequenceIDs, data, n);
			}

	        if (matrixFileInput.get() != null) {
	        	writeDistances(matrixFileInput.get(), distance, taxa);
	        }
		} else {
			taxa = new ArrayList<>();
			distance = readDistances(matrixFileInput.get(), taxa);
			n = distance.length;
		}
		
		double [][] V_T = getPoints(distance);
		
		// resulting points to output
		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
			Log.info("Writing results in " + outputInput.get().getPath());
		}
		for (int i = 0; i < n; i++) {
			out.print(taxa.get(i) + ",");
			for (int j = 0; j < dimensionInput.get(); j++) {
				out.print(V_T[sequenceID[i]][j]);
				out.print(",");
			}
			out.println();
		}
		if (outputInput.get() != null) {
			out.close();
		}
		Log.warning("Done");
	}
	
	
	
	public static Map<String,double[]>  getPoints(Alignment data, int dimension, boolean distanceBased) {
		Alignment2PCA pca = new Alignment2PCA();
		pca.initByName("data", data, "dimension", dimension, "distanceBased", distanceBased);
		int n = data.getTaxonCount();
		List<String> taxa = data.getTaxaNames();

		// remove noninformative sites 
		ConsensusAlignment c = new ConsensusAlignment();
		c.initByName("data", data);
		data = c;
		
		List<Integer> usedSequenceIDs = new ArrayList<>();
		int [] sequenceID = pca.calcSequenceIDs(data, usedSequenceIDs);
		
		// create pairwise distance matrix
		double [][] distance;
		if (distanceBased) {
			distance = pca.calcDistanceMatrix(sequenceID, usedSequenceIDs, data, n);
		} else {
			distance = pca.calcMatrix(sequenceID, usedSequenceIDs, data, n);
		}
		double [][] V_T = pca.getPoints(distance);
		
		Map<String,double[]> map = new HashMap<>();
		
		for (int i = 0; i < n; i++) {
			String taxon = taxa.get(i);
			double [] trait = new double[dimension];
			for (int j = 0; j < dimension; j++) {
				trait[j] = V_T[sequenceID[i]][j];
			}
			map.put(taxon, trait);
		}
		return map;
	}
	
	
	
	// calculate binary matrix from site patterns: 
	// column[i][j*stateCount+k] = 1 means taxon i has at site j a state of k
	private double[][] calcMatrix(int[] sequenceID, List<Integer> usedSequenceIDs, Alignment data, int n) {
		int uniqueSequenceCount = sequenceID[sequenceID.length - 1] + 1;
		int stateCount = data.getMaxStateCount();

		int consensusSiteCount = 0;
		for (int j = 0; j < data.getSiteCount(); j++) {
			int pIndex = data.getPatternIndex(j);
			if (pIndex >= 0) {
				consensusSiteCount++;
			}
		}
		
		int [][] d = new int[uniqueSequenceCount][consensusSiteCount * stateCount];

		consensusSiteCount = 0;
		for (int j = 0; j < data.getSiteCount(); j++) {
			int pIndex = data.getPatternIndex(j);
			if (pIndex >= 0) {
				int [] pattern = data.getPattern(pIndex);
				for (int i = 0; i < uniqueSequenceCount; i++) {
					int state = pattern[usedSequenceIDs.get(i)];
					if (state < stateCount) {
						d[i][consensusSiteCount * stateCount + state] = 1;
//					} else {
//						for (int k = 0; k < stateCount; k++) {
//							d[i][consensusSiteCount * stateCount + k] = 1;
//						}
					}
				}
				consensusSiteCount++;
			}
		}
		
		// Get rid of all uninformative columns (all zero columns, all one columns, etc).		
		int [] length = new int[d[0].length];
		java.util.Arrays.fill(length, -1);
		int columnCount = 0;
		for (int i = 0; i < length.length; i++) {
			int sum = 0;
			for (int j = 0; j < uniqueSequenceCount; j++) {
				sum += d[j][i];
			}
			// only include column that have at least 1 ones and not all ones
			// not quite right, since unique sequences can be represented multiple times
			if (sum > 0 && sum < uniqueSequenceCount) {
				length[i] = sum;
				columnCount++;
			}
		}

		double [][] matrix = new double[uniqueSequenceCount][columnCount];
		int k = 0;
		for (int i = 0; i < length.length; i++) {
			if (length[i] > 0) {
				for (int j = 0; j < uniqueSequenceCount; j++) {
					matrix[j][k] = d[j][i];
				}
				k++;
			}
		}
		
		// Normalise length of vectors to 1
		for (int j = 0; j < uniqueSequenceCount; j++) {
			double sum = 0;
			double [] m = matrix[j];
			for (int i = 0; i < m.length; i++) {
				sum += m[i];
			}
			sum = Math.sqrt(sum);
			for (int i = 0; i < m.length; i++) {
				m[i] /= sum;
			}
		}

		Log.info("Matrix1 has dimension " + d.length + " x " + d[0].length);
		Log.info("Matrix2 has dimension " + matrix.length + " x " + matrix[0].length);
		return matrix;
	}

	
	private int[] calcSequenceIDs(Alignment data, List<Integer> usedSequenceIDs) {
		int [] seqIDs = new int[data.getTaxonCount()];

		Map<String,Integer> sequenceMap = new HashMap<>();
		int k = 0;
		List<String> taxaNames = data.getTaxaNames();
		for (int i = 0; i < data.getTaxonCount(); i++) {
			String taxon = taxaNames.get(i);
			String seq = data.getSequenceAsString(taxon); 
			if (!sequenceMap.containsKey(seq)) {
				usedSequenceIDs.add(i);
				sequenceMap.put(seq, k++);
			}
			seqIDs[i] = sequenceMap.get(seq);
		}
		return seqIDs;
	}

	private double[][] getPoints(double[][] distance) {
        Log.info("Starting SVD " + distance.length + "x" + distance[0].length);
		long start = System.currentTimeMillis();

		// Convert to RealMatrix
        RealMatrix dataMatrix = MatrixUtils.createRealMatrix(distance);

        // Compute the SVD (Singular Value Decomposition)
        org.apache.commons.math3.linear.SingularValueDecomposition svd = new org.apache.commons.math3.linear.SingularValueDecomposition(dataMatrix);

        // Get the first two principal components
        RealMatrix principalComponents = svd.getU().getSubMatrix(0, dataMatrix.getRowDimension() - 1, 0, dimensionInput.get()-1);

		long end = System.currentTimeMillis();
        Log.info("SVD finished in " + (end - start)/1000 + " seconds");
        
        return principalComponents.getData();
	}
	
	private double[][] getPointsX(double[][] distance) {
		// decompose into U x Lambda x V
		long start = System.currentTimeMillis();
		DoubleMatrix2D d = new DenseDoubleMatrix2D(distance);
		SingularValueDecomposition svd = new SingularValueDecomposition(d);
		DoubleMatrix2D V = svd.getV();
		long end = System.currentTimeMillis();
		Log.info("SVD finished in " + (end - start)/1000 + " seconds");
		
		if (debug) {
			DoubleMatrix2D U = svd.getU();
			double [] lambda = svd.getSingularValues();
			Log.info("d=");
			Log.info(d.toString());
			Log.info("U=");
			Log.info(U.toString());
			Log.info("lambda=" + Arrays.toString(lambda));
			Log.info("V=");
			Log.info(V.toString());
		}

		return V.toArray();
	}

	private double[][] readDistances(OutFile matrixFile, List<String> taxa) throws IOException {
        BufferedReader fin = new BufferedReader(new FileReader(matrixFile));
        String str = null;
        str = fin.readLine();
        int n = Integer.parseInt(str);
        double [][] d = new double[n][n];
        for (int i = 0; i < n; i++) {
        	str = fin.readLine();
        	String [] strs = str.split(",");
        	taxa.add(strs[0]);
        	for (int j = 1; j < i+1; j++) {
        		double val = Double.parseDouble(strs[j]);
        		d[i][j-1] = val;
        		d[j-1][i] = val;
        	}
        }
        fin.close();
		return d;
	}

	private void writeDistances(OutFile matrixFile, double [][] d, List<String> taxa) throws IOException {
		PrintStream out = new PrintStream(matrixFile);
		int n = taxa.size();
		out.println(n);
		for (int i = 0; i < n; i++) {
			out.print(taxa.get(i) + ",");
        	for (int j = 0; j < i; j++) {
        		out.print(d[i][j] + ",");
        	}
        	out.println();
		}
		out.close();
	}
	
	private double[][] calcDistanceMatrix(int[] sequenceID, List<Integer> usedSequenceIDs, Alignment data, int nx) {
		int uniqueSequenceCount = sequenceID[sequenceID.length - 1] + 1;
		Log.info("Creating "+uniqueSequenceCount+"x" + uniqueSequenceCount+ " distance matrix");
		long start = System.currentTimeMillis();


		Distance distance = new JukesCantorDistance();
    	((Distance.Base) distance).setPatterns(data);
        final double[][] distance0 = new double[uniqueSequenceCount][uniqueSequenceCount];
        for (int i0 = 0; i0 < uniqueSequenceCount; i0++) {
            distance0[i0][i0] = 0;
            int i = usedSequenceIDs.get(i0);
            for (int j0 = i0 + 1; j0 < uniqueSequenceCount; j0++) {
                int j = usedSequenceIDs.get(j0);
                distance0[i0][j0] = distance.pairwiseDistance(i, j);
                distance0[j0][i0] = distance0[i0][j0];
            }
        }
        
		long end = System.currentTimeMillis();
		Log.info("Distance matrix finished in " + (end - start)/1000 + " seconds");
		return distance0;
	}

	private static Alignment getAlignment(File file) {
        Set<String> importerClasses = Utils.loadService(AlignmentImporter.class);        
        for (String _class: importerClasses) {
        	try {
				AlignmentImporter importer = (AlignmentImporter) BEASTClassLoader.forName(_class).newInstance();
				if (importer.canHandleFile(file)) {
					List<BEASTInterface> os = importer.loadFile(file);
					Alignment a = (Alignment) os.get(0);
					return a;
				}
        	} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }

		throw new IllegalArgumentException("Could not find an importer that can handle the file " + file.getPath());
	}

	public static void main(String[] args) throws Exception {
		new Application(new Alignment2PCA(), "Alignment2PCA", args);
	}
}
