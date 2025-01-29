package targetedbeast.likelihood;

import java.util.Arrays;
import java.util.List;

/**
 * nucleotide implementation of standard likelihood core *
 */
public class RapidLikelihoodCore4 extends RapidLikelihoodCore {

	public RapidLikelihoodCore4() {
		super(4);
	}
	
	int[] matrixOffset;
	int[] catOffset;
	
	@Override
	public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories,
			boolean useAmbiguities) {
		super.initialize(nodeCount, patternCount, matrixCount, integrateCategories, useAmbiguities);
		matrixOffset = new int[nrOfMatrices];
		catOffset = new int[nrOfMatrices];
		for (int l = 0; l < nrOfMatrices; l++) {
			matrixOffset[l] = l * matrixSize;
			catOffset[l] = l * nrOfPatterns * nrOfStates;
		}
		
	}
	
	@Override
	protected void calculateStatesPruning(int[] stateIndex, double[] matrices, double[] intermediatePartials) {

		int v;
		
		for (int l = 0; l < nrOfMatrices; l++) {
			v = l * nrOfPatterns * nrOfStates;
			for (int k = 0; k < nrOfStates; k++) {			
				int state1 = stateIndex[k];
				int w = l * matrixSize;
				for (int i = 0; i < nrOfStates; i++) {
					intermediatePartials[v] = matrices[w + state1];	
					v++;
					w+= nrOfStates;
				}
			}
		}
	}

	@Override
	protected void calculatePartialPruning(double[] partials, double[] matrices, double[] intermediatePartials,
			int[] calcForPatterns) {
		
		double sum;
		int u,v,w;		
		double p1, p2, p3, p4;

		for (int l = 0; l < nrOfMatrices; l++) {
			u = catOffset[l];
			v = u;
			for (int k = 0; k < nrOfStates; k++) {
				w = matrixOffset[l]; // pick the correct matrix index
				p1 = partials[v];
				p2 = partials[v + 1];
				p3 = partials[v + 2];
				p4 = partials[v + 3];
				
				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				intermediatePartials[u] = sum;
				u++;
				
				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				
				intermediatePartials[u] = sum;
				u++;
				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				
				intermediatePartials[u] = sum;
				u++;
				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				intermediatePartials[u] = sum;
				u++;
				
				v += nrOfStates;
			}
		}	
		
		int k_offset;
		for (int k = 0; k < nrOfPatterns; k++) {
			if (calcForPatterns[k]==-1)
				break;
			
			k_offset = (calcForPatterns[k]) * nrOfStates;
			
			u = k_offset;
			v = u;
			w = 0;
			
			p1 = partials[v];
			p2 = partials[v + 1];
			p3 = partials[v + 2];
			p4 = partials[v + 3];

			sum = matrices[w] * p1;
			w++;
			sum += matrices[w] * p2;
			w++;
			sum += matrices[w] * p3;
			w++;
			sum += matrices[w] * p4;
			w++;
			intermediatePartials[u] = sum;
			u++;
			
			sum = matrices[w] * p1;
			w++;
			sum += matrices[w] * p2;
			w++;
			sum += matrices[w] * p3;
			w++;
			sum += matrices[w] * p4;
			w++;
			intermediatePartials[u] = sum;
			u++;
			
			sum = matrices[w] * p1;
			w++;
			sum += matrices[w] * p2;
			w++;
			sum += matrices[w] * p3;
			w++;
			sum += matrices[w] * p4;
			w++;
			intermediatePartials[u] = sum;
			u++;
			
			sum = matrices[w] * p1;
			w++;
			sum += matrices[w] * p2;
			w++;
			sum += matrices[w] * p3;
			w++;
			sum += matrices[w] * p4;
			w++;
			intermediatePartials[u] = sum;
			
					
			for (int l = 1; l < nrOfMatrices; l++) {
				
				u = catOffset[l] + k_offset;
				v = u;
				
				p1 = partials[v];
				p2 = partials[v + 1];
				p3 = partials[v + 2];
				p4 = partials[v + 3];

				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				intermediatePartials[u] = sum;
				u++;
				
				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				intermediatePartials[u] = sum;
				u++;
				
				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				intermediatePartials[u] = sum;
				u++;
				
				sum = matrices[w] * p1;
				w++;
				sum += matrices[w] * p2;
				w++;
				sum += matrices[w] * p3;
				w++;
				sum += matrices[w] * p4;
				w++;
				intermediatePartials[u] = sum;
			}
		}
		

	}

	@Override
	protected void updatePartials(double[] partials, 
			double[] intermediatePartials1, double[] intermediatePartials2,
			int[] mutations1, 
			int[] mutations2,
			int[] calcForPatterns) {
		
		updateConstantSitesPartials(partials, intermediatePartials1, intermediatePartials2, mutations1, mutations2, calcForPatterns);
		
		int v1, v2, u, pre_v1, pre_v2, pre_u, m1, m2, pattern;
		int m = 0;
		// calculate the product of the two intermediate matrices to get the partials at the parent node
		while (m < calcForPatterns.length && calcForPatterns[m] != -1) { 
			pattern = calcForPatterns[m];
			m1 = mutations1[pattern];
			m2 = mutations2[pattern];
			
			pre_v1 = m1 + m1 + m1 + m1;
			pre_v2 = m2 + m2 + m2 + m2;
			pre_u = pattern + pattern + pattern + pattern;
			
			if (m1 == -2) { // m1 and m2 == -2 should never show up as in this case, it shouldn't be treated as a pattern
				v2 = pre_v2;
				u = pre_u;

				partials[u] = intermediatePartials2[v2];
				u++;
				v2++;
				partials[u] = intermediatePartials2[v2];
				u++;
				v2++;
				partials[u] = intermediatePartials2[v2];
				u++;
				v2++;
				partials[u] = intermediatePartials2[v2];
	
				
				for (int l = 1; l < nrOfMatrices; l++) {
					v2 = pre_v2 + catOffset[l];
					u = pre_u + catOffset[l];
	
					partials[u] = intermediatePartials2[v2];
					u++;
					v2++;
					partials[u] = intermediatePartials2[v2];
					u++;
					v2++;
					partials[u] = intermediatePartials2[v2];
					u++;
					v2++;
					partials[u] = intermediatePartials2[v2];
				}
			}else if (m2 == -2) { 
				v1 = pre_v1;
				u = pre_u;
	
				partials[u] = intermediatePartials1[v1];
				u++;
				v1++;
				partials[u] = intermediatePartials1[v1];
				u++;
				v1++;
				partials[u] = intermediatePartials1[v1];
				u++;
				v1++;
				partials[u] = intermediatePartials1[v1];
	
				
				for (int l = 1; l < nrOfMatrices; l++) {
					v1 = pre_v1 + catOffset[l];
					u = pre_u + catOffset[l];
	
					partials[u] = intermediatePartials1[v1];
					u++;
					v1++;
					partials[u] = intermediatePartials1[v1];
					u++;
					v1++;
					partials[u] = intermediatePartials1[v1];
					u++;
					v1++;
					partials[u] = intermediatePartials1[v1];
				}
			}else {
				
					
					
					v1 = pre_v1;
					v2 = pre_v2;
					u = pre_u;
		
					partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
					u++;
					v1++;
					v2++;
					partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
					u++;
					v1++;
					v2++;
					partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
					u++;
					v1++;
					v2++;
					partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
		
					
					for (int l = 1; l < nrOfMatrices; l++) {
						v1 = pre_v1 + catOffset[l];
						v2 = pre_v2 + catOffset[l];
						u = pre_u + catOffset[l];
		
						partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
						u++;
						v1++;
						v2++;
						partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
						u++;
						v1++;
						v2++;
						partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
						u++;
						v1++;
						v2++;
						partials[u] = intermediatePartials1[v1] * intermediatePartials2[v2];
					}
				}
			
			m++;
		}
	}
	
	private void updateConstantSitesPartials(double[] partials, double[] intermediatePartials1,
			double[] intermediatePartials2, int[] mutations1, int[] mutations2, int[] calcForPatterns) {
		// calculate the product for the first nrOfStates
		for (int l = 0; l < nrOfMatrices; l++) {
			int u = catOffset[l];
			for (int k = 0; k < nrOfStates; k++) {
				partials[u] = intermediatePartials1[u] * intermediatePartials2[u];
				u++;
				partials[u] = intermediatePartials1[u] * intermediatePartials2[u];
				u++;
				partials[u] = intermediatePartials1[u] * intermediatePartials2[u];
				u++;
				partials[u] = intermediatePartials1[u] * intermediatePartials2[u];
				u++;
			}
		}	

	}
	
	
	

}
