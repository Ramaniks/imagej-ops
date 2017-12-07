/*
 * #%L
 * ImageJ software for multidimensional image processing and analysis.
 * %%
 * Copyright (C) 2014 - 2017 Board of Regents of the University of
 * Wisconsin-Madison, University of Konstanz and Brian Northan.
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */

package net.imagej.ops.coloc.maxTKendallTau;

import net.imagej.ops.AbstractOpTest;
import net.imglib2.type.numeric.RealType;

import org.junit.Test;

/**
 * Tests {@link net.imagej.ops.Ops.Coloc.MaxTKendallTau}.
 *
 * @author Ellen T Arena
 * @param <V>
 */
public class MaxTKendallTauTest<V extends RealType<V>> extends AbstractOpTest {
	
	// Ranking data, test the rankTransformation() function. There are
	// two cases: 1) no tie breaking, test if this function
	// transforms 2.1 1.2 3.3 4.6 to 2 1 3 4 and 2) some tie breaking, test if
	// this function transforms 2.1 3 3 4.2 to 1 2 3 4 or 1 3 2 4
	@Test
	public void testRankTransformationNoTie() {
		double[][] values = new double[4][2];
		double[] values1 = {2.1, 1.2, 3.3, 4.6};
		double[] values2 = {2.1, 1.2, 3.3, 4.6};
		for (int i = 0; i < 4; i++) {
			values[i][0] = values1[i];
			System.out.println("values1 = " + values1[i] + " ");
			values[i][1] = values2[i];
			System.out.println("values2 = " + values2[i] + " ");
		}
		double[][] rank = MaxTKendallTau.rankTransformation(values, 0.0, 0.0, 4);
//		double[] valuesRank1 = new double[4];
//		double[] valuesRank2 = new double[4];
		for (int i = 0; i < 4; i++) {
//			valuesRank1[i] = rank[i][0];
			System.out.println("valuesRank1 = " + rank[i][0] + " ");
			System.out.println("valuesRank2 = " + rank[i][1] + " ");
//			valuesRank2[i] = rank[i][1];
		}
	}
	@Test
	public void testRankTransformationTie() {

	}
	// kendall tau, we need to test the function calculateKendallTau, we can use
	// some case like (1 2 3 4 and 4 3 2 1) tau=-1 or (1 2 3 4 and 1 2 3 4) tau=1
	@Test
	public void testCalculateKendallTau() {
		
		
	}

	// then we can test the whole class MTKT together. First, we can test one the
	// image which is 1 to 10 and 10 to 1. Second, we can generate some random
	// image and compare the result with results from R function. Third, we can
	// test the real image and compare the result with results from R function.
	@Test
	public void testMTKT() {
		
		
	}

	// Checks calculated pValue for MTKT.
	@Test
	public void testMTKTpValue() {
		
		
	}
}
