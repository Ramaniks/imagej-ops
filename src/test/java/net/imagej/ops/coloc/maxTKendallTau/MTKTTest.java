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

import static org.junit.Assert.assertEquals;

import net.imagej.ops.Ops;
import net.imagej.ops.Ops.Coloc.MaxTKendallTau;
import net.imagej.ops.coloc.ColocalisationTest;
import net.imagej.ops.coloc.ShuffledView;
import net.imagej.ops.coloc.pValue.PValueResult;
import net.imagej.ops.special.function.BinaryFunctionOp;
import net.imagej.ops.special.function.Functions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;

import org.junit.Test;

/**
 * Tests {@link net.imagej.ops.Ops.Coloc.MaxTKendallTau}.
 *
 * @author Ellen T Arena
 * @author Shulei Wang
 * @param <V>
 */
public class MTKTTest<T extends RealType<T>, U extends RealType<U>> extends ColocalisationTest {
	
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
			values[i][1] = values2[i];
		}
		Img<DoubleType> vImage1 = ArrayImgs.doubles(values1, values1.length);
		Img<DoubleType> vImage2 = ArrayImgs.doubles(values2, values2.length);
		long seed = 0x89302341;
		double[][] rank = MTKT.rankTransformation(vImage1, vImage2, 0.0, 0.0, 4, seed);
		double[] expectedRankOrder = {1, 0, 2, 3};
		for (int i = 0; i < 4; i++) {
			assertEquals(expectedRankOrder[i], rank[i][0], 0.0);
			assertEquals(expectedRankOrder[i], rank[i][1], 0.0);
		}
	}
	@Test
	public void testRankTransformationTie() {
		double[][] values = new double[4][2];
		double[] values1 = {2.1, 3.0, 3.0, 4.2};
		double[] values2 = {2.1, 3.0, 3.0, 4.2};
		for (int i = 0; i < 4; i++) {
			values[i][0] = values1[i];
			values[i][1] = values2[i];
		}
		Img<DoubleType> vImage1 = ArrayImgs.doubles(values1, values1.length);
		Img<DoubleType> vImage2 = ArrayImgs.doubles(values2, values2.length);
		long seed = 0x89302341;
		double[][] rank = MTKT.rankTransformation(vImage1, vImage2, 0.0, 0.0, 4, seed);
		double[] expectedRankOrder1 = {0, 1, 2, 3};
		double[] expectedRankOrder2 = {0, 2, 1, 3};
		for (int i = 0; i < 4; i++) {
			// first element
			assertEquals(expectedRankOrder1[0], rank[0][0], 0.0); 
			assertEquals(expectedRankOrder1[0], rank[0][1], 0.0);
			// second element
			if (rank[1][0] == 1.0) {
				assertEquals(expectedRankOrder1[1], rank[1][0], 0.0);
			} else if (rank[1][0] == 2.0) {
				assertEquals(expectedRankOrder2[1], rank[1][0], 0.0);
			}
			if (rank[1][1] == 1.0) {
				assertEquals(expectedRankOrder1[1], rank[1][1], 0.0);
			} else if (rank[1][1] == 2.0) {
				assertEquals(expectedRankOrder2[1], rank[1][1], 0.0);
			}
			// third element
			if (rank[2][0] == 2.0) {
				assertEquals(expectedRankOrder1[2], rank[2][0], 0.0);
			} else if (rank[2][0] == 1.0) {
				assertEquals(expectedRankOrder2[2], rank[2][0], 0.0);
			}
			if (rank[2][1] == 2.0) {
				assertEquals(expectedRankOrder1[2], rank[2][1], 0.0);
			} else if (rank[2][1] == 1.0) {
				assertEquals(expectedRankOrder2[2], rank[2][1], 0.0);
			}
			// fourth element
			assertEquals(expectedRankOrder1[3], rank[3][0], 0.0); 
			assertEquals(expectedRankOrder1[3], rank[3][1], 0.0);
		}
	}

	// then we can test the whole class MTKT together. First, we can test one the
	// image which is 1 to 10 and 10 to 1. Second, we can generate some random
	// image and compare the result with results from R function. Third, we can
	// test the real image and compare the result with results from R function.
	@Test
	public void testMTKTZeroCorr() {
		double result = (Double) ops.run(MTKT.class, zeroCorrelationImageCh1, zeroCorrelationImageCh2);
		assertEquals(2.7109853865740496, result, 0.0);
	}
	@Test
	public void testMTKTPosCorr() {
		double result = (Double) ops.run(MTKT.class, positiveCorrelationImageCh1, positiveCorrelationImageCh2);
		assertEquals(57.04572197953461, result, 0.0);
	}

	// Checks calculated pValue for MTKT.
	@Test
	public void testMTKTpValue() {
		final double mean = 0.2;
		final double spread = 0.1;
		final double[] sigma = new double[] { 3.0, 3.0 };
		Img<FloatType> ch1 = ColocalisationTest.produceMeanBasedNoiseImage(new FloatType(), 24, 24,
			mean, spread, sigma, 0x01234567);
		Img<FloatType> ch2 = ColocalisationTest.produceMeanBasedNoiseImage(new FloatType(), 24, 24,
			mean, spread, sigma, 0x98765432);
		BinaryFunctionOp<RandomAccessibleInterval<FloatType>, RandomAccessibleInterval<FloatType>, Double> op =
			Functions.binary(ops, MaxTKendallTau.class, Double.class, ch1, ch2);
		PValueResult value = (PValueResult) ops.run(Ops.Coloc.PValue.class, new PValueResult(), ch1, ch2, op);
		assertEquals(0.274, value.getPValue(), 0.0);
	}
	@Test
	public void testMTKTpValueZeroCorr() {
		BinaryFunctionOp<Img<UnsignedByteType>, Img<UnsignedByteType>, Double> op =
				Functions.binary(ops, MaxTKendallTau.class, Double.class, zeroCorrelationImageCh1, zeroCorrelationImageCh2);
		final int[] blockSize = new int[zeroCorrelationImageCh1.numDimensions()];
		for (int d = 0; d < blockSize.length; d++) {
			final long size = (long) Math.floor(Math.sqrt(zeroCorrelationImageCh1.dimension(d)));
			blockSize[d] = (int) size;
		}
		RandomAccessibleInterval<UnsignedByteType> cropCh1 = ShuffledView.cropAtMin(zeroCorrelationImageCh1, blockSize);
		RandomAccessibleInterval<UnsignedByteType> cropCh2 = ShuffledView.cropAtMin(zeroCorrelationImageCh2, blockSize);
		PValueResult value = (PValueResult) ops.run(Ops.Coloc.PValue.class, new PValueResult(), cropCh1, cropCh2, op);
		assertEquals(0.609, value.getPValue(), 0.0);
	}
	@Test
	public void testMTKTpValuePosCorr() {
//		BinaryFunctionOp<Img<UnsignedByteType>, Img<UnsignedByteType>, Double> op =
//				Functions.binary(ops, MaxTKendallTau.class, Double.class, positiveCorrelationImageCh1, positiveCorrelationImageCh2);
//		final int[] blockSize = new int[positiveCorrelationImageCh1.numDimensions()];
//		for (int d = 0; d < blockSize.length; d++) {
//			final long size = (long) Math.floor(Math.sqrt(positiveCorrelationImageCh1.dimension(d)));
//			blockSize[d] = (int) size;
//		}
//		RandomAccessibleInterval<UnsignedByteType> cropCh1 = ShuffledView.cropAtMin(positiveCorrelationImageCh1, blockSize);
//		RandomAccessibleInterval<UnsignedByteType> cropCh2 = ShuffledView.cropAtMin(positiveCorrelationImageCh2, blockSize);
//		PValueResult value = (PValueResult) ops.run(Ops.Coloc.PValue.class, new PValueResult(), cropCh1, cropCh2, op);
//		System.out.println("posCorr p value = " + value.getPValue());
////			//assertEquals(0.274, value.getPValue(), 0.0);
	}
}
