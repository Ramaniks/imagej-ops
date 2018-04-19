/*-
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import net.imagej.ops.Contingent;
import net.imagej.ops.Ops;
import net.imagej.ops.coloc.ColocUtil;
import net.imagej.ops.coloc.IntArraySorter;
import net.imagej.ops.coloc.MergeSort;
import net.imagej.ops.special.function.AbstractBinaryFunctionOp;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.IntervalIndexer;
import net.imglib2.util.Intervals;
import net.imglib2.util.IterablePair;
import net.imglib2.util.Pair;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

import org.scijava.plugin.Plugin;
import org.scijava.util.IntArray;

/**
 * This algorithm calculates Maximum Trunctated Kendall Tau (MTKT) from Wang et
 * al. (2017); computes thresholds using Otsu method.
 *
 * @param <T> Type of the first image
 * @param <U> Type of the second image
 * 
 * @author Shulei Wang
 * @author Ellen T Arena
 */
@Plugin(type = Ops.Coloc.MaxTKendallTau.class)
public class MTKT<T extends RealType<T>, U extends RealType<U>>
	extends AbstractBinaryFunctionOp<RandomAccessibleInterval<T>, RandomAccessibleInterval<U>, Double> implements
	Ops.Coloc.MaxTKendallTau, Contingent
{
	@Override
	public Double calculate(final RandomAccessibleInterval<T> image1, final RandomAccessibleInterval<U> image2) {
		// check image sizes
		// TODO: Add these checks to conforms().
		if (Intervals.equalDimensions(image1, image2)) {
			throw new IllegalArgumentException("Image dimensions do not match");
		}
		final long n1 = Intervals.numElements(image1);
		if (n1 > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("Image dimensions too large: " + n1);
		}
		final int n = (int) n1;

		// compute thresholds
		final double thresh1 = threshold(image1);
		final double thresh2 = threshold(image2);

		double[][] rank = rankTransformation(image1, image2, thresh1, thresh2, n);

		double maxtau = calculateMaxKendallTau(rank, thresh1, thresh2, n);

		return maxtau;
	}

	<V extends RealType<V>> double threshold(final RandomAccessibleInterval<V> image) {
		// call Otsu if explicit threshold was not given
		final Histogram1d<V> histogram = ops().image().histogram(Views.iterable(image));
		return ops().threshold().otsu(histogram).getRealDouble();
	}

	static <T extends RealType<T>, U extends RealType<U>> double[][] rankTransformation(final RandomAccessibleInterval<T> image1, final RandomAccessibleInterval<U> image2, final double thres1,
		final double thres2, final int n)
	{
		
/////////////////////////////////////////////////////////////////////////////////////////// Shulei's original version <<below>>
//		final double[][] tempRank = new double[n][2];
//		for (int i = 0; i < n; i++) {
//			tempRank[i][0] = values[i][0];
//			tempRank[i][1] = values[i][1];
//		}
//
//		Arrays.sort(tempRank, new Comparator<double[]>() {
//
//			@Override
//			public int compare(final double[] row1, final double[] row2) {
//				return Double.compare(row1[1], row2[1]);
//			}
//		});
//
//		int start = 0;
//		int end = 0;
//		int rank = 0;
//		while (end < n - 1) {
//			while (Double.compare(tempRank[start][1], tempRank[end][1]) == 0) {
//				end++;
//				if (end >= n) break;
//			}
//			for (int i = start; i < end; i++) {
//				tempRank[i][1] = rank + Math.random();
//			}
//			rank++;
//			start = end;
//		}
//
//		Arrays.sort(tempRank, new Comparator<double[]>() {
//
//			@Override
//			public int compare(final double[] row1, final double[] row2) {
//				return Double.compare(row1[1], row2[1]);
//			}
//		});
//
//		for (int i = 0; i < n; i++) {
//			tempRank[i][1] = i + 1;
//		}
//
//		// second
//		Arrays.sort(tempRank, new Comparator<double[]>() {
//
//			@Override
//			public int compare(final double[] row1, final double[] row2) {
//				return Double.compare(row1[0], row2[0]);
//			}
//		});
//
//		start = 0;
//		end = 0;
//		rank = 0;
//		while (end < n - 1) {
//			while (Double.compare(tempRank[start][0], tempRank[end][0]) == 0) {
//				end++;
//				if (end >= n) break;
//			}
//			for (int i = start; i < end; i++) {
//				tempRank[i][0] = rank + Math.random();
//			}
//			rank++;
//			start = end;
//		}
//
//		Arrays.sort(tempRank, new Comparator<double[]>() {
//
//			@Override
//			public int compare(final double[] row1, final double[] row2) {
//				return Double.compare(row1[0], row2[0]);
//			}
//		});
//
//		for (int i = 0; i < n; i++) {
//			tempRank[i][0] = i + 1;
//		}
//
//		final List<Integer> validIndex = new ArrayList<Integer>();
//		for (int i = 0; i < n; i++) {
//			if (tempRank[i][0] >= thres1 && tempRank[i][1] >= thres2) {
//				validIndex.add(i);
//			}
//		}
//
//		final int rn = validIndex.size();
//		final double[][] finalrank = new double[rn][2];
//		int index = 0;
//		for (final Integer i : validIndex) {
//			finalrank[index][0] = tempRank[i][0];
//			finalrank[index][1] = tempRank[i][1];
//			index++;
//		}
//
//		return finalrank;
/////////////////////////////////////////////////////////////////////////////////////////// Shulei's original version <<above>>		
		
		// FIRST...
		final IntArray rankIndex1 = rankSamples(image1);
		final IntArray rankIndex2 = rankSamples(image2);

		//////////////////////////////// TODO: Confirm dealing with thresholds (same as Shulei's method) <<below>>
		List<Integer> validIndex = new ArrayList<Integer>();
		for (int i = 0; i < n; i++)
		{
			if(rankIndex1.get(i) >= thres1 && rankIndex2.get(i) >= thres2)
			{
				validIndex.add(i);
			}
		}
		int rn=validIndex.size();
		double[][] finalRanks = new double[rn][2];
		int index = 0;
		for( Integer i : validIndex ) {
			finalRanks[index][0] = Math.floor(rankIndex1.get(i));
			finalRanks[index][1] = Math.floor(rankIndex2.get(i));
			index++;
		}
		////////////////////////////////Dealing with thresholds (same as Shulei's method) <<above>>
		return finalRanks;
	}

	private static <V extends RealType<V>> IntArray rankSamples(RandomAccessibleInterval<V> image) {
		final long elementCount = Intervals.numElements(image);
		if (elementCount > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("Image dimensions too large: " + elementCount);
		}
		final int n = (int) elementCount;

		// NB: Initialize rank index in random order, to ensure random tie-breaking.
		final IntArray rankIndex = new IntArray(n);		
		for (int i = 0; i < n; i++) {
			rankIndex.setValue(i, i);
		}
		Collections.shuffle(rankIndex);

		final V a = Util.getTypeFromInterval(image).createVariable();
		final RandomAccess<V> ra = image.randomAccess();
		Collections.sort(rankIndex, new Comparator<Integer>() {
			@Override
			public int compare(final Integer indexA, final Integer indexB) {
				IntervalIndexer.indexToPosition(indexA, image, ra);
				a.set(ra.get());
				IntervalIndexer.indexToPosition(indexB, image, ra);
				final V b = ra.get();
				return a.compareTo(b);
			}
		});

		return rankIndex;
		// TODO: Consider returning Img<UnsignedIntType> or Img<UnsignedLongType>.
		// We could do this now by wrapping rankIndex via ArrayImgs.unsignedInts(rankIndex.getArray(), dims).
		// But better to use UnsignedLongType indices, or even Img<D> where D is N-dimensional coord type?
	}

	static double calculateMaxKendallTau(final double[][] rank,
		final double thresholdRank1, final double thresholdRank2, final int n)
	{
		final int rn = rank.length;
		int an;
		final double step = 1 + 1.0 / Math.log(Math.log(n));
		double tempOff1 = 1;
		double tempOff2;
		List<Integer> activeIndex;
		double sdTau;
		double kendallTau;
		double normalTau;
		double maxNormalTau = Double.MIN_VALUE;

		while (tempOff1 * step + thresholdRank1 < n) {
			tempOff1 *= step;
			tempOff2 = 1;
			while (tempOff2 * step + thresholdRank2 < n) {
				tempOff2 *= step;

				activeIndex = new ArrayList<Integer>();
				for (int i = 0; i < rn; i++) {
					if (rank[i][0] >= n - tempOff1 && rank[i][1] >= n - tempOff2) {
						activeIndex.add(i);
					}
				}
				an = activeIndex.size();
				if (an > 1) {
					kendallTau = calculateKendallTau(rank, activeIndex);
					sdTau = Math.sqrt(2.0 * (2 * an + 5) / 9 / an / (an - 1));
					normalTau = kendallTau / sdTau;
				}
				else {
					normalTau = Double.MIN_VALUE;
				}
				if (normalTau > maxNormalTau) maxNormalTau = normalTau;
			}
		}

		return maxNormalTau;
	}

	static double calculateKendallTau(final double[][] rank,
		final List<Integer> activeIndex)
	{
		final int an = activeIndex.size();
		final double[][] partRank = new double[2][an];
		int indicatr = 0;
		for (final Integer i : activeIndex) {
			partRank[0][indicatr] = rank[i][0];
			partRank[1][indicatr] = rank[i][1];
			indicatr++;
		}
		final double[] partRank1 = partRank[0];
		final double[] partRank2 = partRank[1];

		final int[] index = new int[an];
		for (int i = 0; i < an; i++) {
			index[i] = i;
		}

		IntArraySorter.sort(index, (a, b) -> Double.compare(partRank1[a], partRank1[b]));

		final MergeSort mergeSort = new MergeSort(index, (a, b) -> Double.compare(partRank2[a], partRank2[b]));

		final long n0 = an * (long) (an - 1) / 2;
		final long S = mergeSort.sort();

		return (n0 - 2 * S) / (double) n0;

	}

	@Override
	public boolean conforms() {
		return true;
	}
}
