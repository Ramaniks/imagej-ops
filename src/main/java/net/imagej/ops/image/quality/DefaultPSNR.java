/*
 * #%L
 * ImageJ software for multidimensional image processing and analysis.
 * %%
 * Copyright (C) 2014 - 2018 ImageJ developers.
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

package net.imagej.ops.image.quality;

import net.imagej.ops.Contingent;
import net.imagej.ops.OpService;
import net.imagej.ops.Ops;
import net.imagej.ops.map.Maps;
import net.imagej.ops.special.hybrid.AbstractBinaryHybridCF;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;

import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * Computes peak signal-to-noise ratio (PSNR) between a reference image and a
 * (noisy) test image. The resulting PSNR is expressed in decibel.
 * <p>
 * Computations are based on the definitions of Gonzalez (R.C. Gonzalez and R.E.
 * Woods, "Digital Image Processing," Prentice Hall 2008).
 * </p>
 *
 * @author Stefan Helfrich (University of Konstanz)
 * @param <I> type of input elements
 */
@Plugin(type = Ops.Image.PSNR.class)
public class DefaultPSNR<I extends RealType<I>> extends
	AbstractBinaryHybridCF<IterableInterval<I>, IterableInterval<I>, DoubleType>
	implements Ops.Image.PSNR, Contingent
{

	@Parameter
	private OpService opService;

	@Override
	public void compute(final IterableInterval<I> input1,
		final IterableInterval<I> input2, final DoubleType output)
	{
		final Cursor<I> cursor = input1.cursor();
		final Cursor<I> cursor2 = input2.cursor();
		double max = opService.stats().max(input1).getRealDouble();
		max *= max;

		double denominatorSum = 0d;
		while (cursor.hasNext()) {
			final double r = cursor.next().getRealDouble();
			final double t = cursor2.next().getRealDouble();
			denominatorSum += Math.pow(r - t, 2);
		}

		denominatorSum *= 1d / Intervals.numElements(input1);
		final double psnr = 10 * Math.log10(max / denominatorSum);

		output.setReal(psnr);
	}

	@Override
	public boolean conforms() {
		return Intervals.equalDimensions(in1(), in2()) && //
			Maps.compatible(in1(), in2());
	}

	@Override
	public DoubleType createOutput(final IterableInterval<I> input1,
		final IterableInterval<I> input2)
	{
		return new DoubleType();
	}

}