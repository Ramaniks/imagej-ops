/*
 * #%L
 * SciJava Common shared library for SciJava software.
 * %%
 * Copyright (C) 2009 - 2016 Board of Regents of the University of
 * Wisconsin-Madison, Broad Institute of MIT and Harvard, and Max Planck
 * Institute of Molecular Cell Biology and Genetics.
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

package net.imagej.type;

import java.lang.reflect.Type;

import net.imglib2.img.NativeImg;

import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.type.TypeExtractor;
import org.scijava.type.TypeService;

/**
 * {@link TypeExtractor} plugin which operates on {@link NativeImg} objects.
 *
 * @author Curtis Rueden
 */
@Plugin(type = TypeExtractor.class)
public class NativeImgExtractor implements TypeExtractor<NativeImg<?, ?>> {

	@Parameter
	private TypeService typeService;

	@Override
	public Type typeOf(final NativeImg<?, ?> o, int n) {
		if (n < 0 || n > 1) throw new IndexOutOfBoundsException("" + n);

		if (o.size() == 0) return null;

		if (n == 0) return typeService.typeOf(o.update(null));
		return typeService.typeOf(o.firstElement());
	}

	@Override
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public Class<NativeImg<?, ?>> getRawType() {
		return (Class) NativeImg.class;
	}

}