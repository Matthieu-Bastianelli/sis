/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.sis.internal.storage;

import java.util.List;
import org.opengis.util.GenericName;
import org.apache.sis.coverage.SampleDimension;
import org.apache.sis.coverage.grid.GridExtent;
import org.apache.sis.coverage.grid.GridGeometry;
import org.apache.sis.coverage.grid.GridCoverage;
import org.apache.sis.test.TestCase;
import org.junit.jupiter.api.Test;

import static org.junit.Assert.*;


/**
 * Tests {@link AbstractGridResource} and {@link AbstractGridResource.RangeArgument}.
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @version 1.0
 * @since   1.0
 * @module
 */
public final strictfp class AbstractGridResourceTest  extends TestCase {
    /**
     * A resource performing no operation.
     */
    private final AbstractGridResource resource = new AbstractGridResource((AbstractGridResource) null) {
        @Override public GenericName           getIdentifier()       {throw new UnsupportedOperationException();}
        @Override public GridGeometry          getGridGeometry()     {throw new UnsupportedOperationException();}
        @Override public List<SampleDimension> getSampleDimensions() {throw new UnsupportedOperationException();}
        @Override public GridCoverage read(GridGeometry d, int... r) {throw new UnsupportedOperationException();}
    };

    /**
     * Tests {@link AbstractGridResource.RangeArgument} for data organized in a banded sample model.
     * This is the state when no {@code insert} method is invoked.
     */
    @Test
    public void testRangeArgumentForBandedModel() {
        final AbstractGridResource.RangeArgument r = resource.validateRangeArgument(7, new int[] {4, 6, 2});
        assertEquals("numBands",    3, r.getNumBands());
        assertEquals("first",       4, r.getFirstSpecified());
        assertEquals("source",      2, r.getSourceIndex(0));           // Expect sorted source indices: {2, 4, 6}.
        assertEquals("source",      4, r.getSourceIndex(1));
        assertEquals("source",      6, r.getSourceIndex(2));
        assertEquals("target",      2, r.getTargetIndex(0));           // Expect original positions of sorted indices: {2, 0, 1}.
        assertEquals("target",      0, r.getTargetIndex(1));
        assertEquals("target",      1, r.getTargetIndex(2));
        assertEquals("subsampled",  2, r.getSubsampledIndex(0));       // Expect equivalent to getSourceIndex(i).
        assertEquals("subsampled",  4, r.getSubsampledIndex(1));
        assertEquals("subsampled",  6, r.getSubsampledIndex(2));
        assertEquals("pixelStride", 1, r.getPixelStride());
    }

    /**
     * Tests {@link AbstractGridResource.RangeArgument} for data organized in an interleaved sample model.
     * This is the state when the {@code insert} methods are invoked.
     */
    @Test
    public void testRangeArgumentForInterleavedModel() {
        final AbstractGridResource.RangeArgument r = resource.validateRangeArgument(7, new int[] {4, 6, 2});
        assertEquals(3, r.insertBandDimension(new GridExtent(360, 180), 2).getDimension());
        assertArrayEquals(new int[] {3, 1, 2}, r.insertSubsampling(new int[] {3, 1}, 2));
        assertEquals("numBands",    3, r.getNumBands());
        assertEquals("first",       4, r.getFirstSpecified());
        assertEquals("source",      2, r.getSourceIndex(0));           // Expect sorted source indices: {2, 4, 6}.
        assertEquals("source",      4, r.getSourceIndex(1));
        assertEquals("source",      6, r.getSourceIndex(2));
        assertEquals("target",      2, r.getTargetIndex(0));           // Expect original positions of sorted indices: {2, 0, 1}.
        assertEquals("target",      0, r.getTargetIndex(1));
        assertEquals("target",      1, r.getTargetIndex(2));
        assertEquals("subsampled",  0, r.getSubsampledIndex(0));       // Expect source indices divided by 2 minus 1.
        assertEquals("subsampled",  1, r.getSubsampledIndex(1));
        assertEquals("subsampled",  2, r.getSubsampledIndex(2));
        assertEquals("pixelStride", 3, r.getPixelStride());
    }
}
