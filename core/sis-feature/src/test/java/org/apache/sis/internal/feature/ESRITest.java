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
package org.apache.sis.internal.feature;

import com.esri.core.geometry.Polyline;
import org.apache.sis.util.StringBuilders;
import org.junit.jupiter.api.Test;

import static org.junit.Assert.*;


/**
 * Tests {@link ESRI} implementation.
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @version 1.0
 * @since   1.0
 * @module
 */
public final strictfp class ESRITest extends GeometriesTestCase {
    /**
     * Creates a new test case.
     */
    public ESRITest() {
        super(new ESRI());
    }

    /**
     * Tests {@link ESRI#createPolyline(int, Vector...)}.
     */
    @Test
    @Override
    public void testCreatePolyline() {
        super.testCreatePolyline();
        final Polyline poly = (Polyline) geometry;
        assertEquals("pathCount", 2, poly.getPathCount());
    }

    /**
     * Tests {@link Geometries#tryMergePolylines(Object, Iterator)}.
     */
    @Test
    @Override
    public void testTryMergePolylines() {
        super.testTryMergePolylines();
        final Polyline poly = (Polyline) geometry;
        assertEquals("pathCount", 3, poly.getPathCount());
    }

    /**
     * Verifies that a WKT is equal to the expected one. This method modifies the expected WKT
     * by transforming single geometries into multi-geometries, since the ESRI library formats
     * geometries that way (at least with the objects that we use).
     */
    @Override
    void assertWktEquals(String expected, final String actual) {
        assertTrue(actual.startsWith("MULTI"));
        final StringBuilder b = new StringBuilder(expected.length() + 7).append("MULTI").append(expected);
        StringBuilders.replace(b, "(", "((");
        StringBuilders.replace(b, ")", "))");
        expected = b.toString();
        super.assertWktEquals(expected, actual);
    }
}
