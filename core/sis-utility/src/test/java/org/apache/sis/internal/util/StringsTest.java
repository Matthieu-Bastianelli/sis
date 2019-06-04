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
package org.apache.sis.internal.util;

import org.apache.sis.test.TestCase;
import org.apache.sis.util.Characters;
import org.junit.jupiter.api.Test;

import static org.junit.Assert.*;


/**
 * Tests the {@link Strings} class.
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @version 0.8
 * @since   0.3
 * @module
 */
public final strictfp class StringsTest extends TestCase {
    /**
     * Tests the {@link Strings#toString(Class, Object[])} method.
     */
    @Test
    public void testToString() {
        assertEquals("Number[base=“decimal”, value=20]", Strings.toString(Number.class, "base", "decimal", "value", 20));
    }

    /**
     * Tests the {@link Strings#toUpperCase(String, Characters.Filter)} method.
     */
    @Test
    public void testToUpperCase() {
        final String expected = "WGS84";
        assertSame  (expected, Strings.toUpperCase(expected, Characters.Filter.LETTERS_AND_DIGITS));
        assertEquals(expected, Strings.toUpperCase("WGS 84", Characters.Filter.LETTERS_AND_DIGITS));
        assertEquals(expected, Strings.toUpperCase("wgs 84", Characters.Filter.LETTERS_AND_DIGITS));
    }
}
