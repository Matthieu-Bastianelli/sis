/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */
package org.apache.sis.referencing;

import org.apache.sis.internal.referencing.Formulas;
import org.apache.sis.internal.referencing.Resources;
import org.apache.sis.measure.Units;
import org.apache.sis.referencing.crs.HardCodedCRS;
import org.junit.Test;
import static org.opengis.test.Assert.*;
import org.opengis.referencing.operation.TransformException;

/**
 * Tests {@link GeodesicsOnEllipsoid}.
 *
 *
 */
public class GeodesicOnEllipsoidTest {

    /**
     * Test if the tested ellipsoid is compatible with the one used in the
     * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * article of Charles F F Karney (SRI International)</a>.
     */
    @Test
    public void testGeodesicsOnEllipsoidConstructor() {

        GeodesicsOnEllipsoid testedEarth = new GeodesicsOnEllipsoid(HardCodedCRS.WGS84);

        assertEquals("Inverse Flattening", 298.257223563, testedEarth.ellipsoid.getInverseFlattening(), 10e-9);
        assertEquals("Major Axis", 6378137, testedEarth.ellipsoid.getSemiMajorAxis(), 0.1);
        assertEquals("MinorAxis", 6356752.314245, testedEarth.ellipsoid.getSemiMinorAxis(), 10e-6);
        assertEquals("thirs flattening n", 0.00167922038638370, testedEarth.getThirdFlattening(), 10e-17);
        assertEquals("thirs flattening n", 0.00673949674227643, testedEarth.getSecondEccentricity() * testedEarth.getSecondEccentricity(), 10e-17);

    }

    /**
     * Tests azimuths at poles. Non-regression Test
     * {@link GeodeticCalculatorTest#testAzimuthAtPoles()}
     */
    @Test
    public void testAzimuthAtPoles() {

        GeodesicsOnEllipsoid testedEarth = new GeodesicsOnEllipsoid(HardCodedCRS.WGS84);

        testedEarth.setStartGeographicPoint(90, 30);
        final double tolerance = 0.2;
        testedEarth.setEndGeographicPoint(20, 20);
        assertEquals(-170, testedEarth.getStartingAzimuth(), tolerance);
        testedEarth.setEndGeographicPoint(20, 40);
        assertEquals(170, testedEarth.getStartingAzimuth(), tolerance);
        testedEarth.setEndGeographicPoint(20, 30);
        assertEquals(180, testedEarth.getStartingAzimuth(), tolerance);
        testedEarth.setEndGeographicPoint(-20, 30);
        assertEquals(180, testedEarth.getStartingAzimuth(), tolerance);
        testedEarth.setEndGeographicPoint(-90, 30);
        assertEquals(180, testedEarth.getStartingAzimuth(), tolerance);

        testedEarth.setStartGeographicPoint(90, 0);
        testedEarth.setEndGeographicPoint(20, 20);
        assertEquals(160, testedEarth.getStartingAzimuth(), tolerance);
        testedEarth.setEndGeographicPoint(20, -20);
        assertEquals(-160, testedEarth.getStartingAzimuth(), tolerance);
        testedEarth.setEndGeographicPoint(20, 0);
        assertEquals(180, testedEarth.getStartingAzimuth(), tolerance);
        testedEarth.setEndGeographicPoint(-90, 0);
        assertEquals(180, testedEarth.getStartingAzimuth(), tolerance);

    }

    /**
     * Direct geodesic problem solving test with
     * {@link GeodesicsOnEllipsoid#computeEndPoint()}.
     *
     * This test method is based on the given sample (page 46) of
     * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * article of Charles F F Karney (SRI International)</a>
     *
     * The test {@link #testGeodesicsOnEllipsoidConstructor()} must be passed to
     * ensure that the current test is based on the correct ellipsoid of revolution.
     * 
     *       φ1 = 40°  (Given)
     *       α1 = 30°  (Given)
     *       S12 = geodesicDistance = 10 000 000 m  (Given)
     * 
     * In order to facilitate the possible debugging of {@link GeodesicsOnEllipsoid#computeEndPoint()}
     * Most of the the intermediate variables from Karney's article are checked 
     * in this method. The unchecked parameters are indicated in code's comments.
     * 
     * 
     *       
     *
     */
    @Test
    public void directCalculationOnSampleTest() throws TransformException {
        GeodesicsOnEllipsoid testedEarth = new GeodesicsOnEllipsoid(HardCodedCRS.WGS84);

        testedEarth.setStartGeographicPoint(40, 10);
        assertEquals("Valid Start Point", testedEarth.getStartPoint().getCoordinate()[0], 10.0, Formulas.ANGULAR_TOLERANCE);
        assertEquals("Valid Start Point", testedEarth.getStartPoint().getCoordinate()[1], 40.0, Formulas.ANGULAR_TOLERANCE);

        //Check that the ending point is unreachable without azimuth and geodesic precisions.
        boolean testIndicator = false;
        try {
            testedEarth.getEndPoint();
        } catch (IllegalStateException exception) {
            assertTrue("IllegalStateException", (exception.getMessage() != null)
                    && (exception.getMessage().equals(Resources.format(Resources.Keys.AzimuthAndDistanceNotSet))));
            testIndicator = true;
        }
        assertTrue(testIndicator); //ensure the assertion on IllegaleStateException were raised

        testedEarth.setStartingAzimuth(30);
        assertTrue("Distance Unit Meters", Units.METRE.equals(testedEarth.getDistanceUnit()));
        testedEarth.setGeodesicDistance(10000000);
        try {
            testedEarth.isEndPointComputable();
        } catch (Exception exception) {
            testIndicator = false;
        }
        assertTrue(testIndicator); //test that no exception were thrown by #isisEndPointComputable()

        testedEarth.getEndPoint();

        // Test of start point auxiliary parameters :
        assertEquals("Start point reduce latitude β1", 39.90527714601, Math.toDegrees(testedEarth.auxiliarySpheres.β1), 10e-12);
        assertEquals("Gesodesic's Equatorial Azimuth α0", 22.55394020262, Math.toDegrees(testedEarth.auxiliarySpheres.α0), 10e-12);
        assertEquals("Arc lentgth from Equator of the starting point σ1 on auxiliary sphere", 43.99915364500, Math.toDegrees(testedEarth.auxiliarySpheres.σ1), 10e-12);
        assertEquals("ω1, spherical longitudes of the starting point (from Equatorial point E) on the auxiliary sphere ", 20.32371827837, Math.toDegrees(testedEarth.auxiliarySpheres.ω1), 10e-12);
        
        /**
         * Untested parameters indicated to facilitate potential debug :
         * 
         * k² = 0.00574802962857
         * ε = 0.00143289220416
         * A1 coefficient = 1.00143546236207
         * I1(σ1) = 0.76831538886412
         */
                
        assertEquals("Start point's distance from Equatorial point E along the geodesic s1", 4883990.626232, testedEarth.spherical2EllipsoidalArcLength(testedEarth.auxiliarySpheres.σ1), 10e-7);
        
        /**
         *  s2 = 14 883 990.626 232
         *  τ2 = 133.96266050208 °
         */
        
        assertEquals("Arc lentgth from Equator of the ending point σ2 on auxiliary sphere", 133.92164083038, Math.toDegrees(testedEarth.auxiliarySpheres.σ2), 10e-12);
        assertEquals("End point reduce latitude β2", 41.69771809250, Math.toDegrees(testedEarth.auxiliarySpheres.β2), 10e-12);
        assertEquals("ω2, spherical longitudes of the starting point (from Equatorial point E) on the auxiliary sphere", 158.28412147112, Math.toDegrees(testedEarth.auxiliarySpheres.ω2), 10e-12);
        
        /**
         * A3 coefficient = 0.99928424306
         * I3(σ1) = 0.76773786069
         * I3(σ2) = 2.33534322170
         */
        
        assertEquals("λ01, longitudes of the starting point (from Equatorial point E)", 20.26715038016, Math.toDegrees(testedEarth.spherical2EllipsoidalLongitude(testedEarth.auxiliarySpheres.ω1, testedEarth.auxiliarySpheres.σ1)), 10e-12);
        assertEquals("λ02, longitudes of the ending point (from Equatorial point E)", 158.11205042393, Math.toDegrees(testedEarth.spherical2EllipsoidalLongitude(testedEarth.auxiliarySpheres.ω2, testedEarth.auxiliarySpheres.σ2)), 10e-12);
        
        //==========
        // Solution :
        //==========
        assertEquals("λ12, delta longitudes from the starting to the ending point", 137.84490004377,  Math.toDegrees(testedEarth.λ2 - testedEarth.λ1), 10e-12);
        assertEquals("End point latitude φ2", 41.79331020506, testedEarth.getEndPoint().getCoordinate()[1], 10e-12);
        assertEquals("End point azimuth α2", 149.09016931807, testedEarth.getEndingAzimuth(), 10e-12);

    }
}
