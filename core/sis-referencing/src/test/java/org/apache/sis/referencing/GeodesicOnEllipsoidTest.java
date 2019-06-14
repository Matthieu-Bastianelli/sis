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

import static java.lang.Math.toRadians;
import javax.measure.Unit;
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
     * article ofCharles F. F. Karney (SRI International)</a>
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
     * Test on a sample of the direct problem resolution.
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
        assertTrue(testIndicator); //test tha no exception were thrown by #isisEndPointComputable()
        
        testedEarth.getEndPoint();
        
        assertEquals("End point latitude", 41.79331020506, testedEarth.getEndPoint().getCoordinate()[1], 10e-11);
        assertEquals("End point azimuth", 149.09016931807, testedEarth.getEndingAzimuth(), 10e-11);
        

    }
}
