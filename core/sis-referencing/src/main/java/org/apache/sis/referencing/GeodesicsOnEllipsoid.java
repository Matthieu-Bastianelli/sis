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

import org.apache.sis.internal.referencing.Resources;
import org.opengis.geometry.coordinate.Position;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import static java.lang.Math.*;

/**
 *
 * @version 1.0
 *
 * @see
 * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">Algortihms
 * of geodesics from Charles F. F. Karney (SRI International)</a>
 *
 * @since 1.0
 * @module
 */
final class GeodesicsOnEllipsoid extends GeodeticCalculator {

    final AuxiliarySphereParameters auxiliarySphereParameters = new AuxiliarySphereParameters();

    
    /**
     * Construct a new geodetic calculator for ellipsoid expecting coordinates
     * in the supplis CRS.
     *
     * The {@link GeodeticCalculator} constructor is called.
     *
     * !!!!!!!!!! NOTE A SUPPRIMER: remarque surAPI du constructeur pas
     * comprise!!!!!!!!
     *
     * @param crs the referencing system for the {@link Position} arguments and
     * return values.
     */
    public GeodesicsOnEllipsoid(CoordinateReferenceSystem crs) {
        super(crs);
        if (ellipsoid.getInverseFlattening() == Double.POSITIVE_INFINITY) {
            throw new IllegalStateException(Resources.format(Resources.Keys.AmbiguousEllipsoid_1, 0)); //PAS SUR du message
        }
    }

    /**
     * Computes the end point from the start point, the azimuth and the geodesic
     * distance. This method should be invoked if the end point or ending
     * azimuth is requested while {@link #END_POINT} validity flag is not set.
     *
     * <p>
     * This implementation computes approximations of {@link #φ2}, {@link #λ2}
     * and ∂φ/∂λ derivatives foolowing the method described in :
     * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">Algortihms
     * of geodesics from Charles F. F. Karney (SRI International)</a>
     * </p>
     *
     * @throws IllegalStateException if the azimuth and the distance have not
     * been set.
     */
    @Override
    void computeEndPoint() {
        isEndPointComputable(); // Throws IllegalStateException

        //Given φ1, dλ1, dφ1 and geodesicDistance (s12 in C.F.F Karney's article) :
        
        final double α1 = atan2(dλ1, dφ1);  // Valable pour une ellipse????     // tan(π/2 − θ)  =  1/tan(θ)
        auxiliarySphereParameters.computeFromDirectProblem(ellipsoid.getInverseFlattening(), φ1, α1);
        
        final double e2nd = sqrt( pow(ellipsoid.getSemiMajorAxis(),2) - pow(ellipsoid.getSemiMinorAxis(),2) )/ellipsoid.getSemiMinorAxis();  //e'² = (a² - b²) /b²
        

    }

    /**
     * =========================================================================
     * A METTRE DANS ELLIPSOID.JAVA???
     * =========================================================================
     *
     */
    private static class AuxiliarySphereParameters {

        /**
         * Azimuth at the intercepting point of the equator (0° latitude) and
         * the direction given by an α azimuth at the β reduced latitude.
         */
        double α0;
        /**
         * Reduced latitudes of the starting point (respectively the ending
         * point).
         */
        double β1, β2;

        /**
         * Distance (spherical arc length) between the starting point (resp the
         * ending point) and the intercepting point of the equator (0° latitude)
         * and the 'direction' given by the α1 (resp α2) azimuth.
         */
        double σ1, σ2;

        /**
         * Spherical longitudes on the auxiliary sphere of the starting and
         * ending point.
         * 
         * The spherical longitude of a point on the auxiliary sphere is the angle 
         * formed by the meridian of the intercepting point of the equator (0° latitude) and
         * the direction given by the point given azimuth.
         * 
         */
        double ω1, ω2;

        final void computeFromDirectProblem(final double flattening, final double φ, final double azimuth) {

            β1 = atan((1 - (1 / flattening)) * tan(φ));
            this.computeα0(β1, azimuth);
            this.σ1 = computeσ(β1);
            this.ω1 = computeω(σ1);

        }

        /**
         * Compute the azimuth α0 at the intercepting point of the equator (0°
         * latitude) and the direction given by a reduced latitude and an
         * associated azimuth.
         *
         * @param β reduced latitude
         * @param azimuth the given azimuth.
         * @return α0 value
         */
        private void computeα0(final double β, final double azimuth) {
            this.α0 = atan2(sin(azimuth) * cos(β), hypot(cos(azimuth), sin(azimuth) * sin(β)));
            //indicateur pour contrôle!!
        }

        /**
         * Compute the distance at the equator on the direction given by the
         *
         * @param α the given azimuth.
         * @param β the given reduced latitude
         * @return σ
         */
        private double computeσ(double β) {
            //contrôle!!
            return atan2(sin(β), cos(α0) * cos(β));
        }

        /**
         * Compute the spherical longitude between 
         *
         * @param σ
         * @return ω
         */
        private double computeω(final double σ) {
            //contrôle!!
            return atan2( sin(α0)*sin(σ) , cos(σ) );
        }

    }//end of class AuxiliarySphereTools
    /**
     * =========================================================================
     */

}
