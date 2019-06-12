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

    final AuxiliarySphereParameters auxiliarySphere = new AuxiliarySphereParameters();

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
     * The series' coefficients of Fourier and Taylor considered to compute the
     * geodesic are those presented in C.F.F Karney's article : Eq. 15 and 17
     *
     * @throws IllegalStateException if the azimuth and the distance have not
     * been set.
     */
    @Override
    void computeEndPoint() {
        isEndPointComputable(); // Throws IllegalStateException

        //Given φ1, dλ1, dφ1 and geodesicDistance (s12 in C.F.F Karney's article) :
        final double α1 = atan2(dλ1, dφ1);  // Valable pour une ellipse????     // tan(π/2 − θ)  =  1/tan(θ)
        auxiliarySphere.computeFromDirectProblem(ellipsoid.getInverseFlattening(), φ1, α1); // compute β1, α0, σ1, ω1
        final double s1 = spherical2EllipsoidalArcLength(auxiliarySphere.α0, auxiliarySphere.σ1); //compute s1

    }

    /**
     * Compute s value from σ and α0.
     *
     * @param αo azimuth in the forward direction at the equator.
     * @param σ spherical arc length on the auxiliary sphere between the equatorial
     * point which forward direction following αo azimuth cross the considered point.
     * (σ1 for starting point, σ2 for ending point).
     * 
     * @return the estimation of the associated ellipsoidal arc length s in
     * C.F.F Karney 's article.
     */
    private double spherical2EllipsoidalArcLength(final double αo, final double σ) {

//        //computation of the ellipsoid's second eccentricity :
//        final double e2nd = sqrt(pow(ellipsoid.getSemiMajorAxis(), 2) - pow(ellipsoid.getSemiMinorAxis(), 2)) / ellipsoid.getSemiMinorAxis();  //e'² = (a² - b²) /b²
        //Expansion parameter:
        final double ε = ((sqrt(1 + pow(compute2ndEccentriciy() * cos(αo), 2))) - 1)
                / ((sqrt(1 + pow(compute2ndEccentriciy() * cos(αo), 2))) + 1);

        final double A1 = (((pow(ε, 6) / 256) + (pow(ε, 4) / 64) + (pow(ε, 2) / 4) + 1)) / (1 - ε);

//========================== TO DELETE AFTER TESTS =============================
//        //Equation 7 together with equation 15 :
//        return ((-7 / 2048) * pow(ε, 6) * sin(2 * 6 * σ) //Term with C16
//                - (7 / 1280) * pow(ε, 5) * sin(2 * 5 * σ) //Term with C15
//                + (3 / 512) * pow(ε, 6) * sin(2 * 4 * σ) + (-5 / 512) * pow(ε, 4) * sin(2 * 4 * σ) //Term with C14
//                + (3 / 256) * pow(ε, 5) * sin(2 * 3 * σ) - (1 / 48) * pow(ε, 3) * sin(2 * 3 * σ) //Term with C13
//                - (9 / 2048) * pow(ε, 6) * sin(2 * 2 * σ) + (1 / 32) * pow(ε, 4) * sin(2 * 2 * σ) - (1 / 16) * pow(ε, 2) * sin(2 * 2 * σ) //Term with C12
//                - (1 / 32) * pow(ε, 5) * sin(2 * 1 * σ) + (3 / 16) * pow(ε, 3) * sin(2 * 1 * σ) - (1 / 2) * ε * sin(2 * 1 * σ) //Term with C11
//                + σ) * A1 * ellipsoid.getSemiMinorAxis();
//==============================================================================
        //Equation 7 together with equation 15 using term_C1l(....)
        return (term_C1l((-7 / 2048), ε, 6, 6, σ) + term_C1l((-9 / 2048), ε, 6, 2, σ) + term_C1l((3 / 512), ε, 6, 4, σ) //Termes ε ^ 6
                + term_C1l((-7 / 1280), ε, 5, 5, σ) + term_C1l((3 / 256), ε, 5, 3, σ) + term_C1l((-1 / 32), ε, 5, 1, σ)//Termes ε ^ 5
                + term_C1l((-5 / 512), ε, 4, 4, σ) + term_C1l((1 / 32), ε, 4, 2, σ) //Termes ε ^ 4
                + term_C1l((-1 / 48), ε, 3, 3, σ) + term_C1l((3 / 16), ε, 3, 1, σ) //Termes ε ^ 3
                + term_C1l((-1 / 16), ε, 2, 2, σ) //Termes ε ^ 2
                + term_C1l((-1 / 2), ε, 1, 1, σ) //Termes ε ^ 1
                + σ) * A1 * ellipsoid.getSemiMinorAxis();

    }

    /**
     * Utility method to express the terms of C1l coefficient of equation 15
     * (Used to compute ellipsoidal ArcLength {
     *
     * @see #spherical2EllipsoidalArcLength}) in
     * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * Algortihms of geodesics from Charles F. F. Karney (SRI International)</a>
     *
     * computed_term = fraction * ε ^ exposant * sin(2 * l * σ)
     *
     * @param fraction
     * @param ε
     * @param exposant
     * @param l
     * @param σ
     * @return value of the coefficient of the form 
     * computed_term = fraction * ε ^ exposant * sin(2 * l * σ)
     */
    private double term_C1l(final double fraction, final double ε, final double exposant, final double l, final double σ) {
        return fraction * pow(ε, exposant) * sin(2 * l * σ);
    }

    /**
     * Compute the second eccentricity of the ellipsoid.
     *
     * @return e' = (a²-b²) / b²
     */
    private double compute2ndEccentriciy() {
        return sqrt(pow(ellipsoid.getSemiMajorAxis(), 2) - pow(ellipsoid.getSemiMinorAxis(), 2)) / ellipsoid.getSemiMinorAxis();
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
         * Distance (spherical arc length) between the intercepting E point of
         * the equator (0° latitude) and the starting P point (resp the ending
         * point) along the 'direction' given by the α1 (resp α2) azimuth.
         */
        double σ1, σ2;

        /**
         * Spherical longitudes on the auxiliary sphere of the starting and
         * ending point.
         *
         * The spherical longitude of a point on the auxiliary sphere is the
         * angle formed by the meridian of the intercepting point of the equator
         * (0° latitude) and the direction given by the point given azimuth.
         *
         */
        double ω1, ω2;

        /**
         * compute β1, α0, σ1, ω1 from ellipsoidal 's flattening , starting
         * point latitude and azimuth
         *
         * @param flattening ellipsoidal 's flattening
         * @param φ starting point's latitude
         * @param azimuth starting point's azimuth
         */
        final void computeFromDirectProblem(final double flattening, final double φ, final double azimuth) {
            β1 = atan((1 - (1 / flattening)) * tan(φ));
            this.computeα0(β1, azimuth);
            this.σ1 = computeσ(β1);
            this.ω1 = computeω(σ1);
        }

        /**
         * Compute the azimuth α0 at the intercepting E point of the equator (0°
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
         * Compute the spherical longitude between the 2 Meridians of E and P
         * points. {
         *
         * @see Fig.1 of
         * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
         * Algortihms of geodesics from Charles F. F. Karney (SRI
         * International)</a>
         *
         * @param σ
         * @return ω
         */
        private double computeω(final double σ) {
            //contrôle!!
            return atan2(sin(α0) * sin(σ), cos(σ));
        }

    }//end of class AuxiliarySphereTools
    /**
     * =========================================================================
     */

}
