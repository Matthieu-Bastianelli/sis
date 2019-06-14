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

import org.opengis.geometry.coordinate.Position;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import static java.lang.Math.*;
import org.apache.sis.util.ArgumentChecks;

/**
 *
 * @version 1.0
 *
 * @see
 * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">Algortihms
 * of geodesics from Charles F. F. Karney (SRI International)</a>
 *
 * As in the above article, we note E the equatorial point at which the geodesic
 * crosses the equator in the northward direction.
 *
 * @since 1.0
 * @module
 */
final class GeodesicsOnEllipsoid extends GeodeticCalculator {

    private final double SECOND_ECCENTRICITY;
    private final double THIRD_FLATTENING;
    private double ε;

    private static final int EQUATORIAL_AZIMUTH = 64, AUXILIARY_START_POINT = 128, AUXILIARY_END_POINT = 256, EXPAND_PARAMETER = 512;

    final AuxiliarySpheresParameters auxiliarySpheres;

    /**
     * Construct a new geodetic calculator for ellipsoid expecting coordinates
     * in the supplis CRS.
     *
     * The {@link GeodeticCalculator} constructor is called.
     *
     * !!!!!!!!!! NOTE A SUPPRIMER: remarque surAPI du super.constructeur pas
     * comprise!!!!!!!!
     *
     * @param crs the referencing system for the {@link Position} arguments and
     * return values.
     */
    public GeodesicsOnEllipsoid(CoordinateReferenceSystem crs) {
        super(crs);
        ArgumentChecks.ensureFinite("Inverse flattening of the ellipsoid", ellipsoid.getInverseFlattening());
        ArgumentChecks.ensureStrictlyPositive("Ellipsoid semi minor axis", ellipsoid.getSemiMinorAxis());

        // Compute the second eccentricity of the ellipsoid. e' = (a²-b²) / b²
        SECOND_ECCENTRICITY = sqrt(pow(ellipsoid.getSemiMajorAxis(), 2) - pow(ellipsoid.getSemiMinorAxis(), 2)) / ellipsoid.getSemiMinorAxis();

        // Compute the third flattening of the ellipsoid. n = (a-b) / (a+b)
        THIRD_FLATTENING = (ellipsoid.getSemiMajorAxis() - ellipsoid.getSemiMinorAxis()) / (ellipsoid.getSemiMajorAxis() + ellipsoid.getSemiMinorAxis());

        auxiliarySpheres = new AuxiliarySpheresParameters();
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

        auxiliarySpheres.parametersFromDirectProblem(ellipsoid.getInverseFlattening(), φ1, α1); // compute β1, α0, σ1, ω1

        final double s2 = spherical2EllipsoidalArcLength(auxiliarySpheres.σ1) + this.geodesicDistance; //s2 = s1 + s12
        auxiliarySpheres.resolveEndPoint(s2);

        // λ01 is the longitude angle from the equatorial point E
        // to the staring point. It matchs λ1  of the C.F.F. Karney, 2013 article  
        final double λ01 = spherical2EllipsoidalLongitude(auxiliarySpheres.ω1, auxiliarySpheres.σ1);
        final double λ02 = spherical2EllipsoidalLongitude(auxiliarySpheres.ω2, auxiliarySpheres.σ2);

        this.φ2 = reduced2latitude(auxiliarySpheres.β2); //contrôle!!
        this.λ2 = this.λ1 + λ02 - λ01; //λ2 =  λ1 + Δλ 
        auxiliarySpheres.computeEndPointAzimuth();  // affects dφ2 and dλ2

        validity |= END_POINT;
    }

    private void setExpandParameter() {
        if (isInvalid(AUXILIARY_START_POINT)) {
            throw new IllegalStateException("Auxiliary parameter α0 needed to compute the expand parameter ε");
        }

        ε = ((sqrt(1 + pow(SECOND_ECCENTRICITY * cos(auxiliarySpheres.α0), 2))) - 1)
                / ((sqrt(1 + pow(SECOND_ECCENTRICITY * cos(auxiliarySpheres.α0), 2))) + 1);
    }

    /**
     * Compute σ value from s
     *
     * @param s
     * @return σ
     */
    private double ellipsoidal2SphericalArcLength(final double s) {

        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        final double τ = s / computeBA1Coefficient(); //b * A1 can't be equal to 0 as b and A1 >0 //DOUBLON DE CALCUL BA1 AVEC spherical2EllipsoidalArcLength(....)!!!!

//========================== TO DELETE AFTER TESTS =============================
//        return term_Cil_sin2l(0.5, ε, 1, 1, τ) + term_Cil_sin2l(-0.28125, ε, 3, 1, τ) + term_Cil_sin2l(205/1536, ε, 5, 1, τ) //C'11
//                + term_Cil_sin2l(0.3125, ε, 2, 2, τ) + term_Cil_sin2l(37/96, ε, 4, 2, τ) + term_Cil_sin2l(1335/4096, ε, 6, 2, τ) //C'12
//                + term_Cil_sin2l(29/96, ε, 3, 3, τ) + term_Cil_sin2l( 0.5859375, ε, 5, 3, τ) //C'13
//                + term_Cil_sin2l(539/1536, ε, 4, 4, τ) + term_Cil_sin2l(0.933984375, ε, 6, 4, τ)//C'14
//                + term_Cil_sin2l(3467/7680, ε, 5, 5, τ)//C'15
//                + term_Cil_sin2l(38081/61440, ε, 6, 6, τ)//C'16
//                + τ;
//==============================================================================
        return term_Cil_sin2l(38081.0 / 61440.0, ε, 6, 6, τ) + term_Cil_sin2l(0.933984375, ε, 6, 4, τ) + term_Cil_sin2l(1335.0 / 4096.0, ε, 6, 2, τ) //Termes ε ^ 6
                + term_Cil_sin2l(3467.0 / 7680.0, ε, 5, 5, τ) + term_Cil_sin2l(0.5859375, ε, 5, 3, τ) + term_Cil_sin2l(205.0 / 1536.0, ε, 5, 1, τ) //Termes ε ^ 5
                + term_Cil_sin2l(539.0 / 1536.0, ε, 4, 4, τ) + term_Cil_sin2l(37.0 / 96.0, ε, 4, 2, τ) //Termes ε ^ 4
                + term_Cil_sin2l(29.0 / 96.0, ε, 3, 3, τ) + term_Cil_sin2l(-0.28125, ε, 3, 1, τ) //Termes ε ^ 3
                + term_Cil_sin2l(0.3125, ε, 2, 2, τ) //Termes ε ^ 2
                + term_Cil_sin2l(0.5, ε, 1, 1, τ) //Termes ε ^ 1
                + τ;
    }

    /**
     * Compute s value from σ and α0.
     *
     * @param σ spherical arc length on the auxiliary sphere between the
     * equatorial point which forward direction following αo azimuth cross the
     * considered point. (σ1 for starting point, σ2 for ending point).
     *
     * @return s : the estimation of the associated ellipsoidal arc length in
     * C.F.F Karney 's article.
     */
    private double spherical2EllipsoidalArcLength(final double σ) {

        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        final double bA1 = computeBA1Coefficient(); //b * A1

//========================== TO DELETE AFTER TESTS =============================
        //Equation 7 together with equation 15 :
//        return ((-7 / 2048) * pow(ε, 6) * sin(2 * 6 * σ) //Term with C16
//                - (7 / 1280) * pow(ε, 5) * sin(2 * 5 * σ) //Term with C15
//                + (3 / 512) * pow(ε, 6) * sin(2 * 4 * σ) - (5 / 512) * pow(ε, 4) * sin(2 * 4 * σ) //Term with C14
//                + (3 / 256) * pow(ε, 5) * sin(2 * 3 * σ) - (1 / 48) * pow(ε, 3) * sin(2 * 3 * σ) //Term with C13
//                - (9 / 2048) * pow(ε, 6) * sin(2 * 2 * σ) + (1 / 32) * pow(ε, 4) * sin(2 * 2 * σ) - (1 / 16) * pow(ε, 2) * sin(2 * 2 * σ) //Term with C12
//                - (1 / 32) * pow(ε, 5) * sin(2 * 1 * σ) + (3 / 16) * pow(ε, 3) * sin(2 * 1 * σ) - (1 / 2) * ε * sin(2 * 1 * σ) //Term with C11
//                + σ) * bA1;

//        return ((  (-7 / 2048) * pow(ε, 6)                                                   )* sin(2 * 6 * σ) //Term with C16
//                +( (-7 / 1280) * pow(ε, 5)                                                    )* sin(2 * 5 * σ) //Term with C15
//                +( (3 / 512)   * pow(ε, 6) - (5 / 512)* pow(ε, 4)                             )* sin(2 * 4 * σ) //Term with C14
//                +( (3 / 256)   * pow(ε, 5) - (1 / 48) * pow(ε, 3)                             )* sin(2 * 3 * σ) //Term with C13
//                +( (-9 / 2048) * pow(ε, 6) + (1 / 32) * pow(ε, 4) - (1 / 16) * pow(ε, 2)      )* sin(2 * 2 * σ) //Term with C12
//                +( (-1 / 32)   * pow(ε, 5) + (3 / 16) * pow(ε, 3) - (1 / 2)  * ε              ) * sin(2 * 1 * σ) //Term with C11
//                + σ) * bA1;
//==============================================================================
        // Equation 7 together with equation 15 using term_Cil(....)
        // The same terms than in Eq.18 are considered.
        return (term_Cil_sin2l((-7.0 / 2048.0), ε, 6, 6, σ) + term_Cil_sin2l((-9.0 / 2048.0), ε, 6, 2, σ) + term_Cil_sin2l((3.0 / 512.0), ε, 6, 4, σ) //Termes ε ^ 6
                + term_Cil_sin2l((-7.0 / 1280.0), ε, 5, 5, σ) + term_Cil_sin2l((3.0 / 256.0), ε, 5, 3, σ) + term_Cil_sin2l((-1.0 / 32.0), ε, 5, 1, σ)//Termes ε ^ 5
                + term_Cil_sin2l((-5.0 / 512.0), ε, 4, 4, σ) + term_Cil_sin2l((1.0 / 32.0), ε, 4, 2, σ) //Termes ε ^ 4
                + term_Cil_sin2l((-1.0 / 48.0), ε, 3, 3, σ) + term_Cil_sin2l((3.0 / 16.0), ε, 3, 1, σ) //Termes ε ^ 3
                + term_Cil_sin2l((-1.0 / 16.0), ε, 2, 2, σ) //Termes ε ^ 2
                + term_Cil_sin2l((-1.0 / 2.0), ε, 1, 1, σ) //Termes ε ^ 1
                + σ)
                * bA1;

    }

    /**
     * Compute the longitude angle λ from the equatorial point E and the
     * considered point on the ellipsoid from the associated with ω longitude on
     * the auxiliary sphere.
     *
     * Eq.8 together with Eq. 23 in
     * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * C.F.F karney 2013</a>
     *
     * @param ω longitude on the auxiliary sphere.
     * @param σ spherical arc length beetwen the E Equatorial point with αo
     * forward azimuth and the point of the auxiliary sphere associated with the
     * point of the ellipsoid the angle is computed.
     * @param n third flattening of the ellipsoid
     *
     * @return λ longitude angle from the equatorial point E to the considered
     * point
     */
    private double spherical2EllipsoidalLongitude(final double ω, final double σ) {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }

        final double n = THIRD_FLATTENING; //used for lisibility
        // I3(σ) = A3 * (σ + sum)
        final double A3 = computeA3Coefficient(n);

//========================== TO DELETE AFTER TESTS =============================
//        //Sum from equation 23 :
//        final double sum = term_Cil_sin2l(0.008203125, ε, 5, 5, σ) //C35
//                + term_Cil_sin2l((7-14*n)/512, ε, 4, 4, σ) + term_Cil_sin2l(0.013671875, ε, 5, 4, σ)  //C34
//                + term_Cil_sin2l((5-9*n+5*n*n)/192, ε, 3, 3, σ) + term_Cil_sin2l((3/128 - 5*n/192), ε, 4, 3, σ) + term_Cil_sin2l(0.013671875, ε, 5, 3, σ)//C33
//                + term_Cil_sin2l((2-3*n+n*n)/32, ε, 2, 2, σ) +term_Cil_sin2l((3-2*n-3*n*n)/64, ε, 3, 2, σ) + term_Cil_sin2l((3+n)/128, ε, 4, 2, σ) + term_Cil_sin2l(0.01953125, ε, 5, 2, σ)  //C32
//                + term_Cil_sin2l((1-n)/4, ε, 1, 1, σ) + term_Cil_sin2l((1-n*n)/8, ε, 2, 1, σ) + term_Cil_sin2l((3+3*n-n*n)/64, ε, 3, 1, σ) +  term_Cil_sin2l((5+2*n)/128, ε, 4 , 1, σ) +term_Cil_sin2l(0.0234375, ε, 5, 1, σ)  //C31
//                + σ;
//==============================================================================
        // Computation of the sum in Equation 23 of the C.F.F Karney's article.
        // The same terms than in Eq.25 are considered.
        final double sumC3l = term_Cil_sin2l(0.008203125, ε, 5, 5, σ) + term_Cil_sin2l(0.013671875, ε, 5, 4, σ) + term_Cil_sin2l(0.013671875, ε, 5, 3, σ) + term_Cil_sin2l(0.01953125, ε, 5, 2, σ) + term_Cil_sin2l(0.0234375, ε, 5, 1, σ) // ε⁵ terms
                + term_Cil_sin2l((7 - 14 * n) / 512.0, ε, 4, 4, σ) + term_Cil_sin2l((3.0 / 128.0 - 5 * n / 192.0), ε, 4, 3, σ) + term_Cil_sin2l((3 + n) / 128.0, ε, 4, 2, σ) + term_Cil_sin2l((5 + 2 * n) / 128.0, ε, 4, 1, σ) // ε⁴ terms
                + term_Cil_sin2l((5 - 9 * n + 5 * n * n) / 192, ε, 3, 3, σ) + term_Cil_sin2l((3 - 2 * n - 3 * n * n) / 64, ε, 3, 2, σ) + term_Cil_sin2l((3 + 3 * n - n * n) / 64.0, ε, 3, 1, σ) // ε³ terms
                + term_Cil_sin2l((2 - 3 * n + n * n) / 32, ε, 2, 2, σ) + term_Cil_sin2l((1 - n * n) / 8, ε, 2, 1, σ) // ε² terms
                + term_Cil_sin2l((1 - n) / 4.0, ε, 1, 1, σ) //ε term
                + σ;

        return ω - (1 / ellipsoid.getInverseFlattening()) * sin(auxiliarySpheres.α0) * A3 * sumC3l;
    }

    /**
     * return the latitude on the ellipsoid from the reduced latitude on the
     * associated auxiliary sphere.
     *
     * @param β the reduced latitude on the auxiliary sphere
     * @return φ the latitude on the ellipsoid.
     */
    private double reduced2latitude(final double β) {
        ArgumentChecks.ensureStrictlyPositive("Semi Minor Axis", ellipsoid.getSemiMinorAxis());

        return atan(tan(β) / (1 - 1 / ellipsoid.getInverseFlattening())); // flattening != 1 as semi minor axis b >0 from constructor.
    }

    /**
     * Utility method to express the terms of C1l coefficient of equation 15
     * (Used to compute ellipsoidal ArcLength {
     *
     * @see #spherical2EllipsoidalArcLength}) in
     * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * Algortihms of geodesics from Charles F. F. Karney (SRI International)</a>
     *
     * computed_term = coef * param ^ exposant * sin(2 * l * σ)
     *
     * @param coef real coefficient
     * @param param expansion parameter
     * @param exposant of the expansion parameter for the considered term.
     * @param l number of the fourier serie's term the current term is computed
     * for.
     * @param σ spheric arc length of the computation
     *
     * @return value of the coefficient of the form computed_term = coef * ε ^
     * exposant * sin(2 * l * σ)
     */
    private double term_Cil_sin2l(final double coef, final double param, final double exposant, final double l, final double σ) {
        return coef * pow(param, exposant) * sin(2 * l * σ);
    }

    /**
     * Return the A1 coefficient necessary to process arc length convertions
     * between ellipsoid and the auxiliary spheres.
     * {@link #ellipsoidal2SphericalArcLength(double)} and 
     * {@link #spherical2EllipsoidalArcLength(double) }
     *
     * @param n : third flattening of the ellipsoid
     * @return A3 value from terms suggest in C.F.F. Karney, 2013.
     */
    private double computeBA1Coefficient() {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        if (ε == 1) {
            throw new IllegalStateException("Impossible to compute BA1Coefficient for 1 value of ε expand parameter.");
        }
        return ((((pow(ε, 6) / 256) + (pow(ε, 4) / 64) + (pow(ε, 2) / 4) + 1)) / (1 - ε)) * ellipsoid.getSemiMinorAxis(); //b * A1 , with b the ellipsoid semi minor axis
    }

    /**
     * Return the A3 coefficient necessary to convert the spherical longitude in
     * ellipsoidal longitude.
     * {@link #spherical2EllipsoidalLongitude(double, double)}
     *
     * @param n : third flattening of the ellipsoid
     * @return A3 value from terms suggest in C.F.F. Karney, 2013.
     */
    private double computeA3Coefficient(final double n) {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }

        return (-3 / 128) * pow(ε, 5)
                + ((3 / 64) + (n / 32)) * pow(ε, 4)
                + (1 / 16) * (1 + 3 * n + n * n) * pow(ε, 3)
                + (1 / 8) * (2 + n + 3 * n * n) * pow(ε, 2)
                + 0.5 * (1 - n) * ε
                + 1;
    }

//    /**
//     * Compute the second eccentricity of the ellipsoid.
//     *
//     * @return e' = (a²-b²) / b²
//     */
//    private double secondEccentriciy() {
//        return sqrt(pow(ellipsoid.getSemiMajorAxis(), 2) - pow(ellipsoid.getSemiMinorAxis(), 2)) / ellipsoid.getSemiMinorAxis();
//    }
//
//    /**
//     * Compute the third flattening of the ellipsoid.
//     *
//     * @return (a-b) / (a+b)
//     */
//    private double thirdFlattening() {
//        return (ellipsoid.getSemiMajorAxis() - ellipsoid.getSemiMinorAxis()) / (ellipsoid.getSemiMajorAxis() + ellipsoid.getSemiMinorAxis());
//    }
    //==========================================================================
    //Internal Classe
    //==========================================================================
    /**
     * Internal classe used to handle parameters on used auxiliary spheres.
     *
     */
    private class AuxiliarySpheresParameters {

        /**
         * Azimuth at the intercepting point of the equator (0° latitude) and
         * the direction given by an α azimuth at the β reduced latitude.
         */
        public double α0;
        /**
         * Reduced latitudes of the starting point (respectively the ending
         * point).
         */
        public double β1, β2;

        /**
         * Distance (spherical arc length) between the intercepting E point of
         * the equator (0° latitude) and the starting P point (resp the ending
         * point) along the 'direction' given by the α1 (resp α2) azimuth.
         */
        public double σ1, σ2;

        /**
         * Spherical longitudes on the auxiliary sphere of the starting and
         * ending point.
         *
         * The spherical longitude of a point on the auxiliary sphere is the
         * angle formed by the meridian of the intercepting point of the equator
         * (0° latitude) and the direction given by the point given azimuth.
         *
         */
        public double ω1, ω2;

        /**
         * compute β1, α0, σ1, ω1 from ellipsoidal 's flattening , starting
         * point latitude and azimuth
         *
         * @param inverseFlattening inverse of the ellipsoidal's flattening
         * @param φ starting point's latitude
         * @param azimuth starting point's azimuth
         */
        final void parametersFromDirectProblem(final double inverseFlattening, final double φ, final double azimuth) {
            β1 = computeβ(inverseFlattening, φ);
            this.computeα0(β1, azimuth);
            this.σ1 = computeσ(β1, azimuth);
            this.σ1 = toRadians(43.99915364500);
            this.ω1 = computeω(σ1);
            validity |= AUXILIARY_START_POINT; //validity : bitmask of GeodesicsOnEllipsoid
        }

        /**
         * Compute the parameters on the auxiliary associated with the ending
         * point from the length of the geodesic on the ellipsoid s2 and the
         * starting point.
         *
         * @param s2
         * @return
         */
        final void resolveEndPoint(final double s2) {
            σ2 = ellipsoidal2SphericalArcLength(s2);
            β2 = computeβ(σ2);
            ω2 = computeω(σ2);
            validity |= AUXILIARY_END_POINT; //validity : bitmask of GeodesicsOnEllipsoid
        }

        /**
         * Compute the azimuth α0 at the intercepting E point of the equator (0°
         * latitude) and the direction given by a reduced latitude and an
         * associated azimuth.
         *
         * Eq. 10 in
         * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
         * Algortihms of geodesics from Charles F. F. Karney (SRI
         * International)</a>
         *
         * @param β reduced latitude.
         * @param azimuth the associated azimuth.
         */
        private void computeα0(final double β, final double azimuth) {
            this.α0 = atan2(sin(azimuth) * cos(β), hypot(cos(azimuth), sin(azimuth) * sin(β)));
            validity |= EQUATORIAL_AZIMUTH; //validity : bitmask of GeodesicsOnEllipsoid
        }

        /**
         * Compute distance on the great circle σ between the E point with α0
         * azimuth of the equator and the current point with the reduced
         * latitude β
         *
         * Eq. 11 in
         * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
         * Algortihms of geodesics from Charles F. F. Karney (SRI
         * International)</a>
         * α0 must be known.
         *
         * @param β the given reduced latitude
         * @param azimuth the associated azimuth.
         * @return distance on the great circle (EP fig. 2 of the article)
         */
        private double computeσ(double β, final double azimuth) {
            if (isInvalid(EQUATORIAL_AZIMUTH)) {
                throw new IllegalStateException("Uncomputed auxiliary parameter α0");
            }

            return atan2(sin(β), (cos(azimuth) * cos(β)));
        }

        /**
         * Compute the spherical longitude between the 2 Meridians of E and P
         * points.
         *
         * Eq. 12
         *
         * Fig.1 of
         * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
         * Algortihms of geodesics from Charles F. F. Karney (SRI
         * International)</a>
         *
         * @param σ
         * @return ω
         */
        private double computeω(final double σ) {
            if (isInvalid(EQUATORIAL_AZIMUTH)) {
                throw new IllegalStateException("Uncomputed auxiliary parameter α0");
            }
            return atan2(sin(α0) * sin(σ), cos(σ));
        }

        /**
         * Compute the reduce latitude β from φ
         *
         * Eq. 6
         *
         * @param inverseFlattening inverse of the ellipsoidal's flattening
         * @param φ latitude on the ellipsoid
         * @return β associated reduced latitude
         */
        private double computeβ(final double inverseFlattening, final double φ) {
            if ((inverseFlattening == 0) | (inverseFlattening == 1)) {
                throw new IllegalStateException("Uncomputed auxiliary parameter α0");
            }
            return atan((1 - (1 / inverseFlattening)) * tan(φ1));
        }

        /**
         * Compute the reduce latitude β from σ
         *
         * Eq. 13
         *
         * @param σ
         * @return β
         */
        private double computeβ(final double σ) {
            if (isInvalid(EQUATORIAL_AZIMUTH)) {
                throw new IllegalStateException("Uncomputed auxiliary parameter α0");
            }
            return atan2(cos(α0) * sin(σ), hypot(cos(α0) * cos(σ), sin(α0)));
        }

        /**
         * Compute the latitude α2 from σ2 and extract the
         *
         * Eq. 14 adapted for ending point.
         *
         * @param σ
         * @return α
         */
        private void computeEndPointAzimuth() {
            if (isInvalid(EQUATORIAL_AZIMUTH | AUXILIARY_END_POINT)) {
                throw new IllegalStateException("Unresolved auxiliary sphere's parameters for end point of the geodesic.");
            }
            final double azimuth = atan2(sin(α0), cos(α0) * cos(σ2));

            dφ2 = cos(azimuth);                                 // sin(π/2 − θ)  =  cos(θ)
            dλ2 = sin(azimuth);                                 // cos(π/2 − θ)  =  sin(θ)
            validity |= ENDING_AZIMUTH;
        }

        /**
         * Compute the latitude α2 from σ2 and extract the
         *
         * Eq. 14 adapted for ending point.
         *
         * @param σ
         * @return α
         */
        private void computeStartPointAzimuth() {
            if (isInvalid(EQUATORIAL_AZIMUTH | AUXILIARY_START_POINT)) {
                throw new IllegalStateException("Unresolved auxiliary sphere's parameters for end point of the geodesic.");
            }
            final double azimuth = atan2(sin(α0), cos(α0) * cos(σ1));

            dφ1 = cos(azimuth);       // sin(π/2 − θ)  =  cos(θ)
            dλ1 = sin(azimuth);       // cos(π/2 − θ)  =  sin(θ)
            validity |= STARTING_AZIMUTH;
        }
    }//end of internal class
    //==========================================================================

    public double getThirdFlattening() {
        return THIRD_FLATTENING;
    }

    public double getSecondEccentricity() {
        return SECOND_ECCENTRICITY;
    }

}
