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

    /**
     * HashMap containing the computed powers of the expansion parameter ε 
     */
    private double ε;
//    final private HashMap<Byte,Double> εPowersMap = new HashMap<>();
            
    private static final int DIRECT_AUXILIARY_PARAMETERS = 64, EXPAND_PARAMETER=128;

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
        ArgumentChecks.ensureFinite("inverse flattening of the ellipsoid", ellipsoid.getInverseFlattening());
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
        validity |= auxiliarySphere.parametersFromDirectProblem(ellipsoid.getInverseFlattening(), φ1, α1); // compute β1, α0, σ1, ω1
        
        //------------------------A DEPLACER DANS UNE M2THODE TIERCE--------------------------------
        final double s1 = spherical2EllipsoidalArcLength(auxiliarySphere.α0, auxiliarySphere.σ1); //compute s1

        
        //----------------------------------------------------------------------
        
        // λ01 is the longitude angle from the equatorial point E
        // to the staring point. It matchs λ1  of the C.F.F. Karney, 2013 article  
        final double λ01 = spherical2EllipsoidalLongitude(auxiliarySphere.α0, auxiliarySphere.ω1, auxiliarySphere.σ1);

        //λ02 -> Δλ = λ12 = λ02-λ01 
    }

    private void setExpandParameter(){
        if(isInvalid(DIRECT_AUXILIARY_PARAMETERS)){
            throw new IllegalStateException("Auxiliary parameter α0 needed to compute the expand parameter ε");
        }
        
        ε = ((sqrt(1 + pow(secondEccentriciy() * cos(auxiliarySphere.α0), 2))) - 1) / 
                ((sqrt(1 + pow(secondEccentriciy() * cos(auxiliarySphere.α0), 2))) + 1);
    }
    
    //COMMENT TO DELETE : Classe possible pour ne pas recalculer à chaque fois les puissances d'epsi
//    private double εPowers(byte power){
//        
//        if(εPowersMap.containsKey(power)){
//            return εPowersMap.get(power);
//        }
//       
//        if(power==1){
//            εPowersMap.put(power, 
//                    ((sqrt(1 + pow(secondEccentriciy() * cos(auxiliarySphere.α0), 2))) - 1) / 
//                    ((sqrt(1 + pow(secondEccentriciy() * cos(auxiliarySphere.α0), 2))) + 1));
//        } else {
//            εPowersMap.put(power, 
//                    pow(εPowers(Byte.valueOf(1),power));
//        }
//            
//    }
    
    /**
     * Compute s value from σ and α0.
     *
     * @param αo azimuth in the forward direction at the equator.
     * @param σ spherical arc length on the auxiliary sphere between the
     * equatorial point which forward direction following αo azimuth cross the
     * considered point. (σ1 for starting point, σ2 for ending point).
     *
     * @return the estimation of the associated ellipsoidal arc length s in
     * C.F.F Karney 's article.
     */
    private double spherical2EllipsoidalArcLength(final double αo, final double σ) {

        //Expansion parameter needed to power 6:
        
//            ε = ((sqrt(1 + pow(secondEccentriciy() * cos(αo), 2))) - 1)
//                    / ((sqrt(1 + pow(secondEccentriciy() * cos(αo), 2))) + 1);

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
        // Equation 7 together with equation 15 using term_Cil(....)
        // The same terms than in Eq.18 are considered.
        return (term_Cil_sin2l((-7 / 2048), ε, 6, 6, σ) + term_Cil_sin2l((-9 / 2048), ε, 6, 2, σ) + term_Cil_sin2l((3 / 512), ε, 6, 4, σ) //Termes ε ^ 6
                + term_Cil_sin2l((-7 / 1280), ε, 5, 5, σ) + term_Cil_sin2l((3 / 256), ε, 5, 3, σ) + term_Cil_sin2l((-1 / 32), ε, 5, 1, σ)//Termes ε ^ 5
                + term_Cil_sin2l((-5 / 512), ε, 4, 4, σ) + term_Cil_sin2l((1 / 32), ε, 4, 2, σ) //Termes ε ^ 4
                + term_Cil_sin2l((-1 / 48), ε, 3, 3, σ) + term_Cil_sin2l((3 / 16), ε, 3, 1, σ) //Termes ε ^ 3
                + term_Cil_sin2l((-1 / 16), ε, 2, 2, σ) //Termes ε ^ 2
                + term_Cil_sin2l((-1 / 2), ε, 1, 1, σ) //Termes ε ^ 1
                + σ)
                * A1 * ellipsoid.getSemiMinorAxis();

    }

    /**
     * Utility method to express the terms of C1l coefficient of equation 15
     * (Used to compute ellipsoidal ArcLength {
     *
     * @see #spherical2EllipsoidalArcLength}) in
     * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * Algortihms of geodesics from Charles F. F. Karney (SRI International)</a>
     *
     * computed_term = coef * ε ^ exposant * sin(2 * l * σ)
     *
     * @param coef real coefficient 
     * @param ε expansion parameter
     * @param exposant of the expansion parameter for the considered term.
     * @param l number of the fourier serie's term the current term is computed
     * for.
     * @param σ spheric arc length of the computation
     *
     * @return value of the coefficient of the form 
     * computed_term = coef * ε ^ exposant * sin(2 * l * σ)
     */
    private double term_Cil_sin2l(final double coef, final double ε, final double exposant, final double l, final double σ) {
        return coef * pow(ε, exposant) * sin(2 * l * σ);
    }

    /**
     * Compute the longitude angle λ from the equatorial point E and the
     * considered point on the ellipsoid from the associated with ω longitude on
     * the auxiliary sphere.
     *
     * {
     *
     * @see Eq. 8 together with Eq. 23 in C.F.F karney 2013}
     *
     * @param αo azimuth in the forward direction at the equator.
     * @param ω longitude on the auxiliary sphere.
     * @return λ longitude angle from the equatorial point E to the considered
     * point
     */
    private double spherical2EllipsoidalLongitude(final double αo, final double ω, final double σ) {

        
        final double n = thirdFlattening();
        
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
        final double sumC3l = term_Cil_sin2l(0.008203125, ε, 5, 5, σ) + term_Cil_sin2l(0.013671875, ε, 5, 4, σ)+ term_Cil_sin2l(0.013671875, ε, 5, 3, σ) + term_Cil_sin2l(0.01953125, ε, 5, 2, σ)+term_Cil_sin2l(0.0234375, ε, 5, 1, σ)  // ε⁵ terms
                + term_Cil_sin2l((7-14*n)/512, ε, 4, 4, σ) + term_Cil_sin2l((3/128 - 5*n/192), ε, 4, 3, σ) + term_Cil_sin2l((3+n)/128, ε, 4, 2, σ) +  term_Cil_sin2l((5+2*n)/128, ε, 4 , 1, σ)  // ε⁴ terms
                + term_Cil_sin2l((5-9*n+5*n*n)/192, ε, 3, 3, σ) + term_Cil_sin2l((3-2*n-3*n*n)/64, ε, 3, 2, σ) + term_Cil_sin2l((3+3*n-n*n)/64, ε, 3, 1, σ)  // ε³ terms
                + term_Cil_sin2l((2-3*n+n*n)/32, ε, 2, 2, σ) + term_Cil_sin2l((1-n*n)/8, ε, 2, 1, σ)   // ε² terms
                + term_Cil_sin2l((1-n)/4, ε, 1, 1, σ)     //ε term
                + σ;
        
        return ω - (1 / ellipsoid.getInverseFlattening()) * sin(αo) * A3 * sumC3l;
    }

    /**
     * Return the A3 coefficient necessary to converte the spherical longitude
     * in ellipsoidal longitude. {@link #spherical2EllipsoidalLongitude(double, double)} 
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

    /**
     * Compute the second eccentricity of the ellipsoid.
     *
     * @return e' = (a²-b²) / b²
     */
    private double secondEccentriciy() {
        return sqrt(pow(ellipsoid.getSemiMajorAxis(), 2) - pow(ellipsoid.getSemiMinorAxis(), 2)) / ellipsoid.getSemiMinorAxis();
    }

    /**
     * Compute the third flattening of the ellipsoid.
     *
     * @return (a-b) / (a+b)
     */
    private double thirdFlattening() {
        return (ellipsoid.getSemiMajorAxis() - ellipsoid.getSemiMinorAxis()) / (ellipsoid.getSemiMajorAxis() + ellipsoid.getSemiMinorAxis());
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
         * @param inverseFlattening inverse of the ellipsoidal's flattening
         * @param φ starting point's latitude
         * @param azimuth starting point's azimuth
         */
        final int parametersFromDirectProblem(final double inverseFlattening, final double φ, final double azimuth) {
//            ArgumentChecks.ensureFinite("inverseFlattening", inverseFlattening);    //placed in the constructor
            β1 = atan((1 - (1 / inverseFlattening)) * tan(φ));
            this.computeα0(β1, azimuth);
            this.σ1 = computeσ(β1);
            this.ω1 = computeω(σ1);
            return DIRECT_AUXILIARY_PARAMETERS;
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
