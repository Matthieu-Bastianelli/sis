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
import org.apache.sis.internal.referencing.Resources;
import org.apache.sis.util.ArgumentChecks;

/**
 *Classe extending GeodeticCalculator and which aims to implement the direct and
 * inverse the direct and inverse geodesic problem solving methods described 
 * in {@see
 * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
 * Algortihms of geodesics from Charles F F Karney (SRI International)</a>}.
 *
 * All angles are in radians.
 * 
 * As in the cited article, we note E the equatorial point at which the geodesic
 * crosses the equator in the northward direction.
 *
 * @version 1.0
 * @since 1.0
 * @module
 */
final class GeodesicsOnEllipsoid extends GeodeticCalculator {

    private final double squared_ECCENTRICITY, second_ECCENTRICITY;
    private final double N; // Third flattening of the ellipsoid
    private double ε; // expand parameter used to resolve geodesic problems

    private static final int EQUATORIAL_AZIMUTH = 64, AUXILIARY_START_POINT = 128, AUXILIARY_END_POINT = 256, EXPAND_PARAMETER = 512,
            REDUCED_LATITUDES=1024;

    final AuxiliarySpheresParameters auxiliarySpheres;
    
    /**
     * Error tolerance for λ12 computation during the inverse problem solving
     * In radians.
    */
    private double λ12_tolerance = 5*10e-7; //Less than 0.1°
    /**
     * Maximum number of iterations allowed to approximate the starting point azimuth 
     * when solving the inverse problem.
    */
    private double maxIterationNumber = 20; 

    
    /**
     * Limit for λ variation to consider the starting and ending points as
     * nearly antipodal.
     */
    private double nearly_antipodal_λ12 = PI*(1-1/360);
            

    /**
     * Construct a new geodetic calculator for ellipsoid expecting coordinates
     * in the supplis CRS.
     *
     * The {@link GeodeticCalculator} constructor is called.
     *
     * @param crs the referencing system for the {@link Position} arguments and
     * return values.
     */
    public GeodesicsOnEllipsoid(CoordinateReferenceSystem crs) {
        super(crs);
        ArgumentChecks.ensureFinite("Inverse flattening of the ellipsoid", ellipsoid.getInverseFlattening());
        ArgumentChecks.ensureStrictlyPositive("Ellipsoid semi minor axis", ellipsoid.getSemiMinorAxis());

        squared_ECCENTRICITY = (2 - (1 / ellipsoid.getInverseFlattening())) / ellipsoid.getInverseFlattening();
        // Compute the second eccentricity of the ellipsoid. e' = (a²-b²) / b²
        second_ECCENTRICITY = sqrt(ellipsoid.getSemiMajorAxis()*ellipsoid.getSemiMajorAxis() - ellipsoid.getSemiMinorAxis()*ellipsoid.getSemiMinorAxis()) / ellipsoid.getSemiMinorAxis();

        // Compute the third flattening of the ellipsoid. n = (a-b) / (a+b)
        N = (ellipsoid.getSemiMajorAxis() - ellipsoid.getSemiMinorAxis()) / (ellipsoid.getSemiMajorAxis() + ellipsoid.getSemiMinorAxis());

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

        //Given φ1,  and geodesicDistance (s12 in C.F.F Karney's article) :
        
//        final double α1 = atan2();     // tan(π/2 − θ)  =  1/tan(θ)

        auxiliarySpheres.parametersFromDirectProblem(ellipsoid.getInverseFlattening(), φ1); // compute β1, α0, σ1, ω1

        final double s2 = spherical2EllipsoidalArcLength(auxiliarySpheres.σ1) + this.geodesicDistance; //s2 = s1 + s12
        auxiliarySpheres.resolveEndPoint(s2);

        // λ01 is the longitude angle from the equatorial point E
        // to the staring point. It matchs λ1  of the C.F.F. Karney, 2013 article  
//        final double λ01 = spherical2EllipsoidalLongitude(auxiliarySpheres.ω1, auxiliarySpheres.σ1);
        final double Δλ = spherical2EllipsoidalLongitude(auxiliarySpheres.ω2, auxiliarySpheres.σ2) 
                - spherical2EllipsoidalLongitude(auxiliarySpheres.ω1, auxiliarySpheres.σ1); //Δλ <-> λ12

        this.φ2 = reduced2latitude(auxiliarySpheres.sin_β2, auxiliarySpheres.cos_β2);
        this.λ2 = this.λ1 + Δλ; //λ2 =  λ1 + Δλ 
        auxiliarySpheres.computeEndPointAzimuth();  // affects dφ2 and dλ2

        validity |= END_POINT;
    }

    private void setExpandParameter() {
        if (isInvalid(EQUATORIAL_AZIMUTH)) {
            throw new IllegalStateException("Auxiliary parameter α0 needed to compute the expand parameter ε");
        }

//        ε = (sqrt(1 + SECOND_ECCENTRICITY * auxiliarySpheres.cos_α0 * SECOND_ECCENTRICITY * auxiliarySpheres.cos_α0) - 1)
//                / ((sqrt(1 + SECOND_ECCENTRICITY * auxiliarySpheres.cos_α0*SECOND_ECCENTRICITY * auxiliarySpheres.cos_α0)) + 1);
//        
        final double hypot_with_e2nd_cosa0 = hypot(1, second_ECCENTRICITY * auxiliarySpheres.cos_α0);

        ε = (hypot_with_e2nd_cosa0 - 1) / (hypot_with_e2nd_cosa0 + 1);
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
        final double τ = s / (computeA1Coefficient()*ellipsoid.getSemiMinorAxis()); //b * A1 can't be equal to 0 as b and A1 >0 //DOUBLON DE CALCUL BA1 AVEC spherical2EllipsoidalArcLength(....)!!!!

//========================== C'_1i (τ) terms ===================================
//     double C_11 = (((205. / 1536) * ε * ε + (-9. / 32) ) * ε * ε + (1. / 2) ) * ε   * sin(2 * τ);       //C'11 (τ)
//     double C_12 = (((1335. /4096 )* ε * ε + (-37. /96))* ε * ε + (5. / 16)) * ε * ε * sin(2 * 2 * τ);   //C'12 (τ)
//     double C_13 = ((-75. / 128)* ε * ε + (29. / 96))* ε * ε * ε                     * sin(2 * 3 * τ);   //C'13 (τ)
//     double C_14 = ((-2391. /2560)* ε * ε + (539. / 1536))* ε * ε * ε * ε            * sin(2 * 4 * τ);   //C'14 (τ)
//     double C_15 = (3467. / 7680) * ε * ε * ε  * ε * ε                               * sin(2 * 5 * τ);  // C'15 (τ)
//     double C_16 = (38081. / 61440) * ε * ε * ε  * ε * ε  * ε                        * sin(2 * 6 * τ);  // C'16 (τ)
//==============================================================================

        final double sin_2_4_τ = sin(2 * 4 * τ);
        final double sin_2_3_τ = sin(2 * 3 * τ);
        final double sin_2_2_τ = sin(2 * 2 * τ); 
        final double sin_2_τ = sin(2 * τ);  

        // Equation 20 ; The above terms from Eq.21 are considered.
        // τ + Sum( C'_1i (τ) ) with i={1 , ... , 6}
        return   (((((( 
                  (38081. / 61440) * sin(2 * 6 * τ)  + (-2391. /2560) * sin_2_4_τ  + (1335. /4096 ) * sin_2_2_τ ) * ε  //Terms ε ^ 6
               +  (3467. / 7680)   * sin(2 * 5 * τ)  + (-75. / 128)   * sin_2_3_τ  + (205. / 1536)  * sin_2_τ   ) * ε  //Terms ε ^ 5
               +  (539. / 1536)    * sin_2_4_τ       + (-37. /96)     * sin_2_2_τ                               ) * ε //Terms ε ^ 4
               +  (29. / 96)       * sin_2_3_τ       + (-9. / 32)     * sin_2_τ                                 ) * ε //Terms ε ^ 3
               +  (5. / 16)        * sin_2_2_τ                                                                  ) * ε //Terms ε ^ 2
               +  (1. / 2)         * sin_2_τ                                                                    ) * ε //Terms ε
               +   τ ;
    }

    /**
     * Compute s value from σ and α0.
     * 
     * s = b*I1(σ) 
     *
     * @param σ spherical arc length on the auxiliary sphere between the
     * equatorial point which forward direction following αo azimuth cross the
     * considered point. (σ1 for starting point, σ2 for ending point).
     *
     * @return s : the estimation of the associated ellipsoidal arc length
     */
    double spherical2EllipsoidalArcLength(final double σ) {
        return compute_I1(σ) * ellipsoid.getSemiMinorAxis();
    }
    
    /**
     * See {@link #spherical2EllipsoidalArcLength(double)} :
     * Compute the I1(σ) coefficient needed for spherical to ellipsoidal 
     * Arc Length conversion.
     * 
     * @param σ spherical arc length on the auxiliary sphere.
     * @return I1(σ) coefficient according to Eq 15 in 
     *  <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * C.F.F Karney 's article.</a>
     */
    private double compute_I1(final double σ) {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        final double A1 = computeA1Coefficient();

//=========================== C_1i (σ) terms ===================================
        //Equation 7 together with equation 15 :
//        return ((  (-7. / 2048) * pow(ε, 6)                                                    )* sin(2 * 6 * σ) //Term with C16
//                +( (-7. / 1280) * pow(ε, 5)                                                    )* sin(2 * 5 * σ) //Term with C15
//                +( ((3. / 512)   * pow(ε, 2) - (5. / 512)) * pow(ε, 4)                         )* sin(2 * 4 * σ) //Term with C14
//                +( ((3. / 256)   * pow(ε, 2) - (1. / 48))  * pow(ε, 3)                         )* sin(2 * 3 * σ) //Term with C13
//                +( (((-9. / 2048) * pow(ε, 2) + (1. / 32))  * pow(ε, 2) - (1 / 16)) * pow(ε, 2) * sin(2 * 2 * σ) //Term with C12
//                +( (((-1. / 32)   * pow(ε, 2) + (3. / 16))  * pow(ε, 2) - (1 / 2) ) * ε           ) * sin(2 * 1 * σ) //Term with C11
//                + σ) * bA1;
//==============================================================================

        final double sin_2_4_σ = sin(2 * 4 * σ);
        final double sin_2_3_σ = sin(2 * 3 * σ);
        final double sin_2_2_σ = sin(2 * 2 * σ); 
        final double sin_2_σ = sin(2 * σ); 
        // Equation 7 together with equation 15 
        // The same terms than in Eq.17 are considered.
        // σ + Sum( C_1i (σ) ) with i={1 , ... , 6}
        return (((((((   
                  (-7. / 2048) * sin(2 * 6 * σ) + (-9. / 2048) * sin_2_2_σ  + (3. / 512) * sin_2_4_σ )*ε  //Terms ε ^ 6
                + (-7. / 1280) * sin(2 * 5 * σ) + (3. / 256)   * sin_2_3_σ  + (-1. / 32) * sin_2_σ   )*ε  //Terms ε ^ 5
                + (-5. / 512)  * sin_2_4_σ      + (1. / 32)    * sin_2_2_σ                           )*ε  //Terms ε ^ 4
                + (-1. / 48)   * sin_2_3_σ      + (3. / 16)    * sin_2_σ                             )*ε  //Terms ε ^ 3
                + (-1. / 16)   * sin_2_2_σ                                                           )*ε  //Terms ε ^ 2
                + (-1. / 2)    * sin_2_σ                                                             )*ε  //Terms ε
                +  σ)
                * A1;
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
    double spherical2EllipsoidalLongitude(final double ω, final double σ) {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }

        // I3(σ) = A3 * (σ + sum)
        final double A3 = computeA3Coefficient();

//============================ C_3i (σ) terms ==================================
//     double C_31 = (((((3. / 128) * ε + ((5 + 2 * n) / 128.0) )* ε + ((3 + 3 * n - n * n) / 64.) ) * ε +  ((1 - n * n) / 8.) )* ε + ((1 - n) / 4.) ) * ε  * sin(2 * σ);       //C31 (σ)
//     double C_32 = ((((5. / 256)* ε + ((3 + n) / 128.0))* ε +((3 - 2 * n - 3 * n * n) / 64.))* ε + ((2 - 3 * n + n * n) / 32.)) * ε * ε                   * sin(2 * 2 * σ);   //C32 (σ)
//     double C_33 = (((7. / 512)* ε + (3. / 128 - 5 * n / 192.0) )* ε + ((5 - 9 * n + 5 * n * n) / 192.))* ε * ε * ε                                       * sin(2 * 3 * σ);   //C33 (σ)
//     double C_34 = ((7. / 512) * ε + ((7 - 14 * n) / 512.))* ε * ε * ε * ε                                                                                * sin(2 * 4 * σ);   //C34 (σ)
//     double C_35 = (21. / 2560) * ε * ε * ε * ε * ε                                                                                                       * sin(2 * 5 * σ);  //C35 (σ)
//==============================================================================

        final double sin_2_4_σ = sin(2 * 4 * σ);          //Already computed in #spherical2EllipsoidalArcLength
        final double sin_2_3_σ = sin(2 * 3 * σ);
        final double sin_2_2_σ = sin(2 * 2 * σ); 
        final double sin_2_σ = sin(2 * σ); 

        // Computation of the sum in Equation 23 of the C.F.F Karney's article.
        // The same terms than in Eq.25 are considered.
        // Sum( C_3i (σ) ) with i={1 , ... , 5}
        final double sumC3l = ((((
                   (21. / 2560)                     * sin(2 * 5 * σ) + (7. / 512)                 * sin_2_4_σ + (7. / 512)                  * sin_2_3_σ + (5. / 256)            * sin_2_2_σ + (3. / 128) * sin_2_σ  ) * ε //Terms ε ^ 5
                 + ((7 - 14 * N) / 512.)            * sin_2_4_σ + ((3. / 128 - 5 * N / 192.0)     * sin_2_3_σ + ((3 + N) / 128.0)           * sin_2_2_σ + ((5 + 2 * N) / 128.0) * sin_2_σ  ) * ε //Terms ε ^ 4
                 + ((5 - 9 * N + 5 * N * N) / 192.) * sin_2_3_σ + ((3 - 2 * N - 3 * N * N) / 64.) * sin_2_2_σ + ((3 + 3 * N - N * N) / 64.) * sin_2_σ      ) * ε //Terms ε ^ 3
                 + ((2 - 3 * N + N * N) / 32.)      * sin_2_2_σ + ((1 - N * N) / 8.)              * sin_2_σ     ) * ε //Terms ε ^ 2
                 + ((1 - N) / 4.)                   * sin_2_σ    ) * ε //Terms ε
                 +  σ;                        

        return ω - auxiliarySpheres.sin_α0 * A3 * sumC3l / ellipsoid.getInverseFlattening();
    }

    /**
     * return the latitude on the ellipsoid from the reduced latitude on the
     * associated auxiliary sphere.
     *
     * @param β the reduced latitude on the auxiliary sphere
     * @return φ the latitude on the ellipsoid.
     */
    private double reduced2latitude(final double sin_β, final double cos_β) {
        ArgumentChecks.ensureStrictlyPositive("Semi Minor Axis", ellipsoid.getSemiMinorAxis());

//        return atan(tan(β) / (1 - 1 / ellipsoid.getInverseFlattening())); // flattening != 1 as semi minor axis b >0 from constructor.
        return atan((sin_β/cos_β) / (1 - 1 / ellipsoid.getInverseFlattening())); // flattening != 1 as semi minor axis b >0 from constructor.
    }

    /**
     * See {@link #spherical2EllipsoidalArcLength(double)} :
     * Compute the I1(σ) coefficient needed for spherical to ellipsoidal 
     * Arc Length conversion.
     * 
     * @param σ spherical arc length on the auxiliary sphere.
     * @return I2(σ) coefficient according to Eq 41 together with Eq42 and 43 in 
     *  <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
     * C.F.F Karney 's article.</a>
     */ 
    private double compute_I2(final double σ) {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        final double A2 = computeA2Coefficient();

        final double sin_2_4_σ = sin(2 * 4 * σ); //Doublon de calculs possible avec calculs de I1
        final double sin_2_3_σ = sin(2 * 3 * σ);
        final double sin_2_2_σ = sin(2 * 2 * σ); 
        final double sin_2_σ = sin(2 * σ); 
        
        // Equation 41 together with equation 43 
        // The same terms than in Eq.17 are considered.
        // σ + Sum( C_2i (σ) ) with i={1 , ... , 6}
        return (((((((   
                  (77. / 2048) * sin(2 * 6 * σ) + (35. / 2048) * sin_2_2_σ  + (7. / 512) * sin_2_4_σ )*ε  //Terms ε ^ 6
                + (63. / 1280) * sin(2 * 5 * σ) + (5. / 256)   * sin_2_3_σ  + (1. / 32) * sin_2_σ    )*ε  //Terms ε ^ 5
                + (35. / 512)  * sin_2_4_σ      + (1. / 32)    * sin_2_2_σ                           )*ε  //Terms ε ^ 4
                + (5. / 48)   * sin_2_3_σ      + (1. / 16)    * sin_2_σ                              )*ε  //Terms ε ^ 3
                + (3. / 16)   * sin_2_2_σ                                                            )*ε  //Terms ε ^ 2
                + (1. / 2)    * sin_2_σ                                                              )*ε  //Terms ε
                +  σ)
                * A2;
    }
    
    private double compute_J(final double σ) {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        final double A2 = computeA2Coefficient();
        final double A1 = computeA1Coefficient();

        final double sin_2_6_σ = sin(2 * 6 * σ);
        final double sin_2_5_σ = sin(2 * 5 * σ);
        final double sin_2_4_σ = sin(2 * 4 * σ); //Doublon de calculs possible avec calculs de I1
        final double sin_2_3_σ = sin(2 * 3 * σ);
        final double sin_2_2_σ = sin(2 * 2 * σ); 
        final double sin_2_σ = sin(2 * σ); 
        
        // Equation 41 together with equation 43 
        // The same terms than in Eq.17 are considered.
        // σ + Sum( C_2i (σ) ) with i={1 , ... , 6}
        return (((((((   
                  (((-7. / 2048) * sin_2_6_σ + (-9. / 2048) * sin_2_2_σ  + (3. / 512) * sin_2_4_σ))*A1 - ((77. / 2048) * sin_2_6_σ + (35. / 2048) * sin_2_2_σ  + (7. / 512) * sin_2_4_σ)*A2       )*ε  //Terms ε ^ 6
                +((-7. / 1280)   * sin_2_5_σ + (3. / 256)   * sin_2_3_σ  + (-1. / 32) * sin_2_σ)*A1    - ((63. / 1280) * sin_2_5_σ + (5. / 256)   * sin_2_3_σ  + (1. / 32) * sin_2_σ)*A2          )*ε  //Terms ε ^ 5
                +((-5. / 512)    * sin_2_4_σ      + (1. / 32)    * sin_2_2_σ)*A1                       - ((35. / 512)  * sin_2_4_σ      + (1. / 32)    * sin_2_2_σ)*A2                            )*ε  //Terms ε ^ 4
                +((-1. / 48)     * sin_2_3_σ      + (3. / 16)    * sin_2_σ)*A1                         - ((5. / 48)    * sin_2_3_σ      + (1. / 16)    * sin_2_σ )*A2                             )*ε  //Terms ε ^ 3
                +((-1. / 16)     * sin_2_2_σ)*A1                                                       - ((3. / 16)    * sin_2_2_σ)*A2                                                            )*ε  //Terms ε ^ 2
                +((-1. / 2)      * sin_2_σ)*A1                                                         - ((1. / 2)     * sin_2_σ) *A2                                                             )*ε  //Terms ε
                + σ*A1                                                                                 -  σ*A2);
        
        
    }
    
    /**
     * Return the A1 coefficient necessary to process arc length convertions
     * between ellipsoid and the auxiliary spheres.
     * {@link #ellipsoidal2SphericalArcLength(double)} , 
     * {@link #spherical2EllipsoidalArcLength(double) } 
     * and {@link #compute_I2(double)}
     *
     * @param n : third flattening of the ellipsoid
     * @return A1 value from terms suggest in C.F.F. Karney, 2013.
     */
    private double computeA1Coefficient() {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        if (ε == 1) {
            throw new IllegalStateException("Impossible to compute BA1Coefficient for 1 value of ε expand parameter.");
        }
        return (  ((((1. / 256) *ε*ε  //Terme ε ^ 6
                + (1. / 64))    *ε*ε  //Terme ε ^ 4
                + (1. / 4))     *ε*ε  //Terme ε ^ 2
                + 1)                
                / (1 - ε) ); //A1
    }

    /**
     * Return the A2 coefficient necessary to compute I2 coefficient.
     * {@link #compute_I2(double)}
     *
     * @return A2 value from terms suggest in C.F.F. Karney, 2013, Eq 42.
     */
    private double computeA2Coefficient() {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }
        if (ε == 1) {
            throw new IllegalStateException("Impossible to compute BA2Coefficient for 1 value of ε expand parameter.");
        }
        return (  ((((25. / 256) *ε*ε  //Terme ε ^ 6
                + (9. / 64))    *ε*ε  //Terme ε ^ 4
                + (1. / 4))     *ε*ε  //Terme ε ^ 2
                + 1)                
                / (1 - ε) ); //A1
    }
    
    /**
     * Return the A3 coefficient necessary to convert the spherical longitude in
     * ellipsoidal longitude.
     * {@link #spherical2EllipsoidalLongitude(double, double)}
     *
     * @param n : third flattening of the ellipsoid
     * @return A3 value from terms suggest in C.F.F. Karney, 2013.
     */
    private double computeA3Coefficient() {
        if (isInvalid(EXPAND_PARAMETER)) {
            setExpandParameter();
        }

        
        return (((( (-3. / 128)                         * ε
                -   ((3. / 64) + (N / 32.))           ) * ε
                -   (1. / 16) * (1. + 3. * N + N * N) ) * ε
                -   (1. / 8) * (2. + N + 3. * N * N)  ) * ε
                -   (1. / 2) * (1. - N)               ) * ε
                +   1.;
    }

    public double getThirdFlattening() {
        return N;
    }

    public double getSecondEccentricity() {
        return second_ECCENTRICITY;
    }
    
    
    /**
     * Computes the geodesic distance and azimuths from the start point and end point.
     * This method should be invoked if the distance or an azimuth is requested while
     * {@link #STARTING_AZIMUTH}, {@link #ENDING_AZIMUTH} or {@link #GEODESIC_DISTANCE}
     * validity flag is not set.
     *
     *The problem solving considers the following canonical configuration :
     * Canonical configuration for problem solving :
     * 0 >= φ1 ; -φ1 >= φ2 >= φ1 ; PI >= λ12 (= λ2 - λ1) >= 0
     * Thus the computed distance is always the shortest distance. Distances for
     * values of λ12 in ]PI , 2PI[ aren't computed.
     * 
     * 
     * @throws IllegalStateException if the distance or azimuth has not been set.
     */
    @Override
    void computeDistance() {
        if (isInvalid(START_POINT | END_POINT)) {
            throw new IllegalStateException(Resources.format(
                    Resources.Keys.StartOrEndPointNotSet_1, Integer.signum(validity & START_POINT)));
        }
        
            //************
            //Not nearly antipodial points -> What's considered so?
            //************
        
        validity &= ~(EQUATORIAL_AZIMUTH | AUXILIARY_START_POINT | AUXILIARY_END_POINT | EXPAND_PARAMETER);
        
        
        // Variables used to indicate if coordinate changes were performed to reach 
        // The considered condition for the inverse problem solving.
        boolean inversionIndex = false;
        boolean oppositeIndex = false;
        boolean translatedEndPointIndex = false;
        
        final double λ12 = λ2 - λ1;
        double α1Approximation;
//        double dφ1_approximation;
//        double dλ1_approximation;

        //Particular cases :
        if ((λ12 == 0) || (λ12 == PI)) {
//            dφ1 = cos(λ12);
//            dλ1 = sin(λ12);
            α1Approximation = λ12;
        } else if ((φ1 == 0) && (φ2 == 0) && (λ12 <= (1 - 1 / ellipsoid.getInverseFlattening()) * PI)) {
//            double α1 = PI / 2;
//            dφ1 = cos(α1);
//            dλ1 = sin(α1);
            α1Approximation = PI / 2;
        } else {

            
            //------------------------------------------------------------------
            // Get the canonical configuration for problem solving :
            //------------------------------------------------------------------
            // φ1 <= 0 ; φ1 <= φ2 <= -φ1 ; 0 <= λ12 (= λ2 - λ1) <= PI
            //-φ1 >= φ2 >= φ1
            if (abs(φ2) > abs(φ1)) {
                startAndEndInversion();
                inversionIndex = true;
            }
            // φ1 <= 0
            if (φ1 > 0) {
                φ1 = -φ1;
                φ2 = -φ2;
                oppositeIndex = true;
            }
            // 0 <= λ12 <= PI
//        double λ12 = λ2 - λ1;
            if ((λ12 < 0) || (λ12 > PI)) {

                // Method used to reach PI >= λ12 (= λ2 - λ1) >= 0 recquired configuration.
                // As consequence of revolution of the ellipsoid this modification of the
                // end point preserves geodesic distance.
                λ2 = λ1 + abs(λ12);
                translatedEndPointIndex = true;
            }
            
            //------------------------------------------------------------------
            // At this state, it si considered :
            //----------------------------------
            // α1 in [0,PI] and λ12 in [0 , PI]
            // α2 restricted in [0,PI/2] as it correspond to the first intersection 
            // of the geodesic with the φ2 latitude 
            // Compute the  associated with the problem
            
            if (λ12>nearly_antipodal_λ12) {
                α1Approximation = startingAzimuthNearlyAntipodal();       //not valid yet
            } else {
                α1Approximation = startingAzimuthFirstApproximation(); //reduced latitudes are computed
            }
        
        }
        
        double λ12_approx = Double.MIN_VALUE;
        double derivative_λ12_on_dα1_approx = Double.MIN_VALUE;

        
        int count=0;

        // While loop used to improve the azimtuh approximation 
        while (λ12 - λ12_approx > λ12_tolerance) {
            if(count > maxIterationNumber){
                throw new ArithmeticException("Maximum iteration reached when trying to solve the inverse geodesic problem. Try to raise the λ12_tolerance value OR the maxIterationNumber.");
            } else if (count > 0) { // this bloc isn't executed the first time it's reached.

                if ((abs(dφ1) == 0)
                        && (auxiliarySpheres.cos_β1 == auxiliarySpheres.cos_β2)
                        && (auxiliarySpheres.sin_β1 == auxiliarySpheres.sin_β2)) {

                    derivative_λ12_on_dα1_approx = signum(PI / 2 - α1Approximation)
                            + sqrt(1 - squared_ECCENTRICITY * auxiliarySpheres.cos_β1 * auxiliarySpheres.cos_β1) / auxiliarySpheres.sin_β1;

                } else {
                    derivative_λ12_on_dα1_approx = auxiliarySpheres.derivativeDλ12Dα1(); //λ12'(α1Approximation)
                }

                if (derivative_λ12_on_dα1_approx == 0) {
                    break; // derivative_λ12_on_dα1_approx = 0 means no improvement.
                }
                α1Approximation = α1Approximation - (λ12_approx - λ12) / derivative_λ12_on_dα1_approx;
            }
            
            //Hybrid geodesic problem solving ,  α1 is approximated, computes :
            //(σ1,σ2,ω1,ω2, associated sinus and cosinus) and end point 's azimuth (dλ2, dφ2)
            λ12_approx = auxiliarySpheres.hybridInverseProblemSolving(α1Approximation);  //λ12(α1Approximation)
            count++;
        }

        geodesicDistance = spherical2EllipsoidalArcLength(auxiliarySpheres.σ2) - spherical2EllipsoidalArcLength(auxiliarySpheres.σ1);

        //********************************************
        // Return to the initial configuration :
        //********************************************
        //Ellipsoid of revolution.
           
        if(translatedEndPointIndex){
                λ2 = λ1 - abs(λ12);
                dλ1 = - dλ1; // dφ1=dφ1;
                dλ2=  - dλ2;  //  dφ2=dφ2;
//                translatedEndPointIndex = false;
        }   
        if(oppositeIndex){
                φ1 = -φ1;
                φ2 = -φ2;
                                
                dφ1=-dφ1;   // dλ1 = dλ1; 
                dφ2=-dφ2;   // dλ2= dλ2;
//              oppositeIndex = false;
        }   
        if(inversionIndex){
                startAndEndInversion();
                final double dφ1_inter = dφ1;
                final double dλ1_inter = dλ1;
                dφ1 = - dφ2;
                dλ1 = - dλ2;
                dφ1 = - dφ1_inter;
                dλ1 = - dλ1_inter;
//                inversionIndex = false;
        }


    }
    
    
    
    /**
     * Return an approximation of the azimuth at the starting point.
     * 
     * According to C.F.F Karney's article, it is provide assuming a relationship
     * between the variation of longitude (λ12) on the ellipsoid and on the auxiliary 
     * sphere (ω12):
     *          ω12 = (λ2 - λ1) / ωOverbar;
     * 
     * And solving for the great circle on the auxiliary sphere.
     * 
     * @return Approximation of starting point's azimuth
     */
    double startingAzimuthFirstApproximation() {
        if (isInvalid(REDUCED_LATITUDES)) {
            auxiliarySpheres.computeReducedLatitudes();
//             throw new IllegalStateException("uncomputed βi reduced latitudes of starting and ending point.");
        }

        final double cos_β1 = auxiliarySpheres.cos_β1;
        final double sin_β1 = auxiliarySpheres.sin_β1;
        final double cos_β2 = auxiliarySpheres.cos_β2;
        final double sin_β2 = auxiliarySpheres.sin_β2;

//        final double e2 = (2 - (1 / ellipsoid.getInverseFlattening())) / ellipsoid.getInverseFlattening()  //e2 = f (2-f);
        double ωOverbar = sqrt(1 - (squared_ECCENTRICITY * (cos_β1 + cos_β2) * (cos_β1 + cos_β2) / 4));

        double ω12 = (λ2 - λ1) / ωOverbar; // λ12 = (λ2-λ1)

        final double sin_ω12 = sin(ω12);
        final double cos_ω12 = cos(ω12);

        return atan2(cos_β2 * sin_ω12, (cos_β1 * sin_β2 - sin_β1 * cos_β2 * cos_ω12));      // table3 example :77.043 533 542 96651 value found with  != 77.043 533 542 37 expected from the article 
//            double α2Approx = atan2(cos_β1*sin_ω12 , (- sin_β1*cos_β2 + cos_β1*sin_β2*cos_ω12) );     // table3 example :77.043 508 44 6 33922 value found with  != 77.043 508 44 9 13 expected from the article 
//            
//            double σ12 = atan2( hypot(cos_β1*sin_β2 - sin_β1*cos_β2*cos_ω12, cos_β2*sin_ω12 ), sin_β1*sin_β2 + cos_β1*cos_β2*cos_ω12 );

    }
    double startingAzimuthNearlyAntipodal() {
        if (isInvalid(REDUCED_LATITUDES)) {
            auxiliarySpheres.computeReducedLatitudes();
//             throw new IllegalStateException("uncomputed βi reduced latitudes of starting and ending point.");
        }
        
            throw new IllegalStateException("unfinished implementation.");

            //*********************************************************
            //TODO FINISH IMPLEMENTING COMPUTING μ and particular case.
            // then to choose the calling condition in {@link #computeDistance()} l.539
            //*********************************************************
//        final double cos_β1 = auxiliarySpheres.cos_β1;
//        final double β1 = acos(cos_β1);
//        final double β2 = acos(auxiliarySpheres.cos_β2);
//        final double Δ = cos_β1*PI/ellipsoid.getInverseFlattening(); //correspond to Δ/(cos_β1*a) in (CFF Karney,2013)
//
//        final double x = (dλ2-dλ1 - PI) / Δ; //Eq 53
//        final double y = (β2+β1) / (Δ*cos_β1); //Eq 53
//        final double μ = ; //Eq 55  Appendice B from CFF Karney 2011 give a method of resolution.
//        
//        if(y==0){
//            throw new IllegalStateException("undefined yet result for 0 value of y.");
////            return atan2(-x, signum(?)*sqrt(max(0.1-x*x))); //Eq57 limit y -> 0 du cas générale nécessite l'expression de μ pour bonne compréhension. Quel signe et que signifie le sqrt(max(0.1-x²))?
//        }else{
//            return atan2(-x/(1+μ), y/μ);//Eq56 
//        
//        }
            //*********************************************************

    }
    
    
    
    private void startAndEndInversion(){
        final double φ1_inter = φ1;
        final double λ1_inter = λ1;
        φ1 = φ2;
        φ2 = φ1_inter;
        λ1 = λ2;
        λ2 = λ1_inter;
    }

    /**
     * Set the error tolerance for λ12 computation during the inverse probleme
     * solving.
     * 
     * Degrees input are converted in radians.
     * 
     * @param λ12_tolerance in degrees
     */
    public void setλ12_tolerance(double λ12_tolerance) {
        this.λ12_tolerance = λ12_tolerance;
    }

    /**
     * Return the error tolerance for λ12 computation during the inverse probleme
     * solving.
     * 
     * @return λ12_tolerance in degrees
     */
    public double getλ12_tolerance() {
        return toDegrees(λ12_tolerance);
    }

    /**
     * Returns the maximum number of iterations allowed to approximate the
     * starting point azimuth when solving the inverse problem. 
     * 
     * @return maxIterationNumber
     */
    public double getMaxIterationNumber() {
        return maxIterationNumber;
    }

    /**
     * Defines the maximum number of iterations allowed to approximate the
     * starting point azimuth when solving the inverse problem. 
     * @param maxIterationNumber 
     */
    public void setMaxIterationNumber(double maxIterationNumber) {
        this.maxIterationNumber = maxIterationNumber;
    }

    /**
     * Defines nearly_antipodal_λ12, the limit for λ variation to consider the
     * starting and ending points as nearly antipodal.
     * 
     * @param nearly_antipodal_λ12 
     */
    public void setNearly_antipodal_λ12(double nearly_antipodal_λ12) {
        this.nearly_antipodal_λ12 = nearly_antipodal_λ12;
    }

    /**
     * Return the limit value for λ variation to consider the
     * starting and ending points as nearly antipodal.
     * 
     * @return nearly_antipodal_λ12
     */
    public double getNearly_antipodal_λ12() {
        return nearly_antipodal_λ12;
    }
    
        
        
        
    
     
             
             
    
    //==========================================================================
    //Internal Classe
    //==========================================================================
    /**
     * Internal classe used to handle parameters on used auxiliary spheres.
     *
     */
    class AuxiliarySpheresParameters {

        /**
         * Azimuth at the intercepting point of the equator (0° latitude) and
         * the direction given by an α azimuth at the β reduced latitude.
         * (In radians)
         */
        double α0;
        double sin_α0, cos_α0;
        /**
         * Reduced latitudes of the starting point (respectively the ending
         * point).
         * (In radians)
         */
//        double β1, β2;
        double sin_β1, cos_β1;
        double sin_β2, cos_β2;

        /**
         * Spherical arc length between the intercepting E point of
         * the equator (0° latitude) and the starting P point (resp the ending
         * point) along the 'direction' given by the α1 (resp α2) azimuth.
         * (In radians)
         */
        double σ1, σ2;
        double sin_σ1, cos_σ1, sin_σ2, cos_σ2;

        /**
         * Spherical longitudes on the auxiliary sphere of the starting and
         * ending point.
         *
         * The spherical longitude of a point on the auxiliary sphere is the
         * angle formed by the meridian of the intercepting point of the equator
         * (0° latitude) and the direction given by the point given azimuth.
         *(In radians)
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
        final void parametersFromDirectProblem(final double inverseFlattening, final double φ) {
            if (isInvalid(START_POINT | STARTING_AZIMUTH)) {
                throw new IllegalStateException(Resources.format(Resources.Keys.StartOrEndPointNotSet_1, 0));
            }
            if ((inverseFlattening == 0) | (inverseFlattening == 1)) {
                throw new IllegalStateException("Unvalide ellipsoid : 0 or 1 value of inverse flattening.");
            }
            double tan_β1 =  (1 - (1 / inverseFlattening)) * tan(φ1);
            
//            double β1 = atan(tan_β1);
                    
//            double sin_azimuth = sin(azimuth); -> dλ1
//            double cos_azimuth = cos(azimuth); -> dφ1
//            sin_β1 = sin(β1);
//            cos_β1 = cos(β1);
            cos_β1 = cos(atan(tan_β1));
            sin_β1 = tan_β1 * cos_β1;

            α0 = atan2(dλ1 * cos_β1, hypot(dφ1, dλ1 * sin_β1));
//            sin_α0 = sin(α0);
            sin_α0 = dλ1 * cos_β1; // as sin_α0 = sin_α1*cos_β1 =sin_α2*cos_β2  (Eq.5 from C.F.F Karney)
            cos_α0 = cos(α0);        
            
            validity |= EQUATORIAL_AZIMUTH; //validity : bitmask of GeodesicsOnEllipsoid
            
            //-----------
            //compute σ1
            //-----------
            σ1 = computeσi(sin_β1, cos_β1, dλ1, dφ1); //compute σ1
            sin_σ1 = sin(σ1);
            cos_σ1 = cos(σ1);
//            if ( ((dλ1 == 1)&&(dφ1 == 0) ) || ((sin_β1 == 1)&&(cos_β1 == 0) )){  // <=> ((α1 = PI/2) || (β1 == PI / 2))
//                σ1 = 0;
//                sin_σ1 = 0;
//                cos_σ1 = 1;
//            } else{
//                σ1 = atan2(sin_β1, (dφ1 * cos_β1));
//                sin_σ1 = sin(σ1);
//                cos_σ1 = cos(σ1);
//            }
            //-----------
            
            ω1 = computeωi(sin_σ1, cos_σ1);  //<=>  atan2(sin_α0 * sin_σ1, cos_σ1);
            
            validity |= AUXILIARY_START_POINT; //validity : bitmask of GeodesicsOnEllipsoid
            
            
            //--------------------
            //alternative method :
            //--------------------
//            β1 = computeβ(inverseFlattening, φ);
//            this.computeα0(β1, azimuth);
//            σ1 = computeσi(sin_β1, cos_β1, dλ1, dφ1); //compute σ1
//            sin_σ1 = sin(σ1);
//            cos_σ1 = cos(σ1);
//            this.ω1 = computeωi(sin_σ1, cos_σ1);
//            validity |= AUXILIARY_START_POINT; //validity : bitmask of GeodesicsOnEllipsoid
            //--------------------
            
        }

        /**
         * Compute the parameters on the auxiliary associated with the ending
         * point from the length of the geodesic on the ellipsoid s2 and the
         * starting point.
         *
         * @param s2 : ending point's distance from Equatorial point E along the geodesic
         * @return
         */
        final void resolveEndPoint(final double s2) {
            if (isInvalid(EQUATORIAL_AZIMUTH)) {
                throw new IllegalStateException("Uncomputed auxiliary parameter α0");
            }
            σ2 = ellipsoidal2SphericalArcLength(s2);
            sin_σ2 = sin(σ2);
            cos_σ2 = cos(σ2);
            
            final double β2 = atan2(cos_α0 * sin_σ2, hypot(cos_α0 * cos_σ2, sin_α0));
            sin_β2 = sin(β2);
            cos_β2 = cos(β2);
            ω2 =computeωi(sin_σ2, cos_σ2); //<=> ω2 = atan2(sin_α0 * sin_σ2, cos_σ2);
            
            validity |= (AUXILIARY_END_POINT|REDUCED_LATITUDES); //validity : bitmask of GeodesicsOnEllipsoid
            
            //--------------------
            //alternative method :
            //--------------------
//            σ2 = ellipsoidal2SphericalArcLength(s2);
//            β2 = computeβ(σ2);
//            ω2 = computeω(σ2);
//            validity |= AUXILIARY_END_POINT; //validity : bitmask of GeodesicsOnEllipsoid
            //--------------------
        }
         
        /**
         * Compute the longitude variation, λ12, from the starting azimuth
         * approximation in input. The computation follows the hybrid inverse
         * problem solving described in part 4 of 
         * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
         * Algortihms of geodesics from Charles F. F. Karney (SRI
         * International)</a>
         * 
         * Starting and ending point must be known.
         * 
         * This method computes the auxiliary sphere's values associated with
         * those parameters.
         *
         * @param α1Approximation : approximation of the starting azimuth
         * @return associated longitude variation.
         */
        private double hybridInverseProblemSolving(final double α1Approximation) {
            if (isInvalid(REDUCED_LATITUDES)) {
                auxiliarySpheres.computeReducedLatitudes();
            }
            
            double sin_α1_approx = sin(α1Approximation); // dλ1 approximation 
            double cos_α1_approx = cos(α1Approximation); // dφ1 approximation 

            //Compute α0, cos_α0, sin_α0 associated with the given α1 Approximation
            
            //************* Mieux? *************************
//            computeα0fromβ1(cos_α1_approx, sin_α1_approx);
            //**********************************************

            this.α0 = atan2(cos_α1_approx * cos_β1, hypot(cos_α1_approx, sin_α1_approx * sin_β1));
            sin_α0 = sin_α1_approx * cos_β1;
            cos_α0 = cos(α0);
            validity|= EQUATORIAL_AZIMUTH;

            double sin_α2_approx = sin_α0 / cos_β2; // as sin_α0 = sin_α1*cos_β1 =sin_α2*cos_β2  (Eq.5 from C.F.F Karney)
            double cos_α2_approx = sqrt(cos_α1_approx * cos_α1_approx * cos_β1 * cos_β1 + (cos_β2 * cos_β2 - cos_β1 * cos_β1)) / cos_β2; //(Eq.45 from C.F.F Karney)
            
            σ1 = computeσi(sin_β1, cos_β1, sin_α1_approx, cos_α1_approx); //compute σ1
            sin_σ1 = sin(σ1);
            cos_σ1 = cos(σ1);
            ω1 = computeωi(sin_σ1, cos_σ1);  //compute ω1
            σ2 = computeσi(sin_β2, cos_β2, sin_α2_approx, cos_α2_approx); //compute σ2 
            sin_σ2 = sin(σ2);
            cos_σ2 = cos(σ2);
            ω2 = computeωi(sin_σ2, cos_σ2);  //compute ω2
            
            dλ1 = sin_α1_approx;
            dφ1 = cos_α1_approx;
            dλ2 = sin_α2_approx;
            dφ2 = cos_α2_approx;
                    
            validity |= (AUXILIARY_START_POINT | AUXILIARY_END_POINT | STARTING_AZIMUTH |ENDING_AZIMUTH);
            
            return spherical2EllipsoidalLongitude(ω2, σ2) - spherical2EllipsoidalLongitude(ω1, σ1);

        }
        
        /**
         * Compute the derivative value dλ12(α1)/dα1 for Newton's method application 
         * for the computation of the start point azimuth α1.
         * 
         * Computation based on equations 38 of C.F.F Karney article. 
         *  
         * ---------------------------------------------------------------------
         * Important : 
         * -----------
         * EQUATORIAL_AZIMUTH, AUXILIARY_START_POINT, AUXILIARY_END_POINT, 
         * and ENDING_AZIMUTH  must be valid as :
         * A first approximation of the starting point's azimuth α1 (dλ1, dφ1)
         * must have been done and used to compute auxiliary parameters 
         * (σ1,σ2,ω1,ω2, associated sinus and cosinus) and 
         * end point 's azimuth (dλ2, dφ2). {@link #hybridInverseProblemSolving(double)} 
         * ---------------------------------------------------------------------
         * 
         * @return derivative value dλ12(α1)/dα1 from parameter
         */
        private double derivativeDλ12Dα1(){
            if(isInvalid(AUXILIARY_START_POINT | AUXILIARY_END_POINT | EQUATORIAL_AZIMUTH | ENDING_AZIMUTH)){
                throw new IllegalStateException("Reduced latitudes and equatorial azimuth must be provided in order to compute derivative value : dλ12(α1)/dα1");
                
            }
            if ((abs(dφ1) ==0)&& (cos_β1==cos_β2)&&(sin_β1==sin_β2) ){  ///particular case which should be carried after azimuth approximation. -> voir equation 47 p.48
                throw new IllegalStateException("Particular case where starting azimuth = pi/2 and equal starting and ending reduced latitude β1,β2. This case should be handle before calling the current method.");
            }
            
            final double k = second_ECCENTRICITY * cos_α0;
            
            return (  hypot(1, k*sin_σ2)*cos_σ1*sin_σ2
                    - hypot(1, k*sin_σ1)*sin_σ1*cos_σ2
                    - cos_σ1*cos_σ2*(compute_J(σ2)- compute_J(σ1))   // -> m12/b
                    ) * ellipsoid.getSemiMinorAxis() / (ellipsoid.getSemiMajorAxis()*dφ2*cos_β2);
            
//            final double I1_1 = compute_I1(σ1);      
//            final double I1_2 = compute_I1(σ2);
//            
//            final double I2_1 = compute_I2(σ1);      
//            final double I2_2 = compute_I2(σ2);            
//            
                            //            final double J_2 = I1_2 - I2_2;
                            //            final double J_1= I1_1 - I2_1;
            
                            //            final double m12_on_b = hypot(1, k*sin_σ2)*cos_σ1*sin_σ2
                            //                                  - hypot(1, k*sin_σ1)*sin_σ1*cos_σ2
                            //                                  - cos_σ1*cos_σ2*(I1_2 - I2_2 - I1_1 + I2_1);
            
//            return (  hypot(1, k*sin_σ2)*cos_σ1*sin_σ2
//                    - hypot(1, k*sin_σ1)*sin_σ1*cos_σ2
//                    - cos_σ1*cos_σ2*(I1_2 - I2_2 - I1_1 + I2_1)   // -> m12/b
//                    ) * ellipsoid.getSemiMinorAxis() / (ellipsoid.getSemiMajorAxis()*dφ2*cos_β2);
        }
        
               
        /**
         * Compute the reduced latitudes parameters associated with the starting
         * and ending latitudes.
         *
         */
        final void computeReducedLatitudes() {
            if (isInvalid(START_POINT | END_POINT)) {
                throw new IllegalStateException(Resources.format(
                        Resources.Keys.StartOrEndPointNotSet_1, Integer.signum(validity & START_POINT)));
            }
//            β1 = computeβ(ellipsoid.getInverseFlattening(), φ1);
            double β1 = atan((1 - (1 / ellipsoid.getInverseFlattening())) * tan(φ1));
            cos_β1 = cos(β1);
            sin_β1 = sin(β1);
//            β2 = computeβ(ellipsoid.getInverseFlattening(), φ2);
            double β2 = atan((1 - (1 / ellipsoid.getInverseFlattening())) * tan(φ2));
            cos_β2 = cos(β2);
            sin_β2 = sin(β2);
            validity |= REDUCED_LATITUDES;
        }
        
//        /**
//         * Compute the azimuth α0 (and associated sin and cos values)at the
//         * intercepting E point of the equator (0° latitude) from starting point
//         * reduced latitude and an input azimuth.
//         * 
//         * dφ1 and dλ1 aren't directly used as it can be computed for starting
//         * azimuth approximations. 
//         *
//         * @param cos_α1_azimuth
//         * @param sin_α1_azimuth
//         */
//        private void computeα0fromβ1(final double cos_α1_azimuth, final double sin_α1_azimuth) {
//            if (isInvalid(AUXILIARY_START_POINT)) {
//                throw new IllegalStateException("Unresolved auxiliary sphere's parameters for end point of the geodesic.");
//            }
//            this.α0 = atan2(sin_α1_azimuth * cos_β1, hypot(cos_α1_azimuth, sin_α1_azimuth * sin_β1));
////            sin_α0 = sin(α0); 
//            sin_α0 = sin_α1_azimuth * cos_β1;
//            cos_α0 = cos(α0);        
////            validity |= EQUATORIAL_AZIMUTH; //validity : bitmask of GeodesicsOnEllipsoid  
//        }
        
         /**
         * Compute the azimuth of the ending point's forward direction α2.
         * Valid equatorial azimuth and ending point are necessary. 
         *
         * Eq. 14 adapted for ending point.
         *
         */
        private void computeEndPointAzimuth() {
            if (isInvalid(EQUATORIAL_AZIMUTH | AUXILIARY_END_POINT)) {
                throw new IllegalStateException("Unresolved auxiliary sphere's parameters for end point of the geodesic.");
            }
            final double azimuth = atan2(sin_α0, cos_α0 * cos_σ2);
            dφ2 = cos(azimuth);                                 // sin(π/2 − θ)  =  cos(θ)
            dλ2 = sin(azimuth);                                 // cos(π/2 − θ)  =  sin(θ)
            validity |= ENDING_AZIMUTH;
        }

        /**
         * Compute the azimuth of the starting point's forward direction α1.
         * Valid equatorial azimuth and starting point are necessary. 
         *
         * Eq. 14 adapted for starting point.
         *
         */
        private void computeStartPointAzimuth() {
            if (isInvalid(EQUATORIAL_AZIMUTH | AUXILIARY_START_POINT)) {
                throw new IllegalStateException("Unresolved auxiliary sphere's parameters for end point of the geodesic.");
            }
//            final double azimuth = atan2(sin(α0), cos(α0) * cos(σ1));
            final double azimuth = atan2(sin_α0, cos_α0 * cos_σ1);

            dφ1 = cos(azimuth);       // sin(π/2 − θ)  =  cos(θ)
            dλ1 = sin(azimuth);       // cos(π/2 − θ)  =  sin(θ)
            validity |= STARTING_AZIMUTH;
        }

        
        //**********************************************************************
        // currently unused methods allowing to compute parameters individually
        //**********************************************************************
        
        
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
            double sin_α = sin(azimuth);
            double cos_β = cos(β);
            this.α0 = atan2(sin_α * cos_β, hypot(cos(azimuth), sin_α * sin(β)));
            sin_α0 = sin_α*cos_β; //Eq.5 from CFF Karney
//            sin_α0 = sin(α0); 
            cos_α0 = cos(α0);        
            validity |= EQUATORIAL_AZIMUTH; //validity : bitmask of GeodesicsOnEllipsoid
        }
        
        
//        private void computeα0fromβ2(final double cos_α2_azimuth, final double sin_α2_azimuth) {
//            this.α0 = atan2(sin_α2_azimuth * cos_β1, hypot(cos_α2_azimuth, sin_α2_azimuth * sin_β1));
//            sin_α0 = sin(α0); 
//            cos_α0 = cos(α0);        
//            validity |= EQUATORIAL_AZIMUTH; //validity : bitmask of GeodesicsOnEllipsoid
//        }

        /**
         * Compute arc length on the great circle σ, sin_σ and cos_σ between the
         * E point with α0 azimuth of the equator and the current point with the
         * reduced latitude β
         *
         * Eq. 11 in
         * <a href="https://link.springer.com/content/pdf/10.1007%2Fs00190-012-0578-z.pdf">
         * Algortihms of geodesics from Charles F. F. Karney (SRI
         * International)</a>
         * α0 must be known.
         *
         * @param sin_βi sinus of the associated reduce latitude. 
         * @param cos_βi cosinus of the associated reduce latitude.
         * @param sin_αi sinus of the associated azimuth. 
         * @param cos_αi cosinus of the associated azimuth.
         * @return σi
         */
        private double computeσi(double sin_βi, final double cos_βi, final double sin_αi, final double cos_αi) {
            if ( ((sin_αi == 1)&&(cos_αi == 0) ) || ((sin_βi == 1)&&(cos_βi == 0) )){  // <=> ((α1 = PI/2) || (β1 == PI / 2)) 
                //Note of the author Table2 from the article.
                return 0;
            } else{
                return atan2(sin_βi, (cos_αi * cos_βi));
            }
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
         * @param sin_σi
         * @param cos_σi
         * @return ωi
         */
        private double computeωi(final double sin_σi, final double cos_σi) {
            if (isInvalid(EQUATORIAL_AZIMUTH)) {
                throw new IllegalStateException("Uncomputed auxiliary parameter α0");
            }
            return atan2(sin_α0 * sin_σi, cos_σi);
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
            return atan((1 - (1 / inverseFlattening)) * tan(φ));
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
            return atan2(cos_α0 * sin(σ), hypot(cos_α0 * cos(σ), sin_α0));
        }
    }//end of internal class
    //==========================================================================
}
