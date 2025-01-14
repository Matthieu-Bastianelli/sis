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
package org.apache.sis.internal.referencing.provider;

import javax.xml.bind.annotation.XmlTransient;
import org.opengis.parameter.ParameterDescriptor;
import org.opengis.parameter.ParameterDescriptorGroup;
import org.opengis.parameter.ParameterNotFoundException;
import org.apache.sis.metadata.iso.citation.Citations;
import org.apache.sis.internal.util.Constants;
import org.apache.sis.parameter.Parameters;
import org.apache.sis.referencing.operation.projection.NormalizedProjection;


/**
 * The provider for <cite>"sinusoidal equal-area"</cite> projection.
 * This is a pseudo-cylindrical (or "false cylindrical") projection.
 * This projection method has no associated EPSG code.
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @version 1.0
 *
 * @see <a href="https://en.wikipedia.org/wiki/Sinusoidal_projection">Sinusoidal projection on Wikipedia</a>
 * @see <a href="http://geotiff.maptools.org/proj_list/sinusoidal.html">GeoTIFF parameters for Sinusoidal</a>
 *
 * @since 1.0
 * @module
 */
@XmlTransient
public final class Sinusoidal extends MapProjection {
    /**
     * For cross-version compatibility.
     */
    private static final long serialVersionUID = -3236247448683326299L;

    /**
     * Name of this projection.
     */
    public static final String NAME = "Sinusoidal";

    /**
     * The operation parameter descriptor for the <cite>Longitude of projection center</cite> (λ₀) parameter value.
     * Valid values range is [-180 … 180]° and default value is 0°.
     */
    public static final ParameterDescriptor<Double> CENTRAL_MERIDIAN = ESRI.CENTRAL_MERIDIAN;

    /**
     * The operation parameter descriptor for the <cite>False easting</cite> (FE) parameter value.
     * Valid values range is unrestricted and default value is 0 metre.
     */
    public static final ParameterDescriptor<Double> FALSE_EASTING = ESRI.FALSE_EASTING;

    /**
     * The operation parameter descriptor for the <cite>False northing</cite> (FN) parameter value.
     * Valid values range is unrestricted and default value is 0 metre.
     */
    public static final ParameterDescriptor<Double> FALSE_NORTHING = ESRI.FALSE_NORTHING;

    /**
     * The group of all parameters expected by this coordinate operation.
     */
    private static final ParameterDescriptorGroup PARAMETERS;
    static {
        PARAMETERS = builder().setCodeSpace(Citations.OGC, Constants.OGC)
                .addName      (NAME)
                .addName      ("Sanson-Flamsteed")
                .addName      (Citations.GEOTIFF,  "CT_Sinusoidal")
                .addIdentifier(Citations.GEOTIFF,  "24")
                .addName      (Citations.PROJ4,    "sinu")
                .createGroupForMapProjection(CENTRAL_MERIDIAN, FALSE_EASTING, FALSE_NORTHING);
    }

    /**
     * Constructs a new provider.
     */
    public Sinusoidal() {
        super(PARAMETERS);
    }

    /**
     * {@inheritDoc}
     *
     * @return the map projection created from the given parameter values.
     */
    @Override
    protected NormalizedProjection createProjection(final Parameters parameters) throws ParameterNotFoundException {
        return new org.apache.sis.referencing.operation.projection.Sinusoidal(this, parameters);
    }
}
