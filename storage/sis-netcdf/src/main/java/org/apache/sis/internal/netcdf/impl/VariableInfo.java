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
package org.apache.sis.internal.netcdf.impl;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.Collection;
import java.util.Collections;
import java.util.regex.Matcher;
import java.io.IOException;
import javax.measure.Unit;
import ucar.nc2.constants.CF;
import ucar.nc2.constants.CDM;
import ucar.nc2.constants._Coordinate;
import org.apache.sis.coverage.grid.GridExtent;
import org.apache.sis.internal.netcdf.Decoder;
import org.apache.sis.internal.netcdf.DataType;
import org.apache.sis.internal.netcdf.Dimension;
import org.apache.sis.internal.netcdf.Grid;
import org.apache.sis.internal.netcdf.Variable;
import org.apache.sis.internal.netcdf.Resources;
import org.apache.sis.internal.storage.io.ChannelDataInput;
import org.apache.sis.internal.storage.io.HyperRectangleReader;
import org.apache.sis.internal.storage.io.Region;
import org.apache.sis.internal.util.StandardDateFormat;
import org.apache.sis.internal.util.UnmodifiableArrayList;
import org.apache.sis.storage.DataStoreException;
import org.apache.sis.storage.DataStoreContentException;
import org.apache.sis.storage.netcdf.AttributeNames;
import org.apache.sis.util.CharSequences;
import org.apache.sis.util.ArraysExt;
import org.apache.sis.util.Classes;
import org.apache.sis.util.Numbers;
import org.apache.sis.measure.Units;
import org.apache.sis.math.Vector;


/**
 * Description of a variable found in a netCDF file.
 * The natural ordering of {@code VariableInfo} is the order in which the variables appear in the stream of bytes
 * that make the netCDF file. Reading variables in natural order reduces the amount of channel seek operations.
 *
 * @author  Johann Sorel (Geomatys)
 * @author  Martin Desruisseaux (Geomatys)
 * @version 1.0
 * @since   0.3
 * @module
 */
final class VariableInfo extends Variable implements Comparable<VariableInfo> {
    /**
     * The names of attributes where to look for the description to be returned by {@link #getDescription()}.
     * We use the same attributes than the one documented in the {@link ucar.nc2.Variable#getDescription()} javadoc.
     */
    private static final String[] DESCRIPTION_ATTRIBUTES = {
        CDM.LONG_NAME,
        CDM.DESCRIPTION,
        CDM.TITLE,
        CF.STANDARD_NAME
    };

    /**
     * Helper class for reading a sub-area with a subsampling,
     * or {@code null} if {@code dataType} is not a supported type.
     */
    private final HyperRectangleReader reader;

    /**
     * The variable name.
     *
     * @see #getName()
     */
    private final String name;

    /**
     * The dimensions of this variable, in the order they appear in netCDF file. When iterating over the values stored in
     * this variable (a flattened one-dimensional sequence of values), index in the domain of {@code dimensions[length-1]}
     * varies faster, followed by index in the domain of {@code dimensions[length-2]}, <i>etc.</i>
     *
     * @see #getGridDimensions()
     * @see GridInfo#domain
     */
    final DimensionInfo[] dimensions;

    /**
     * The offset (in bytes) to apply for moving to the next record in variable data, or -1 if unknown.
     * A record is a <var>n-1</var> dimensional cube spanning all dimensions except the unlimited one.
     * This concept is relevant only for {@linkplain #isUnlimited() unlimited variables}.
     * For other kind of variables, this field is ignored and its value can be arbitrary.
     *
     * <p>On construction, this field is initialized to the number of bytes in this variable, omitting
     * the unlimited dimension (if any), or 0 if unknown. After {@link #complete(VariableInfo[])} has
     * been invoked, this field is modified as below:</p>
     *
     * <ul>
     *   <li>For {@linkplain #isUnlimited() unlimited variables}, this is the number of bytes to
     *       skip <strong>after</strong> the data has been read for reading the next record.</li>
     *   <li>For other variables, this field is ignored. In current implementation it is left to
     *       the size of this variable, and kept for information purpose only.</li>
     * </ul>
     */
    private long offsetToNextRecord;

    /**
     * The attributes associates to the variable, or an empty map if none.
     * Values can be:
     *
     * <ul>
     *   <li>{@link String} if the attribute contains a single textual value.</li>
     *   <li>{@link Number} if the attribute contains a single numerical value.</li>
     *   <li>{@link Vector} if the attribute contains many numerical values.</li>
     * </ul>
     *
     * If the value is a {@code String}, then leading and trailing spaces and control characters
     * should be trimmed by {@link String#trim()}.
     */
    private final Map<String,Object> attributes;

    /**
     * The netCDF type of data, or {@code null} if unknown.
     */
    private final DataType dataType;

    /**
     * The grid geometry associated to this variable, computed by {@link ChannelDecoder#getGrids()} when first needed.
     * May stay {@code null} if the variable is not a data cube. We do not need disambiguation between the case where
     * the grid has not yet been computed and the case where the computation has been done with {@code null} result,
     * because {@link #getGrid(Adjustment)} should be invoked only once per variable.
     *
     * @see #getGrid(Adjustment)
     */
    GridInfo grid;

    /**
     * {@code true} if this variable seems to be a coordinate system axis, as determined by comparing its name
     * with the name of all dimensions in the netCDF file. This information is computed at construction time
     * because requested more than once.
     */
    boolean isCoordinateSystemAxis;

    /**
     * The values of the whole variable, or {@code null} if not yet read. This vector should be assigned only
     * for relatively small variables, or for variables that are critical to the use of other variables
     * (for example the values in coordinate system axes).
     */
    private transient Vector values;

    /**
     * The {@code flag_meanings} values (used for enumeration values), or {@code null} if this variable is not
     * an enumeration.
     *
     * @see #isEnumeration()
     * @see #meaning(int)
     *
     * @todo Need to be consistent with {@code VariableWrapper}. We could move this field to {@link FeaturesInfo},
     *       or provides the same functionality in {@code VariableWrapper}. Whatever solution is chosen,
     *       {@code RasterResource.createEnumeration(…)} needs to use the mechanism common to both implementations.
     */
    private final String[] meanings;

    /**
     * Creates a new variable.
     *
     * @param  decoder     the netCDF file where this variable is stored.
     * @param  input       the channel together with a buffer for reading the variable data.
     * @param  name        the variable name.
     * @param  dimensions  the dimensions of this variable.
     * @param  attributes  the attributes associates to the variable, or an empty map if none.
     * @param  dataType    the netCDF type of data, or {@code null} if unknown.
     * @param  size        the variable size. May be inaccurate and ignored.
     * @param  offset      the offset where the variable data begins in the netCDF file.
     * @throws ArithmeticException if the variable size exceeds {@link Long#MAX_VALUE}.
     * @throws DataStoreContentException if a logical error is detected.
     */
    VariableInfo(final Decoder            decoder,
                 final ChannelDataInput   input,
                 final String             name,
                 final DimensionInfo[]    dimensions,
                 final Map<String,Object> attributes,
                       DataType           dataType,
                 final int                size,
                 final long               offset) throws DataStoreContentException
    {
        super(decoder);
        this.name       = name;
        this.dimensions = dimensions;
        this.attributes = attributes;
        final Object isUnsigned = getAttributeValue(CDM.UNSIGNED, "_unsigned");
        if (isUnsigned instanceof String) {
            dataType = dataType.unsigned(Boolean.valueOf((String) isUnsigned));
        }
        this.dataType = dataType;
        /*
         * The 'size' value is provided in the netCDF files, but doesn't need to be stored since it
         * is redundant with the dimension lengths and is not large enough for big variables anyway.
         * Instead we compute the length ourselves, excluding the unlimited dimension.
         */
        if (dataType != null && (offsetToNextRecord = dataType.size()) != 0) {
            for (int i=0; i<dimensions.length; i++) {
                final DimensionInfo dim = dimensions[i];
                if (!dim.isUnlimited) {
                    offsetToNextRecord = Math.multiplyExact(offsetToNextRecord, dim.length());
                } else if (i != 0) {
                    // Unlimited dimension, if any, must be first in a netCDF 3 classic format.
                    throw new DataStoreContentException(getLocale(), Decoder.FORMAT_NAME, input.filename, null);
                }
            }
            reader = new HyperRectangleReader(dataType.number, input, offset);
        } else {
            reader = null;
        }
        /*
         * If the value that we computed ourselves does not match the value declared in the netCDF file,
         * maybe for some reason the writer used a different layout.  For example maybe it inserted some
         * additional padding.
         */
        if (size != -1) {                           // Maximal unsigned value, means possible overflow.
            final long expected = paddedSize();
            final long actual = Integer.toUnsignedLong(size);
            if (actual != expected) {
                if (expected != 0) {
                    warning(ChannelDecoder.class, "readVariables",          // Caller of this constructor.
                            Resources.Keys.MismatchedVariableSize_3, getFilename(), name, actual - expected);
                }
                if (actual > offsetToNextRecord) {
                    offsetToNextRecord = actual;
                }
            }
        }
        /*
         * If the "_CoordinateAliasForDimension" attribute is defined, then its value will be used
         * instead of the variable name when determining if the variable is a coordinate system axis.
         * "_CoordinateVariableAlias" seems to be a legacy attribute name for the same purpose.
         */
        if (dimensions.length == 1) {
            Object value = getAttributeValue(_Coordinate.AliasForDimension, "_coordinatealiasfordimension");
            if (value == null) {
                value = getAttributeValue("_CoordinateVariableAlias", "_coordinatevariablealias");
                if (value == null) {
                    value = name;
                }
            }
            isCoordinateSystemAxis = dimensions[0].name.equals(value);
        }
        /*
         * Verify if this variable is an enumeration. If yes, we remove the attributes that define the
         * enumeration since those attributes may be verbose and "pollute" the variable definition.
         */
        if (!attributes.isEmpty()) {    // For avoiding UnsupportedOperationException if unmodifiable map.
            final Object flags = attributes.remove(AttributeNames.FLAG_MEANINGS);
            if (flags != null) {
                meanings = (String[]) CharSequences.split(flags.toString(), ' ');
                return;
            }
        }
        meanings = null;
    }

    /**
     * Returns the value of {@link #offsetToNextRecord} with padding applied. This value is true
     * only if this method is invoked <strong>before</strong> {@link #complete(VariableInfo[])}.
     */
    private long paddedSize() {
        return Math.addExact(offsetToNextRecord, Integer.BYTES - 1) & ~(Integer.BYTES - 1);
    }

    /**
     * Performs the final adjustment of the {@link #offsetToNextRecord} field of all the given variables.
     * This method applies padding except for the special case documented in netCDF specification:
     * <cite>"In the special case when there is only one {@linkplain #isUnlimited() record variable}
     * and it is of type character, byte, or short, no padding is used between record slabs,
     * so records after the first record do not necessarily start on four-byte boundaries"</cite>.
     *
     * <p>After padding has been applied, this method set the {@link #offsetToNextRecord} of all unlimited
     * variables to the number of bytes to skip before reading the next record.</p>
     *
     * @throws ArithmeticException if the stride between two records exceeds {@link Long#MAX_VALUE}.
     */
    static void complete(final VariableInfo[] variables) {
        final Set<CharSequence> referencedAsAxis = new HashSet<>();
        final VariableInfo[] unlimited = new VariableInfo[variables.length];
        int     count        = 0;               // Number of valid elements in the 'unlimited' array.
        long    recordStride = 0;               // Sum of the size of all variables having a unlimited dimension.
        boolean isUnknown    = false;           // True if 'total' is actually unknown.
        for (final VariableInfo variable : variables) {
            // Opportunistically store names of all axes listed in "coordinates" attributes of all variables.
            referencedAsAxis.addAll(Arrays.asList(variable.getCoordinateVariables()));
            if (variable.isUnlimited()) {
                final long paddedSize = variable.paddedSize();
                unlimited[count++] = variable;
                isUnknown |= (paddedSize == 0);
                recordStride = Math.addExact(recordStride, paddedSize);
            }
        }
        if (isUnknown) {
            // If at least one unlimited variable has unknown size, the whole record stride is unknown.
            for (int i=0; i<count; i++) {
                unlimited[i].offsetToNextRecord = -1;
            }
        } else if (count == 1) {
            unlimited[0].offsetToNextRecord = 0;        // Special case cited in method javadoc.
        } else for (int i=0; i<count; i++) {
            unlimited[i].offsetToNextRecord = recordStride - unlimited[i].offsetToNextRecord;
        }
        /*
         * If some variables have a "coordinates" attribute listing names of variables used as axes,
         * mark those variables as axes. We perform this check as a complement for the check done in
         * the constructor in order to detect two-dimensional axes. Example:
         *
         *    dimensions:
         *        Y = 471 ;
         *        X = 720 ;
         *    variables:
         *        float lon(Y, X) ;                       // Not detected as coordinate axis by the constructor.
         *            lon:long_name = "longitude" ;
         *            lon:units = "degree_east" ;
         *        float lat(Y, X) ;                       // Not detected as coordinate axis by the constructor.
         *            lat:long_name = "latitude" ;
         *            lat:units = "degree_north" ;
         *        short temperature(Y, X) ;
         *            temperature:long_name = "Potential Temperature" ;
         *            temperature:coordinates = "lon lat" ;
         */
        if (!referencedAsAxis.isEmpty()) {
            for (final VariableInfo variable : variables) {
                if (referencedAsAxis.remove(variable.name)) {
                    variable.isCoordinateSystemAxis = true;
                    if (referencedAsAxis.isEmpty()) break;
                }
            }
        }
    }

    /**
     * Returns the name of the netCDF file containing this variable, or {@code null} if unknown.
     */
    @Override
    public String getFilename() {
        if (reader != null) {
            final String filename = reader.filename();
            if (filename != null) return filename;
        }
        return super.getFilename();
    }

    /**
     * Returns the name of this variable.
     */
    @Override
    public String getName() {
        return name;
    }

    /**
     * Returns the description of this variable, or {@code null} if none.
     * This method searches for the first attribute named {@code "long_name"},
     * {@code "description"}, {@code "title"} or {@code "standard_name"}.
     */
    @Override
    public String getDescription() {
        for (final String attributeName : DESCRIPTION_ATTRIBUTES) {
            final String value = getAttributeAsString(attributeName);
            if (value != null) {
                return value;
            }
        }
        return null;
    }

    /**
     * Returns the unit of measurement as a string, or {@code null} if none.
     */
    @Override
    protected String getUnitsString() {
        return getAttributeAsString(CDM.UNITS);
    }

    /**
     * Parses the given unit symbol and set the {@link #epoch} if the parsed unit is a temporal unit.
     * This method is called by {@link #getUnit()}.
     */
    @Override
    protected Unit<?> parseUnit(String symbols) {
        final Matcher parts = TIME_UNIT_PATTERN.matcher(symbols);
        if (parts.matches()) {
            /*
             * If we enter in this block, the unit is of the form "days since 1970-01-01 00:00:00".
             * The TIME_PATTERN splits the string in two parts, "days" and "1970-01-01 00:00:00".
             * The parse method will replace the space between date and time by 'T' letter.
             */
            epoch = StandardDateFormat.parseInstantUTC(parts.group(2));
            symbols = parts.group(1);
        }
        return Units.valueOf(symbols);
    }

    /**
     * Returns the type of data, or {@code UNKNOWN} if the data type is unknown to this method.
     * If this variable has a {@code "_Unsigned = true"} attribute, then the returned data type
     * will be a unsigned variant.
     */
    @Override
    public DataType getDataType() {
        return dataType;
    }

    /**
     * Returns {@code true} if this variable is an enumeration.
     */
    final boolean isEnumeration() {
        return meanings != null;
    }

    /**
     * Returns whether this variable can grow. A variable is unlimited if at least one of its dimension is unlimited.
     * In netCDF 3 classic format, only the first dimension can be unlimited.
     */
    @Override
    protected boolean isUnlimited() {
        // The constructor verified that only the first dimension is unlimited.
        return (dimensions.length != 0) && dimensions[0].isUnlimited;
    }

    /**
     * Returns {@code true} if this variable seems to be a coordinate system axis,
     * determined by comparing its name with the name of all dimensions in the netCDF file.
     * Also determined by inspection of {@code "coordinates"} attribute on other variables.
     */
    @Override
    protected boolean isCoordinateSystemAxis() {
        return isCoordinateSystemAxis;
    }

    /**
     * Returns the value of the {@code "_CoordinateAxisType"} attribute, or {@code null} if none.
     */
    final String getAxisType() {
        final Object value = getAttributeValue(_Coordinate.AxisType, "_coordinateaxistype");
        return (value instanceof String) ? (String) value : null;
    }

    /**
     * Returns the names of variables to use as axes for this variable, or an empty array if none.
     */
    final CharSequence[] getCoordinateVariables() {
        return CharSequences.split(getAttributeAsString(CF.COORDINATES), ' ');
    }

    /**
     * Returns a builder for the grid geometry of this variable, or {@code null} if this variable is not a data cube.
     * The grid geometry builders are opportunistically cached in {@code VariableInfo} instances after they have been
     * computed by {@link ChannelDecoder#getGrids()}. This method delegates to the super-class method only if the grid
     * requires more analysis than the one performed by {@link ChannelDecoder}.
     *
     * @see ChannelDecoder#getGrids()
     */
    @Override
    protected Grid getGrid(final Adjustment adjustment) throws IOException, DataStoreException {
        if (grid == null) {
            decoder.getGrids();                               // Force calculation of grid geometries if not already done.
            if (grid == null) {                               // May have been computed as a side-effect of decoder.getGrids().
                grid = (GridInfo) super.getGrid(adjustment);  // Non-null if grid dimensions are different than this variable.
            }
        }
        return grid;
    }

    /**
     * Returns the dimensions of this variable in the order they are declared in the netCDF file.
     * The dimensions are those of the grid, not the dimensions (or axes) of the coordinate system.
     * In ISO 19123 terminology, the {@linkplain Dimension#length() dimension lengths} give the upper
     * corner of the grid envelope plus one. The lower corner is always (0, 0, …, 0).
     */
    @Override
    public List<Dimension> getGridDimensions() {
        return UnmodifiableArrayList.wrap(dimensions);
    }

    /**
     * Returns the names of all attributes associated to this variable.
     *
     * @return names of all attributes associated to this variable.
     */
    @Override
    public Collection<String> getAttributeNames() {
        return Collections.unmodifiableSet(attributes.keySet());
    }

    /**
     * Returns the type of the attribute of the given name, or {@code null} if the given attribute is not found.
     * If the attribute contains more than one value, then this method returns a {@code Vector.class} subtype.
     *
     * @param  attributeName  the name of the attribute for which to get the type, preferably in lowercase.
     * @return type of the given attribute, or {@code null} if the attribute does not exist.
     */
    @Override
    public Class<?> getAttributeType(final String attributeName) {
        return Classes.getClass(getAttributeValue(attributeName));
    }

    /**
     * Returns the value of the given attribute, or {@code null} if none.
     * This method should be invoked only for hard-coded names that mix lower-case and upper-case letters.
     *
     * @param  attributeName  name of attribute to search, in the expected case.
     * @param  lowerCase      the all lower-case variant of {@code attributeName}.
     * @return variable attribute value of the given name, or {@code null} if none.
     */
    private Object getAttributeValue(final String attributeName, final String lowerCase) {
        Object value = getAttributeValue(attributeName);
        if (value == null) {
            value = attributes.get(lowerCase);
        }
        return value;
    }

    /**
     * Returns the value of the given attribute, or {@code null} if none.
     * This method does not search the lower-case variant of the given name because the argument given to this method
     * is usually a hard-coded value from {@link CF} or {@link CDM} conventions, which are already in lower-cases.
     *
     * <p>All other {@code getAttributeFoo(…)} methods in this class or parent class ultimately invoke this method.
     * This provides a single point to override if the functionality needs to be extended.</p>
     *
     * @param  attributeName  name of attribute to search, in the expected case.
     * @return variable attribute value of the given name, or {@code null} if none.
     */
    @Override
    protected Object getAttributeValue(final String attributeName) {
        return attributes.get(attributeName);
    }

    /**
     * Sets the values in this variable. The values are normally read from the netCDF file by the {@link #read()} method,
     * but this {@code setValues(Object)} method may also be invoked if we want to overwrite those values.
     *
     * @param  array  the values as an array of primitive type (for example {@code float[]}.
     */
    final void setValues(final Object array) {
        Vector data = createDecimalVector(array, dataType.isUnsigned);
        /*
         * This method is usually invoked with vector of increasing or decreasing values. Set a tolerance threshold to the
         * precision of greatest (in magnitude) number, provided that this precision is not larger than increment. If values
         * are not sorted in increasing or decreasing order, the tolerance computed below will be smaller than it could be.
         * This is okay since it will cause more conservative compression (i.e. it does not increase the risk of data loss).
         */
        double tolerance = 0;
        if (Numbers.isFloat(data.getElementType())) {
            final int n = data.size() - 1;
            if (n >= 0) {
                double first = data.doubleValue(0);
                double last  = data.doubleValue(n);
                double inc   = Math.abs((last - first) / n);
                if (!Double.isNaN(inc)) {
                    double ulp = Math.ulp(Math.max(Math.abs(first), Math.abs(last)));
                    tolerance = Math.min(inc, ulp);
                }
            }
        }
        values = data.compress(tolerance);
        values = SHARED_VECTORS.unique(values);
    }

    /**
     * Reads all the data for this variable and returns them as an array of a Java primitive type.
     * Multi-dimensional variables are flattened as a one-dimensional array (wrapped in a vector).
     * Fill values/missing values are replaced by NaN if {@link #hasRealValues()} is {@code true}.
     * The vector is cached and returned as-is in all future invocation of this method.
     *
     * @throws ArithmeticException if the size of the variable exceeds {@link Integer#MAX_VALUE}, or other overflow occurs.
     */
    @Override
    @SuppressWarnings("ReturnOfCollectionOrArrayField")
    public Vector read() throws IOException, DataStoreContentException {
        if (values == null) {
            if (reader == null) {
                throw new DataStoreContentException(unknownType());
            }
            final int    dimension   = dimensions.length;
            final long[] lower       = new long[dimension];
            final long[] upper       = new long[dimension];
            final int [] subsampling = new int [dimension];
            for (int i=0; i<dimension; i++) {
                upper[i] = dimensions[(dimension - 1) - i].length();
                subsampling[i] = 1;
            }
            final Region region = new Region(upper, lower, upper, subsampling);
            applyUnlimitedDimensionStride(region);
            Object array = reader.read(region);
            replaceNaN(array);
            /*
             * If we can convert a double[] array to a float[] array, we should do that before
             * to invoke 'setValues(array)' - we can not rely on data.compress(tolerance). The
             * reason is because we assume that float[] arrays are accurate in base 10 even if
             * the data were originally stored as doubles. The Vector class does not make such
             * assumption since it is specific to what we observe with netCDF files. To enable
             * this assumption, we need to convert to float[] before createDecimalVector(…).
             */
            if (array instanceof double[]) {
                final float[] copy = ArraysExt.copyAsFloatsIfLossless((double[]) array);
                if (copy != null) array = copy;
            }
            setValues(array);
        }
        return values;
    }

    /**
     * If this variable uses the unlimited dimension, we have to skip the records of all other unlimited variables
     * before to reach the next record of this variable.  Current implementation can do that only if the number of
     * bytes to skip is a multiple of the data type size. It should be the case most of the time because variables
     * in netCDF files have a 4 bytes padding. It may not work however if the variable uses {@code long} or
     * {@code double} type.
     */
    private void applyUnlimitedDimensionStride(final Region region) throws DataStoreContentException {
        if (isUnlimited()) {
            final int dataSize = reader.dataSize();
            if (offsetToNextRecord < 0 || (offsetToNextRecord % dataSize) != 0) {
                throw new DataStoreContentException(resources()
                        .getString(Resources.Keys.CanNotComputeVariablePosition_2, getFilename(), name));
            }
            region.increaseStride(dimensions.length - 1, offsetToNextRecord / dataSize);
        }
    }

    /**
     * Reads a subsampled sub-area of the variable.
     * Multi-dimensional variables are flattened as a one-dimensional array (wrapped in a vector).
     * Array elements are in "natural" order (inverse of netCDF order).
     *
     * @param  area         indices of cell values to read along each dimension, in "natural" order.
     * @param  subsampling  subsampling along each dimension. 1 means no subsampling.
     * @return the data as an array of a Java primitive type.
     * @throws ArithmeticException if the size of the region to read exceeds {@link Integer#MAX_VALUE}, or other overflow occurs.
     */
    @Override
    public Vector read(final GridExtent area, final int[] subsampling) throws IOException, DataStoreException {
        if (reader == null) {
            throw new DataStoreContentException(unknownType());
        }
        if (values != null) {
            throw new DataStoreException();     // TODO: create a view.
        }
        /*
         * NetCDF sorts datas in reverse dimension order. Example:
         *
         * DIMENSIONS:
         *   time: 3
         *   lat : 2
         *   lon : 4
         *
         * VARIABLES:
         *   temperature (time,lat,lon)
         *
         * DATA INDICES:
         *   (0,0,0) (0,0,1) (0,0,2) (0,0,3)
         *   (0,1,0) (0,1,1) (0,1,2) (0,1,3)
         *   (1,0,0) (1,0,1) (1,0,2) (1,0,3)
         *   (1,1,0) (1,1,1) (1,1,2) (1,1,3)
         *   (2,0,0) (2,0,1) (2,0,2) (2,0,3)
         *   (2,1,0) (2,1,1) (2,1,2) (2,1,3)
         */
        final int dimension = dimensions.length;
        final long[] size  = new long[dimension];
        final long[] lower = new long[dimension];
        final long[] upper = new long[dimension];
        for (int i=0; i<dimension; i++) {
            lower[i] = area.getLow(i);
            upper[i] = Math.incrementExact(area.getHigh(i));
            size [i] = dimensions[(dimension - 1) - i].length();
        }
        final Region region = new Region(size, lower, upper, subsampling);
        applyUnlimitedDimensionStride(region);
        final Object array = reader.read(region);
        replaceNaN(array);
        return Vector.create(array, dataType.isUnsigned);
    }

    /**
     * Returns a coordinate for this two-dimensional grid coordinate axis.
     * This is (indirectly) a callback method for {@link Grid#getAxes(Decoder)}.
     *
     * @throws ArithmeticException if the axis size exceeds {@link Integer#MAX_VALUE}, or other overflow occurs.
     */
    @Override
    protected double coordinateForAxis(final int j, final int i) throws IOException, DataStoreException {
        assert j >= 0 && j < dimensions[0].length : j;
        assert i >= 0 && i < dimensions[1].length : i;
        final long n = dimensions[1].length();
        return read().doubleValue(Math.toIntExact(i + n*j));
    }

    /**
     * Returns the meaning of the given ordinal value, or {@code null} if none.
     * Callers must have verified that {@link #isEnumeration()} returned {@code true}
     * before to invoke this method
     *
     * @param  ordinal  the ordinal of the enumeration for which to get the value.
     * @return the value associated to the given ordinal, or {@code null} if none.
     */
    final String meaning(final int ordinal) {
        return (ordinal >= 0 && ordinal < meanings.length) ? meanings[ordinal] : null;
    }

    /**
     * Returns the error message for an unknown data type.
     */
    private String unknownType() {
        return resources().getString(Resources.Keys.UnsupportedDataType_3, getFilename(), name, dataType);
    }

    /**
     * Returns -1 if this variable is located before the other variable in the stream of bytes that make
     * the netCDF file, or +1 if it is located after.
     */
    @Override
    public int compareTo(final VariableInfo other) {
        int c = Long.compare(reader.origin, other.reader.origin);
        if (c == 0) c = name.compareTo(other.name);                 // Should not happen, but we are paranoiac.
        return c;
    }
}
