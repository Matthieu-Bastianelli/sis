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
package org.apache.sis.internal.coverage;

import java.nio.Buffer;
import java.nio.ReadOnlyBufferException;
import java.awt.Point;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferShort;
import java.awt.image.DataBufferUShort;
import java.awt.image.DataBufferInt;
import java.awt.image.DataBufferFloat;
import java.awt.image.DataBufferDouble;
import java.awt.image.SampleModel;
import java.awt.image.BandedSampleModel;
import java.awt.image.ComponentSampleModel;
import java.awt.image.PixelInterleavedSampleModel;
import java.awt.image.RasterFormatException;
import java.awt.image.WritableRaster;
import org.apache.sis.util.ArraysExt;
import org.apache.sis.util.Numbers;
import org.apache.sis.util.Static;
import org.apache.sis.util.Workaround;
import org.apache.sis.util.ArgumentChecks;
import org.apache.sis.internal.feature.Resources;


/**
 * Creates rasters from given properties.
 *
 * @author  Martin Desruisseaux (IRD, Geomatys)
 * @version 1.0
 * @since   1.0
 * @module
 */
public final class RasterFactory extends Static {
    /**
     * Do not allow instantiation of this class.
     */
    private RasterFactory() {
    }

    /**
     * Wraps the given data buffer in a raster.
     * The sample model type is selected according the number of bands and the pixel stride.
     * The number of bands is determined by {@code bandOffsets.length}, which should be one of followings:
     *
     * <ul>
     *   <li>For banded sample model, all {@code bandOffsets} can be zero.</li>
     *   <li>For interleaved sample model ({@code buffer.getNumBanks()} = 1), each band needs a different offset.
     *       They may be 0, 1, 2, 3….</li>
     * </ul>
     *
     * @param  buffer          buffer that contains the sample values.
     * @param  width           raster width in pixels.
     * @param  height          raster height in pixels.
     * @param  pixelStride     number of data elements between two samples for the same band on the same line.
     * @param  scanlineStride  number of data elements between a given sample and the corresponding sample in the same column of the next line.
     * @param  bankIndices     bank indices for each band, or {@code null} for 0, 1, 2, 3….
     * @param  bandOffsets     number of data elements from the first element of the bank to the first sample of the band.
     * @param  location        the upper-left corner of the raster, or {@code null} for (0,0).
     * @return a raster built from given properties.
     * @throws NullPointerException if {@code buffer} is {@code null}.
     * @throws RasterFormatException if the width or height is less than or equal to zero, or if there is an integer overflow.
     *
     * @see WritableRaster#createInterleavedRaster(DataBuffer, int, int, int, int, int[], Point)
     * @see WritableRaster#createBandedRaster(DataBuffer, int, int, int, int[], int[], Point)
     */
    @SuppressWarnings("fallthrough")
    public static WritableRaster createRaster(final DataBuffer buffer,
            final int width, final int height, final int pixelStride, final int scanlineStride,
            int[] bankIndices, final int[] bandOffsets, final Point location)
    {
        /*
         * We do not verify the argument validity. Since this class is internal, caller should have done verification
         * itself. Furthermore those arguments are verified by WritableRaster constructors anyway.
         */
        final int dataType = buffer.getDataType();
        /*
         * This SampleModel variable is a workaround for WritableRaster static methods not supporting all data types.
         * If 'dataType' is unsupported, then we create a SampleModel ourselves in the 'switch' statements below and
         * use it for creating a WritableRaster at the end of this method. This variable, together with the 'switch'
         * statements, may be removed in a future SIS version if all types become supported by the JDK.
         */
        @Workaround(library = "JDK", version = "10")
        final SampleModel model;
        if (buffer.getNumBanks() == 1 && (bankIndices == null || bankIndices[0] == 0)) {
            /*
             * Sample data are stored for all bands in a single bank of the DataBuffer, in an interleaved fashion.
             * Each sample of a pixel occupies one data element of the DataBuffer, with a different offset since
             * the buffer beginning. The number of bands is inferred from bandOffsets.length.
             */
            switch (dataType) {
                case DataBuffer.TYPE_BYTE:
                case DataBuffer.TYPE_USHORT: {
                    // 'scanlineStride' and 'pixelStride' really interchanged in that method signature.
                    return WritableRaster.createInterleavedRaster(buffer, width, height, scanlineStride, pixelStride, bandOffsets, location);
                }
                case DataBuffer.TYPE_INT: {
                    if (bandOffsets.length == 1 && pixelStride == 1) {
                        /*
                         * From JDK javadoc: "To create a 1-band Raster of type TYPE_INT, use createPackedRaster()".
                         * However this would require the creation of a PackedColorModel subclass. For SIS purposes,
                         * it is easier to create a banded sample model.
                         */
                        return WritableRaster.createBandedRaster(buffer, width, width, scanlineStride, new int[1], bandOffsets, location);
                    }
                    // else fallthrough.
                }
                default: {
                    model = new PixelInterleavedSampleModel(dataType, width, height, pixelStride, scanlineStride, bandOffsets);
                    break;
                }
            }
        } else {
            /*
             * Sample data are stored in different banks (arrays) for each band. If all pixels are consecutive (pixelStride = 1),
             * we have the classical banded sample model. Otherwise the type is not well identified; neither interleaved or banded.
             */
            if (bankIndices == null) {
                bankIndices = ArraysExt.range(0, bandOffsets.length);
            }
            if (pixelStride == 1) {
                switch (dataType) {
                    case DataBuffer.TYPE_BYTE:
                    case DataBuffer.TYPE_USHORT:
                    case DataBuffer.TYPE_INT: {
                        // This constructor supports only above-cited types.
                        return WritableRaster.createBandedRaster(buffer, width, width, scanlineStride, bankIndices, bandOffsets, location);
                    }
                    default: {
                        model = new BandedSampleModel(dataType, width, height, scanlineStride, bankIndices, bandOffsets);
                        break;
                    }
                }
            } else {
                model = new ComponentSampleModel(dataType, width, height, pixelStride, scanlineStride, bankIndices, bandOffsets);
            }
        }
        return WritableRaster.createWritableRaster(model, buffer, location);
    }

    /**
     * Returns the {@link DataBuffer} constant for the given type. The given {@code sample} class
     * should be a primitive type such as {@link Float#TYPE}. Wrappers class are also accepted.
     *
     * @param  sample    the primitive type or its wrapper class. May be {@code null}.
     * @param  unsigned  whether the type should be considered unsigned.
     * @return the {@link DataBuffer} type, or {@link DataBuffer#TYPE_UNDEFINED}.
     */
    public static int getType(final Class<?> sample, final boolean unsigned) {
        switch (Numbers.getEnumConstant(sample)) {
            case Numbers.BYTE:    if (unsigned) return DataBuffer.TYPE_BYTE; else break;
            case Numbers.SHORT:   return unsigned ? DataBuffer.TYPE_USHORT : DataBuffer.TYPE_SHORT;
            case Numbers.INTEGER: if (!unsigned) return DataBuffer.TYPE_INT; else break;
            case Numbers.FLOAT:   return DataBuffer.TYPE_FLOAT;
            case Numbers.DOUBLE:  return DataBuffer.TYPE_DOUBLE;
        }
        return DataBuffer.TYPE_UNDEFINED;
    }

    /**
     * Wraps the backing arrays of given NIO buffers into Java2D buffers.
     * This method wraps the underlying array of primitive types; data are not copied.
     * For each buffer, the data starts at {@linkplain Buffer#position() buffer position}
     * and ends at {@linkplain Buffer#limit() limit}.
     *
     * @param  dataType  type of buffer to create as one of {@link DataBuffer} constants.
     * @param  data      the data, one for each band.
     * @return buffer of the given type, or {@code null} if {@code dataType} is unrecognized.
     * @throws UnsupportedOperationException if a buffer is not backed by an accessible array.
     * @throws ReadOnlyBufferException if a buffer is backed by an array but is read-only.
     * @throws ArrayStoreException if the type of a backing array is not {@code dataType}.
     * @throws ArithmeticException if a buffer position overflows the 32 bits integer capacity.
     * @throws RasterFormatException if buffers do not have the same amount of remaining values.
     */
    public static DataBuffer wrap(final int dataType, final Buffer... data) {
        final int numBands = data.length;
        final Object[] arrays;
        switch (dataType) {
            case DataBuffer.TYPE_USHORT: // fall through
            case DataBuffer.TYPE_SHORT:  arrays = new short [numBands][]; break;
            case DataBuffer.TYPE_INT:    arrays = new int   [numBands][]; break;
            case DataBuffer.TYPE_BYTE:   arrays = new byte  [numBands][]; break;
            case DataBuffer.TYPE_FLOAT:  arrays = new float [numBands][]; break;
            case DataBuffer.TYPE_DOUBLE: arrays = new double[numBands][]; break;
            default: return null;
        }
        final int[] offsets = new int[numBands];
        int length = 0;
        for (int i=0; i<numBands; i++) {
            final Buffer buffer = data[i];
            ArgumentChecks.ensureNonNullElement("data", i, buffer);
            arrays [i] = buffer.array();
            offsets[i] = Math.addExact(buffer.arrayOffset(), buffer.position());
            final int r = buffer.remaining();
            if (i == 0) length = r;
            else if (length != r) {
                throw new RasterFormatException(Resources.format(Resources.Keys.MismatchedBandSize));
            }
        }
        switch (dataType) {
            case DataBuffer.TYPE_BYTE:   return new DataBufferByte  (  (byte[][]) arrays, length, offsets);
            case DataBuffer.TYPE_SHORT:  return new DataBufferShort ( (short[][]) arrays, length, offsets);
            case DataBuffer.TYPE_USHORT: return new DataBufferUShort( (short[][]) arrays, length, offsets);
            case DataBuffer.TYPE_INT:    return new DataBufferInt   (   (int[][]) arrays, length, offsets);
            case DataBuffer.TYPE_FLOAT:  return new DataBufferFloat ( (float[][]) arrays, length, offsets);
            case DataBuffer.TYPE_DOUBLE: return new DataBufferDouble((double[][]) arrays, length, offsets);
            default: return null;
        }
    }
}
