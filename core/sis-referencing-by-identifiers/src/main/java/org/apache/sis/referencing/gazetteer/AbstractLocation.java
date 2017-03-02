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
package org.apache.sis.referencing.gazetteer;

import java.util.Collection;
import java.util.Collections;
import org.opengis.metadata.citation.Party;
import org.opengis.metadata.extent.TemporalExtent;
import org.opengis.metadata.extent.GeographicExtent;
import org.opengis.metadata.extent.GeographicBoundingBox;
import org.opengis.geometry.coordinate.Position;
import org.opengis.referencing.gazetteer.Location;
import org.opengis.referencing.gazetteer.LocationType;
import org.opengis.util.InternationalString;
import org.apache.sis.util.iso.Types;
import org.apache.sis.metadata.iso.extent.Extents;


/**
 * Identifiable geographic place. A geographic place may be identified by a name (for example “Eiffel Tower”),
 * by a postcode, or any other method specified by the {@linkplain ModifiableLocationType#getIdentifications()
 * location type identifications}.
 *
 * <p>ISO 19112 describes the following properties as mandatory, but Apache SIS relaxes this restriction by
 * providing default values (possibly {@code null}) in most cases:</p>
 * <ul>
 *   <li><b>geographic identifier</b> (the value, for example a name or code)</li>
 *   <li><b>geographic extent</b> (the position of the identified thing)</li>
 *   <li><b>administrator</b> (who is responsible for this identifier)</li>
 *   <li><b>location type</b> (which specifies the nature of the identifier and its associated geographic location)</li>
 * </ul>
 *
 * The following properties are optional:
 * <ul>
 *   <li><b>temporal extent</b></li>
 *   <li><b>alternative geographic identifier</b></li>
 *   <li><b>position</b> (mandatory if the geographic identifier contains insufficient information to identify location)</li>
 *   <li><b>parent location instance</b></li>
 *   <li><b>child location instance</b></li>
 * </ul>
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @since   0.8
 * @version 0.8
 * @module
 */
public abstract class AbstractLocation implements Location {
    /**
     * The geographic identifier, or {@code null} if unspecified.
     *
     * @see #getGeographicIdentifier()
     */
    private final CharSequence identifier;

    /**
     * Creates a new location for the given geographic identifier.
     * This constructor accepts a {@code null} argument, but this is not recommended.
     *
     * @param identifier  the geographic identifier to be returned by {@link #getGeographicIdentifier()}.
     */
    protected AbstractLocation(final CharSequence identifier) {
        this.identifier = identifier;
    }

    /**
     * Returns a unique identifier for the location instance. The methods of identifying locations is specified
     * by the {@linkplain ModifiableLocationType#getIdentifications() location type identifications}.
     *
     * <div class="note"><b>Examples:</b>
     * if {@link LocationType#getIdentifications()} contain “name”, then geographic identifiers may be country
     * names like “Japan” or “France”, or places like “Eiffel Tower”. If location type identifications contain
     * “code”, then geographic identifiers may be “SW1P 3AD” postcode.
     * </div>
     *
     * In order to ensure that a geographic identifier is unique within a wider geographic domain,
     * the geographic identifier may need to include an identifier of an instance of a parent location type,
     * for example “Paris, Texas”.
     *
     * @return unique identifier for the location instance.
     *
     * @see ModifiableLocationType#getIdentifications()
     */
    @Override
    public InternationalString getGeographicIdentifier() {
        return Types.toInternationalString(identifier);
    }

    /**
     * Returns other identifier(s) for the location instance.
     * The default implementation returns an empty set.
     *
     * @return other identifier(s) for the location instance, or an empty collection if none.
     */
    @Override
    public Collection<? extends InternationalString> getAlternativeGeographicIdentifiers() {
        return Collections.emptySet();
    }

    /**
     * Returns the date of creation of this version of the location instance.
     * The default implementation returns {@code null}.
     *
     * @return date of creation of this version of the location instance, or {@code null} if none.
     */
    @Override
    public TemporalExtent getTemporalExtent() {
        return null;
    }

    /**
     * Returns a description of the location instance. This properties is mandatory according ISO 19112,
     * but Apache SIS nevertheless allows {@code null} value. If non-null, SIS implementations typically
     * provide instances of {@linkplain org.apache.sis.metadata.iso.extent.DefaultGeographicBoundingBox
     * geographic bounding boxes}.
     *
     * @return description of the location instance, or {@code null} if none.
     *
     * @see org.apache.sis.metadata.iso.extent.DefaultGeographicBoundingBox
     * @see org.apache.sis.metadata.iso.extent.DefaultBoundingPolygon
     */
    @Override
    public GeographicExtent getGeographicExtent() {
        return null;
    }

    /**
     * Returns coordinates of a representative point for the location instance.
     * The default implementation returns the centroid of the {@linkplain #getGeographicExtent() geographic extent}
     * if that extent is geographic bounding box.
     *
     * @return coordinates of a representative point for the location instance, or {@code null} if none.
     */
    @Override
    public Position getPosition() {
        final GeographicExtent extent = getGeographicExtent();
        return (extent instanceof GeographicBoundingBox) ? Extents.centroid((GeographicBoundingBox) extent) : null;
    }

    /**
     * Returns a description of the nature of this geographic identifier.
     *
     * @return the nature of the identifier and its associated geographic location.
     */
    @Override
    public abstract LocationType getLocationType();

    /**
     * Returns the organization responsible for defining the characteristics of the location instance.
     * The default implementation returns the {@linkplain ModifiableLocationType#getOwner() owner}.
     *
     * @return organization responsible for defining the characteristics of the location instance, or {@code null}.
     *
     * @see ModifiableLocationType#getOwner()
     * @see ReferencingByIdentifiers#getOverallOwner()
     */
    @Override
    public Party getAdministrator() {
        final LocationType type = getLocationType();
        return (type != null) ? type.getOwner() : null;
    }

    /**
     * Returns location instances of a different location type, for which this location instance is a sub-division.
     * The default implementation returns an empty list.
     *
     * @return parent locations, or an empty collection if none.
     *
     * @see ModifiableLocationType#getParents()
     */
    @Override
    public Collection<? extends Location> getParents() {
        return Collections.emptyList();
    }

    /**
     * Returns location instances of a different location type which subdivides this location instance.
     * The default implementation returns an empty list.
     *
     * @return child locations, or an empty collection if none.
     *
     * @see ModifiableLocationType#getChildren()
     */
    @Override
    public Collection<? extends Location> getChildren() {
        return Collections.emptyList();
    }
}