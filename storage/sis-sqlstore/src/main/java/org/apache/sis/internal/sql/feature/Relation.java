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
package org.apache.sis.internal.sql.feature;

import java.util.Map;
import java.util.LinkedHashMap;
import java.util.Collection;
import java.util.Objects;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.DatabaseMetaData;
import org.apache.sis.internal.util.CollectionsExt;
import org.apache.sis.internal.metadata.sql.Reflection;
import org.apache.sis.storage.DataStoreContentException;
import org.apache.sis.storage.DataStoreException;
import org.apache.sis.util.collection.TreeTable;
import org.apache.sis.util.Debug;


/**
 * Description of a relation between two tables, as defined by foreigner keys.
 * Each {@link Table} may contain an arbitrary amount of relations.
 * Relations are defined in two directions:
 *
 * <ul>
 *   <li>{@link Direction#IMPORT}: primary keys of <em>other</em> tables are referenced by the foreigner keys
 *       of the table containing this {@code Relation}.</li>
 *   <li>{@link Direction#EXPORT}: foreigner keys of <em>other</em> tables are referencing the primary keys
 *       of the table containing this {@code Relation}.</li>
 * </ul>
 *
 * Instances of this class are created from the results of {@link DatabaseMetaData#getImportedKeys​ getImportedKeys​}
 * or {@link DatabaseMetaData#getExportedKeys​ getExportedKeys​} with {@code (catalog, schema, table)} parameters.
 *
 * @author  Johann Sorel (Geomatys)
 * @author  Martin Desruisseaux (Geomatys)
 * @version 1.0
 * @since   1.0
 * @module
 */
final class Relation extends TableReference {
    /**
     * Whether another table is <em>using</em> or is <em>used by</em> the table containing the {@link Relation}.
     */
    enum Direction {
        /**
         * Primary keys of other tables are referenced by the foreigner keys of the table containing the {@code Relation}.
         * In other words, the table containing {@code Relation} is <em>using</em> the {@link Relation#table}.
         *
         * @see DatabaseMetaData#getImportedKeys(String, String, String)
         */
        IMPORT(Reflection.PKTABLE_CAT, Reflection.PKTABLE_SCHEM, Reflection.PKTABLE_NAME,
               Reflection.PKCOLUMN_NAME, Reflection.FKCOLUMN_NAME),

        /**
         * Foreigner keys of other tables are referencing the primary keys of the table containing the {@code Relation}.
         * In other words, the table containing {@code Relation} is <em>used by</em> {@link Relation#table}.
         *
         * @see DatabaseMetaData#getExportedKeys(String, String, String)
         */
        EXPORT(Reflection.FKTABLE_CAT, Reflection.FKTABLE_SCHEM, Reflection.FKTABLE_NAME,
               Reflection.FKCOLUMN_NAME, Reflection.PKCOLUMN_NAME);

        /*
         * Note: another possible type of relation is the one provided by getCrossReference(…).
         * Inconvenient is that it requires to know the tables on both side of the relation.
         * But advantage is that it work with any set of columns having unique values
         * (not necessarily the primary key).
         */

        /**
         * The database {@link Reflection} key to use for fetching the name of other table column.
         * That column is part of a primary key if the direction is {@link #IMPORT}, or part of a
         * foreigner key if the direction is {@link #EXPORT}.
         */
        final String catalog, schema, table, column;

        /**
         * The database {@link Reflection} key to use for fetching the name of the column in the table
         * containing the {@code Relation}. That column is part of a foreigner key if the direction is
         * {@link #IMPORT}, or part of a primary key if the direction is {@link #EXPORT}.
         */
        final String containerColumn;

        /**
         * Creates a new {@code Direction} enumeration value.
         */
        private Direction(final String catalog, final String schema, final String table,
                          final String column, final String containerColumn)
        {
            this.catalog         = catalog;
            this.schema          = schema;
            this.table           = table;
            this.column          = column;
            this.containerColumn = containerColumn;
        }
    }

    /**
     * The columns of the other table that constitute a primary or foreigner key. Keys are the columns of the
     * other table and values are columns of the table containing this {@code Relation}.
     */
    private final Map<String,String> columns;

    /**
     * Whether entries in foreigner table will be deleted if the primary keys in the referenced table is deleted.
     * This is used as a hint for detecting what may be the "main" table in a relation.
     */
    final boolean cascadeOnDelete;

    /**
     * The other table which is related to the table containing this relation.
     * This is set during {@link Table} construction and should not be modified after that point.
     */
    private Table related;

    /**
     * The name of the feature property where the association to {@link #related} table will be stored.
     */
    private String propertyName;

    /**
     * Creates a new relation for an imported key. The given {@code ResultSet} must be positioned
     * on the first row of {@code DatabaseMetaData.getImportedKeys​(catalog, schema, table)} result,
     * and the result must be sorted in the order of the given keys:
     *
     * <ol>
     *   <li>{@link Direction#catalog}</li>
     *   <li>{@link Direction#schema}</li>
     *   <li>{@link Direction#table}</li>
     * </ol>
     *
     * Note that JDBC specification ensures this order if {@link Direction#IMPORT} is used with the result of
     * {@code getImportedKeys​} and {@link Direction#EXPORT} is used with the result of {@code getExportedKeys​}.
     *
     * <p>After construction, the {@code ResultSet} will be positioned on the first row of the next relation,
     * or be closed if the last row has been reached. This constructor always moves the given result set by at
     * least one row, unless an exception occurs.</p>
     */
    Relation(final Direction dir, final ResultSet reflect) throws SQLException, DataStoreContentException {
        super(reflect.getString(dir.catalog),
              reflect.getString(dir.schema),
              reflect.getString(dir.table),
              reflect.getString(Reflection.FK_NAME));

        final Map<String,String> m = new LinkedHashMap<>();
        boolean cascade = false;
        do {
            final String column = reflect.getString(dir.column);
            if (m.put(column, reflect.getString(dir.containerColumn)) != null) {
                throw new DataStoreContentException(Resources.format(Resources.Keys.DuplicatedEntity_2, "Column", column));
            }
            if (!cascade) {
                cascade = reflect.getInt(Reflection.DELETE_RULE) == DatabaseMetaData.importedKeyCascade;
            }
            if (!reflect.next()) {
                reflect.close();
                break;
            }
        } while (table.equals(reflect.getString(dir.table)) &&                  // Table name is mandatory.
                 Objects.equals(schema,  reflect.getString(dir.schema)) &&      // Schema and catalog may be null.
                 Objects.equals(catalog, reflect.getString(dir.catalog)));

        columns = CollectionsExt.compact(m);
        cascadeOnDelete = cascade;
    }

    /**
     * Invoked after construction for setting the target table.
     * Shall be invoked exactly once.
     */
    final void setRelatedTable(final Table table, final String property) throws DataStoreException {
        if (related != null) {
            throw new DataStoreException(Resources.format(Resources.Keys.InternalError));     // Should never happen.
        }
        related = table;
        propertyName = property;
    }

    /**
     * Returns the other table which is related to the table containing this relation.
     * See {@link Direction} for more details about the relation.
     */
    final Table getRelatedTable() {
        return related;
    }

    /**
     * Returns the name of the feature property where the association to {@link #related} table will be stored.
     */
    final String getPropertyName() {
        return propertyName;
    }

    /**
     * Returns the primary keys of the table related by the foreigner keys.
     */
    final Collection<String> getPrimaryKeys() {
        return columns.keySet();
    }

    /**
     * Returns the foreigner keys of the table that contains this relation. This method returns only
     * the foreigner keys known to this relation; this is not necessarily all the table foreigner keys.
     * Some columns may be used in more than one relation.
     */
    final Collection<String> getForeignerKeys() {
        return columns.values();
    }

    /**
     * Creates a tree representation of this relation for debugging purpose.
     *
     * @param  parent  the parent node where to add the tree representation.
     * @param  arrow   the symbol to use for relating the columns of two tables in a foreigner key.
     */
    @Debug
    void appendTo(final TreeTable.Node parent, final String arrow) {
        final TreeTable.Node node = newChild(parent, freeText);
        for (final Map.Entry<String,String> e : columns.entrySet()) {
            newChild(node, e.getValue() + arrow + e.getKey());
        }
    }

    /**
     * Formats a graphical representation of this relation for debugging purpose. This representation can
     * be printed to the {@linkplain System#out standard output stream} (for example) if the output device
     * uses a monospaced font and supports Unicode.
     */
    @Override
    public String toString() {
        return toString(this, (n) -> appendTo(n, " — "));
    }
}