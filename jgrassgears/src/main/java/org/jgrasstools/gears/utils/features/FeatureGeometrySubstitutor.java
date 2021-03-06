/*
 * JGrass - Free Open Source Java GIS http://www.jgrass.org 
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Library General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Library General Public License
 * along with this library; if not, write to the Free Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.jgrasstools.gears.utils.features;

import java.util.List;

import org.geotools.factory.FactoryRegistryException;
import org.geotools.feature.SchemaException;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.AttributeDescriptor;
import org.opengis.feature.type.FeatureType;
import org.opengis.feature.type.GeometryDescriptor;

import com.vividsolutions.jts.geom.Geometry;

/**
 * Utility to clone features by modify only the geometry.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class FeatureGeometrySubstitutor {
    private SimpleFeatureType newFeatureType;

    public FeatureGeometrySubstitutor( SimpleFeatureType oldFeatureType ) throws Exception {
        this(oldFeatureType, null);
    }

    /**
     * @param oldFeatureType the {@link FeatureType} of the existing features.
     * @throws FactoryRegistryException 
     * @throws SchemaException
     */
    public FeatureGeometrySubstitutor( SimpleFeatureType oldFeatureType, Class< ? > newGeometryType ) throws Exception {

        List<AttributeDescriptor> oldAttributeDescriptors = oldFeatureType.getAttributeDescriptors();

        // create the feature type
        SimpleFeatureTypeBuilder b = new SimpleFeatureTypeBuilder();
        b.setName(oldFeatureType.getName());
        b.setCRS(oldFeatureType.getCoordinateReferenceSystem());

        if (newGeometryType == null) {
            b.addAll(oldAttributeDescriptors);
        } else {
            for( AttributeDescriptor attributeDescriptor : oldAttributeDescriptors ) {
                if (attributeDescriptor instanceof GeometryDescriptor) {
                    b.add("the_geom", newGeometryType);
                } else {
                    b.add(attributeDescriptor);
                }
            }
        }
        newFeatureType = b.buildFeatureType();
    }

    public SimpleFeatureType getNewFeatureType() {
        return newFeatureType;
    }

    /**
     * @param oldFeature the feature from which to clone the existing attributes from.
     * @param newGeometry new geometry to insert.
     * @return the new created feature, as merged from the old feature plus the new attributes.
     */
    public SimpleFeature substituteGeometry( SimpleFeature oldFeature, Geometry newGeometry ) {
        return substituteGeometry(oldFeature, newGeometry, null);
    }

    public SimpleFeature substituteGeometry( SimpleFeature oldFeature, Geometry newGeometry, String id ) {
        Object[] attributes = oldFeature.getAttributes().toArray();
        Object[] newAttributes = new Object[attributes.length];
        System.arraycopy(attributes, 0, newAttributes, 0, attributes.length);
        newAttributes[0] = newGeometry;
        // create the feature
        SimpleFeatureBuilder builder = new SimpleFeatureBuilder(newFeatureType);
        builder.addAll(newAttributes);
        SimpleFeature f = builder.buildFeature(id);
        return f;
    }

}