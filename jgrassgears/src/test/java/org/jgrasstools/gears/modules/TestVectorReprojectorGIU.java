/*
 * This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.jgrasstools.gears.modules;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.jgrasstools.gears.io.vectorreader.OmsVectorReader;
import org.jgrasstools.gears.io.vectorwriter.OmsVectorWriter;
import org.jgrasstools.gears.modules.v.vectorreprojector.OmsVectorReprojector;
import org.jgrasstools.gears.utils.HMTestCase;
import org.jgrasstools.gears.utils.HMTestMaps;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

/**
 * Test for the {@link OmsVectorReprojector}.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestVectorReprojectorGIU extends HMTestCase {
    private double delta = 0.0000001;

    public void testVectorReprojector() throws Exception {

    	 OmsVectorReader reader = new OmsVectorReader();
         reader.file ="/Users/giuseppeformetta/Desktop/UNICAL/GiuseppeLUGLIO/F/Frane.shp";
         reader.process();
         SimpleFeatureCollection readFC = reader.outVector;


        OmsVectorReprojector reprojector = new OmsVectorReprojector();
        reprojector.inVector = readFC;
        reprojector.pCode = "EPSG:3004";
        reprojector.pm = pm;
        reprojector.process();

        CoordinateReferenceSystem sourceCRS = HMTestMaps.crs;
        CoordinateReferenceSystem targetCRS = CRS.decode("EPSG:4326");

        OmsVectorWriter writer = new OmsVectorWriter();
        writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/GiuseppeLUGLIO/F/FraneREP.shp";
        writer.inVector = reprojector.outVector;
        writer.process();


    }
}
