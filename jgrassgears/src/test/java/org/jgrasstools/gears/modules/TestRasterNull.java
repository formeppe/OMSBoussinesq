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

import static java.lang.Double.NaN;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.modules.r.rasternull.OmsRasterNull;
import org.jgrasstools.gears.utils.HMTestCase;
import org.jgrasstools.gears.utils.HMTestMaps;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
/**
 * Test {@link OmsRasterNull}.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestRasterNull extends HMTestCase {

    public void testRasterNull() throws Exception {
        double[][] data = new double[][]{//
        /*    */{5, 5, 5, 5, 5, 5, 5, 5, 5, 5}, //
                {5, 5, 6, 6, 6, 6, 6, 6, 6, 5}, //
                {5, 7, 6, 6, 6, 6, 6, 7, 7, 5}, //
                {5, 5, 5, 7, 6, 6, 6, 6, 5, 5}, //
                {5, 3, 4, 5, 5, 5, 5, 5, 5, 5}, //
                {5, 2, 3, 3, 4, 4, 4, 3, 3, 5}, //
                {5, 4, 4, 4, 4, 4, 5, 4, 4, 5}, //
                {5, 5, 5, 5, 5, 5, 5, 5, 5, 5}};

        HashMap<String, Double> envelopeParams = HMTestMaps.envelopeParams;
        CoordinateReferenceSystem crs = HMTestMaps.crs;
        GridCoverage2D inCoverage = CoverageUtilities.buildCoverage("data", data, envelopeParams, crs, true);

        OmsRasterNull transformer = new OmsRasterNull();
        transformer.inRaster = inCoverage;
        transformer.pValue = 4.0;
        transformer.process();
        GridCoverage2D outCoverage = transformer.outRaster;

        double[][] expected1 = new double[][]{//
        /*    */{5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}, //
                {5.0, 5.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 5.0}, //
                {5.0, 7.0, 6.0, 6.0, 6.0, 6.0, 6.0, 7.0, 7.0, 5.0}, //
                {5.0, 5.0, 5.0, 7.0, 6.0, 6.0, 6.0, 6.0, 5.0, 5.0}, //
                {5.0, 3.0, NaN, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}, //
                {5.0, 2.0, 3.0, 3.0, NaN, NaN, NaN, 3.0, 3.0, 5.0}, //
                {5.0, NaN, NaN, NaN, NaN, NaN, 5.0, NaN, NaN, 5.0}, //
                {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0} //
        };
        checkMatrixEqual(outCoverage.getRenderedImage(), expected1);

        transformer = new OmsRasterNull();
        transformer.inRaster = inCoverage;
        transformer.pValue = 4.0;
        transformer.pNull = 1.0;
        transformer.process();
        outCoverage = transformer.outRaster;

        double[][] expected2 = new double[][]{//
        /*    */{5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}, //
                {5.0, 5.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 5.0}, //
                {5.0, 7.0, 6.0, 6.0, 6.0, 6.0, 6.0, 7.0, 7.0, 5.0}, //
                {5.0, 5.0, 5.0, 7.0, 6.0, 6.0, 6.0, 6.0, 5.0, 5.0}, //
                {5.0, 3.0, 1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}, //
                {5.0, 2.0, 3.0, 3.0, 1.0, 1.0, 1.0, 3.0, 3.0, 5.0}, //
                {5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 1.0, 1.0, 5.0}, //
                {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0} //
        };
        checkMatrixEqual(outCoverage.getRenderedImage(), expected2);
    }
}
