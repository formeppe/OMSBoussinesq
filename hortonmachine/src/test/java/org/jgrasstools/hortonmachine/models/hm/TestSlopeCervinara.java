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
package org.jgrasstools.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.geomorphology.slope.OmsSlope;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
/**
 * Tests the {@link OmsSlope} module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestSlopeCervinara extends HMTestCase {

    public void testSlope() throws Exception {
    	OmsRasterReader reader = new OmsRasterReader();
		reader.file = "/Users/giuseppeformetta/Desktop/Software/l32632/Cervinaralab/cell/dtmCervinaraLab";
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D pit = reader.outRaster;

    	 reader = new OmsRasterReader();
		reader.file = "/Users/giuseppeformetta/Desktop/Software/l32632/Cervinaralab/cell/flow";
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D flow = reader.outRaster;
        OmsSlope slope = new OmsSlope();
        slope.inPit = pit;
        slope.inFlow = flow;
        slope.pm = pm;

        slope.process();

        GridCoverage2D slopeCoverage = slope.outSlope;
        checkMatrixEqual(slopeCoverage.getRenderedImage(), HMTestMaps.slopeData, 0.01);
    }

}
