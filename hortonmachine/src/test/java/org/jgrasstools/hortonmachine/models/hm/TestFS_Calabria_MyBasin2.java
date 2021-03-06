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
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.demmanipulation.pitfiller.OmsPitfiller;
import org.jgrasstools.hortonmachine.modules.geomorphology.ab.OmsAb;
import org.jgrasstools.hortonmachine.modules.geomorphology.flow.OmsFlowDirections;
import org.jgrasstools.hortonmachine.modules.geomorphology.slope.OmsSlope;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.fs.OmsFs;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test ab.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestFS_Calabria_MyBasin2 extends HMTestCase {

	public void testAb() throws Exception {

		OmsRasterReader reader = new OmsRasterReader();
		//reader.file = "/Users/giuseppeformetta/Desktop/UNICAL/DEM5m/loc32633/mapset2/cell/mybasin_2_SLOPE";
		reader.file = "/Users/giuseppeformetta/Desktop/UNICAL/DEM5m/slope_mybasin4.asc";
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slope = reader.outRaster;

		OmsFs fs = new OmsFs();
//		fs.pathToSoilMecFile = "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin2/soil/soilMec0001MOD.txt";
//		fs.pathToSoilFile = "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin2/soil/soil0001.txt";
		fs.pathToSoilMecFile = "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin4/soil/soilMec0001MOD.txt";
		fs.pathToSoilFile = "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin4/soil/soil0001.txt";
		fs.nLayers = 8;
		fs.doUseGEOtopOutputs = true;
		fs.inSlope = slope;
		fs.pCohesion = 1;
		fs.pathToGEOtopOutputMaps="/Users/giuseppeformetta/Desktop/UNICAL/DEM5m/loc32633/M3/cell";
		fs.headerPsiGEOtop="pressureLIQ";
		fs.headerThetaGEOtop="thetaliq";
		fs.process();

	}

}
