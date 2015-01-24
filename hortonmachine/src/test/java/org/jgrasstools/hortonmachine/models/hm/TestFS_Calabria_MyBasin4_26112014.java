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
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.fs.OmsFs16042014;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test ab.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestFS_Calabria_MyBasin4_26112014 extends HMTestCase {

	public void testAb() throws Exception {

		for (int i = 600; i < 1450; i=i+200) {
			OmsFs16042014 fs = new OmsFs16042014();
			// fs.pathToSoilMecFile =
			// "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin2/soil/soilMec0001MOD.txt";
			// fs.pathToSoilFile =
			// "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin2/soil/soil0001.txt";
			fs.pathToSoilMecFile = "/Users/giuseppeformetta/Desktop/UNICAL/Submitted_IWL/LAST_SIMU_26112014/UnicalMyBasin4/soil/soilMec0001.txt";
			fs.pathToSoilFile = "/Users/giuseppeformetta/Desktop/UNICAL/Submitted_IWL/LAST_SIMU_26112014/UnicalMyBasin4/soil/soil0001.txt";
			fs.nLayers = 8;
			fs.time = i;
			fs.pathToSlopeAngleMap = "/Users/giuseppeformetta/Desktop/UNICAL/Submitted_IWL/LAST_SIMU_26112014/UnicalMyBasin4/input_maps/slope_mybasin_4.asc";
			fs.pathToGEOtopOutputMaps = "/Users/giuseppeformetta/Desktop/UNICAL/Submitted_IWL/LAST_SIMU_26112014/UnicalMyBasin4/output_maps";
			fs.headerPsiGEOtop = "pressureLIQ";
			fs.headerThetaGEOtop = "thetaliq";
			fs.process();
		}

	}

}
