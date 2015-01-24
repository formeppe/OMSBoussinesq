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
public class TestFS_LU_26112014 extends HMTestCase {

	public void testAb() throws Exception {

		for (int i = 1; i < 57; i++) {
			OmsFs16042014 fs = new OmsFs16042014();
			// fs.pathToSoilMecFile =
			// "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin2/soil/soilMec0001MOD.txt";
			// fs.pathToSoilFile =
			// "/Users/giuseppeformetta/GEOtop_211212/oms3.prj.geotop/data/UnicalMyBasin2/soil/soil0001.txt";
			/////////////CONCAVE//////////////////////
			fs.pathToSoilMecFile = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/geotop_SOIL/concave/soil/soilMec0001.txt";
			fs.pathToSoilFile = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/geotop_SOIL/concave/soil/soil0001.txt";
			fs.nLayers = 40;
			fs.time = i;
			fs.pathToSlopeAngleMap = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/geotop_SOIL/concave/input_maps/slope.asc";
			fs.pathToGEOtopOutputMaps = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/geotop_SOIL/concave/output_maps";
			fs.headerPsiGEOtop = "pressureLIQ";
			fs.headerThetaGEOtop = "thetaliq";
			//////////////////////////////////////////
			
			/////////////PLANAR//////////////////////

//			fs.pathToSoilMecFile = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/planar/soil/soilMec0001.txt";
//			fs.pathToSoilFile = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/planar/soil/soil0001.txt";
//			fs.nLayers = 40;
//			fs.time = i;
//			fs.pathToSlopeAngleMap = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/planar/input_maps/slope.asc";
//			fs.pathToGEOtopOutputMaps = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/planar/output_maps";
//			fs.headerPsiGEOtop = "pressureLIQ";
//			fs.headerThetaGEOtop = "thetaliq";
			
			/////////////CONVEX//////////////////////

//			
//			fs.pathToSoilMecFile = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/convex/soil/soilMec0001.txt";
//			fs.pathToSoilFile = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/convex/soil/soil0001.txt";
//			fs.nLayers = 40;
//			fs.time = i;
//			fs.pathToSlopeAngleMap = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/convex/input_maps/slope.asc";
//			fs.pathToGEOtopOutputMaps = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/WORKLU/convex/output_maps";
//			fs.headerPsiGEOtop = "pressureLIQ";
//			fs.headerThetaGEOtop = "thetaliq";
		fs.process();
		}

	}

}
