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
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.fs.OmsFs01072014;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test ab.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestFS_01072014 extends HMTestCase {

	public void testAb() throws Exception {

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = "/Users/giuseppeformetta/Desktop/UNICAL/sito7slope.asc";
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slope = reader.outRaster;

		OmsFs01072014 fs = new OmsFs01072014();
		fs.pathToSoilMecFile = "/Users/giuseppeformetta/Desktop/simOKDRAKE/UnicalSITO7/soil/soilMec0001.txt";
		fs.pathToSoilFile = "/Users/giuseppeformetta/Desktop/simOKDRAKE/UnicalSITO7/soil/soil0001.txt";
		fs.nLayers = 15;
		fs.inSlope = slope;
		fs.pathToGEOtopOutputMaps="/Users/giuseppeformetta/Desktop/simOKDRAKE/UnicalSITO7/testFolder";
		fs.headerPsiGEOtop="pressureLIQ";
		fs.headerThetaGEOtop="thetaliq";
		fs.process();

	}

}
