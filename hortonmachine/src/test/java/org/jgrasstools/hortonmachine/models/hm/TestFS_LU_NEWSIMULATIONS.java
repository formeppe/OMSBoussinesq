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


import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.fs.OmsFs16042014;
import org.jgrasstools.hortonmachine.utils.HMTestCase;


/**
 * Test ab.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestFS_LU_NEWSIMULATIONS extends HMTestCase {

	public void testAb() throws Exception {

		//String pathToSimulationOutput = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I3/DRY/";
		//String pathToSimulationOutput = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I3/MEDIUM/";
		String pathToSimulationOutput = "/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I3/WET/";

		for (int i = 1; i < 58; i++) {
			OmsFs16042014 fs = new OmsFs16042014();
			fs.pathToSoilMecFile = pathToSimulationOutput
					+ "concave/soil/soilMec0001.txt";
			fs.pathToSoilFile = pathToSimulationOutput
					+ "concave/soil/soil0001.txt";
			fs.nLayers = 30;
			fs.time = i;
			fs.pathToSlopeAngleMap = pathToSimulationOutput
					+ "concave/input_maps/slope.asc";
			fs.pathToGEOtopOutputMaps = pathToSimulationOutput
					+ "concave/output_maps";
			fs.headerPsiGEOtop = "pressureLIQ";
			fs.headerThetaGEOtop = "thetaliq";
			fs.process();
		}
		for (int i = 1; i < 58; i++) {
			OmsFs16042014 fs = new OmsFs16042014();
			fs.pathToSoilMecFile = pathToSimulationOutput
					+ "convex/soil/soilMec0001.txt";
			fs.pathToSoilFile = pathToSimulationOutput
					+ "convex/soil/soil0001.txt";
			fs.nLayers = 30;
			fs.time = i;
			fs.pathToSlopeAngleMap = pathToSimulationOutput
					+ "convex/input_maps/slope.asc";
			fs.pathToGEOtopOutputMaps = pathToSimulationOutput
					+ "convex/output_maps";
			fs.headerPsiGEOtop = "pressureLIQ";
			fs.headerThetaGEOtop = "thetaliq";
			fs.process();
		}
		for (int i = 1; i < 58; i++) {
			OmsFs16042014 fs = new OmsFs16042014();
			fs.pathToSoilMecFile = pathToSimulationOutput
					+ "planar/soil/soilMec0001.txt";
			fs.pathToSoilFile = pathToSimulationOutput
					+ "planar/soil/soil0001.txt";
			fs.nLayers = 30;
			fs.time = i;
			fs.pathToSlopeAngleMap = pathToSimulationOutput
					+ "planar/input_maps/slope.asc";
			fs.pathToGEOtopOutputMaps = pathToSimulationOutput
					+ "planar/output_maps";
			fs.headerPsiGEOtop = "pressureLIQ";
			fs.headerThetaGEOtop = "thetaliq";
			fs.process();
		}
		// ////////////////////////////////////////


	}

}
