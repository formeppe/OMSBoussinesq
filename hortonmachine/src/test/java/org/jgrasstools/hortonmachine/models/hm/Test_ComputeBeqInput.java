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
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.boussinesq.ComputeBeqInput;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.DoGridMesh;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.shalstab.OmsShalstab;
import org.jgrasstools.hortonmachine.modules.statistics.PS;
import org.jgrasstools.hortonmachine.modules.statistics.ROC;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test the {@link OmsShalstab} module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class Test_ComputeBeqInput extends HMTestCase {
	public void testShalstab() throws Exception {

		// String pathToSlope =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_slope";
//		String pathToSlope = "/Users/giuseppeformetta/Desktop/testMesh/SmallExample.asc";
		String pathToSlope = "/Users/giuseppeformetta/Desktop/testMesh/pitLW100.asc";

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		ComputeBeqInput c = new ComputeBeqInput();
		c.inDtm = slopeCoverage;
		c.pPath = "/Users/giuseppeformetta/Desktop/BEQ";

		c.valueEtaInitialCondVtoAddtoDem = 0.1;
		c.valueCV = 10.0;
		c.valueMV = 0.3;
		c.valueEtaDricheletV = -9999;
		c.valuePorosityV = 0.4;
		c.valueSourceV = 0.1;
		c.process();

		// OmsRasterWriter writer = new OmsRasterWriter();
		// writer.inRaster = g.outShalstab;
		// writer.file =
		// "/Users/giuseppeformetta/Desktop/testMesh/pitLWNUMBERED.asc";
		// writer.process();

	}

}
