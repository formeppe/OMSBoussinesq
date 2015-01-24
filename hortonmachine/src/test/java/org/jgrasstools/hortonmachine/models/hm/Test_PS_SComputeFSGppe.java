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
public class Test_PS_SComputeFSGppe extends HMTestCase {
	public void testShalstab() throws Exception {

		// String pathToSlope =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_slope";
		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_slope_Cal";

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		// String pathToTCA =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca";
		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca_Cal";
		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		// String pathToSoilDepthMap =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/SpessoreDiSuolo2";
		String pathToSoilDepthMap = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/SpessoreDiSuolo_Cal";

		reader = new OmsRasterReader();
		reader.file = pathToSoilDepthMap;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D sDepthCoverage = reader.outRaster;
		OmsRasterReader reader2 = new OmsRasterReader();
		reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/MeasuredForRoc";
		reader2.fileNovalue = -9999.0;
		reader2.geodataNovalue = Double.NaN;
		reader2.process();
		GridCoverage2D M2 = reader2.outRaster;

		PS p = new PS();
		p.slope = slopeCoverage;
		p.sdepth = sDepthCoverage;
		p.tca = abCoverage;
		p.inMeasuredForRoc = M2;

		// shalstab.pTrasmissivity = xx[0][numpart];
		// shalstab.pCohesion =xx[1][numpart] ;
		// shalstab.pRho = xx[2][numpart];
		// shalstab.pTgphi = xx[3][numpart];
		// shalstab.pQ = xx[4][numpart];
		double[] pmin = { 5.0, 0.0, 1.6, 0.1, 50, 1.0, 0.1 };
		double[] pmax = { 100.0, 1000, 2.5, 0.9, 200, 5.0, 3.0 };

		p.ModelName = "ComputeFS";
		p.ModelName = "ComputeFSconGamma";
		// p.ModelName = "ComputeFSROC";
		// p.ModelName = "ComputeFS_Success_Index";
		// p.ModelName = "ComputeFS_Critical_Success_Index";
		pmin = new double[] { 5.0, 0.0, 10.6, 16.0, 0.2, 50.0, .5, 0.1 };
		pmax = new double[] { 100.0, 15.0, 18.0, 25.0, 0.9, 200, 3.0, 2.8 };
		p.ParRange_minn = pmin;
		p.ParRange_maxn = pmax;
		p.kmax = 150000;
		p.p = 20;
		p.parameters = 8;

		p.process();

	}

}
