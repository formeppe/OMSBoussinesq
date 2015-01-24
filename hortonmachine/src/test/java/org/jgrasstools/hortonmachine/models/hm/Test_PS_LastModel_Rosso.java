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
public class Test_PS_LastModel_Rosso extends HMTestCase {
	public void testShalstab() throws Exception {

		// String pathToSlope =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_slope";
		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/slope_CAL_30";

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		// String pathToTCA =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca";
		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/tca_CAL_30";
		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		String pathToDem = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/pit_CAL_30";
		reader = new OmsRasterReader();
		reader.file = pathToDem;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D dem = reader.outRaster;

		OmsRasterReader reader2 = new OmsRasterReader();
		reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/LAST_FRANE_CUT_CAL";
		reader2.fileNovalue = -9999.0;
		reader2.geodataNovalue = Double.NaN;
		reader2.process();
		GridCoverage2D M2 = reader2.outRaster;

		PS p = new PS();
		p.inDem = dem;
		p.slope = slopeCoverage;
		p.tca = abCoverage;
		p.ModelName = "LASTMODEL_ROSSO_Par_CONST";

//		shalstab.pSdepth = xx[numpart][5];
//		shalstab.pTrasmissivity = xx[numpart][0];
//		shalstab.pCohesion = xx[numpart][1];
//		shalstab.pSoilPorosity = 0.5;
//		shalstab.pTgphi = xx[numpart][2];
//		shalstab.pQ = xx[numpart][3];
//		shalstab.pDegreeOfSat =0.5;
//		shalstab.pGs =  xx[numpart][4];
//		shalstab.pDuration = xx[numpart][6];

		double[] pmin = new double[] { 10.0, 100.0,  0.2, 50.0, 2.0, 1.0,0.1 };
		double[] pmax = new double[] { 150.0, 50000.0,  0.8, 300,2.8, 3.0, 2.0 };
		// double[] pmin = new double[] { 5.0, 0.0, 0.1, 0.1, 50.0, 0.1, 0.5};
		// double[] pmax = new double[] { 150.0, 80000.0, 0.6, 0.8, 550, 0.90,
		// 4.0 };
		//
		p.ParRange_minn = pmin;
		p.ParRange_maxn = pmax;
		p.kmax = 30;
		p.p = 50;
		p.parameters = 8;
		p.inMeasuredForRoc = M2;
		p.pPath = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/ResultPSO";

		p.process();

	}

}
