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

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeFS_GSE;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.shalstab.OmsShalstab;
import org.jgrasstools.hortonmachine.modules.statistics.ROC;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test the {@link OmsShalstab} module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class Test_M2_HSS_LUGLIO2014 extends HMTestCase {
	public void testShalstab() throws Exception {

		String period = "CAL";
		String model = "M2";
		// type= DISTR or UNIF
		String type = "UNIF";
		String index = "HSS";
		String size = "30";
		String all = model + "_" + type + "_" + period + "_" + index + "_"
				+ size;

		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset"
				+ size + "m/cell/slope_" + period + "_" + size;

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		// String pathToTCA =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca";
		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset"
				+ size + "m/cell/tca_" + period + "_" + size;

		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		String pathToDem = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset"
				+ size + "m/cell/pit_" + period + "_" + size;
		reader = new OmsRasterReader();
		reader.file = pathToDem;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D dem = reader.outRaster;

		int nsample = 200;
		 String pPath =
		 "/Users/giuseppeformetta/Desktop/Sensitivity/M2_HSS_Fi";
		 double par[] = UNIFORM(0.2, 0.7, nsample);

//		 String pPath =
//		 "/Users/giuseppeformetta/Desktop/Sensitivity/M2_HSS_Q";
//		 double par[] = UNIFORM(10, 200, nsample);
		//
//		 String pPath =
//		 "/Users/giuseppeformetta/Desktop/Sensitivity/M2_HSS_SD";
//		 double par[] = UNIFORM(0.8, 5.0, nsample);

//		 String pPath =
//		 "/Users/giuseppeformetta/Desktop/Sensitivity/M2_HSS_T";
//		 double par[] = UNIFORM(10, 150.0, nsample);
		//
//		String pPath = "/Users/giuseppeformetta/Desktop/Sensitivity/M2_HSS_C";
//		double par[] = UNIFORM(0, 50000.0, nsample);
		FileWriter Rstatfile = new FileWriter(pPath);
		PrintWriter errestat = new PrintWriter(Rstatfile);
		for (int i = 0; i < par.length; i++) {
			OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
			shalstab.inSlope = slopeCoverage;
			shalstab.pThresholdFS = 1.0;
			shalstab.inTca = abCoverage;
			shalstab.inDem = dem;
			shalstab.pTrasmissivity = 33.22048342651098;
			shalstab.pCohesion =  49626.62822943622;
			shalstab.pSoilPorosity = 0.5;
			shalstab.pTgphi =par[i];// 0.7914637826086255;
			shalstab.pQ = 293.4415367897272;
			shalstab.pDegreeOfSat = 0.5;
			shalstab.pSdepth = 1.676043690869428;
			shalstab.pGs = 2.3991479957129926;

			shalstab.pMode = 0;
			// Per ROSSO
			// shalstab.pMode = 1;
			// shalstab.pDuration = xx[numpart][8];
			shalstab.process();

			GridCoverage2D classiCoverage = shalstab.outShalstab;
			OmsRasterWriter writer = new OmsRasterWriter();
			writer.inRaster = classiCoverage;
			writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset"
					+ size + "m/cell/Classi_" + all;
			writer.process();

			OmsRasterReader reader2 = new OmsRasterReader();
			reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset"
					+ size + "m/cell/LAST_FRANE_CUT_" + period;
			reader2.fileNovalue = -9999.0;
			reader2.geodataNovalue = Double.NaN;
			reader2.process();
			GridCoverage2D M2 = reader2.outRaster;

			ROC ab = new ROC();
			ab.flagFalseModelled = 1;
			ab.flagTrueModelled = 0;
			ab.flagNegativeMeas = 1;
			ab.flagPositiveMeas = 0;
			ab.inMeasuredMap = M2;
			ab.inModelledMap = shalstab.outShalstab;
			ab.pPath = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/RESULTS/"
					+ all;

			ab.process();

			classiCoverage = ab.outClass;
			writer = new OmsRasterWriter();
			writer.inRaster = classiCoverage;
			writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset"
					+ size + "m/cell/Classi_ROC_" + all;
			writer.process();

			errestat.println(par[i] + "," + ab.outHss);

		}
		Rstatfile.close();

	}

	public static double[] UNIFORM(double xmin, double xmax, int nsample) {
		// Latin Hypercube sampling
		// double[][] LHSresult=new double [1][1];
		double[] s = new double[nsample];

		double[] ran = new double[nsample];
		for (int row = 0; row < ran.length; row++) {

			s[row] = ((xmax - xmin)) * Math.random() + xmin;
		}

		return s;
	}
}