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
public class TestShalstabGppe_CAL_LUGLIO2014 extends HMTestCase {
	public void testShalstab() throws Exception {

		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_SLOPE_VAL";

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		// String pathToTCA =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca";
		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_TCA_VAL";
		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		int nsample = 200;

		double par[] = UNIFORM(0.1, 0.8, nsample);

//		String pathResult = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/ChanginFi.csv";
//		FileWriter Rstatfile = new FileWriter(pathResult);
//		PrintWriter errestat = new PrintWriter(Rstatfile);

		// PARAMETRI OK
		// shalstab.pTrasmissivity = 10;
		// shalstab.pCohesion = 580.0;
		// shalstab.pRho = 1.61;
		// shalstab.pTgphi =0.46;
		// shalstab.pQ = 420;

		//for (int i = 0; i < nsample; i++) {
			OmsShalstab shalstab = new OmsShalstab();
			shalstab.inSlope = slopeCoverage;
			shalstab.inTca = abCoverage;
			shalstab.pSdepth = 4.897128626197997;
			shalstab.pTrasmissivity =5.81563822532721;
			shalstab.pCohesion = 0;
			shalstab.pRho = 1.6922920263162284;
			shalstab.pTgphi = 0.5568543372089303;
			shalstab.pQ = 436.0466317593652;

			shalstab.process();

			GridCoverage2D classiCoverage = shalstab.outShalstab2Classes;
			OmsRasterWriter writer = new OmsRasterWriter();
			writer.inRaster = classiCoverage;
			writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_Val_Shalstab_omogeneo";
			writer.process();

			OmsRasterReader reader2 = new OmsRasterReader();
			reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_FraneMisurate_VAL";
			reader2.fileNovalue = -9999.0;
			reader2.geodataNovalue = Double.NaN;
			reader2.process();
			GridCoverage2D M2 = reader2.outRaster;

			ROC ab = new ROC();
			ab.pPath="/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/RisultatiLuglio2014/shalstab_Val";
			ab.flagFalseModelled = 1;
			ab.flagTrueModelled = 0;
			ab.flagNegativeMeas = 1;
			ab.flagPositiveMeas = 0;
			ab.inMeasuredMap = M2;
			ab.inModelledMap = shalstab.outShalstab2Classes;
			ab.process();
			
			classiCoverage = ab.outClass;
			writer = new OmsRasterWriter();
			writer.inRaster = classiCoverage;
			writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_Roc_Val_Shalstab_omogeneo";
			writer.process();

//			errestat.println(par[i] + "," + ab.outAccuracy + ","
//					+ ab.outPrecision + "," + ab.outTruePositiveRate + ","
//					+ ab.outFalsePositiveRate);
//
//		}
//		Rstatfile.close();

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