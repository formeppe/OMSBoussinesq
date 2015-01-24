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
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeSoilThickness;
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
public class TestLASTMODELGppeROSSOCAL extends HMTestCase {
	public void testShalstab() throws Exception {

		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_slope_Cal";
		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca_Cal";
		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		String pathToDem = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/a3Dem_CAL";
		reader = new OmsRasterReader();
		reader.file = pathToDem;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D dem = reader.outRaster;

		String pathToSoilDepthMap = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/SpessoreDiSuolo_Cal";
		reader = new OmsRasterReader();
		reader.file = pathToSoilDepthMap;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D sDepthCoverage = reader.outRaster;

		int nsample = 200;

		double par[] = UNIFORM(0.1, 0.8, nsample);

		// 3.1596642919584963 815.752254485915 1.6177897547435263
		// 0.3384032596142335 115.04247705630965 g_best_value=
		// 0.5602893043737108

		// String pathResult =
		// "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/ChanginFi.csv";
		// FileWriter Rstatfile = new FileWriter(pathResult);
		// PrintWriter errestat = new PrintWriter(Rstatfile);

		// PARAMETRI OK
		// shalstab.pTrasmissivity = 10;
		// shalstab.pCohesion = 580.0;
		// shalstab.pRho = 1.61;
		// shalstab.pTgphi =0.46;
		// shalstab.pQ = 420;

		// for (int i = 0; i < nsample; i++) {

		OmsRasterReader reader2 = new OmsRasterReader();
		reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/MeasuredForRoc";
		reader2.fileNovalue = -9999.0;
		reader2.geodataNovalue = Double.NaN;
		reader2.process();
		GridCoverage2D M2 = reader2.outRaster;

		OmsComputeSoilThickness s = new OmsComputeSoilThickness();
		s.doP1 = true;
		s.p1 = 2.81;
		s.p2 = 2.29;
		s.inSlope = slopeCoverage;
		s.process();

		GridCoverage2D classiCoverage = s.outSoilDepth;
		OmsRasterWriter writer = new OmsRasterWriter();
		writer.inRaster = classiCoverage;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
		writer.process();

		OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
		shalstab.inSlope = slopeCoverage;
		shalstab.pThresholdFS = 1.0;
		shalstab.inTca = abCoverage;
		shalstab.pGs = 2.5;

		shalstab.inSdepth = s.outSoilDepth;
		shalstab.inDem = dem;
		shalstab.pTrasmissivity = 69.0;
		shalstab.pCohesion = 930.0;
		shalstab.pSoilPorosity = 0.83;
		shalstab.pTgphi = 0.6;
		shalstab.pQ = 198.0;
		shalstab.pDegreeOfSat = 0.84;
		shalstab.pMode = 1;
		shalstab.pDuration=2.0;
		shalstab.process();

		classiCoverage = shalstab.outShalstab;
		writer = new OmsRasterWriter();
		writer.inRaster = classiCoverage;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_SI_P1";
		writer.process();

		ROC ab = new ROC();
		ab.flagFalseModelled = 1;
		ab.flagTrueModelled = 0;
		ab.flagNegativeMeas = 1;
		ab.flagPositiveMeas = 0;
		ab.inMeasuredMap = M2;
		ab.inModelledMap = shalstab.outShalstab;
		ab.process();

		classiCoverage = ab.outClass;
		writer = new OmsRasterWriter();
		writer.inRaster = classiCoverage;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/ClassiROC_Cal";
		writer.process();

		// errestat.println(par[i] + "," + ab.outAccuracy + ","
		// + ab.outPrecision + "," + ab.outTruePositiveRate + ","
		// + ab.outFalsePositiveRate);
		//
		// }
		// Rstatfile.close();

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