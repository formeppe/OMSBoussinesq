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
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeFS;
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
public class TestComputeFSGppe extends HMTestCase {
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

		String pathToSoilDepthMap = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/SpessoreDiSuolo_Cal";
		reader = new OmsRasterReader();
		reader.file = pathToSoilDepthMap;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();

		String pathResult = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/NEWROC3.csv";
		FileWriter Rstatfile = new FileWriter(pathResult);
		PrintWriter errestat = new PrintWriter(Rstatfile);
		// Dall curva roc[10.050682549183866, 609.614198619692,
		// 1.6151833097103352, 0.21055910067190148, 199.99350685013553]
		
		// Dall ottimizzazione di SUccessIndex= 10.073383125273232
		// 292.2189102335744 1.829297572753146 0.2895850461058376
		// 199.85830906375955 g_best_value= 0.40346709351310905
		
		// Dall ottimizzazione Math.pow(((1 - ab.outTruePositiveRate)* (1 -
		// ab.outTruePositiveRate) + (ab.outFalsePositiveRate *
		// ab.outFalsePositiveRate)),0.5); 10.183164574555931 ;75.4882148347961;1.6952745365050945;0.3954484275197825 ; 191.9506397533307 ;
		
		//Dall ottimizzazione di CSI 10.337769723306726  621.5265191848322  1.600202688940084  0.2559842822384563  194.22105422117775  g_best_value= 0.8025049166752924

		
		
		int nsample = 1000;
		double par[] = UNIFORM(0.1, 2.5, nsample);
		for (int i = 0; i < nsample; i++) {

			GridCoverage2D sDepthCoverage = reader.outRaster;

			OmsComputeFS shalstab = new OmsComputeFS();
			shalstab.inSlope = slopeCoverage;
			shalstab.inTca = abCoverage;
			shalstab.inSdepth = sDepthCoverage;
			shalstab.pTrasmissivity = 10.337769723306726;//10.073383125273232;//10.050682549183866;// 10.183164574555931 ;
			shalstab.pCohesion =621.5265191848322 ;// 292.2189102335744;//609.614198619692;// 75.4882148347961;
			shalstab.pRho = 1.600202688940084;// 1.829297572753146;//1.6151833097103352;// 1.6952745365050945;
			shalstab.pTgphi = 0.2559842822384563;//0.2895850461058376;//0.21055910067190148;// 0.3954484275197825 ;//0.46;
			shalstab.pQ = 194.22105422117775;//199.85830906375955 ;//199.99350685013553;// 191.9506397533307 ;
			shalstab.pThresholdFS = par[i];
			shalstab.pm = pm;

			shalstab.process();

			GridCoverage2D classiCoverage = shalstab.outShalstab;
			OmsRasterWriter writer = new OmsRasterWriter();
			writer.inRaster = classiCoverage;
			writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_Compute_FS_CAL";
			writer.process();

			OmsRasterReader reader2 = new OmsRasterReader();
			reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/MeasuredForRoc_Cal";
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
			ab.inModelledMap = classiCoverage;
			ab.process();

			classiCoverage = ab.outClass;
			writer = new OmsRasterWriter();
			writer.inRaster = classiCoverage;
			writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/ClassiROC_ComputeFS_CAL";
			writer.process();

			errestat.println(par[i] + "," + ab.outAccuracy + ","
					+ ab.outPrecision + "," + ab.outTruePositiveRate + ","
					+ ab.outFalsePositiveRate);

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