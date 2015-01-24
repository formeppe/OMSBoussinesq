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
import org.jgrasstools.gears.modules.r.rangelookup.OmsRangeLookup;
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
public class TestROC_SHALSTAB_FEB2014 extends HMTestCase {
	public void testShalstab() throws Exception {

		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_SLOPE_CAL";
		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_TCA_CAL";
		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		String pathToDem = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_PIT_CAL";
		reader = new OmsRasterReader();
		reader.file = pathToDem;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D dem = reader.outRaster;

		OmsRasterReader reader2 = new OmsRasterReader();
		reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/NewForRoc_CAL_OK";
		reader2.fileNovalue = -9999.0;
		reader2.geodataNovalue = Double.NaN;
		reader2.process();
		GridCoverage2D M2 = reader2.outRaster;

		String pathResult = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/ChanginFi.csv";
		FileWriter Rstatfile = new FileWriter(pathResult);
		PrintWriter errestat = new PrintWriter(Rstatfile);

		// PARAMETRI OK
		// shalstab.pTrasmissivity = 10;
		// shalstab.pCohesion = 580.0;
		// shalstab.pRho = 1.61;
		// shalstab.pTgphi =0.46;
		// shalstab.pQ = 420;
		int nsample = 10;

		double par1[] = UNIFORM(0.6, 0.85, nsample);
		double par2[] = UNIFORM(0.5, 0.75, nsample);
		double par3[] = UNIFORM(0.3, 0.65, nsample);
		
		OmsRasterReader reader33 = new OmsRasterReader();
		reader33.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Tess_A3_CAL";
		reader33.fileNovalue = -9999.0;
		reader33.geodataNovalue = Double.NaN;
		reader33.process();
		GridCoverage2D tessitura = reader33.outRaster;
		
		for (int i = 0; i < nsample; i++) {
			for (int j = 0; j < nsample; j++) {

				for (int k = 0; k < nsample; k++) {

					OmsRangeLookup range = new OmsRangeLookup();

					range = new OmsRangeLookup();
					range.pm = pm;
					range.inRaster = tessitura;
					range.pRanges = "[0 1],(1.0 2],(2.0 3]";
					range.pClasses = par1[i] + "," + par2[j] + "," + par3[k];
					range.process();

					GridCoverage2D classiCoverage = range.outRaster;
					OmsRasterWriter writer = new OmsRasterWriter();
					writer.inRaster = classiCoverage;
					writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/tanphi";
					writer.process();

					
					OmsShalstab shalstab = new OmsShalstab();
					shalstab.inSlope = slopeCoverage;
					shalstab.inTca = abCoverage;
					shalstab.pSdepth = 3.32;
					shalstab.pTrasmissivity = 7.81;
					shalstab.pCohesion = 0;
					shalstab.pRho = 1.75;
					//shalstab.pTgphi = 0.57;
					shalstab.inTgphi=range.outRaster;
					shalstab.pQ = 550;
					shalstab.pm = pm;

					shalstab.process();

					ROC ab = new ROC();
					ab.flagFalseModelled = 1;
					ab.flagTrueModelled = 0;
					ab.flagNegativeMeas = 1;
					ab.flagPositiveMeas = 0;
					ab.inMeasuredMap = M2;
					ab.inModelledMap = shalstab.outShalstab2Classes;
					ab.process();

					errestat.println(par1[i] + "," + par2[j] + "," + par3[k]
							+ "," + ab.outTruePositiveRate + ","
							+ ab.outFalsePositiveRate);

				}
			}
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