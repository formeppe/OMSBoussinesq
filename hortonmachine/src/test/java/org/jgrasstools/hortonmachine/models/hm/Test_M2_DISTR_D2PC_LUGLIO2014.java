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
public class Test_M2_DISTR_D2PC_LUGLIO2014 extends HMTestCase {
	public void testShalstab() throws Exception {

		String period = "VAL";
		String model = "M2";
		// type= DISTR or UNIF
		String type = "DISTR";
		String index = "D2PC";
		String all = model + "_" + type + "_" + period + "_" + index;

		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_SLOPE_"
				+ period;

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		// String pathToTCA =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca";
		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_TCA_"
				+ period;
		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		OmsComputeSoilThickness s = new OmsComputeSoilThickness();
		s.doP1 = true;
		s.p1 = 2.1908907498688155;
		s.p2 = 0.4952446225846152;
		s.inSlope = slopeCoverage;
		s.process();

		GridCoverage2D classiCoverage = s.outSoilDepth;
		OmsRasterWriter writer = new OmsRasterWriter();
		writer.inRaster = classiCoverage;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
		writer.process();

		OmsRasterReader reader33 = new OmsRasterReader();
		reader33.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/NewA3TESS_"
				+ period;
		reader33.fileNovalue = -9999.0;
		reader33.geodataNovalue = Double.NaN;
		reader33.process();
		GridCoverage2D tessitura = reader33.outRaster;

		OmsRangeLookup range = new OmsRangeLookup();
		range.pm = pm;
		range.inRaster = tessitura;
		range.pRanges = "[0 1],(1.0 2],(2.0 3]";
		range.pClasses = 9683.892896806026 + "," + 9683.892896806026 * 2.0
				+ "," + 9683.892896806026 * 3.5;
		range.process();
		GridCoverage2D cohesion = range.outRaster;

		writer.inRaster = cohesion;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/cohesion";
		writer.process();

		range = new OmsRangeLookup();
		range.pm = pm;
		range.inRaster = tessitura;
		range.pRanges = "[0 1],(1.0 2],(2.0 3]";
		range.pClasses = 0.5963442503336484 + "," + 0.5963442503336484 * 0.85
				+ "," + 0.5963442503336484 * 0.7;
		range.process();
		GridCoverage2D frictionAngle = range.outRaster;
		writer.inRaster = frictionAngle;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/frictionAngle";
		writer.process();

		range = new OmsRangeLookup();
		range.pm = pm;
		range.inRaster = tessitura;
		range.pRanges = "[0 1],(1.0 2],(2.0 3]";
		range.pClasses = 96.90476182316314 + "," + 96.90476182316314 * 0.001
				+ "," + 96.90476182316314 * 0.0000001;
		range.process();
		GridCoverage2D conductivity = range.outRaster;

		writer.inRaster = conductivity;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/conductivity";
		writer.process();

		String pathToDem = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_PIT_"
				+ period;
		reader = new OmsRasterReader();
		reader.file = pathToDem;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D dem = reader.outRaster;

		reader33 = new OmsRasterReader();
		reader33.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/VoidRatioMap_A3_"
				+ period;
		reader33.fileNovalue = -9999.0;
		reader33.geodataNovalue = Double.NaN;
		reader33.process();
		GridCoverage2D porosityMap = reader33.outRaster;

		OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
		shalstab.inSlope = slopeCoverage;
		shalstab.pThresholdFS = 1.0;
		shalstab.inTca = abCoverage;
		shalstab.inDem = dem;

		shalstab.inSdepth = s.outSoilDepth;
		shalstab.inTrasmissivity = conductivity;
		shalstab.computeTransmissivity = true;
		shalstab.inCohesion = cohesion;
		// shalstab.pCohesion = 0.0;
		shalstab.inTgphi = frictionAngle;
		shalstab.inSoilPorosity = porosityMap;
		shalstab.pQ = 468.26125807289833;
		shalstab.pDegreeOfSat = 0.3475473814662574;
		shalstab.pGs = 2.65;

		shalstab.pMode = 0;
		// Per ROSSO
		// shalstab.pMode = 1;
		// shalstab.pDuration = xx[numpart][8];
		shalstab.process();

		classiCoverage = shalstab.outShalstab;
		writer = new OmsRasterWriter();
		writer.inRaster = classiCoverage;
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_"
				+ all;
		writer.process();

		OmsRasterReader reader2 = new OmsRasterReader();
		reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_FraneMisurate_"
				+ period;
		reader2.fileNovalue = -9999.0;
		reader2.geodataNovalue = Double.NaN;
		reader2.process();
		GridCoverage2D M2 = reader2.outRaster;

		ROC ab = new ROC();
		ab.pPath = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/RisultatiLuglio2014/"
				+ all;
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
		writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_ROC_"
				+ all;
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