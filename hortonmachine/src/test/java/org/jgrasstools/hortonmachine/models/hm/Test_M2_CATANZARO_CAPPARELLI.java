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
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeFS_GSE_PROBABILISTIC;
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
public class Test_M2_CATANZARO_CAPPARELLI extends HMTestCase {
	public void testShalstab() throws Exception {

		String period = "VAL";
		String model = "M2";
		// type= DISTR or UNIF
		String type = "UNIF";
		String index = "ACC";
		String size = "30";
		String all = model + "_" + type + "_" + period + "_" + index + "_"
				+ size;

		String pathToSlope = "/Users/giuseppeformetta/Documents/catanzaro/m/cell/slope";

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		// String pathToTCA =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca";
		String pathToTCA = "/Users/giuseppeformetta/Documents/catanzaro/m/cell/tca";
		reader = new OmsRasterReader();
		reader.file = pathToTCA;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D abCoverage = reader.outRaster;

		String pathToDem = "/Users/giuseppeformetta/Documents/catanzaro/m/cell/pit";
		reader = new OmsRasterReader();
		reader.file = pathToDem;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D dem = reader.outRaster;

		// String pathToSoilDepthMap =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/SpessoreDiSuolo_Cal";
		//
		// reader = new OmsRasterReader();
		// reader.file = pathToSoilDepthMap;
		// reader.fileNovalue = -9999.0;
		// reader.geodataNovalue = Double.NaN;
		// reader.process();
		// GridCoverage2D sDepthCoverage = reader.outRaster;

		OmsComputeSoilThickness c = new OmsComputeSoilThickness();
		c.doP1 = true;
		c.inSlope = slopeCoverage;
		c.p1 = 2.0;
		c.p2 = 1.8;
		c.process();
		GridCoverage2D classiCoverage = c.outSoilDepth;
		OmsRasterWriter writer = new OmsRasterWriter();
		writer.inRaster = classiCoverage;
		writer.file = "/Users/giuseppeformetta/Documents/catanzaro/m/cell/ModelSoilThickness";
		writer.process();

		OmsComputeFS_GSE_PROBABILISTIC shalstab = new OmsComputeFS_GSE_PROBABILISTIC();
		shalstab.inSlope = slopeCoverage;
		shalstab.doProbabilistic = true;
		shalstab.inParVectMin = new double[] {0.4,0,2.3,0.3,0.2,50,80};
		shalstab.inParVectMax = new double[] {0.75,40000,2.7,0.3,0.2,150,150};

		shalstab.pThresholdFS = 0.9;
		shalstab.inTca = abCoverage;
		shalstab.inDem = dem;
		shalstab.pTrasmissivity = 50;
		shalstab.pCohesion = 0;
		shalstab.pSoilPorosity = 0.4;
		shalstab.pTgphi = 0.5317094;
		shalstab.pQ = 200;
		shalstab.pDegreeOfSat = 0.1;
		shalstab.inSdepth = c.outSoilDepth;
		// shalstab.pSdepth = 3;

		shalstab.pGs = 2.2;

		shalstab.pMode = 0;

		shalstab.process();

		classiCoverage = shalstab.outShalstab;
		writer = new OmsRasterWriter();
		writer.inRaster = classiCoverage;
		writer.file = "/Users/giuseppeformetta/Documents/catanzaro/m/cell/ModelResults";
		writer.process();

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