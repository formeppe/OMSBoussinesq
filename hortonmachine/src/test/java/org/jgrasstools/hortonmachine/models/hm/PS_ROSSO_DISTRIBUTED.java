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
import org.jgrasstools.gears.modules.r.rangelookup.OmsRangeLookup;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.Hydrometers;
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
public class PS_ROSSO_DISTRIBUTED extends HMTestCase {
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
		reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/FINALECAL";
		reader2.fileNovalue = -9999.0;
		reader2.geodataNovalue = Double.NaN;
		reader2.process();
		GridCoverage2D M2 = reader2.outRaster;

		OmsRasterReader reader33 = new OmsRasterReader();
		reader33.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/porosity_CAL";
		reader33.fileNovalue = -9999.0;
		reader33.geodataNovalue = Double.NaN;
		reader33.process();
		GridCoverage2D porosityMap = reader33.outRaster;

		OmsRasterReader reader333 = new OmsRasterReader();
		reader333.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_k_m_al_g_CAL";
		reader333.fileNovalue = -9999.0;
		reader333.geodataNovalue = Double.NaN;
		reader333.process();
		GridCoverage2D ConductivityMap = reader333.outRaster;

		OmsRasterReader reader3333 = new OmsRasterReader();

		reader3333.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/cohesion_CAL";
		reader3333.fileNovalue = -9999.0;
		reader3333.geodataNovalue = Double.NaN;
		reader3333.process();
		GridCoverage2D CohesionMap = reader3333.outRaster;
		OmsRasterReader reader33333 = new OmsRasterReader();

		reader33333.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/tanFI_CAL";
		reader33333.fileNovalue = -9999.0;
		reader33333.geodataNovalue = Double.NaN;
		reader33333.process();
		GridCoverage2D FrictionAngleMap = reader33333.outRaster;

		// OmsRangeLookup range = new OmsRangeLookup();
		// range.pm = pm;
		// range.inRaster = tessitura;
		// range.pRanges = "[0 1],(1.0 2],(2.0 3]";
		// range.pClasses = "22000,65000,80000";
		// range.process();
		// GridCoverage2D out = range.outRaster;
		//
		// OmsRasterWriter writer = new OmsRasterWriter();
		// writer.inRaster = out;
		// writer.file =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/CoesionMap_A3_CAL";
		// writer.process();
		//
		// range = new OmsRangeLookup();
		// range.pm = pm;
		// range.inRaster = tessitura;
		// range.pRanges = "[0 1],(1.0 2],(2.0 3]";
		// range.pClasses = "0.5609821,0.5093498,0.4545605";
		// range.process();
		// out = range.outRaster;
		//
		// writer = new OmsRasterWriter();
		// writer.inRaster = out;
		// writer.file =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/FrictionAngleMap_A3_CAL";
		// writer.process();
		//
		// range = new OmsRangeLookup();
		// range.pm = pm;
		// range.inRaster = tessitura;
		// range.pRanges = "[0 1],(1.0 2],(2.0 3]";
		// range.pClasses = "6.0,3.0,0.2";
		// range.process();
		// out = range.outRaster;
		//
		// writer = new OmsRasterWriter();
		// writer.inRaster = out;
		// writer.file =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/HydraulicCondictivityMap_A3_CAL";
		// writer.process();
		// //
		//
		// range = new OmsRangeLookup();
		// range.pm = pm;
		// range.inRaster = tessitura;
		// range.pRanges = "[0 1],(1.0 2],(2.0 3]";
		// range.pClasses = "0.41,0.38,0.49";
		// range.process();
		// out = range.outRaster;
		//
		// writer = new OmsRasterWriter();
		// writer.inRaster = out;
		// writer.file =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/VoidRatioMap_A3_CAL";
		// writer.process();
		//

		PS p = new PS();
		p.slope = slopeCoverage;
		p.tca = abCoverage;
		p.inMeasuredForRoc = M2;
		p.hydraulicCondictivityMap = ConductivityMap;
		p.porosityMap = porosityMap;
		p.cohesionMap = CohesionMap;
		p.frictionAngleMap = FrictionAngleMap;
		p.inDem = dem;

		// shalstab.pQ = xx[numpart][0];
		// shalstab.pDegreeOfSat = xx[numpart][1];
		// shalstab.pGs = xx[numpart][2];

		double[] pmin = { 5, 10, 0.1, 0.1, 0, 0, 0, 0.5, 0.3, 0.3,
				0.5, 0.1 };
		double[] pmax = { 100.0, 500, 0.6, 3.0, 23000, 75000, 90000, 0.8, 0.7,
				0.6, 4.0, 2.0 };

		p.ParRange_minn = pmin;
		p.ParRange_maxn = pmax;

		p.ModelName = "ROSSODISTRIBUTED";

		// shalstab.pTrasmissivity = xx[numpart][0];
		// shalstab.pQ = xx[numpart][1];
		// shalstab.pDegreeOfSat = xx[numpart][2];

		p.kmax = 150000;
		p.p = 20;
		p.parameters = 8;

		p.process();

	}
}
