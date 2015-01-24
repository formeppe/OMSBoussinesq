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
public class Test_PS_ShalstabGppe_DISTRIBUTED extends HMTestCase {
	public void testShalstab() throws Exception {

		// String pathToSlope =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_slope";
		String pathToSlope = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_SLOPE_CAL";

		OmsRasterReader reader = new OmsRasterReader();
		reader.file = pathToSlope;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slopeCoverage = reader.outRaster;

		// String pathToTCA =
		// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_tca";
		String pathToTCA = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/last_TCA_CAL";
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

		PS p = new PS();
		p.slope = slopeCoverage;
		p.sdepth = sDepthCoverage;
		p.porosityMap = porosityMap;
		p.cohesionMap = CohesionMap;
		p.hydraulicCondictivityMap = ConductivityMap;
		p.frictionAngleMap = FrictionAngleMap;
		p.tca = abCoverage;

		// OmsShalstab shalstab = new OmsShalstab();
		// shalstab.inSlope = slope;
		// shalstab.inTca = tca;
		// shalstab.inSdepth = s.outSoilDepth;
		// shalstab.pTrasmissivity = xx[numpart][0];
		// shalstab.pCohesion = 0;
		// shalstab.pRho = xx[numpart][1];
		// shalstab.inTgphi = frictionAngle;
		// shalstab.pQ = xx[numpart][2];
		// shalstab.pm = pm;

		double[] pmin = { 5, 1.6, 10, 0.45, 1.5, 0.1 };
		double[] pmax = { 100, 2.5, 600, 0.6, 5.0, 1.5 };

		p.ParRange_minn = pmin;
		p.ParRange_maxn = pmax;

		p.ModelName = "SHALSTABDISTRIBUTED";

		p.kmax = 150000;
		p.p = 20;
		p.parameters = 5;

		p.process();

	}

}
