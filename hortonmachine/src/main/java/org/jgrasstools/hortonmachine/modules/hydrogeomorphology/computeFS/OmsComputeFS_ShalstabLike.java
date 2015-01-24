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
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS;

import static java.lang.Math.atan;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_AUTHORCONTACTS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_AUTHORNAMES;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_KEYWORDS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_LABEL;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_LICENSE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_NAME;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_STATUS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inCohesion_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inQ_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inRho_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inSdepth_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inSlope_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inTca_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inTgphi_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inTrasmissivity_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_outQcrit_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_outShalstab_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pCohesion_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pQ_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pRho_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pRock_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pSdepth_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pTgphi_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pTrasmissivity_DESCRIPTION;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.util.HashMap;

import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import oms3.annotations.Unit;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.coverage.ConstantRandomIter;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.geomorphology.ab.OmsAb;
import org.jgrasstools.hortonmachine.modules.geomorphology.curvatures.OmsCurvatures;
import org.jgrasstools.hortonmachine.modules.geomorphology.tca.OmsTca;

@Description(OMSSHALSTAB_DESCRIPTION)
@Author(name = OMSSHALSTAB_AUTHORNAMES, contact = OMSSHALSTAB_AUTHORCONTACTS)
@Keywords(OMSSHALSTAB_KEYWORDS)
@Label(OMSSHALSTAB_LABEL)
@Name(OMSSHALSTAB_NAME)
@Status(OMSSHALSTAB_STATUS)
@License(OMSSHALSTAB_LICENSE)
public class OmsComputeFS_ShalstabLike extends JGTModel {

	@Description(OMSSHALSTAB_inSlope_DESCRIPTION)
	@In
	public GridCoverage2D inSlope = null;

	@Description(OMSSHALSTAB_inTca_DESCRIPTION)
	@In
	public GridCoverage2D inDem = null;

	@Description(OMSSHALSTAB_inTca_DESCRIPTION)
	@In
	public GridCoverage2D inTca = null;

	@Description(OMSSHALSTAB_inTrasmissivity_DESCRIPTION)
	@Unit("m^2/day")
	@In
	public GridCoverage2D inTrasmissivity = null;

	@Description(OMSSHALSTAB_pTrasmissivity_DESCRIPTION)
	@Unit("m^2/day")
	@In
	public double pTrasmissivity = -1.0;

	@Description(OMSSHALSTAB_pTrasmissivity_DESCRIPTION)
	@Unit("m^2/day")
	@In
	public double pThresholdFS = -1.0;

	@Description(OMSSHALSTAB_inTgphi_DESCRIPTION)
	@In
	public GridCoverage2D inTgphi = null;

	@Description(OMSSHALSTAB_pTgphi_DESCRIPTION)
	@In
	public double pTgphi = -1.0;

	@Description(OMSSHALSTAB_inCohesion_DESCRIPTION)
	@Unit("Pa")
	@In
	public GridCoverage2D inCohesion = null;

	@Description(OMSSHALSTAB_inCohesion_DESCRIPTION)
	@Unit("Pa")
	@In
	public GridCoverage2D inRhoS = null;

	@Description(OMSSHALSTAB_inCohesion_DESCRIPTION)
	@Unit("Pa")
	@In
	public GridCoverage2D inGammaS = null;

	@Description(OMSSHALSTAB_pCohesion_DESCRIPTION)
	@Unit("Pa")
	@In
	public double pCohesion = -1.0;

	@Description(OMSSHALSTAB_inSdepth_DESCRIPTION)
	@Unit("m")
	@In
	public GridCoverage2D inSdepth = null;

	@Description(OMSSHALSTAB_pSdepth_DESCRIPTION)
	@Unit("m")
	@In
	public double pSdepth = -1.0;

	@Description(OMSSHALSTAB_inQ_DESCRIPTION)
	@Unit("mm/day")
	@In
	public GridCoverage2D inQ = null;

	@Description(OMSSHALSTAB_pQ_DESCRIPTION)
	@Unit("mm/day")
	@In
	public double pQ = -1.0;

	@Description(OMSSHALSTAB_pQ_DESCRIPTION)
	@Unit("mm/day")
	@In
	boolean computeTransmissivity = false;

	@Description(OMSSHALSTAB_inRho_DESCRIPTION)
	@In
	public GridCoverage2D inRho = null;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	public double pRho = -1.0;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	public double prhoS = -1.0;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	public double pGammaP = -1.0;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	public double pRhoW = -1.0;

	@Description(OMSSHALSTAB_pRock_DESCRIPTION)
	@In
	public double pRock = -9999.0;

	@Description(OMSSHALSTAB_outQcrit_DESCRIPTION)
	@Out
	public GridCoverage2D outQcrit = null;

	@Description(OMSSHALSTAB_outShalstab_DESCRIPTION)
	@Out
	public GridCoverage2D outShalstab = null;

	@Description(OMSSHALSTAB_outShalstab_DESCRIPTION)
	@Out
	public GridCoverage2D outShalstab2Classes = null;

	public final double EPS = 0.01;

	/**
	 * Value to be given to pixels if <code>h_s < eps</code>.
	 */
	private static final double ROCK = 8888.0;

	@Execute
	public void process() throws Exception {
		if (!concatOr(outShalstab == null, doReset)) {
			return;
		}
		checkNull(inSlope, inDem, inTca);

		if (pRock == -9999.0)
			pRock = 5.67;

		RenderedImage slopeRI = inSlope.getRenderedImage();

		OmsCurvatures c = new OmsCurvatures();
		c.inElev = inDem;
		c.process();

		OmsAb ab = new OmsAb();
		ab.inPlan = c.outPlan;
		ab.inTca = inTca;
		ab.process();

		RenderedImage abRI = ab.outAb.getRenderedImage();

		RandomIter trasmissivityIter = null;
		if (inTrasmissivity != null) {
			RenderedImage trasmissivityRI = inTrasmissivity.getRenderedImage();
			trasmissivityIter = RandomIterFactory.create(trasmissivityRI, null);
		} else {
			trasmissivityIter = new ConstantRandomIter(pTrasmissivity);
		}

		RandomIter tghiIter = null;
		if (inTgphi != null) {
			RenderedImage tgphiRI = inTgphi.getRenderedImage();
			tghiIter = RandomIterFactory.create(tgphiRI, null);
		} else {
			tghiIter = new ConstantRandomIter(pTgphi);
		}

		RandomIter cohesionIter = null;
		if (inCohesion != null) {
			RenderedImage cohesionRI = inCohesion.getRenderedImage();
			cohesionIter = RandomIterFactory.create(cohesionRI, null);
		} else {
			cohesionIter = new ConstantRandomIter(pCohesion);
		}

		RandomIter hsIter = null;
		if (inSdepth != null) {
			RenderedImage hsRI = inSdepth.getRenderedImage();
			hsIter = RandomIterFactory.create(hsRI, null);
		} else {
			hsIter = new ConstantRandomIter(pSdepth);
		}

		RandomIter qIter = null;
		if (inQ != null) {
			RenderedImage qRI = inQ.getRenderedImage();
			qIter = RandomIterFactory.create(qRI, null);
		} else {
			qIter = new ConstantRandomIter(pQ);
		}

		RandomIter rhoIter = null;
		if (inRho != null) {
			RenderedImage rhoRI = inRho.getRenderedImage();
			rhoIter = RandomIterFactory.create(rhoRI, null);
		} else {
			rhoIter = new ConstantRandomIter(pRho);
		}

		RandomIter rhoSIter = null;
		if (inRhoS != null) {
			RenderedImage rhoSRI = inRhoS.getRenderedImage();
			rhoSIter = RandomIterFactory.create(rhoSRI, null);
		} else {
			rhoSIter = new ConstantRandomIter(prhoS);
		}

		RandomIter rhoWIter = null;
		rhoWIter = new ConstantRandomIter(998.0);

		qcrit(slopeRI, abRI, trasmissivityIter, tghiIter, cohesionIter, hsIter,
				qIter, rhoIter, rhoSIter, rhoWIter);
	}

	/**
	 * Calculates the trasmissivity in every pixel of the map.
	 */
	private void qcrit(RenderedImage slope, RenderedImage ab,
			RandomIter trasmissivityRI, RandomIter frictionRI,
			RandomIter cohesionRI, RandomIter hsIter, RandomIter effectiveRI,
			RandomIter densityRI, RandomIter rhoSIter, RandomIter rhoWIter) {
		HashMap<String, Double> regionMap = CoverageUtilities
				.getRegionParamsFromGridCoverage(inSlope);
		int cols = regionMap.get(CoverageUtilities.COLS).intValue();
		int rows = regionMap.get(CoverageUtilities.ROWS).intValue();

		RandomIter slopeRI = RandomIterFactory.create(slope, null);
		RandomIter abRI = RandomIterFactory.create(ab, null);

		WritableRaster classiWR = CoverageUtilities.createDoubleWritableRaster(
				cols, rows, null, null, null);
		WritableRandomIter classiIter = RandomIterFactory.createWritable(
				classiWR, null);

		WritableRaster classi2WR = CoverageUtilities
				.createDoubleWritableRaster(cols, rows, null, null, null);
		WritableRandomIter classi2Iter = RandomIterFactory.createWritable(
				classi2WR, null);

		/*
		 * build the class matrix 1=inc inst 2=inc stab 3=stab 4=instab 0 inst 1
		 * stab rock=presence of rock
		 */
		pm.beginTask("Creating stability map...", rows);

		for (int j = 0; j < rows; j++) {
			pm.worked(1);
			for (int i = 0; i < cols; i++) {

				double slopeValue = slopeRI.getSampleDouble(i, j, 0);
				double abValue = abRI.getSampleDouble(i, j, 0);// VERIFICARE
																// UNITA DI
																// MISURA
																// METTERE M2
				double tangPhiValue = frictionRI.getSampleDouble(i, j, 0);
				double cohValue = cohesionRI.getSampleDouble(i, j, 0);
				double rhoValue = densityRI.getSampleDouble(i, j, 0);
				double hsValue = hsIter.getSampleDouble(i, j, 0);
				double rhoS = rhoSIter.getSampleDouble(i, j, 0);
				double rainfall_m_giorno = effectiveRI.getSampleDouble(i, j, 0) / 1000;
				double rhoW = rhoWIter.getSampleDouble(i, j, 0);
				double transmissivity_m2_g = trasmissivityRI.getSampleDouble(i,
						j, 0);
				if (computeTransmissivity) {
					transmissivity_m2_g = trasmissivityRI.getSampleDouble(i, j,
							0) * hsValue;
				}
				if (!isNovalue(slopeValue) && !isNovalue(abValue)
						&& !isNovalue(tangPhiValue) && !isNovalue(cohValue)
						&& !isNovalue(rhoValue) && !isNovalue(hsValue)
						&& !isNovalue(rhoS) && !isNovalue(rhoW)
						&& !isNovalue(transmissivity_m2_g)) {
					if (hsValue <= EPS || slopeValue > pRock) {
						classiIter.setSample(i, j, 0, ROCK);
						classi2Iter.setSample(i, j, 0, ROCK);
					} else {
						// if(slopeValue==0){slopeValue=0.001;}

						// VERICARE abvalue unita di misura
						double minvalue = Math.min(
								((rainfall_m_giorno / transmissivity_m2_g)
										* abValue / Math.sin(Math
										.atan(slopeValue))), 1.0);
						// System.out.println(minvalue+" "+abValue/
						// Math.sin(Math.atan(slopeValue)));
						double cosquadroalpha = Math.pow(
								(Math.cos(Math.atan(slopeValue))), 2);
						double fsnumC = cohValue
								/ (rhoS * 9.8 * hsValue
										* Math.sin(Math.atan(slopeValue)) * Math
											.cos(Math.atan(slopeValue)));
						double fsnumPhi = ((1 - (rhoW / rhoS) * minvalue) * cosquadroalpha)
								/ (Math.sin(Math.atan(slopeValue)) * Math
										.cos(Math.atan(slopeValue)))
								* tangPhiValue;
						double fs = fsnumC + fsnumPhi;
						if (fs >= pThresholdFS) {
							classiIter.setSample(i, j, 0, 1);

						} else {
							classiIter.setSample(i, j, 0, 0);

						}
					}
				} else {
					classiIter.setSample(i, j, 0, doubleNovalue);
				}
			}
		}
		pm.done();

		outShalstab = CoverageUtilities.buildCoverage("classi", classiWR,
				regionMap, inSlope.getCoordinateReferenceSystem());

	}
}
