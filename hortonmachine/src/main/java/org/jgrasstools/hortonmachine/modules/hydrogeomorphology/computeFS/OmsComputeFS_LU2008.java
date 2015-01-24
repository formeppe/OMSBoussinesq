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

@Description(OMSSHALSTAB_DESCRIPTION)
@Author(name = OMSSHALSTAB_AUTHORNAMES, contact = OMSSHALSTAB_AUTHORCONTACTS)
@Keywords(OMSSHALSTAB_KEYWORDS)
@Label(OMSSHALSTAB_LABEL)
@Name(OMSSHALSTAB_NAME)
@Status(OMSSHALSTAB_STATUS)
@License(OMSSHALSTAB_LICENSE)
public class OmsComputeFS_LU2008 extends JGTModel {

	@Description(OMSSHALSTAB_inSlope_DESCRIPTION)
	@In
	public GridCoverage2D inSlope = null;

	@Description(OMSSHALSTAB_inSlope_DESCRIPTION)
	@In
	public int pMode = 0;

	@Description(OMSSHALSTAB_inTca_DESCRIPTION)
	@In
	public GridCoverage2D inTca = null;

	@Description(OMSSHALSTAB_inTrasmissivity_DESCRIPTION)
	@Unit("-")
	@In
	public GridCoverage2D innVG = null;

	@Description(OMSSHALSTAB_pTrasmissivity_DESCRIPTION)
	@Unit("-")
	@In
	public double pnVG = -1.0;

	@Description(OMSSHALSTAB_inTca_DESCRIPTION)
	@In
	public GridCoverage2D inDem = null;

	@Description(OMSSHALSTAB_pTrasmissivity_DESCRIPTION)
	@Unit("m^2/day")
	@In
	public double pThresholdFS = -1.0;
	

	@Description(OMSSHALSTAB_pTrasmissivity_DESCRIPTION)
	@Unit("m")
	@In
	public double pHwt = -1.0;

	@Description(OMSSHALSTAB_inTgphi_DESCRIPTION)
	@In
	public GridCoverage2D inTgphi = null;

	@Description(OMSSHALSTAB_pTgphi_DESCRIPTION)
	@In
	public double pTgphi = -1.0;

	@Description(OMSSHALSTAB_inCohesion_DESCRIPTION)
	@Unit("kPa")
	@In
	public GridCoverage2D inCohesion = null;

	@Description(OMSSHALSTAB_pCohesion_DESCRIPTION)
	@Unit("kPa")
	@In
	public double pCohesion = -1.0;

	@Description(OMSSHALSTAB_inCohesion_DESCRIPTION)
	@Unit("1/KPa")
	@In
	public GridCoverage2D inAlphaVG = null;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	@Unit("1/KPa")
	public double palphaVG = -1.0;

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

	@Description(OMSSHALSTAB_inRho_DESCRIPTION)
	@In
	@Unit("kN/m3")
	public GridCoverage2D inGammaSoil = null;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	@Unit("kN/m3")
	public double pGammaSoil = -1.0;

	@Description(OMSSHALSTAB_inRho_DESCRIPTION)
	@In
	public GridCoverage2D inDuration = null;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	public double pDuration = -1.0;

	@Description(OMSSHALSTAB_inRho_DESCRIPTION)
	@In
	@Unit("mm/g")
	public GridCoverage2D inKs = null;

	@Description(OMSSHALSTAB_pRho_DESCRIPTION)
	@In
	@Unit("mm/g")
	public double pKs = -1.0;

	@Description(OMSSHALSTAB_pRock_DESCRIPTION)
	@In
	public double pRock = -9999.0;

	@Description(OMSSHALSTAB_outShalstab_DESCRIPTION)
	@Out
	public GridCoverage2D outShalstab = null;

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
		checkNull(inSlope, inTca);

		if (pRock == -9999.0)
			pRock = 5.67;

		OmsCurvatures c = new OmsCurvatures();
		c.inElev = inDem;
		c.process();

		OmsAb ab = new OmsAb();
		ab.inPlan = c.outPlan;
		ab.inTca = inTca;
		ab.process();

		RenderedImage abRI = ab.outAb.getRenderedImage();

		RenderedImage slopeRI = inSlope.getRenderedImage();

		RandomIter nVGIter = null;
		if (innVG != null) {
			RenderedImage nVGRI = innVG.getRenderedImage();
			nVGIter = RandomIterFactory.create(nVGRI, null);
		} else {
			nVGIter = new ConstantRandomIter(pnVG);
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

		RandomIter gammaSoilIter = null;
		if (inGammaSoil != null) {
			RenderedImage gammaSoilRI = inGammaSoil.getRenderedImage();
			gammaSoilIter = RandomIterFactory.create(gammaSoilRI, null);
		} else {
			gammaSoilIter = new ConstantRandomIter(pGammaSoil);
		}

		RandomIter alphaVGIter = null;
		if (inAlphaVG != null) {
			RenderedImage alphaVGRI = inAlphaVG.getRenderedImage();
			alphaVGIter = RandomIterFactory.create(alphaVGRI, null);
		} else {
			alphaVGIter = new ConstantRandomIter(palphaVG);
		}

		RandomIter ksSIter = null;
		if (inKs != null) {
			RenderedImage ksRI = inKs.getRenderedImage();
			ksSIter = RandomIterFactory.create(ksRI, null);
		} else {
			ksSIter = new ConstantRandomIter(pKs);
		}

		RandomIter gammaWIter = null;
		gammaWIter = new ConstantRandomIter(9.789);// KN/m3

		qcrit(slopeRI, abRI, nVGIter, tghiIter, cohesionIter, hsIter, qIter,
				gammaSoilIter, alphaVGIter, ksSIter, gammaWIter);
	}

	/**
	 * Calculates the trasmissivity in every pixel of the map.
	 */
	private void qcrit(RenderedImage slope, RenderedImage ab, RandomIter nVGRI,
			RandomIter frictionRI, RandomIter cohesionRI, RandomIter hsIter,
			RandomIter effectiveRI, RandomIter gammaSoilIter,
			RandomIter alphaVGIter, RandomIter ksIter, RandomIter gammaWIter) {

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

		pm.beginTask("Creating stability map...", rows);

		for (int j = 0; j < rows; j++) {
			pm.worked(1);
			for (int i = 0; i < cols; i++) {

				double slopeValue = slopeRI.getSampleDouble(i, j, 0);
				// no unit
				double abValue = abRI.getSampleDouble(i, j, 0);
				// m
				double tangPhiValue = frictionRI.getSampleDouble(i, j, 0);
				// no units
				double cohValue = cohesionRI.getSampleDouble(i, j, 0);
				// KPa
				double gammaSoilValue = gammaSoilIter.getSampleDouble(i, j, 0);
				// KN/m3
				double ksValue = ksIter.getSampleDouble(i, j, 0);
				// mm/g
				double hsValue = hsIter.getSampleDouble(i, j, 0);
				// m
				double alphaVGValue = alphaVGIter.getSampleDouble(i, j, 0);
				// no KPa-1
				double rainfall_mm_giorno = effectiveRI
						.getSampleDouble(i, j, 0);
				// mm/g
				double gammaW = gammaWIter.getSampleDouble(i, j, 0);
				// 9.8 kN/m3
				double nValue = nVGRI.getSampleDouble(i, j, 0);
				// -

				if (!isNovalue(slopeValue) && !isNovalue(abValue)
						&& !isNovalue(tangPhiValue) && !isNovalue(cohValue)
						&& !isNovalue(alphaVGValue) && !isNovalue(hsValue)
						&& !isNovalue(ksValue) && !isNovalue(gammaW)
						&& !isNovalue(gammaSoilValue) && !isNovalue(nValue)) {
					if (hsValue <= EPS || slopeValue > pRock) {
						classiIter.setSample(i, j, 0, ROCK);
					} else {
						double alphaSlope = Math.atan(slopeValue);
						double fsC = 0;
						double fsPhi = 0;
						double fsSigmaS = 0;
						double zzz = pHwt-hsValue;

						fsC = 2
								* cohValue
								/ (gammaSoilValue * zzz * Math
										.sin(2.0 * alphaSlope));
						fsPhi = tangPhiValue / Math.tan(alphaSlope);
						double sigmaS = 0;
						double uA_uW = 0;
//						nValue = 4.75;
//						alphaVGValue = 0.08;
//						ksValue=10e-6;
//						rainfall_mm_giorno=1.5*1e-7;
						//while (zzz > 0) {

							//zzz = zzz - 0.1;
							// double zzz = 5.0 - hsValue;

							double log = (1.0 + (-rainfall_mm_giorno) / ksValue)
									* Math.exp(-gammaW * alphaVGValue * zzz)
									- (-rainfall_mm_giorno) / ksValue;
							uA_uW = -(1 / alphaVGValue) * Math.log(log);

							if (uA_uW <= 0.0) {
								sigmaS = -uA_uW;
							} else {
								double ter1 = alphaVGValue * uA_uW;
								double ter2 = Math.pow(ter1, nValue);
								double ter3 = Math.pow(1 + ter2, (nValue - 1.0)
										/ nValue);
								sigmaS = -uA_uW / (ter3);
							}
							//System.out.println(zzz + "  " + sigmaS);
						//}
						fsSigmaS = (sigmaS / (gammaSoilValue * hsValue))
								* (Math.tan(alphaSlope) + (1 / (Math
										.tan(alphaSlope)))) * tangPhiValue;

						double fs = fsC + fsPhi - fsSigmaS;
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
