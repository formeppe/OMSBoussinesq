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
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inSlope_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_inTca_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_pRock_DESCRIPTION;

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
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;

@Description(OMSSHALSTAB_DESCRIPTION)
@Author(name = OMSSHALSTAB_AUTHORNAMES, contact = OMSSHALSTAB_AUTHORCONTACTS)
@Keywords(OMSSHALSTAB_KEYWORDS)
@Label(OMSSHALSTAB_LABEL)
@Name(OMSSHALSTAB_NAME)
@Status(OMSSHALSTAB_STATUS)
@License(OMSSHALSTAB_LICENSE)
public class OmsComputeSoilThickness extends JGTModel {

	@Description(OMSSHALSTAB_inSlope_DESCRIPTION)
	@In
	public GridCoverage2D inSlope = null;

	@Description(OMSSHALSTAB_inTca_DESCRIPTION)
	@In
	public GridCoverage2D inTca = null;

	@Description("Pamaterisazion 1: on the slope")
	@In
	public boolean doP1 = false;

	@Description("Pamaterisazion 2: on the tca")
	@In
	public boolean doP2 = false;

	@Description("Parameter 1")
	@Unit("-")
	@In
	public double p1 = -1.0;

	@Description("Parameter 2")
	@Unit("-")
	@In
	public double p2 = -1.0;

	@Description(OMSSHALSTAB_pRock_DESCRIPTION)
	@In
	public double pRock = -9999.0;

	@Description("Soil Map")
	@Out
	public GridCoverage2D outSoilDepth = null;

	public final double EPS = 0.01;

	/**
	 * Value to be given to pixels if <code>h_s < eps</code>.
	 */

	@Execute
	public void process() throws Exception {
		if (!concatOr(outSoilDepth == null, doReset)) {
			return;
		}
		if (doP1 == true) {
			checkNull(inSlope);
		}
		if (doP2 == true) {
			checkNull(inTca);
		}

		if (pRock == -9999.0)
			pRock = 5.67;

		RenderedImage abRI = null;
		RenderedImage slopeRI = null;
		if (doP1 == true) {

			slopeRI = inSlope.getRenderedImage();
		}
		if (doP2 == true) {
			abRI = inTca.getRenderedImage();
		}
		qcrit(slopeRI, abRI);
	}

	/**
	 * Calculates the trasmissivity in every pixel of the map.
	 */
	private void qcrit(RenderedImage slope, RenderedImage ab) {
		HashMap<String, Double> regionMap = CoverageUtilities
				.getRegionParamsFromGridCoverage(inSlope);
		int cols = regionMap.get(CoverageUtilities.COLS).intValue();
		int rows = regionMap.get(CoverageUtilities.ROWS).intValue();
		RandomIter slopeRI = null;
		if (doP1 == true) {
			slopeRI = RandomIterFactory.create(slope, null);
		}
		RandomIter abRI = null;
		if (doP2 == true) {
			abRI = RandomIterFactory.create(ab, null);
		}
		WritableRaster classiWR = CoverageUtilities.createDoubleWritableRaster(
				cols, rows, null, null, null);
		WritableRandomIter classiIter = RandomIterFactory.createWritable(
				classiWR, null);

		pm.beginTask("Creating soil depth map...", rows);
		for (int j = 0; j < rows; j++) {
			pm.worked(1);
			for (int i = 0; i < cols; i++) {
				double slopeValue = 0;
				if (doP1 == true) {
					slopeValue = slopeRI.getSampleDouble(i, j, 0);
				}
				double abValue = 0;
				if (doP2 == true) {

					abValue = abRI.getSampleDouble(i, j, 0);
				}
				if (!isNovalue(slopeValue) || !isNovalue(abValue)) {
					if (doP1) {
						double hsValue = p1 - p2 * slopeValue;
//						if (hsValue > 5.0) {
//							hsValue = 5;
//						}
						if (hsValue < EPS) {
							hsValue = EPS;
						}
						classiIter.setSample(i, j, 0, hsValue);

					} else {
						if (doP2) {
							double hsValue = p1 + p2 * Math.log(abValue);
//							if (hsValue > 5.0) {
//								hsValue = 5;
//							}
							if (hsValue < EPS) {
								hsValue = EPS;
							}
							classiIter.setSample(i, j, 0, hsValue);
						} else {
							System.out
									.println("No parameterisation rule was choosen");
						}

					}
				}

				else {
					classiIter.setSample(i, j, 0, doubleNovalue);
				}
			}
		}
		pm.done();

		outSoilDepth = CoverageUtilities.buildCoverage("classi", classiWR,
				regionMap, inSlope.getCoordinateReferenceSystem());

	}
}
