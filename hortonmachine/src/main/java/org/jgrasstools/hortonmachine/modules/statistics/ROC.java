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
package org.jgrasstools.hortonmachine.modules.statistics;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;
import static org.jgrasstools.gears.libs.modules.ModelsEngine.go_downstream;
import static org.jgrasstools.gears.libs.modules.ModelsEngine.isSourcePixel;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_outShalstab_DESCRIPTION;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

import oms3.annotations.Author;
import oms3.annotations.Documentation;
import oms3.annotations.Label;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;

import org.geotools.coverage.grid.GridCoverage2D;
import org.hamcrest.core.IsNull;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.i18n.HortonMessageHandler;

@Description("Marks all the outlets of the considered region on the drainage directions map with the conventional value 10.")
@Documentation("Markoutlets.html")
@Author(name = "Antonello Andrea, Erica Ghesla, Cozzini Andrea, Franceschi Silvia, Pisoni Silvano, Rigon Riccardo", contact = "http://www.hydrologis.com, http://www.ing.unitn.it/dica/hp/?user=rigon")
@Keywords("Outlets, Dem, Raster, FlowDirections, DrainDir")
@Label(JGTConstants.DEMMANIPULATION)
@Name("markoutlets")
@Status(Status.TESTED)
@License("General Public License Version 3 (GPLv3)")
public class ROC extends JGTModel {
	@Description("The map of flowdirections.")
	@In
	public GridCoverage2D inModelledMap = null;

	@Description("The map of flowdirections.")
	@In
	public double flagPositiveMeas = -999;

	@Description("The map of flowdirections.")
	@In
	public double flagNegativeMeas = -999;

	@Description("The map of flowdirections.")
	@In
	public double flagTrueModelled = -999;

	@Description("The map of flowdirections.")
	@In
	public double flagFalseModelled = -999;

	@Description("The map of flowdirections.")
	@In
	public GridCoverage2D inConditioningMap = null;

	@Description("The map of flowdirections.")
	@In
	public GridCoverage2D inMeasuredMap = null;

	@Description("The map of flowdirections.")
	@In
	public double threshold = -999;

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outTruePositiveRate = -999;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outFalsePositiveRate = -999;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outAccuracy = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outSuccessIndex = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outD2PerfectClassification = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outTrueSkillStatistic = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outHss = -999.0;
	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outOddsRatio = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outOddsRatioSkillScore = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outAverageIndex = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outCriticalSuccessIndex = -999.0;
	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outProbFalseDetection = -999.0;
	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outEquitableSuccesIndex = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outBiasScore = -999.0;
	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outFalseAlarmRatio = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public double outPrecision = -999.0;

	@Description("The map of the flowdirections with outlet marked.")
	@Out
	public String pPath = null;

	@Description(OMSSHALSTAB_outShalstab_DESCRIPTION)
	@Out
	public GridCoverage2D outClass = null;

	private HortonMessageHandler msg = HortonMessageHandler.getInstance();

	@Execute
	public void process() throws IOException {

		checkNull(inModelledMap);
		checkNull(inMeasuredMap);

		HashMap<String, Double> regionMap = CoverageUtilities
				.getRegionParamsFromGridCoverage(inMeasuredMap);
		int nCols = regionMap.get(CoverageUtilities.COLS).intValue();
		int nRows = regionMap.get(CoverageUtilities.ROWS).intValue();

		WritableRaster classiWR = CoverageUtilities.createDoubleWritableRaster(
				nCols, nRows, null, null, null);
		WritableRandomIter classiIter = RandomIterFactory.createWritable(
				classiWR, null);

		RenderedImage modelledRI = inModelledMap.getRenderedImage();
		RenderedImage measuredRI = inMeasuredMap.getRenderedImage();

		RandomIter conditioningIter = null;
		RenderedImage conditioningRI = null;
		if (inConditioningMap != null) {
			conditioningRI = inConditioningMap.getRenderedImage();
			conditioningIter = RandomIterFactory.create(conditioningRI, null);
		}
		//
		// WritableRaster flowWR = CoverageUtilities.createDoubleWritableRaster(
		// nCols, nRows, null, null, doubleNovalue);
		// CoverageUtilities.setNovalueBorder(flowWR);

		RandomIter modelledIter = RandomIterFactory.create(modelledRI, null);
		RandomIter measuredIter = RandomIterFactory.create(measuredRI, null);

		pm.beginTask(msg.message("markoutlets.working"), 2 * nRows); //$NON-NLS-1$

		int[] punto = new int[2];

		double tp = 0, fp = 0, tn = 0, fn = 0, contaNA = 0;
		for (int i = 0; i < nRows; i++) {
			for (int j = 0; j < nCols; j++) {
				punto[0] = j;
				punto[1] = i;

				double measuredSample = measuredIter.getSampleDouble(punto[0],
						punto[1], 0);
				double modelledSample = modelledIter.getSampleDouble(punto[0],
						punto[1], 0);

				// 1=tp; 2=fn; 3=tn;4=fp
				if (inConditioningMap != null) {
					double condtioningSample = conditioningIter
							.getSampleDouble(punto[0], punto[1], 0);

					if (!isNovalue(measuredSample)
							&& !isNovalue(condtioningSample)
							&& !isNovalue(modelledSample)) {

						if (condtioningSample > threshold) {

							if (measuredSample == flagPositiveMeas
									&& modelledSample == flagTrueModelled) {
								tp += 1;
								classiIter.setSample(punto[0], punto[1], 0, 1);
							}
							if (measuredSample == flagPositiveMeas
									&& modelledSample == flagFalseModelled) {
								fn += 1;
								classiIter.setSample(punto[0], punto[1], 0, 2);

							}
							if (measuredSample == flagNegativeMeas
									&& modelledSample == flagTrueModelled) {
								tn += 1;
								classiIter.setSample(punto[0], punto[1], 0, 3);

							}
							if (measuredSample == flagNegativeMeas
									&& modelledSample == flagFalseModelled) {
								fp += 1;
								classiIter.setSample(punto[0], punto[1], 0, 4);

							}

						}

					} else {

						classiIter.setSample(punto[0], punto[1], 0,
								doubleNovalue);

					}

				} else {

					if (!isNovalue(measuredSample)

					&& !isNovalue(modelledSample)) {

						if (measuredSample == flagPositiveMeas
								&& modelledSample == flagTrueModelled) {
							tp += 1;

							classiIter.setSample(punto[0], punto[1], 0, 1);

						}
						if (measuredSample == flagPositiveMeas
								&& modelledSample == flagFalseModelled) {
							fn += 1;
							classiIter.setSample(punto[0], punto[1], 0, 2);

						}
						if (measuredSample == flagNegativeMeas
								&& modelledSample == flagFalseModelled) {
							tn += 1;
							classiIter.setSample(punto[0], punto[1], 0, 3);

						}
						if (measuredSample == flagNegativeMeas
								&& modelledSample == flagTrueModelled) {
							fp += 1;
							classiIter.setSample(punto[0], punto[1], 0, 4);

						}

					} else {
						contaNA += 1;
						classiIter.setSample(punto[0], punto[1], 0,
								doubleNovalue);

					}

				}
			}

			pm.worked(1);
		}
		pm.done();

		outClass = CoverageUtilities.buildCoverage("classi", classiWR,
				regionMap, inModelledMap.getCoordinateReferenceSystem());

		// ////////////// 1 ////////////////////
		if (tp > 1 && tn > 1) {
			// 0-1, optimal=1
			outTruePositiveRate = tp / (tp + fn);
			// /////////////////////////////////////

			// ////////////// 2 ////////////////////
			// 0-1, optimal=0
			outFalsePositiveRate = fp / (fp + tn);
			// /////////////////////////////////////

			// ////////////// 3 ////////////////////
			// 0-1, optimal=1
			outAccuracy = (tp + tn) / (tp + fn + fp + tn);

			// /////////////////////////////////////

			// ////////////// 3 ////////////////////
			outPrecision = (tp) / (tp + fp);
			// /////////////////////////////////////

			// ////////////// 4 ////////////////////
			// 0-oo, optimal=1
			outBiasScore = (tp + fp) / (tp + fn);
			// /////////////////////////////////////

			// ////////////// 5 ////////////////////
			// 0-1, optimal=0
			outFalseAlarmRatio = fp / (tp + fp);
			// /////////////////////////////////////

			// ////////////// 6 ////////////////////
			// 0-1, optimal=1
			outCriticalSuccessIndex = tp / (tp + fn + fp);
			// /////////////////////////////////////

			// ////////////// 7 ////////////////////
			// -1/3-1, optimal=1
			double R = (tp + fn) * (tp + fp) / (tp + fn + fp + tn);
			outEquitableSuccesIndex = (tp - R) / (tp + fp + fn - R);
			// /////////////////////////////////////

			// ////////////// 8 ////////////////////
			// 0-1, optimal=1
			outSuccessIndex = 0.5 * (((tp) / (tp + fn)) + ((tn) / (fp + tn)));
			// /////////////////////////////////////

			// ////////////// 9 ////////////////////
			// 0-1, optimal=0
			outD2PerfectClassification = Math.pow((1 - outTruePositiveRate)
					* (1 - outTruePositiveRate) + outFalsePositiveRate
					* outFalsePositiveRate, 0.5);

			// ////////////// 10 ////////////////////
			//
			//0-1
			outAverageIndex = ((tp / (tp + fn)) + (tp / (tp + fp))
					+ (tn / (fp + tn)) + (tn / (fn + tn))) / 4.0;

			// ////////////// 11 ////////////////////
			// -1 and 1, optimal=1
			outTrueSkillStatistic = outTruePositiveRate - outFalsePositiveRate;

			// ////////////// 12 ////////////////////
			// -oo and 1, optimal=1
			outHss = (2 * (tp * tn) - (fp * fn))
					/ ((tp + fn) * (fn + tn) + (tp + fp) * (fp + tn));

			// ////////////// 13 ////////////////////
			// 0 and +oo, optimal=1
			outOddsRatio = (tp * tn) / (fp * fn);

			// ////////////// 14 ////////////////////
			// -1 and +1, optimal=1
			outOddsRatioSkillScore = (outOddsRatio - 1) / (outOddsRatio + 1);

			// outProbFalseDetection = fp / (tn + fp);
		} else {
			// 0-1, optimal=1
			outTruePositiveRate = 0;
			// /////////////////////////////////////

			// ////////////// 2 ////////////////////
			// 0-1, optimal=0
			outFalsePositiveRate = 1;
			// /////////////////////////////////////

			// ////////////// 3 ////////////////////
			// 0-1, optimal=1
			outAccuracy = 0;

			// /////////////////////////////////////

			// ////////////// 3 ////////////////////
			outPrecision = (tp) / (tp + fp);
			// /////////////////////////////////////

			// ////////////// 4 ////////////////////
			// 0-oo, optimal=1
			outBiasScore = 1000;
			// /////////////////////////////////////

			// ////////////// 5 ////////////////////
			// 0-1, optimal=1
			outFalseAlarmRatio = 0;
			// /////////////////////////////////////

			// ////////////// 6 ////////////////////
			// 0-1, optimal=1
			outCriticalSuccessIndex = 0;
			// /////////////////////////////////////

			// ////////////// 7 ////////////////////
			// -1/3-1, optimal=1
			double R = (tp + fn) * (tp + fp) / (tp + fn + fp + tn);
			outEquitableSuccesIndex = -1;
			;
			// /////////////////////////////////////

			// ////////////// 8 ////////////////////
			// 0-1, optimal=1
			outSuccessIndex = 0;
			// /////////////////////////////////////

			// ////////////// 9 ////////////////////
			// 0-1, optimal=0
			outD2PerfectClassification = 1.0;

			// ////////////// 10 ////////////////////
			//
			outAverageIndex = ((tp / (tp + fn)) + (tp / (tp + fp))
					+ (tn / (fp + tn)) + (tn / (fn + tn))) / 4.0;

			// ////////////// 11 ////////////////////
			// -1 and 1, optimal=1
			outTrueSkillStatistic = -1;

			// ////////////// 12 ////////////////////
			// -oo and 1, optimal=1
			outHss = -10;

			// ////////////// 13 ////////////////////
			// 0 and +oo, optimal=1
			outOddsRatio = 1000;

			// ////////////// 14 ////////////////////
			// -1 and +1, optimal=1
			outOddsRatioSkillScore = -1;
		}
		if (pPath != null && pPath.length() > 0) {
			FileWriter Rstatfile = new FileWriter(pPath);
			PrintWriter errestat = new PrintWriter(Rstatfile);

			errestat.println("outEquitableSuccesIndex,Accuracy,Precision,outTruePositiveRate,outFalsePositiveRate,outBiasScore,outProbFalseDetection,outCriticalSuccessIndex,outSuccessIndex,outD2PerfectClassification,outAverageIndex,outTrueSkillStatistic,outHss,outOddsRatio,outOddsRatioSkillScore,tp,fn,tn,fp");
			errestat.print(outEquitableSuccesIndex + "," + outAccuracy + ","
					+ outPrecision + "," + outTruePositiveRate + ","
					+ outFalsePositiveRate + "," + outBiasScore + ","
					+ outProbFalseDetection + "," + outCriticalSuccessIndex
					+ "," + outSuccessIndex + "," + outD2PerfectClassification
					+ "," + outAverageIndex + "," + outTrueSkillStatistic + ","
					+ outHss + "," + outOddsRatio + ","
					+ outOddsRatioSkillScore + "," + tp + "," + fn + "," + tn
					+ "," + fp);

			Rstatfile.close();
		}

		System.out.println("tp= " + tp);
		System.out.println("fn= " + fn);
		System.out.println("tn= " + tn);
		System.out.println("fp= " + fp);
		System.out.println("Na= " + contaNA);
		System.out.println("Accuracy= " + outAccuracy + "   Precision= "
				+ outPrecision + "  outTruePositiveRate=" + outTruePositiveRate
				+ "  outFalsePositiveRate=" + outFalsePositiveRate
				+ " outBiasScore=" + outBiasScore + " outProbFalseDetection="
				+ outProbFalseDetection + " outCriticalSuccessIndex="
				+ outCriticalSuccessIndex + " outSuccessIndex="
				+ outSuccessIndex);

	}

}
