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
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.fs;

import static java.lang.Math.atan;
import java.io.BufferedReader;
import java.io.File;
import org.apache.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
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
import umontreal.iro.lecuyer.probdist.NormalDist;
import umontreal.iro.lecuyer.probdist.ContinuousDistribution;
import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.coverage.ConstantRandomIter;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;

@Description(OMSSHALSTAB_DESCRIPTION)
@Author(name = OMSSHALSTAB_AUTHORNAMES, contact = OMSSHALSTAB_AUTHORCONTACTS)
@Keywords(OMSSHALSTAB_KEYWORDS)
@Label(OMSSHALSTAB_LABEL)
@Name(OMSSHALSTAB_NAME)
@Status(OMSSHALSTAB_STATUS)
@License(OMSSHALSTAB_LICENSE)
public class OmsFs19032014 extends JGTModel {

	@Description("Map of slope (tg alpha)")
	@In
	public GridCoverage2D inSlope = null;

	@Description("N layer")
	@In
	public int nLayers = -1;

	@Description(OMSSHALSTAB_pRock_DESCRIPTION)
	@In
	public double pRock = -9999.0;

	@Description("Map of in TreeW")
	@Unit("KPa")
	@In
	public GridCoverage2D inTreeW = null;

	@Description("Value of in TreeW")
	@Unit("KPa")
	@In
	public double pTreeW = -1.0;

	@Description("Friction angle map tg(phi)")
	@In
	public GridCoverage2D inTgphi = null;

	@Description("Friction angle value tg(phi)")
	@In
	public double pTgphi = -1.0;

	@Description("Cohesion map")
	@Unit("KPa")
	@In
	public GridCoverage2D inCohesion = null;

	@Description("Cohesion value")
	@Unit("KPa")
	@In
	public double pCohesion = -1.0;

	@Description("Soil dept map")
	@Unit("m")
	@In
	public GridCoverage2D inSdepth = null;

	@Description("Use GEOtop data")
	@Unit("m")
	@In
	public boolean doUseGEOtopOutputs = false;

	@Description("Soil depth value")
	@Unit("m")
	@In
	public double pSdepth = -1.0;

	@Description("soil specific weight map ()")
	@Unit("KN m-2")
	@In
	public GridCoverage2D inGammaS = null;

	@Description("soil specific weight map ()")
	@Unit("KN m-2")
	@In
	public int time = 0;

	@Description("Path to the soil param file")
	@In
	public String pathToSoilFile = null;

	@Description("Path to the soil complete param file")
	@In
	public String pathToSoilCompleteFile = null;

	@Description("Path to geotopmaps")
	@In
	public String pathToGEOtopOutputMaps = null;

	@Description("Header psi GEOtop output maps")
	@In
	public String headerPsiGEOtop = null;

	@Description("Header theta GEOtop output maps")
	@In
	public String headerThetaGEOtop = null;

	@Description("Path to the soil mec param file")
	@In
	public String pathToSoilMecFile = null;

	@Description("Soil specific weight value")
	@In
	public double pGammaS = -1.0;

	@Description("Out FS map")
	@Out
	public GridCoverage2D outFs = null;

	public final double EPS = 0.01;
	public final double GAMMAW_NM3 = 9810;

	/**
	 * Value to be given to pixels if <code>h_s < eps</code>.
	 */
	private static final double ROCK = 8888.0;
	public GridCoverage2D inTheta = null;
	public GridCoverage2D inPsi = null;

	@Execute
	public void process() throws Exception {
		if (pRock == -9999.0)
			pRock = 5.67;
		if (!concatOr(outFs == null, doReset)) {
			return;
		}

		File folder = new File(pathToGEOtopOutputMaps);
		File[] listOfFiles = folder.listFiles();
		File[] thetaMaps = new File[listOfFiles.length];
		File[] psiMaps = new File[listOfFiles.length];
		int countP = 0;
		int countT = 0;
		int countTime = 0;
		boolean entrato = false;
		int vecttime[] = new int[listOfFiles.length];
		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].getName().startsWith(headerPsiGEOtop)) {
				psiMaps[countP] = listOfFiles[i];
				countP++;
			}
			if (listOfFiles[i].getName().startsWith(headerThetaGEOtop)) {
				thetaMaps[countT] = listOfFiles[i];
				countT++;
			}
		}
		for (int i = 0; i < countT; i++) {

			// int a =
			// Integer.parseInt(thetaMaps[i].getName().split("L")[1].split("N")[1].split(".asc")[0]);
			int a = Integer.parseInt(thetaMaps[i].getName().split("L")[1]
					.split("N")[1]);
			if (countTime > 0) {
				if (a == vecttime[0]) {
					break;
				}
			}
			vecttime[countTime] = a;
			countTime++;

		}

		// checkNull(inSlope);
		//
		// RenderedImage slopeRI = inSlope.getRenderedImage();
		//
		// RandomIter tghiIter = null;
		// if (inTgphi != null) {
		// RenderedImage tgphiRI = inTgphi.getRenderedImage();
		// tghiIter = RandomIterFactory.create(tgphiRI, null);
		// } else {
		// tghiIter = new ConstantRandomIter(pTgphi);
		// }
		//
		RandomIter treeWIter = null;
		if (inTreeW != null) {
			RenderedImage treeWRI = inTreeW.getRenderedImage();
			treeWIter = RandomIterFactory.create(treeWRI, null);
		} else {
			treeWIter = new ConstantRandomIter(pTreeW);
		}
		//
		// RandomIter cohesionIter = null;
		// if (inCohesion != null) {
		// RenderedImage cohesionRI = inCohesion.getRenderedImage();
		// cohesionIter = RandomIterFactory.create(cohesionRI, null);
		// } else {
		// cohesionIter = new ConstantRandomIter(pCohesion);
		// }
		//
		// RandomIter hsIter = null;
		// if (inSdepth != null) {
		// RenderedImage hsRI = inSdepth.getRenderedImage();
		// hsIter = RandomIterFactory.create(hsRI, null);
		// } else {
		// hsIter = new ConstantRandomIter(pSdepth);
		// }
		//
		// RandomIter gammaSIter = null;
		// if (inGammaS != null) {
		// RenderedImage rhoRI = inGammaS.getRenderedImage();
		// gammaSIter = RandomIterFactory.create(rhoRI, null);
		// } else {
		// gammaSIter = new ConstantRandomIter(pGammaS);
		// }
		// 0=Dz, 1=tetaR, 2=tetaS, 3=C', 4= c_root, 5=phi', 6=gammas

		double[][] dataComplete = new double[nLayers][11];
		double[][] dataSoil = new double[nLayers][9];
		double[][] dataMecSoil = new double[nLayers][5];

		if (doUseGEOtopOutputs) {
			File file = new File(pathToSoilFile);
			int row = 0;
			int col = 0;
			BufferedReader bufRdr = new BufferedReader(new FileReader(file));
			String line = null;
			// 0=Dz, 1=Kh, 2=Kv, 3=res, 4=fc, 5=sat, 6=a, 7=n, 8=SS
			// read each line of text file
			int skip = 0;
			while ((line = bufRdr.readLine()) != null && row < dataSoil.length) {
				StringTokenizer st = new StringTokenizer(line, ",");
				if (skip > 0) {
					while (st.hasMoreTokens()) {
						// get next token and store it in the array
						dataSoil[row][col] = Double.parseDouble(st.nextToken());
						System.out.print(dataSoil[row][col] + " ");
						col++;
					}
					col = 0;
					row++;
					System.out.println(" ");
				} else {
					skip = 100;
				}

			}

			System.out.println();

			bufRdr.close();

			file = new File(pathToSoilMecFile);
			bufRdr = new BufferedReader(new FileReader(file));
			row = 0;
			col = 0;
			skip = 0;
			bufRdr = new BufferedReader(new FileReader(file));
			line = null;
			// 0=Dz, 1=c', 2=c_root, 3=phi', 4=gammas
			// read each line of text file
			while ((line = bufRdr.readLine()) != null
					&& row < dataMecSoil.length) {
				if (skip > 0) {
					StringTokenizer st = new StringTokenizer(line, ",");
					while (st.hasMoreTokens()) {
						// get next token and store it in the array
						dataMecSoil[row][col] = Double.parseDouble(st
								.nextToken());
						System.out.print(dataMecSoil[row][col] + " ");
						col++;
					}

					col = 0;
					row++;
					System.out.println();
				} else {
					skip = 100;
				}
			}
			bufRdr.close();
			// fill dataComplete
			// 0=Dz, 1=tetaR, 2=tetaS, 3=c', 4=c_root, 5=phi', 6=gammas, 7=a,
			// 8=n
			for (int i = 0; i < nLayers; i++) {
				dataComplete[i][0] = dataSoil[i][0];
				dataComplete[i][1] = dataSoil[i][3];
				dataComplete[i][2] = dataSoil[i][5];
				dataComplete[i][7] = dataSoil[i][6];
				dataComplete[i][8] = dataSoil[i][7];
				dataComplete[i][3] = dataMecSoil[i][1];
				dataComplete[i][4] = dataMecSoil[i][2];
				dataComplete[i][5] = dataMecSoil[i][3];
				dataComplete[i][6] = dataMecSoil[i][4];
				dataComplete[i][9] = dataSoil[i][1];
				dataComplete[i][10] = dataSoil[i][2];

			}

		} else {
			File file = new File(pathToSoilCompleteFile);
			int row = 0;
			int col = 0;
			BufferedReader bufRdr = new BufferedReader(new FileReader(file));
			String line = null;
			// 0=Dz, 1=tetaR, 2=tetaS, 3=c', 4=c_root, 5=phi', 6=rho_s 1800
			// kg/m3, 7=a,
			// 8=n
			// read each line of text file
			int skip = 0;
			while ((line = bufRdr.readLine()) != null
					&& row < dataComplete.length) {
				StringTokenizer st = new StringTokenizer(line, ",");
				if (skip > 0) {
					while (st.hasMoreTokens()) {
						// get next token and store it in the array
						dataComplete[row][col] = Double.parseDouble(st
								.nextToken());
						System.out.print(dataComplete[row][col] + " ");
						col++;
					}
					col = 0;
					row++;
					System.out.println(" ");
				} else {
					skip = 100;
				}

			}

			System.out.println();

			bufRdr.close();
		}

		// per ogni layer: creo una mappa di coesione, angolo di attrito tetaR e
		// tataS una mappa di spessore di suolo pari al layer o la creo?
		// leggo o creo una mappa di ua o psi e calcolo FS

		// for()
		for (int ii = 0; ii < countTime; ii++) {
			time = vecttime[ii];

			for (int i = 0; i < nLayers; i++) {
				String indiceLayer = null;
				String indicelayerTheta = null;
				String indicetime = null;
				int k = i + 1;
				if (time <= 9) {
					indicetime = "000" + time;
				} else if (time <= 99) {
					indicetime = "00" + time;

				} else if (time <= 999) {
					indicetime = "0" + time;
				} else {
					indicetime = "" + time;
				}

				if (i <= 9) {
					indiceLayer = "000" + i;

				} else if (i <= 99) {
					indiceLayer = "00" + i;

				} else if (i <= 999) {
					indiceLayer = "0" + i;
				} else {
					indiceLayer = "" + i;
				}

				if (k <= 9) {
					indicelayerTheta = "000" + k;

				} else if (i <= 99) {
					indicelayerTheta = "00" + k;

				} else if (i <= 999) {
					indicelayerTheta = "0" + k;

				} else {
					indicelayerTheta = "" + k;
				}

				// String pathToPsi = pathToGEOtopOutputMaps + "/"
				// + headerPsiGEOtop + "L" + indiceLayer + "N"
				// + indicetime+".asc";
				// String pathToTheta = pathToGEOtopOutputMaps + "/"
				// + headerThetaGEOtop + "L" + (indicelayerTheta) + "N"
				// + indicetime+".asc";
				// System.out.println("indiceLayer=" + indiceLayer
				// + "  indicelayerTheta=" + indicelayerTheta
				// + " indicetime" + indicetime + "  pathToPsi="
				// + pathToPsi + "  pathToTheta=" + pathToTheta);
				String pathToPsi = pathToGEOtopOutputMaps + "/"
						+ headerPsiGEOtop + "L" + indiceLayer + "N"
						+ indicetime;
				String pathToTheta = pathToGEOtopOutputMaps + "/"
						+ headerThetaGEOtop + "L" + (indicelayerTheta) + "N"
						+ indicetime;
				System.out.println("indiceLayer=" + indiceLayer
						+ "  indicelayerTheta=" + indicelayerTheta
						+ " indicetime" + indicetime + "  pathToPsi="
						+ pathToPsi + "  pathToTheta=" + pathToTheta);

				OmsRasterReader reader_PSI = new OmsRasterReader();
				reader_PSI = new OmsRasterReader();
				reader_PSI.file = pathToPsi;
				reader_PSI.fileNovalue = -9999.0;
				reader_PSI.geodataNovalue = Double.NaN;
				reader_PSI.process();
				inPsi = reader_PSI.outRaster;

				OmsRasterReader reader_THETA = new OmsRasterReader();

				reader_THETA = new OmsRasterReader();
				reader_THETA.file = pathToTheta;
				reader_THETA.fileNovalue = -9999.0;
				reader_THETA.geodataNovalue = Double.NaN;
				reader_THETA.process();
				inTheta = reader_THETA.outRaster;
				// calcolo Li: manca la W_t
				double L_i = 0;
				for (int c = 0; c <= i; c++) {
					// z e in mm e lo trasformo in m
					double deltazeta_in_m = dataComplete[c][0] / 1000;
					// dataComplete[c][6] e il soil specific weight in KN/m3
					L_i += (deltazeta_in_m * dataComplete[c][6]);
					// dataComplete[c][6] e in KN/m3 e deltaz_in_m e in metri.
					// Risultato: kN/m2 a questo ci devo sommare il three
				}
				// L_i e in KN/m2 gammaw e in N/m3 quindi lo divido epr 1000
				// diverso da silvia. L_i e solo gamma_i per H_i
				// L_i = L_i / (GAMMAW_KNM3 / 1000);
				// ora L_i e in m
				System.out.println();
				RandomIter cohesionIter = null;
				// dataComplete[i][3] + dataComplete[i][4] sono in KN/m2
				double tatalCohesion = dataComplete[i][3] + dataComplete[i][4];
				cohesionIter = new ConstantRandomIter(tatalCohesion);

				RandomIter phiIter = null;
				double phi = dataComplete[i][5];
				phiIter = new ConstantRandomIter(phi);

				RandomIter tetaSatIter = null;
				double tetasat = dataComplete[i][2];
				tetaSatIter = new ConstantRandomIter(tetasat);

				RandomIter tetaResIter = null;
				double tetares = dataComplete[i][1];
				tetaResIter = new ConstantRandomIter(tetares);

				RandomIter khIter = null;
				double kh = dataComplete[i][9];
				khIter = new ConstantRandomIter(kh);

				RandomIter kvIter = null;
				double kv = dataComplete[i][10];
				kvIter = new ConstantRandomIter(kv);

				checkNull(inSlope);

				RenderedImage slopeRI = inSlope.getRenderedImage();
				RenderedImage inThetaRI = null;
				RenderedImage psiRI = null;
				if (doUseGEOtopOutputs) {
					checkNull(inPsi);

					psiRI = inPsi.getRenderedImage();
					checkNull(inTheta);

					inThetaRI = inTheta.getRenderedImage();

				}
				computeFS(slopeRI, psiRI, inThetaRI, cohesionIter, phiIter,
						tetaResIter, tetaSatIter, treeWIter, khIter, kvIter,
						L_i, i, dataComplete);
				String pathToFS = pathToGEOtopOutputMaps + "/" + "FS_N" + "L"
						+ indiceLayer + "N" + indicetime;
				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = outFs;
				writer.file = pathToFS;
				writer.process();

			}
		}
	}

	private void computeFS(RenderedImage slope, RenderedImage psi,
			RenderedImage inTheta, RandomIter cohesionIter,
			RandomIter tgPhiIter, RandomIter tetaResIter,
			RandomIter tetaSatIter, RandomIter treeWIter, RandomIter khIter,
			RandomIter kvIter, double L_i, int iLayer, double[][] dataComplete) {
		HashMap<String, Double> regionMap = CoverageUtilities
				.getRegionParamsFromGridCoverage(inSlope);
		int cols = regionMap.get(CoverageUtilities.COLS).intValue();
		int rows = regionMap.get(CoverageUtilities.ROWS).intValue();

		RandomIter slopeRI = RandomIterFactory.create(slope, null);
		RandomIter psiRI = null;
		RandomIter thetaRI = null;
		if (doUseGEOtopOutputs) {
			psiRI = RandomIterFactory.create(psi, null);
			thetaRI = RandomIterFactory.create(inTheta, null);
		}
		WritableRaster fsWR = CoverageUtilities.createDoubleWritableRaster(
				cols, rows, null, null, null);
		WritableRandomIter fsIter = RandomIterFactory
				.createWritable(fsWR, null);

		pm.beginTask("Creating fs map...", rows);
		for (int j = 0; j < rows; j++) {
			pm.worked(1);
			for (int i = 0; i < cols; i++) {
				double slopeValue = slopeRI.getSampleDouble(i, j, 0);
				// if (slopeValue == 0) {
				// slopeValue = 0.00001;
				// }
				double phiValue = tgPhiIter.getSampleDouble(i, j, 0);
				double cohValue = cohesionIter.getSampleDouble(i, j, 0);
				double thetaValue = 0;
				if (doUseGEOtopOutputs) {
					thetaValue = thetaRI.getSampleDouble(i, j, 0);
				}
				double thetaResValue = tetaResIter.getSampleDouble(i, j, 0);
				double thetaSatValue = tetaSatIter.getSampleDouble(i, j, 0);
				double wtreeValue = treeWIter.getSampleDouble(i, j, 0);
				double khValue = khIter.getSampleDouble(i, j, 0);
				double kvValue = kvIter.getSampleDouble(i, j, 0);
				if (!isNovalue(slopeValue) && !isNovalue(phiValue)
						&& !isNovalue(cohValue) && !isNovalue(thetaValue)
						&& !isNovalue(thetaResValue)
						&& !isNovalue(thetaSatValue) && !isNovalue(wtreeValue)
						&& !isNovalue(khValue) && !isNovalue(kvValue)) {
					if (wtreeValue == -1) {
						wtreeValue = 0;
					}

					L_i = L_i + wtreeValue;
					if (slopeValue > pRock) {
						fsIter.setSample(i, j, 0, ROCK);
					} else {
						double psiValue = 0;
						if (doUseGEOtopOutputs) {
							if (!isNovalue(psiRI.getSampleDouble(i, j, 0))) {
								psiValue = psiRI.getSampleDouble(i, j, 0);
							} else {
								fsIter.setSample(i, j, 0, doubleNovalue);
								break;
							}
						} else {
							double zzz = 0;
							for (int l = 0; l < iLayer; l++) {
								zzz += dataComplete[l][0];
							}
							double alpha = 0.4;
							double q = -1.5 * 10e-7;
							psiValue = -(1 / alpha) * Math.log(1 + q / kvValue)
									* Math.exp(-GAMMAW_NM3 * alpha * zzz);

						}

						// 0=Dz, 1=tetaR, 2=tetaS, 3=c', 4=c_root, 5=phi',
						// 6=gammas, 7=a,
						double cohe = 0;
						for (int layer = 0; layer < iLayer; layer++) {
							double c = 2
									* (dataComplete[layer][3] + dataComplete[layer][4])
									/ (Math.sin(2 * Math.atan(slopeValue))
											* (dataComplete[layer][0] / 1000) * dataComplete[layer][6]);
							cohe = cohe + c;
						}

						double ter1 = 2 * cohValue
								/ (L_i * Math.sin(2 * Math.atan(slopeValue)));
						double primo = ((thetaValue - thetaResValue) / (thetaSatValue - thetaResValue));
						double uw = primo * (psiValue / 1000) * (GAMMAW_NM3);

						double ter2 = 2
								* uw
								* Math.tan(phiValue * Math.PI / 180)
								/ ((L_i * 1000) * sin(2 * Math.atan(slopeValue)));

						double ter3 = Math.tan(phiValue * Math.PI / 180)
								/ Math.tan(Math.atan(slopeValue));

						double fs = ter1 + ter3 - ter2;
						double fsModified = Math.min(fs, 10);
						fsIter.setSample(i, j, 0, fsModified);
						System.out.print(i + "   " + j + "  " + fsModified
								+ "  statistico= ");
						boolean doStatistic = true;
						if (doStatistic) {

							double M_i = L_i / (GAMMAW_NM3 / 1000);
							double A = M_i * (GAMMAW_NM3 / 1000)
									* sin(2 * Math.atan(slopeValue)) / 2;
							double D = M_i
									* Math.tan(Math.atan(slopeValue))
									/ (M_i - (primo * (psiValue / 1000) / (Math
											.cos(Math.atan(slopeValue)) * Math
											.cos(Math.atan(slopeValue)))));

							D = (A / (GAMMAW_NM3 / 1000))
									/ (((Math.cos(Math.atan(slopeValue)) * Math
											.cos(Math.atan(slopeValue)))) * M_i - primo
											* psiValue / 1000);
							double meanFS = (Math.tan(phiValue * Math.PI / 180) / D)
									+ cohValue / A;

							double varPhi = Math.tan(5 * Math.PI / 180);
							double varC = 0;
							// double varFS = varPhi / (D * D) + varC / (A * A);
							double varFS = Math.pow(varPhi / D, 2)
									+ Math.pow(varC / A, 2);
							double fsss = (umontreal.iro.lecuyer.probdist.NormalDist
									.cdf(meanFS, Math.pow(varFS, 0.5), 1)) / 0.5;

							System.out.println(fsss + " mean=" + meanFS
									+ " variance=" + Math.pow(varFS, 0.5)
									+ "  A=" + A + "  D=" + D + "  psivalue="
									+ psiValue + " theta=" + thetaValue);

						}

					}

				} else {
					fsIter.setSample(i, j, 0, doubleNovalue);
				}

			}
		}

		outFs = CoverageUtilities.buildCoverage("fs", fsWR, regionMap,
				inSlope.getCoordinateReferenceSystem());

	}
	/**
	 * Calculates the FS for each pixel.
	 */
	// private void fs_computation(RenderedImage slope, RandomIter frictionRI,
	// RandomIter cohesionRI, RandomIter hsIter, RandomIter effectiveRI,
	// RandomIter densityRI) {
	// HashMap<String, Double> regionMap = CoverageUtilities
	// .getRegionParamsFromGridCoverage(inSlope);
	// int cols = regionMap.get(CoverageUtilities.COLS).intValue();
	// int rows = regionMap.get(CoverageUtilities.ROWS).intValue();
	//
	// RandomIter slopeRI = RandomIterFactory.create(slope, null);
	//
	// WritableRaster qcritWR = CoverageUtilities.createDoubleWritableRaster(
	// cols, rows, null, null, null);
	// WritableRandomIter qcritIter = RandomIterFactory.createWritable(
	// qcritWR, null);
	// WritableRaster classiWR = CoverageUtilities.createDoubleWritableRaster(
	// cols, rows, null, null, null);
	// WritableRandomIter classiIter = RandomIterFactory.createWritable(
	// classiWR, null);
	//
	// pm.beginTask("Creating qcrit map...", rows);
	// for (int j = 0; j < rows; j++) {
	// pm.worked(1);
	// for (int i = 0; i < cols; i++) {
	// double slopeValue = slopeRI.getSampleDouble(i, j, 0);
	// double tanPhiValue = frictionRI.getSampleDouble(i, j, 0);
	// double cohValue = cohesionRI.getSampleDouble(i, j, 0);
	// double rhoValue = densityRI.getSampleDouble(i, j, 0);
	// double hsValue = hsIter.getSampleDouble(i, j, 0);
	//
	// if (!isNovalue(slopeValue) && !isNovalue(tanPhiValue)
	// && !isNovalue(cohValue) && !isNovalue(rhoValue)
	// && !isNovalue(hsValue)) {
	// if (hsValue <= EPS || slopeValue > pRock) {
	// qcritIter.setSample(i, j, 0, ROCK);
	// } else {
	// double checkUnstable = tanPhiValue + cohValue
	// / (9810.0 * rhoValue * hsValue)
	// * (1 + pow(slopeValue, 2));
	// if (slopeValue >= checkUnstable) {
	// /*
	// * uncond unstable
	// */
	// qcritIter.setSample(i, j, 0, 5);
	// } else {
	// double checkStable = tanPhiValue
	// * (1 - 1 / rhoValue) + cohValue
	// / (9810 * rhoValue * hsValue)
	// * (1 + pow(slopeValue, 2));
	// if (slopeValue < checkStable) {
	// /*
	// * uncond. stable
	// */
	// qcritIter.setSample(i, j, 0, 0);
	// } else {
	// double qCrit = trasmissivityRI.getSampleDouble(
	// i, j, 0)
	// * sin(atan(slopeValue))
	// / abRI.getSampleDouble(i, j, 0)
	// * rhoValue
	// * (1 - slopeValue / tanPhiValue + cohValue
	// / (9810 * rhoValue * hsValue * tanPhiValue)
	// * (1 + pow(slopeValue, 2)))
	// * 1000;
	// qcritIter.setSample(i, j, 0, qCrit);
	// /*
	// * see the Qcrit (critical effective
	// * precipitation) that leads the slope to
	// * instability (see article of Montgomery et Al,
	// * Hydrological Processes, 12, 943-955, 1998)
	// */
	// double value = qcritIter.getSampleDouble(i, j,
	// 0);
	// if (value > 0 && value < 50)
	// qcritIter.setSample(i, j, 0, 1);
	// if (value >= 50 && value < 100)
	// qcritIter.setSample(i, j, 0, 2);
	// if (value >= 100 && value < 200)
	// qcritIter.setSample(i, j, 0, 3);
	// if (value >= 200)
	// qcritIter.setSample(i, j, 0, 4);
	// }
	// }
	// }
	// } else {
	// qcritIter.setSample(i, j, 0, doubleNovalue);
	// }
	// }
	// }
	// pm.done();
	//
	// /*
	// * build the class matrix 1=inc inst 2=inc stab 3=stab 4=instab
	// * rock=presence of rock
	// */
	// pm.beginTask("Creating stability map...", rows);
	// double Tq = 0;
	// for (int j = 0; j < rows; j++) {
	// pm.worked(1);
	// for (int i = 0; i < cols; i++) {
	// Tq = trasmissivityRI.getSampleDouble(i, j, 0)
	// / (effectiveRI.getSampleDouble(i, j, 0) / 1000.0);
	// double slopeValue = slopeRI.getSampleDouble(i, j, 0);
	// double abValue = abRI.getSampleDouble(i, j, 0);
	// double tangPhiValue = frictionRI.getSampleDouble(i, j, 0);
	// double cohValue = cohesionRI.getSampleDouble(i, j, 0);
	// double rhoValue = densityRI.getSampleDouble(i, j, 0);
	// double hsValue = hsIter.getSampleDouble(i, j, 0);
	//
	// if (!isNovalue(slopeValue) && !isNovalue(abValue)
	// && !isNovalue(tangPhiValue) && !isNovalue(cohValue)
	// && !isNovalue(rhoValue) && !isNovalue(hsValue)) {
	// if (hsValue <= EPS || slopeValue > pRock) {
	// classiIter.setSample(i, j, 0, ROCK);
	// } else {
	// double checkUncondUnstable = tangPhiValue + cohValue
	// / (9810 * rhoValue * hsValue)
	// * (1 + pow(slopeValue, 2));
	// double checkUncondStable = tangPhiValue
	// * (1 - 1 / rhoValue) + cohValue
	// / (9810 * rhoValue * hsValue)
	// * (1 + pow(slopeValue, 2));
	// double checkStable = Tq
	// * sin(atan(slopeValue))
	// * rhoValue
	// * (1 - slopeValue / tangPhiValue + cohValue
	// / (9810 * rhoValue * hsValue * tangPhiValue)
	// * (1 + pow(slopeValue, 2)));
	// if (slopeValue >= checkUncondUnstable) {
	// classiIter.setSample(i, j, 0, 1);
	// } else if (slopeValue < checkUncondStable) {
	// classiIter.setSample(i, j, 0, 2);
	// } else if (abValue < checkStable
	// && classiIter.getSampleDouble(i, j, 0) != 1
	// && classiIter.getSampleDouble(i, j, 0) != 2) {
	// classiIter.setSample(i, j, 0, 3);
	// } else {
	// classiIter.setSample(i, j, 0, 4);
	// }
	// }
	// } else {
	// classiIter.setSample(i, j, 0, doubleNovalue);
	// }
	// }
	// }
	// pm.done();

	// outQcrit = CoverageUtilities.buildCoverage("qcrit", qcritWR,
	// regionMap,
	// inSlope.getCoordinateReferenceSystem());
	// outShalstab = CoverageUtilities.buildCoverage("classi", classiWR,
	// regionMap, inSlope.getCoordinateReferenceSystem());

	// }

}
