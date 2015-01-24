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
import java.io.FileWriter;

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
public class OmsFs16042014 extends JGTModel {

	@Description("Map of slope (alpha)")
	@In
	public String pathToSlopeAngleMap = null;

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

	@Description("number of layers")
	@Unit("-")
	@In
	public int nLayers = 0;

	@Description("Path to the soil param file")
	@In
	public String pathToSoilFile = null;

	@Description("Path to the soil param file")
	@In
	public String headerPsiGEOtop = null;

	@Description("Path to the soil param file")
	@In
	public String headerThetaGEOtop = null;

	@Description("Path to the soil param file")
	@In
	public int time = 0;

	@Description("Path to geotopmaps")
	@In
	public String pathToGEOtopOutputMaps = null;

	@Description("Path to the soil mec param file")
	@In
	public String pathToSoilMecFile = null;

	@Description("Path to the soil mec param file")
	@In
	public String pathToSoilStatFile = null;

	@Description("Out FS map")
	@Out
	public GridCoverage2D outFs = null;

	@Description("Out FS map")
	@Out
	public GridCoverage2D outSigmaS = null;

	public final double EPS = 0.01;
	public final double GAMMAW_NM3 = 9810;

	/**
	 * Value to be given to pixels if <code>h_s < eps</code>.
	 */
	private static final double ROCK = 8888.0;
	public GridCoverage2D inTheta = null;
	public GridCoverage2D inPsi = null;
	public GridCoverage2D inSlope = null;
	public double matSoilWeigth[] = null;
	double[][] dataComplete = null;

	@Execute
	public void process() throws Exception {
		if (pRock == -9999.0)
			pRock = 5.67;
		if (!concatOr(outFs == null, doReset)) {
			return;
		}

		RandomIter treeWIter = null;
		if (inTreeW != null) {
			RenderedImage treeWRI = inTreeW.getRenderedImage();
			treeWIter = RandomIterFactory.create(treeWRI, null);
		} else {
			treeWIter = new ConstantRandomIter(pTreeW);
		}

		dataComplete = new double[nLayers][11];
		double[][] dataSoil = new double[nLayers][9];
		double[][] dataMecSoil = new double[nLayers][5];

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
					if (col > 8) {
						break;
					}
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
		while ((line = bufRdr.readLine()) != null && row < dataMecSoil.length) {
			if (skip > 0) {
				StringTokenizer st = new StringTokenizer(line, ",");
				while (st.hasMoreTokens()) {
					// get next token and store it in the array
					dataMecSoil[row][col] = Double.parseDouble(st.nextToken());
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
		// 8=n 9=kh 10=kv; 11=MeanC'; 12=MeanCroot; 13=MeanPhi; 14=MeanPsi;
		// 15=MeanTheta;
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

		OmsRasterReader reader_SLOPE = new OmsRasterReader();

		reader_SLOPE = new OmsRasterReader();
		reader_SLOPE.file = pathToSlopeAngleMap;
		reader_SLOPE.fileNovalue = -9999.0;
		reader_SLOPE.geodataNovalue = Double.NaN;
		reader_SLOPE.process();
		inSlope = reader_SLOPE.outRaster;
		HashMap<String, Double> regionMap = CoverageUtilities
				.getRegionParamsFromGridCoverage(inSlope);
		int cols = regionMap.get(CoverageUtilities.COLS).intValue();
		int rows = regionMap.get(CoverageUtilities.ROWS).intValue();
		matSoilWeigth = new double[cols * rows];

		for (int i = 0; i < nLayers; i++) {
			int k = i + 1;

			String indiceLayer = null;
			String indicetime = null;
			String indicelayerTheta = null;

			if (time <= 9) {
				indicetime = "000" + time;
			} else if (time <= 99) {
				indicetime = "00" + time;

			} else if (time <= 999) {
				indicetime = "0" + time;
			} else {
				indicetime = "" + time;
			}

			if (k <= 9) {
				indiceLayer = "000" + k;

			} else if (i <= 99) {
				indiceLayer = "00" + k;

			} else if (i <= 999) {
				indiceLayer = "0" + k;
			} else {
				indiceLayer = "" + k;
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

			String pathToPsi = pathToGEOtopOutputMaps + "/" + headerPsiGEOtop
					+ "L" + indiceLayer + "N" + indicetime + ".asc";
			String pathToTheta = pathToGEOtopOutputMaps + "/"
					+ headerThetaGEOtop + "L" + indicelayerTheta + "N"
					+ indicetime + ".asc";
			try {
				String filename = pathToPsi;
				FileWriter fw = new FileWriter(filename, true); // the true will
																// append the
																// new data
				fw.write("\n");// appends the string to the file
				fw.close();
			} catch (IOException ioe) {
				System.err.println("IOException: " + ioe.getMessage());
			}

			try {
				String filename = pathToTheta;
				FileWriter fw = new FileWriter(filename, true); // the true will
																// append the
																// new data
				fw.write("\n");// appends the string to the file
				fw.close();
			} catch (IOException ioe) {
				System.err.println("IOException: " + ioe.getMessage());
			}
			System.out.println(pathToPsi);

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

			RenderedImage slopeRI = inSlope.getRenderedImage();
			RenderedImage inThetaRI = null;
			RenderedImage psiRI = null;

			checkNull(inPsi);

			psiRI = inPsi.getRenderedImage();
			checkNull(inTheta);

			inThetaRI = inTheta.getRenderedImage();

			computeFS(slopeRI, psiRI, inThetaRI, cohesionIter, phiIter,
					tetaResIter, tetaSatIter, treeWIter, khIter, kvIter, i,
					dataComplete);
			String pathToFS = pathToGEOtopOutputMaps + "/" + "FS_N" + "L" + k
					+ "N" + indicetime + ".asc";
			OmsRasterWriter writer = new OmsRasterWriter();
			writer.inRaster = outFs;
			writer.file = pathToFS;
			writer.process();

			String pathToSigmaS = pathToGEOtopOutputMaps + "/" + "sigmaS_N"
					+ "L" + k + "N" + indicetime + ".asc";
			writer = new OmsRasterWriter();
			writer.inRaster = outSigmaS;
			writer.file = pathToSigmaS;
			writer.process();

		}
	}

	private void computeFS(RenderedImage slope, RenderedImage psi,
			RenderedImage inTheta, RandomIter cohesionIter,
			RandomIter tgPhiIter, RandomIter tetaResIter,
			RandomIter tetaSatIter, RandomIter treeWIter, RandomIter khIter,
			RandomIter kvIter, int iLayer, double[][] dataComplete) {
		HashMap<String, Double> regionMap = CoverageUtilities
				.getRegionParamsFromGridCoverage(inSlope);
		int cols = regionMap.get(CoverageUtilities.COLS).intValue();
		int rows = regionMap.get(CoverageUtilities.ROWS).intValue();

		RandomIter slopeRI = RandomIterFactory.create(slope, null);
		RandomIter psiRI = null;
		RandomIter thetaRI = null;

		psiRI = RandomIterFactory.create(psi, null);
		thetaRI = RandomIterFactory.create(inTheta, null);

		WritableRaster fsWR = CoverageUtilities.createDoubleWritableRaster(
				cols, rows, null, null, null);
		WritableRandomIter fsIter = RandomIterFactory
				.createWritable(fsWR, null);
		WritableRaster sigmaSWR = CoverageUtilities.createDoubleWritableRaster(
				cols, rows, null, null, null);
		WritableRandomIter sigmaSIter = RandomIterFactory.createWritable(
				sigmaSWR, null);

		pm.beginTask("Creating fs map...", rows);

		int conta = 0;
		for (int j = 0; j < rows; j++) {
			pm.worked(1);
			for (int i = 0; i < cols; i++) {

				double slopeValue = slopeRI.getSampleDouble(i, j, 0);
				slopeValue = Math.max(slopeValue, 0.001);

				double phiValue = tgPhiIter.getSampleDouble(i, j, 0);
				double cohValue = cohesionIter.getSampleDouble(i, j, 0);
				double thetaValue = 0;
				thetaValue = thetaRI.getSampleDouble(i, j, 0);

				double thetaResValue = tetaResIter.getSampleDouble(i, j, 0);
				double thetaSatValue = tetaSatIter.getSampleDouble(i, j, 0);
				double wtreeValue = treeWIter.getSampleDouble(i, j, 0);
				double khValue = khIter.getSampleDouble(i, j, 0);
				double kvValue = kvIter.getSampleDouble(i, j, 0);
				double psiValue = psiRI.getSampleDouble(i, j, 0);
				if (!isNovalue(slopeValue) && !isNovalue(phiValue)
						&& !isNovalue(cohValue) && !isNovalue(thetaValue)
						&& !isNovalue(thetaResValue)
						&& !isNovalue(thetaSatValue) && !isNovalue(wtreeValue)
						&& !isNovalue(khValue) && !isNovalue(kvValue)
						&& !isNovalue(psiValue)) {
					if (wtreeValue == -1) {
						wtreeValue = 0;
					}
					matSoilWeigth[conta] = matSoilWeigth[conta]
							+ (1 - 0.25)
							* ((dataComplete[iLayer][0] / 1000) * dataComplete[iLayer][6])
							+ thetaValue * (dataComplete[iLayer][0] / 1000)
							* 10;
//					matSoilWeigth[conta] = matSoilWeigth[conta]
//							+ 
//							 ((dataComplete[iLayer][0] / 1000) * dataComplete[iLayer][6]);
					double L_i = matSoilWeigth[conta];
					conta++;
					L_i = L_i + wtreeValue;
					if (Math.tan(slopeValue * Math.PI / 180) > pRock) {
						fsIter.setSample(i, j, 0, ROCK);
					} else {

						double[] fsValueV = computeFSvalue((slopeValue),
								(phiValue), (cohValue), (thetaValue),
								(thetaResValue), (thetaSatValue), (wtreeValue),
								(khValue), (kvValue), (psiValue), iLayer, L_i);
						double fsValue = fsValueV[0];
						double var = 0;

						// ///////////////////////////COESIONE/////////////////////////
						double epsilon = cohValue / 100;
						double varc = cohValue / 10;
						double cohValuePepsV[] = computeFSvalue((slopeValue),
								(phiValue), (cohValue + epsilon), (thetaValue),
								(thetaResValue), (thetaSatValue), (wtreeValue),
								(khValue), (kvValue), (psiValue), iLayer, L_i);
						double cohValuePeps = cohValuePepsV[0];

						double[] cohValueMepsV = computeFSvalue((slopeValue),
								(phiValue), (cohValue - epsilon), (thetaValue),
								(thetaResValue), (thetaSatValue), (wtreeValue),
								(khValue), (kvValue), (psiValue), iLayer, L_i);
						double cohValueMeps = cohValueMepsV[0];
						double derivativeCohe = (cohValuePeps - cohValueMeps)
								/ (2 * epsilon);
						double varCohe = varc * varc * derivativeCohe
								* derivativeCohe;
						// ///////////////////////////COESIONE/////////////////////////

						// ///////////////////////////ATTRITO/////////////////////////

						epsilon = phiValue / 100;
						double varp = phiValue / 10;
						double[] phiValuePepsV = computeFSvalue((slopeValue),
								(phiValue + epsilon), (cohValue), (thetaValue),
								(thetaResValue), (thetaSatValue), (wtreeValue),
								(khValue), (kvValue), (psiValue), iLayer, L_i);
						double phiValuePeps = phiValuePepsV[0];
						double[] phiValueMepsV = computeFSvalue((slopeValue),
								(phiValue - epsilon), (cohValue), (thetaValue),
								(thetaResValue), (thetaSatValue), (wtreeValue),
								(khValue), (kvValue), (psiValue), iLayer, L_i);
						double phiValueMeps = phiValueMepsV[0];
						double derivativePhi = (phiValuePeps - phiValueMeps)
								/ (2 * epsilon);
						double varPhi = varp * varp * derivativePhi
								* derivativePhi;
						// ///////////////////////////ATTRITO/////////////////////////

						// ///////////////////////////TETA/////////////////////////
						epsilon = thetaValue / 100;
						double vart = thetaValue / 10;

						double[] thetaValuePepsV = computeFSvalue((slopeValue),
								(phiValue + epsilon), (cohValue),
								(thetaValue + epsilon), (thetaResValue),
								(thetaSatValue), (wtreeValue), (khValue),
								(kvValue), (psiValue), iLayer, L_i);
						double thetaValuePeps = thetaValuePepsV[0];
						double[] thetaValueMepsV = computeFSvalue((slopeValue),
								(phiValue), (cohValue), (thetaValue - epsilon),
								(thetaResValue), (thetaSatValue), (wtreeValue),
								(khValue), (kvValue), (psiValue), iLayer, L_i);
						double thetaValueMeps = thetaValueMepsV[0];
						double derivativeTheta = (thetaValuePeps - thetaValueMeps)
								/ (2 * epsilon);
						double varTheta = vart * vart * derivativeTheta
								* derivativeTheta;
						// ///////////////////////////TETA/////////////////////////

						// //
						// ///////////////////////////PSI/////////////////////////
						// epsilon = Math.abs(psiValue) / 100;
						// double varps = Math.abs(psiValue) / 10;
						//
						// double[] psiValuePepsV = computeFSvalue((slopeValue),
						// (phiValue), (cohValue), (thetaValue),
						// (thetaResValue), (thetaSatValue), (wtreeValue),
						// (khValue), (kvValue), (psiValue + epsilon),
						// iLayer, L_i);
						// double psiValuePeps = psiValuePepsV[0];
						// double[] psiValueMepsV = computeFSvalue((slopeValue),
						// (phiValue), (cohValue), (thetaValue),
						// (thetaResValue), (thetaSatValue), (wtreeValue),
						// (khValue), (kvValue), (psiValue - epsilon),
						// iLayer, L_i);
						// double psiValueMeps = psiValueMepsV[0];
						// double derivativePsi = (psiValuePeps - psiValueMeps)
						// / (2 * epsilon);
						// double varPsi = varps * varps * derivativePsi
						// * derivativePsi;
						// ///////////////////////////PSI/////////////////////////

						// var = Math.pow((varPsi + varTheta + varCohe +
						// varPhi),0.5);
						var = Math.pow((varTheta + varCohe + varPhi), 0.5);
						// System.out.println(fsValue+" "+var);
						var = Math.max(var, 0.00001);
						// double p =1-
						// umontreal.iro.lecuyer.probdist.NormalDist.cdf(fsValue,
						// var, fsValue / var);

						// double p = 1-
						// umontreal.iro.lecuyer.probdist.NormalDist
						// .cdf(fsValue, var, fsValue / var);
						//
						// p = 1- umontreal.iro.lecuyer.probdist.NormalDist
						// .cdf(fsValue, var, 1);
						// p=fsValue / var;
						// fsIter.setSample(i, j, 0,p);
						double beta = (fsValue - 1) / var;
						double pf = 1 - umontreal.iro.lecuyer.probdist.NormalDist
								.cdf(0, 1, beta);
						fsIter.setSample(i, j, 0, pf);

						sigmaSIter.setSample(i, j, 0, fsValueV[1]);

						// System.out.println(p);

					}

				} else {
					matSoilWeigth[conta] = 0;
					conta++;

					fsIter.setSample(i, j, 0, doubleNovalue);
					sigmaSIter.setSample(i, j, 0, doubleNovalue);

				}

			}
		}

		outFs = CoverageUtilities.buildCoverage("fs", fsWR, regionMap,
				inSlope.getCoordinateReferenceSystem());
		outSigmaS = CoverageUtilities.buildCoverage("sigmaS", sigmaSWR,
				regionMap, inSlope.getCoordinateReferenceSystem());

	}

	public double[] computeFSvalue(double slopeValue, double phiValue,
			double cohValue, double thetaValue, double thetaResValue,
			double thetaSatValue, double wtreeValue, double khValue,
			double kvValue, double psiValue, int iLayer, double L_i) {
		// 0=Dz, 1=tetaR, 2=tetaS, 3=c', 4=c_root, 5=phi',
		// 6=gammas, 7=a,

		// public double computeFSvalue(double slopeValue, double phiValue,
		// double cohValue, double thetaValue, double thetaResValue,
		// double thetaSatValue, double wtreeValue, double khValue,
		// double kvValue, double psiValue, int iLayer, double L_i) {
		// // 0=Dz, 1=tetaR, 2=tetaS, 3=c', 4=c_root, 5=phi',
		// // 6=gammas, 7=a,

		// double ter1 = 2 * cohValue
		// / (L_i * Math.sin(2 * Math.atan(slopeValue)));
		// if (thetaValue < thetaResValue) {
		// thetaValue = thetaResValue;
		// }
		// if (thetaValue > thetaSatValue) {
		// thetaValue = thetaSatValue;
		// }
		// double primo = ((thetaValue - thetaResValue) / (thetaSatValue -
		// thetaResValue));
		// double uw = primo * (psiValue / 1000) * (GAMMAW_NM3);
		//
		// double ter2 = 2 * uw * Math.tan(phiValue * Math.PI / 180)
		// / ((L_i * 1000) * sin(2 * Math.atan(slopeValue)));
		//
		// double ter3 = Math.tan(phiValue * Math.PI / 180)
		// / Math.tan(Math.atan(slopeValue));
		//
		// double fs = ter1 + ter3 - ter2;
		// double fsModified = Math.min(fs, 10);

		// return fsModified;

		// double beta = slopeValue * Math.PI / 180;
		// double ter2 = 2 * cohValue / (L_i * Math.sin(2 * beta));
		// double ter1 = Math.tan(phiValue * Math.PI / 180) / Math.tan(beta);
		//
		// double ua = 0;
		// double uwKN_m2 = psiValue * 0.01;
		// double sigmaS = 0;
		//
		// if ((ua - uwKN_m2) <= 0) {
		// sigmaS = -(ua - uwKN_m2);
		// } else {
		// double primo = ((thetaValue - thetaResValue) / (thetaSatValue -
		// thetaResValue));
		// sigmaS = -primo * (ua - uwKN_m2);
		// }
		//
		// double ter3 = (sigmaS / (L_i))
		// * (Math.tan(beta) + (1 / Math.tan(beta)))
		// * Math.tan(phiValue * Math.PI / 180);
		//
		// double fs = ter1 + ter2 - ter3;
		// double fsModified = Math.min(fs, fs);
		// double r[] = new double[] { fsModified, sigmaS };

		double beta = slopeValue * Math.PI / 180;
		
		double ter1 = Math.tan(phiValue * Math.PI / 180) / Math.tan(beta);
		double ter2 = 2*cohValue / (L_i * Math.sin(2*beta));


		double uwKN_m2 = psiValue * 0.01;
		double ua=0;
		double sigmaS=0;
		if(ua-uwKN_m2<=0){
			sigmaS=-(ua-uwKN_m2);
		}else{
			double primo = ((thetaValue - thetaResValue) / (thetaSatValue - thetaResValue));
			sigmaS=-primo*(ua-uwKN_m2);
		}
		
		double ter3 =(sigmaS*(Math.tan(beta)+1/(Math.tan(beta)))*Math.tan(phiValue * Math.PI / 180))/L_i;

		double fs = ter1 + ter2 - ter3;
		double fsModified = Math.min(fs, fs);
		double r[] = new double[] { fsModified, sigmaS };

		// return fsModified;
		return r;

	}

}
