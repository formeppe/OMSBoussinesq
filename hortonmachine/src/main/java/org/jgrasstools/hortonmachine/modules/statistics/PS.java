package org.jgrasstools.hortonmachine.modules.statistics;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Map.Entry;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;
import oms3.annotations.Unit;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.DummyProgressMonitor;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.modules.r.rangelookup.OmsRangeLookup;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeFS;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeFS_GSE;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeFS_LU2008;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeFS_ShalstabLike;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.computeFS.OmsComputeSoilThickness;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.shalstab.OmsShalstab;
import org.joda.time.DateTime;

import umontreal.iro.lecuyer.rng.MRG32k3a;
import umontreal.iro.lecuyer.rng.RandomPermutation;
import umontreal.iro.lecuyer.rng.RandomStream;

//import corejava.Format;

//import sun.tools.asm.Cover;
//import umontreal.iro.lecuyer.probdist.*;
//import umontreal.iro.lecuyer.stat.*;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;

public class PS extends JGTModel {

	@Description("The map of the skyview factor.")
	@In
	public GridCoverage2D inskyview = null;

	@Description("The vector of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public GridCoverage2D inNet = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public GridCoverage2D inTopIndex = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public GridCoverage2D inDrainDir = null;

	@Description("The progress monitor.")
	@In
	public String END_DATE = null;

	@Description("The map of the dem.")
	@In
	public GridCoverage2D dem = null;

	@Description("The map of the insolation for February.")
	@In
	public GridCoverage2D insolationFebruary = null;

	@Description("The map of the insolation for February.")
	@In
	public GridCoverage2D insolationJanuary = null;

	@Description("The map of the insolation for March.")
	@In
	public GridCoverage2D insolationMarch = null;

	@Description("The map of the insolation for March.")
	@In
	public double[] inVar = null;

	@Description("The map of the insolation for March.")
	@In
	public double[] inNp = null;

	@Description("The map of the insolation for March.")
	@In
	public double[] inDistance = null;
	@Description("The map of the insolation for March.")
	@In
	public String inmodelname = null;

	@Description("The map of the insolation for April.")
	@In
	public GridCoverage2D insolationApril = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D insolationMay = null;

	@Description("The map of the insolation for June.")
	@In
	public GridCoverage2D insolationJune = null;

	@Description("The map of the insolation for June.")
	@In
	public GridCoverage2D inMeasuredForRoc = null;
	// ////////////////////////////////////////////
	@Description("The map of the skyview factor.")
	@In
	public GridCoverage2D inSkyview = null;
	@Description("The first day of the simulation.")
	@In
	public String tStartDate = null;
	@Description("The time step in minutes of the measurement.")
	@In
	public int inTimestep;
	@Description("The last day of the simulation.")
	@In
	public String tEndDate = null;

	@Description("The map of the dem.")
	@In
	public GridCoverage2D inDem = null;

	@Description("The first day of the simulation.")
	@In
	public String time = null;

	@Description("The Tmelt.")
	@In
	@Unit("C")
	public double pTmelt;
	@Description("Path to the temperature file in input.")
	@In
	public String pPathToRain = null;

	@Description("Path to the humidity file in input.")
	@In
	public String pPathToTemp = null;
	@Description("Path to the temperature file in input.")
	@In
	public String pPathToCondIniI = null;

	@Description("Path to the humidity file in input.")
	@In
	public String pPathToCondIniL = null;

	@Description("The pr of lmax.")
	@In
	@Unit("C")
	public double pR;
	@Description("The path to direct output.")
	@In
	public String pPathSWE;
	@Description("The path to direct output.")
	@In
	public String pPathMelting;
	@Description("The map of the insolation for February.")
	@In
	public GridCoverage2D inInsFeb = null;

	@Description("The map of the insolation for March.")
	@In
	public GridCoverage2D inInsMar = null;

	@Description("The map of the insolation for April.")
	@In
	public GridCoverage2D inInsApr = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D inInsMay = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D slope = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D hydraulicCondictivityMap = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D porosityMap = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D cohesionMap = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D frictionAngleMap = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D tca = null;

	@Description("The map of the insolation for May.")
	@In
	public GridCoverage2D sdepth = null;

	@Description("The map of the insolation for June.")
	@In
	public GridCoverage2D inInsJun = null;
	@Description("The map of the insolation for June.")
	@In
	public double threshold = 0.001;

	@Description("Do you want raster map as output.")
	@In
	public boolean doRaster = true;

	@Description("The temperature hashmap.")
	@In
	@Unit(" C ")
	public HashMap<Integer, double[]> inTemp;

	@Description("The rainfall.")
	@In
	@Unit("mm")
	public HashMap<Integer, double[]> inRainfall;

	@Description("The inital condition I.")
	@In
	@Unit("mm")
	public HashMap<Integer, double[]> inCondIniI;
	@Description("The inital condition L")
	@In
	@Unit("mm")
	public HashMap<Integer, double[]> inCondIniL;
	// ////////////////////////////////////////////////////

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new DummyProgressMonitor();

	@Description("")
	@In
	public int Warmup;

	@Description("Potential ETP vector.")
	@In
	public double[] Extra_PET;

	@Description("Precipitations vector.")
	@In
	public double[] Extra_Precip;

	@Description("The limit of the time series considered in computation.")
	@In
	public double Extra_MaxT = doubleNovalue;

	@Description("The network data.")
	@In
	public SimpleFeatureCollection net;

	@Description("The network data.")
	@In
	public SimpleFeatureCollection hydrome;

	@Description("The network data.")
	@In
	public SimpleFeatureCollection basin;

	@Description("Define the measured streamflow data.")
	@In
	public double[] Measurement_MeasData;

	@Description("Area of the basin.")
	@In
	public double Area;

	@Description("The first value of mesured discharge")
	@In
	@Unit("m3/s")
	public double Q0 = doubleNovalue;

	@Description("kmax")
	@In
	public int kmax;

	@Description("Rain input path")
	@In
	public String pathToRain;

	@Description("Runoff input path")
	@In
	public String pathToRunoff;

	@Description("ModelName")
	@In
	public String ModelName = null;

	@Description("Give the parameter ranges (minimum values).")
	@In
	public double[] ParRange_minn;

	@Description("Give the parameter ranges (maximum values).")
	@In
	public double[] ParRange_maxn;

	@Description("p")
	@In
	public int p;
	@Description("parameters")
	@In
	public int parameters;

	@Description("parameters")
	@In
	public String pPath = null;

	@Description("gbest")
	@Out
	public double outOpt;

	@Description("optimal parameters")
	@Out
	public double[] outOptSet;
	@Description("optimal parameters")
	@Out
	public double[] observedvect;

	private static final String DBL_FMT = "##.##########";
	public int int_Extra_MaxT;
	public double[] costvect;
	public double[][] p_best;
	public double[] g_best;
	public double g_best_value;
	public double[] p_best_value;
	public double[] hist_g_best_value;
	public double observed[];
	public double modelled[];

	public FileWriter Rstatfile;
	public PrintWriter errestat;

	@Execute
	public void process() throws Exception {
		int_Extra_MaxT = (int) Extra_MaxT;
		if (pPath != null) {
			Rstatfile = new FileWriter(pPath);
			errestat = new PrintWriter(Rstatfile);
		}
		double xX[][] = new double[parameters][p];

		// xX = LHSU(ParRange_minn, ParRange_maxn, p);
		xX = LHSU(ParRange_minn, ParRange_maxn, p);
		double x[][] = ReflectBounds(xX, ParRange_maxn, ParRange_minn);
		// for (int i = 0; i < x.length; i++) {
		// for (int j = 0; j < x[0].length;j++){
		// System.out.print(x[i][j]+ " ");
		// }
		// System.out.println();
		//
		// }
		double VelRange_minn[] = new double[ParRange_minn.length];

		double VelRange_max[] = new double[ParRange_minn.length];
		// ipotizzo che la velocita iniziale sia inizializzata compresa tra
		// 1/100 dei valori max e min dei parametri
		for (int i = 0; i < ParRange_maxn.length; i++) {
			VelRange_minn[i] = ParRange_minn[i];
			VelRange_max[i] = ParRange_maxn[i];
		}

		double vel[][] = new double[parameters][p];
		vel = LHSU(VelRange_minn, VelRange_max, p);

		// calculate the cost of each particle
		double[] costvectold = ComputeCostFunction(x);

		int kkk = 0;
		p_best = x;
		p_best_value = costvectold;
		double min = Math.abs(costvectold[0]);
		int posmin = 0;
		for (int i = 1; i < costvectold.length; i++) {
			if (Math.abs(costvectold[i]) < min) {
				min = Math.abs(costvectold[i]);
				posmin = i;
			}
		}
		g_best = new double[x[0].length];
		g_best_value = min;
		for (int i = 0; i < x[0].length; i++) {
			g_best[i] = x[posmin][i];
		}

		hist_g_best_value = new double[kmax];
		boolean fermati = false;
		while (kkk < (kmax - 1) || !fermati) {
			double[][] x_old = x;
			double[][] velnew = Compute_velocity(x_old, vel);
			vel = velnew;
			x = Compute_particle(x_old, velnew);
			costvect = ComputeCostFunction(x);
			p_best = Compute_pBest(x, costvect);

			g_best = Compute_gBest(x, costvect);
			hist_g_best_value[kkk] = g_best_value;
			for (int jjj = 0; jjj < g_best.length; jjj++) {

				System.out.print(g_best[jjj] + "  ");

			}
			System.out.println("g_best_value= " + g_best_value);

			if (kkk > 1000) {
				int sum = 0;
				for (int c = 0; c < 550; c++) {
					if (Math.abs(hist_g_best_value[kkk - c]
							- hist_g_best_value[kkk - c - 1]) < threshold) {
						sum = sum + 1;
					}
				}

				if (sum > 30) {
					fermati = true;
					break;
				}
			}
			if (kkk > kmax - 2) {
				break;
			}
			// aggiunto per NASH ALBAN
			// if (ModelName.equals("PeakFlow")) {
			// if (g_best_value < 0.05) {
			// break;
			// }
			// }
			costvectold = costvect;
			// System.out.println("ite="+kkk);
			// System.out.println("fermati="+fermati);
			kkk++;

		}
		for (int i = 0; i < g_best.length; i++) {
			System.out.println(g_best[i]);

		}
		System.out.println(g_best_value);

		outOpt = g_best_value;
		outOptSet = g_best;
	}

	public boolean StoppingCondition(double vett[], double s) {

		boolean result = false;
		int cnt = 0;
		for (int i = 0; i < vett.length; i++) {
			if (vett[i] < s) {
				cnt += 1;
			}
		}
		if (cnt == vett.length) {
			result = true;
		}

		return result;
	}

	public double[][] LHSU(double[] xmin, double[] xmax, int nsample) {
		RandomStream randomstream = new MRG32k3a();

		// Latin Hypercube sampling
		// double[][] LHSresult=new double [1][1];
		int nvar = xmin.length;
		double[][] ran = new double[nsample][nvar];
		for (int row = 0; row < ran.length; row++) {
			for (int col = 0; col < ran[0].length; col++) {
				ran[row][col] = Math.random();
			}
		}

		double[] rr;
		double[] P = new double[nsample];
		double[][] s = new double[nsample][nvar];
		for (int j = 0; j < nvar; j++) {
			rr = new double[(nsample)];
			RandomPermutation.init(rr, (nsample));
			// RandomStream rs=new MRG32k3a();
			RandomStream rs = randomstream;
			RandomPermutation.shuffle(rr, rs);
			for (int kk = 0; kk < nsample; kk++) {
				P[kk] = (rr[kk] - ran[kk][j]) / nsample;
				s[kk][j] = xmin[j] + P[kk] * (xmax[j] - xmin[j]);
			}

		}
		// System.out.println("ciao");

		return s;

	}

	public double[][] UNIFORM(double[] xmin, double[] xmax, int nsample) {
		// Latin Hypercube sampling
		// double[][] LHSresult=new double [1][1];
		int nvar = xmin.length;
		double[][] s = new double[nsample][nvar];

		double[][] ran = new double[nsample][nvar];
		for (int row = 0; row < ran.length; row++) {

			for (int col = 0; col < ran[0].length; col++) {

				s[row][col] = ((xmax[col] - xmin[col])) * Math.random()
						+ xmin[col];
			}
		}

		return s;

	}

	public double[] ComputeCostFunction(double xx[][]) throws Exception {
		double[] res = new double[xx.length];

		//

		// ///////////////////////////////////////////////////////////////////////////////////////////////////////
		// ///////////////////////////////////////F I N E A D I G
		// E//////////////////////////////////////////////////
		// ///////////////////////////////////////////////////////////////////////////////////////////////////////
		if (ModelName.equals("SHALSTAB")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {


				OmsShalstab shalstab = new OmsShalstab();
				shalstab.inSlope = slope;
				shalstab.inTca = tca;
				shalstab.pSdepth =  xx[numpart][0];
				shalstab.pTrasmissivity = xx[numpart][1];
				shalstab.pCohesion = 0;
				shalstab.pRho = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];
				shalstab.pm = pm;

				shalstab.process();

				OmsRasterReader reader2 = new OmsRasterReader();
				reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/LAST_FRANE_CUT_CAL";
				reader2.fileNovalue = -9999.0;
				reader2.geodataNovalue = Double.NaN;
				reader2.process();
				GridCoverage2D M2 = reader2.outRaster;

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = M2;
				ab.inModelledMap = shalstab.outShalstab2Classes;
				ab.process();

				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = shalstab.outShalstab2Classes;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/Modeled_CAL_30";
				writer.process();
				// 1:res[numpart] = (1.0 - ab.outOddsRatioSkillScore);
				// 2: res[numpart] = 1.0 - ab.outTruePositiveRate;
				// 3: 
				res[numpart] = 1.0 - ab.outAccuracy;
				// 4: res[numpart] = Math.abs(1.0 - ab.outBiasScore);
				// 5: res[numpart] = (1.0 - ab.outCriticalSuccessIndex);
				// 6:res[numpart] = (1.0 - ab.outEquitableSuccesIndex);
				
				// 7:res[numpart] = (1.0 - ab.outSuccessIndex);
				// 8: res[numpart] = (ab.outD2PerfectClassification);
				// 9:res[numpart] = (1.-ab.outTrueSkillStatistic);
				// 10:res[numpart] = (1. - ab.outHss);
				// 11: res[numpart] = (1. - ab.outAverageIndex);
				if (pPath != null) {

					errestat.append(xx[numpart][0] + " " + xx[numpart][1] + " "
							+ xx[numpart][2] + " " + xx[numpart][3] + " "
							+ xx[numpart][4] + " " + res[numpart] + " "
							+ ab.outOddsRatioSkillScore + " " + ab.outBiasScore
							+ " " + ab.outCriticalSuccessIndex + " "
							+ ab.outEquitableSuccesIndex + " "
							+ ab.outSuccessIndex + " "
							+ ab.outD2PerfectClassification + " "
							+ ab.outTrueSkillStatistic + " " + ab.outHss + " "
							+ ab.outAverageIndex);
					errestat.println();
					errestat.flush();

				}
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("SHALSTABDISTRIBUTED")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {
				OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				s.doP1 = true;
				s.p1 = xx[numpart][4];
				s.p2 = xx[numpart][5];
				s.inSlope = slope;
				s.process();

				GridCoverage2D classiCoverage = s.outSoilDepth;
				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				writer.process();

				OmsRasterReader reader33 = new OmsRasterReader();
				reader33.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/NewA3TESS_CAL";
				reader33.fileNovalue = -9999.0;
				reader33.geodataNovalue = Double.NaN;
				reader33.process();
				GridCoverage2D tessitura = reader33.outRaster;

				OmsRangeLookup range = new OmsRangeLookup();

				range = new OmsRangeLookup();
				range.pm = pm;
				range.inRaster = tessitura;
				range.pRanges = "[0 1],(1.0 2],(2.0 3]";
				range.pClasses = xx[numpart][3] + "," + xx[numpart][3] * 0.85
						+ "," + xx[numpart][3] * 0.7;
				range.process();
				GridCoverage2D frictionAngle = range.outRaster;
				writer.inRaster = frictionAngle;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/frictionAngle";
				writer.process();

				OmsShalstab shalstab = new OmsShalstab();
				shalstab.inSlope = slope;
				shalstab.inTca = tca;
				shalstab.inSdepth = s.outSoilDepth;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = 0;
				shalstab.pRho = xx[numpart][1];
				shalstab.inTgphi = frictionAngle;
				shalstab.pQ = xx[numpart][2];
				shalstab.pm = pm;

				shalstab.process();

				OmsRasterReader reader2 = new OmsRasterReader();
				reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/new_A3_FraneMisurate_CAL";
				reader2.fileNovalue = -9999.0;
				reader2.geodataNovalue = Double.NaN;
				reader2.process();
				GridCoverage2D M2 = reader2.outRaster;

				writer = new OmsRasterWriter();
				writer.inRaster = shalstab.outShalstab2Classes;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_PROVA";
				writer.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = M2;
				ab.inModelledMap = shalstab.outShalstab2Classes;
				ab.process();

				// 1: res[numpart] = 1.0 - ab.outOddsRatioSkillScore;
				// 2: res[numpart] = 1.0 - ab.outTruePositiveRate;
				// 3: res[numpart] = 1.0 - ab.outAccuracy;
				// 4: res[numpart] = Math.abs(1.0 - ab.outBiasScore);
				// 5: res[numpart] = (1.0 - ab.outCriticalSuccessIndex);
				// 6: res[numpart] = (1.0 - ab.outEquitableSuccesIndex);
				// 7: res[numpart] = (1.0 - ab.outSuccessIndex);
				// 8:res[numpart] = (ab.outD2PerfectClassification);
				// 9: res[numpart] = (1.-ab.outTrueSkillStatistic);

				// 10: res[numpart] = (1-ab.outHss);
				// 11
				res[numpart] = (1. - ab.outAverageIndex);
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("SHALSTABPARVAR")) {
			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsShalstab shalstab = new OmsShalstab();
				shalstab.inSlope = slope;
				shalstab.inTca = tca;
				shalstab.pSdepth = xx[numpart][0];
				shalstab.pTrasmissivity = xx[numpart][1];
				shalstab.pCohesion = 0;
				shalstab.pRho = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];
				shalstab.pm = pm;

				shalstab.process();

				OmsRasterReader reader2 = new OmsRasterReader();
				reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/MeasuredForRoc_Cal";
				reader2.fileNovalue = -9999.0;
				reader2.geodataNovalue = Double.NaN;
				reader2.process();
				GridCoverage2D M2 = reader2.outRaster;

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = M2;
				ab.inModelledMap = shalstab.outShalstab2Classes;
				ab.process();

				// 1: res[numpart] = 1.0 - ab.outOddsRatioSkillScore;
				// 2: res[numpart] = 1.0 - ab.outTruePositiveRate;
				// 3: res[numpart] = 1.0 - ab.outAccuracy;
				// 4: res[numpart] = Math.abs(1.0 - ab.outBiasScore);
				// 5: res[numpart] = (1.0 - ab.outCriticalSuccessIndex);
				// 6: res[numpart] = (1.0 - ab.outEquitableSuccesIndex);
				// 7: res[numpart] = (1.0 - ab.outSuccessIndex);
				// 8:res[numpart] = (ab.outD2PerfectClassification);
				// 9: res[numpart] = (1-ab.outTrueSkillStatistic);
				// 10:
				res[numpart] = (1 - ab.outHss);

				System.out.println(res[numpart]);

			}
		}

		if (ModelName.equals("ComputeFS")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				s.doP1 = true;
				s.p1 = xx[numpart][5];
				s.p2 = xx[numpart][6];
				s.inSlope = slope;
				s.process();

				OmsComputeFS shalstab = new OmsComputeFS();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;

				shalstab.inSdepth = s.outSoilDepth;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pRho = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];

				shalstab.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = Math
						.pow(((1 - ab.outTruePositiveRate)
								* (1 - ab.outTruePositiveRate) + (ab.outFalsePositiveRate * ab.outFalsePositiveRate)),
								0.5);
				System.out.println(res[numpart]);

			}
		}

		if (ModelName.equals("ComputeFSconGamma")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				s.doP1 = true;
				s.p1 = xx[numpart][6];
				s.p2 = xx[numpart][7];
				s.inSlope = slope;
				s.process();

				GridCoverage2D classiCoverage = s.outSoilDepth;
				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				writer.process();

				OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;

				shalstab.inSdepth = s.outSoilDepth;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pSoilPorosity = xx[numpart][2];
				// shalstab.pGammaS = xx[numpart][3];
				shalstab.pTgphi = xx[numpart][4];
				shalstab.pQ = xx[numpart][5];

				shalstab.process();

				classiCoverage = shalstab.outShalstab;
				writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_SI_P1";
				writer.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = Math
						.pow(((1 - ab.outTruePositiveRate)
								* (1 - ab.outTruePositiveRate) + (ab.outFalsePositiveRate * ab.outFalsePositiveRate)),
								0.5);
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("SHALSTABLIKE")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				s.doP1 = true;
				s.p1 = xx[numpart][5];
				s.p2 = xx[numpart][6];
				s.inSlope = slope;
				s.process();

				GridCoverage2D classiCoverage = s.outSoilDepth;
				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				writer.process();

				OmsComputeFS_ShalstabLike shalstab = new OmsComputeFS_ShalstabLike();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;

				shalstab.inSdepth = s.outSoilDepth;
				shalstab.inDem = inDem;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.prhoS = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];

				shalstab.process();

				classiCoverage = shalstab.outShalstab;
				writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_SI_P1";
				writer.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = Math
						.pow(((1 - ab.outTruePositiveRate)
								* (1 - ab.outTruePositiveRate) + (ab.outFalsePositiveRate * ab.outFalsePositiveRate)),
								0.5);
				System.out.println(res[numpart]);

			}
		}

		if (ModelName.equals("LASTMODEL_SHALSTAB_Par_CONST")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				// OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				// s.doP1 = true;
				// s.p1 = xx[numpart][7];
				// s.p2 = xx[numpart][8];
				// s.inSlope = slope;
				// s.process();
				//
				// GridCoverage2D classiCoverage = s.outSoilDepth;
				// OmsRasterWriter writer = new OmsRasterWriter();
				// writer.inRaster = classiCoverage;
				// writer.file =
				// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				// writer.process();

				OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;
				shalstab.inDem = inDem;

				// shalstab.inSdepth = s.outSoilDepth;
				shalstab.pSdepth = xx[numpart][5];
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pSoilPorosity = 0.5;
				shalstab.pTgphi = xx[numpart][2];
				shalstab.pQ = xx[numpart][3];
				shalstab.pDegreeOfSat =0.5;
				shalstab.pGs =  xx[numpart][4];

				shalstab.pMode = 0;
				// Per ROSSO
				// shalstab.pMode = 1;
				// shalstab.pDuration = xx[numpart][8];
				shalstab.process();

				OmsRasterReader reader2 = new OmsRasterReader();
				reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/LAST_FRANE_CUT_CAL";
				reader2.fileNovalue = -9999.0;
				reader2.geodataNovalue = Double.NaN;
				reader2.process();
				GridCoverage2D M2 = reader2.outRaster;

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = M2;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = shalstab.outShalstab;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/Modeled_CAL_30";
				writer.process();

				// 1:res[numpart] = 1.0 - ab.outOddsRatioSkillScore;
				// 2: res[numpart] = 1.0 - ab.outTruePositiveRate;
				// 3: res[numpart] = 1.0 - ab.outAccuracy;
				// 4: res[numpart] = Math.abs(1.0 - ab.outBiasScore);
				// 5:res[numpart] = (1.0 - ab.outCriticalSuccessIndex);
				// 6:res[numpart] = (1.0 - ab.outEquitableSuccesIndex);
				// 7:res[numpart] = (1.0 - ab.outSuccessIndex);
				// 8:res[numpart] = (ab.outD2PerfectClassification);

				// 9:res[numpart] = (1.-ab.outTrueSkillStatistic);
				// 10:res[numpart] = (1 - ab.outHss);
				// 11:
				res[numpart] = (1. - ab.outAverageIndex);
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("LASTMODELDISTRIBUTED")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				s.doP1 = true;
				s.p1 = xx[numpart][5];
				s.p2 = xx[numpart][6];
				s.inSlope = slope;
				s.process();

				GridCoverage2D classiCoverage = s.outSoilDepth;
				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				writer.process();

				OmsRasterReader reader33 = new OmsRasterReader();
				reader33.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/NewA3TESS_CAL";
				reader33.fileNovalue = -9999.0;
				reader33.geodataNovalue = Double.NaN;
				reader33.process();
				GridCoverage2D tessitura = reader33.outRaster;

				OmsRangeLookup range = new OmsRangeLookup();
				range.pm = pm;
				range.inRaster = tessitura;
				range.pRanges = "[0 1],(1.0 2],(2.0 3]";
				range.pClasses = xx[numpart][2] + "," + xx[numpart][2] * 2.0
						+ "," + xx[numpart][2] * 3.5;
				range.process();
				GridCoverage2D cohesion = range.outRaster;

				writer.inRaster = cohesion;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/cohesion";
				writer.process();

				range = new OmsRangeLookup();
				range.pm = pm;
				range.inRaster = tessitura;
				range.pRanges = "[0 1],(1.0 2],(2.0 3]";
				range.pClasses = xx[numpart][3] + "," + xx[numpart][3] * 0.85
						+ "," + xx[numpart][3] * 0.7;
				range.process();
				GridCoverage2D frictionAngle = range.outRaster;
				writer.inRaster = frictionAngle;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/frictionAngle";
				writer.process();

				range = new OmsRangeLookup();
				range.pm = pm;
				range.inRaster = tessitura;
				range.pRanges = "[0 1],(1.0 2],(2.0 3]";
				range.pClasses = xx[numpart][4] + "," + xx[numpart][4] * 0.001
						+ "," + xx[numpart][4] * 0.0000001;
				range.process();
				GridCoverage2D conductivity = range.outRaster;

				writer.inRaster = conductivity;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/conductivity";
				writer.process();

				OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;
				shalstab.inDem = inDem;

				shalstab.inSdepth = s.outSoilDepth;
				shalstab.inTrasmissivity = conductivity;
				shalstab.computeTransmissivity = true;
				shalstab.inCohesion = cohesion;
				// shalstab.pCohesion = 0.0;
				shalstab.inTgphi = frictionAngle;
				shalstab.inSoilPorosity = porosityMap;
				shalstab.pQ = xx[numpart][0];
				shalstab.pDegreeOfSat = xx[numpart][1];
				shalstab.pGs = 2.65;

				shalstab.pMode = 0;
				// Per ROSSO
				// shalstab.pMode = 1;
				// shalstab.pDuration = xx[numpart][8];
				shalstab.process();

				classiCoverage = shalstab.outShalstab;
				writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_PROVA";
				writer.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				// 1: res[numpart] = 1.0 - ab.outOddsRatioSkillScore;
				// 2: res[numpart] = 1.0 - ab.outTruePositiveRate;
				// 3:
				res[numpart] = 1.0 - ab.outAccuracy;
				// 4: res[numpart] = Math.abs(1.0 - ab.outBiasScore);
				// 5: res[numpart] = (1.0 - ab.outCriticalSuccessIndex);
				// 6: res[numpart] = (1.0 - ab.outEquitableSuccesIndex);
				// 7: res[numpart] = (1.0 - ab.outSuccessIndex);
				// 8:res[numpart] = (ab.outD2PerfectClassification);
				// 9:res[numpart] = (1 - ab.outTrueSkillStatistic);
				// 10:res[numpart] = (1 - ab.outHss);
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("ROSSODISTRIBUTED")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				s.doP1 = true;
				s.p1 = xx[numpart][10];
				s.p2 = xx[numpart][11];
				s.inSlope = slope;
				s.process();

				GridCoverage2D classiCoverage = s.outSoilDepth;
				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				writer.process();

				OmsRasterReader reader33 = new OmsRasterReader();
				reader33.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Tess_A3_CAL";
				reader33.fileNovalue = -9999.0;
				reader33.geodataNovalue = Double.NaN;
				reader33.process();
				GridCoverage2D tessitura = reader33.outRaster;

				OmsRangeLookup range = new OmsRangeLookup();
				range.pm = pm;
				range.inRaster = tessitura;
				range.pRanges = "[0 1],(1.0 2],(2.0 3]";
				range.pClasses = xx[numpart][4] + "," + xx[numpart][5] + ","
						+ xx[numpart][6];
				range.process();
				GridCoverage2D cohesion = range.outRaster;

				range = new OmsRangeLookup();
				range.pm = pm;
				range.inRaster = tessitura;
				range.pRanges = "[0 1],(1.0 2],(2.0 3]";
				range.pClasses = xx[numpart][7] + "," + xx[numpart][8] + ","
						+ xx[numpart][9];
				range.process();
				GridCoverage2D frictionAngle = range.outRaster;
				//
				// range = new OmsRangeLookup();
				// range.pm = pm;
				// range.inRaster = tessitura;
				// range.pRanges = "[0 1],(1.0 2],(2.0 3]";
				// range.pClasses = xx[numpart][4] + "," + xx[numpart][4] * 2
				// + "," + xx[numpart][4] * 30;
				// range.process();
				// GridCoverage2D conductivity = range.outRaster;

				OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;
				shalstab.inDem = inDem;

				shalstab.inSdepth = s.outSoilDepth;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.computeTransmissivity = false;
				shalstab.inCohesion = cohesion;
				// shalstab.pCohesion = 0.0;
				shalstab.inTgphi = frictionAngle;
				shalstab.inSoilPorosity = porosityMap;
				shalstab.pQ = xx[numpart][1];
				shalstab.pDegreeOfSat = xx[numpart][2];
				shalstab.pGs = 2.65;

				// shalstab.pMode = 0;
				// Per ROSSO
				shalstab.pMode = 1;
				shalstab.pDuration = xx[numpart][3];
				shalstab.process();

				// classiCoverage = shalstab.outShalstab;
				// writer = new OmsRasterWriter();
				// writer.inRaster = classiCoverage;
				// writer.file =
				// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_PROVA";
				// writer.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				// 1: res[numpart] = 1.0 - ab.outOddsRatioSkillScore;
				// 2: res[numpart] = 1.0 - ab.outTruePositiveRate;
				// 3: res[numpart] = 1.0 - ab.outAccuracy;
				// 4: res[numpart] = Math.abs(1.0 - ab.outBiasScore);
				// 5: res[numpart] = (1.0 - ab.outCriticalSuccessIndex);
				// 6: res[numpart] = (1.0 - ab.outEquitableSuccesIndex);
				// 7: 
				res[numpart] = (1.0 - ab.outSuccessIndex);
				// 8:res[numpart] = (ab.outD2PerfectClassification);
				// 9:res[numpart] = (1 - ab.outTrueSkillStatistic);
				// 10:res[numpart] = (1 - ab.outHss);
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("LASTMODEL_ROSSO_Par_CONST")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				// OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				// s.doP1 = true;
				// s.p1 = xx[numpart][7];
				// s.p2 = xx[numpart][8];
				// s.inSlope = slope;
				// s.process();
				//
				// GridCoverage2D classiCoverage = s.outSoilDepth;
				// OmsRasterWriter writer = new OmsRasterWriter();
				// writer.inRaster = classiCoverage;
				// writer.file =
				// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				// writer.process();

				OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;
				shalstab.inDem = inDem;

				// shalstab.inSdepth = s.outSoilDepth;
				shalstab.pSdepth = xx[numpart][5];
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pSoilPorosity = 0.5;
				shalstab.pTgphi = xx[numpart][2];
				shalstab.pQ = xx[numpart][3];
				shalstab.pDegreeOfSat =0.5;
				shalstab.pGs =  xx[numpart][4];
				shalstab.pDuration = xx[numpart][6];

				// shalstab.pMode = 0;
				// Per ROSSO
				shalstab.pMode = 1;

				shalstab.process();

				OmsRasterReader reader2 = new OmsRasterReader();
				reader2.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/LAST_FRANE_CUT_CAL";
				reader2.fileNovalue = -9999.0;
				reader2.geodataNovalue = Double.NaN;
				reader2.process();
				GridCoverage2D M2 = reader2.outRaster;

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = M2;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = shalstab.outShalstab;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/SHALSTAB_SENSITIVITY/LastSimulations/Loc_32633_cal/Mapset30m/cell/Modeled_CAL_30";
				writer.process();

				// 1:res[numpart] = 1.0 - ab.outOddsRatioSkillScore;
				// 2: res[numpart] = 1.0 - ab.outTruePositiveRate;
				// 3: res[numpart] = 1.0 - ab.outAccuracy;
				// 4: res[numpart] = Math.abs(1.0 - ab.outBiasScore);
				// 5:res[numpart] = (1.0 - ab.outCriticalSuccessIndex);
				// 6:res[numpart] = (1.0 - ab.outEquitableSuccesIndex);
				// 7:res[numpart] = (1.0 - ab.outSuccessIndex);
				// 8:
				res[numpart] = (ab.outD2PerfectClassification);
				// 9: res[numpart] = (1-ab.outTrueSkillStatistic);
				// 10:res[numpart] = (1 - ab.outHss);
				// 11:res[numpart] = (1. - ab.outAverageIndex);
				System.out.println(res[numpart]);

			}
		}

		if (ModelName.equals("LASTMODEL")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				s.doP1 = true;
				s.p1 = xx[numpart][6];
				s.p2 = xx[numpart][7];
				s.inSlope = slope;
				s.process();

				GridCoverage2D classiCoverage = s.outSoilDepth;
				OmsRasterWriter writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				writer.process();

				OmsComputeFS_GSE shalstab = new OmsComputeFS_GSE();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;
				shalstab.pGs = 2.5;

				shalstab.inSdepth = s.outSoilDepth;
				shalstab.inDem = inDem;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pSoilPorosity = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];
				shalstab.pDegreeOfSat = xx[numpart][5];
				shalstab.pMode = 0;
				// Per ROSSO
				shalstab.pMode = 1;
				shalstab.pDuration = xx[numpart][8];
				shalstab.process();

				classiCoverage = shalstab.outShalstab;
				writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_SI_P1";
				writer.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = Math
						.pow(((1 - ab.outTruePositiveRate)
								* (1 - ab.outTruePositiveRate) + (ab.outFalsePositiveRate * ab.outFalsePositiveRate)),
								0.5);
				System.out.println(res[numpart]);

			}
		}

		if (ModelName.equals("LU2008MODEL")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {
				//
				// OmsComputeSoilThickness s = new OmsComputeSoilThickness();
				// s.doP1 = true;
				// s.p1 = xx[numpart][6];
				// s.p2 = xx[numpart][7];
				// s.inSlope = slope;
				// s.process();

				// GridCoverage2D classiCoverage = s.outSoilDepth;
				OmsRasterWriter writer = new OmsRasterWriter();
				// writer.inRaster = classiCoverage;
				// writer.file =
				// "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Soil";
				// writer.process();
				OmsComputeFS_LU2008 shalstab = new OmsComputeFS_LU2008();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;

				// shalstab.inSdepth = s.outSoilDepth;
				shalstab.inDem = inDem;
				shalstab.pGammaSoil = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pnVG = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pKs = xx[numpart][4];
				shalstab.pQ = xx[numpart][5];

				shalstab.pSdepth = xx[numpart][6];
				shalstab.pHwt = xx[numpart][7];
				shalstab.palphaVG = xx[numpart][8];
				shalstab.process();

				GridCoverage2D classiCoverage = shalstab.outShalstab;
				writer = new OmsRasterWriter();
				writer.inRaster = classiCoverage;
				writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_SI_P1";
				writer.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = Math
						.pow(((1 - ab.outTruePositiveRate)
								* (1 - ab.outTruePositiveRate) + (ab.outFalsePositiveRate * ab.outFalsePositiveRate)),
								0.5);
				System.out.println(res[numpart]);

			}
		}

		if (ModelName.equals("ComputeFSContantsoil")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeFS shalstab = new OmsComputeFS();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;

				shalstab.pSdepth = 1.50;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pRho = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];
				// shalstab.p_SoilAlpha=xx[numpart][5];
				// shalstab.p_SoilBeta=xx[numpart][6];
				shalstab.process();

				shalstab.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = Math
						.pow(((1 - ab.outTruePositiveRate)
								* (1 - ab.outTruePositiveRate) + (ab.outFalsePositiveRate * ab.outFalsePositiveRate)),
								0.5);
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("ComputeFS_Success_Index")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeFS shalstab = new OmsComputeFS();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;

				shalstab.inSdepth = sdepth;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pRho = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];
				shalstab.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = 1 - ab.outSuccessIndex;
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("ComputeFS_Critical_Success_Index")) {

			for (int numpart = 0; numpart < xx.length; numpart++) {

				OmsComputeFS shalstab = new OmsComputeFS();
				shalstab.inSlope = slope;
				shalstab.pThresholdFS = 1.0;
				shalstab.inTca = tca;
				shalstab.inSdepth = sdepth;
				shalstab.pTrasmissivity = xx[numpart][0];
				shalstab.pCohesion = xx[numpart][1];
				shalstab.pRho = xx[numpart][2];
				shalstab.pTgphi = xx[numpart][3];
				shalstab.pQ = xx[numpart][4];

				shalstab.process();

				ROC ab = new ROC();
				ab.flagFalseModelled = 1;
				ab.flagTrueModelled = 0;
				ab.flagNegativeMeas = 1;
				ab.flagPositiveMeas = 0;
				ab.inMeasuredMap = inMeasuredForRoc;
				ab.inModelledMap = shalstab.outShalstab;
				ab.process();

				res[numpart] = 1 - ab.outCriticalSuccessIndex;
				System.out.println(res[numpart]);

			}
		}
		if (ModelName.equals("ComputeFSROC")) {
			double a = 0.01;
			double b = 4.5;
			double nk = (b - a) / 100;

			for (int numpart = 0; numpart < xx.length; numpart++) {
				int ccc = 0;
				double[] vetFp = new double[100];
				double[] vetTp = new double[100];
				double fsThreshold = a;
				while (ccc < 100) {

					OmsComputeFS shalstab = new OmsComputeFS();
					shalstab.inSlope = slope;
					shalstab.pThresholdFS = fsThreshold;
					shalstab.inTca = tca;
					shalstab.inSdepth = sdepth;
					shalstab.pTrasmissivity = xx[numpart][0];
					shalstab.pCohesion = xx[numpart][1];
					shalstab.pRho = xx[numpart][2];
					shalstab.pTgphi = xx[numpart][3];
					shalstab.pQ = xx[numpart][4];

					shalstab.process();
					GridCoverage2D classiCoverage = shalstab.outShalstab;
					OmsRasterWriter writer = new OmsRasterWriter();
					writer.inRaster = classiCoverage;
					writer.file = "/Users/giuseppeformetta/Desktop/UNICAL/MM2/m/cell/Classi_Compute_FS_CAL";
					writer.process();

					ROC ab = new ROC();
					ab.flagFalseModelled = 1;
					ab.flagTrueModelled = 0;
					ab.flagNegativeMeas = 1;
					ab.flagPositiveMeas = 0;
					ab.inMeasuredMap = inMeasuredForRoc;
					ab.inModelledMap = shalstab.outShalstab;
					ab.process();

					vetFp[ccc] = ab.outFalsePositiveRate;
					vetTp[ccc] = ab.outTruePositiveRate;
					ccc++;
					fsThreshold += nk;
				}
				double sum = 0;
				for (int i = 1; i < vetFp.length; i++) {
					sum += (vetFp[i] - vetFp[i - 1])
							* (vetTp[i] + vetTp[i - 1]) * 0.5;

				}

				// res[numpart] = Math
				// .pow(((1 - ab.outTruePositiveRate)
				// * (1 - ab.outTruePositiveRate) + (ab.outFalsePositiveRate *
				// ab.outFalsePositiveRate)),
				// 0.5);
				res[numpart] = 1 - sum;
				System.out.println(res[numpart]);

			}
		}
		return res;

	}

	public double[] Compute_gBest(double xx[][], double[] vettcostnew) {
		double re[] = g_best;
		int pos = 0;
		double min = g_best_value;
		for (int i = 0; i < vettcostnew.length; i++) {
			if (Math.abs(vettcostnew[i]) <= min) {
				g_best_value = Math.abs(vettcostnew[i]);
				min = Math.abs(vettcostnew[i]);
				pos = i;
				for (int ii = 0; ii < xx[0].length; ii++) {
					re[ii] = xx[pos][ii];
				}
			}
		}

		// System.out.println("minimo="+g_best_value);
		return re;
	}

	public double[][] Compute_particle(double pos[][], double[][] vel) {

		double xnew[][] = new double[pos.length][pos[0].length];
		for (int i = 0; i < vel.length; i++) {
			for (int j = 0; j < vel[0].length; j++) {
				xnew[i][j] = pos[i][j] + vel[i][j];
			}

		}
		double[][] xneww = ReflectBounds(xnew, ParRange_maxn, ParRange_minn);
		return xneww;
	}

	public double[][] Compute_velocity(double pos[][], double[][] vel) {

		double velnew[][] = new double[pos.length][pos[0].length];
		for (int i = 0; i < vel.length; i++) {
			for (int j = 0; j < vel[0].length; j++) {
				double c1 = 1.5;
				double r1 = Math.random();
				double c2 = 2.5;
				double r2 = Math.random();
				double inertia = 0.5;
				velnew[i][j] = inertia * vel[i][j] + c1 * r1
						* (p_best[i][j] - pos[i][j]) + c2 * r2
						* (g_best[j] - pos[i][j]);
			}

		}
		return velnew;
	}

	public double[][] Compute_pBest(double currentpos[][], double[] currentbest) {
		double pos_best[][] = p_best;
		for (int i = 0; i < currentbest.length; i++) {
			// per tutti
			if (Math.abs(currentbest[i]) < Math.abs(p_best_value[i])) {
				// per nash
				// if(Math.abs(currentbest[i])>Math.abs(p_best_value[i])){
				p_best_value[i] = Math.abs(currentbest[i]);
				for (int j = 0; j < currentpos[0].length; j++) {
					pos_best[i][j] = currentpos[i][j];
				}
			}
		}
		return pos_best;
	}

	public double[][] ReflectBounds(double[][] neww, double[] ParRange_maxnn,
			double[] ParRange_minnn) {
		// ParRange_maxnn=new double []{15000,50,0.99,0.3,0.9};
		// ParRange_minnn=new double []{1,0.01,0.01,0.0001,0.01};
		// Checks the bounds of the parameters
		// First determine the size of new
		// int nmbOfIndivs=neww.length;
		// int Dim=neww[0].length;
		double[][] y = neww;
		for (int row = 0; row < neww.length; row++) {
			for (int col = 0; col < neww[0].length; col++) {
				if (y[row][col] < ParRange_minnn[col]) {
					y[row][col] = 2 * ParRange_minnn[col] - y[row][col];
				}
				if (y[row][col] > ParRange_maxnn[col]) {
					y[row][col] = 2 * ParRange_maxnn[col] - y[row][col];
				}
			}
		}
		for (int row = 0; row < neww.length; row++) {
			for (int col = 0; col < neww[0].length; col++) {
				if (y[row][col] < ParRange_minn[col]) {
					y[row][col] = ParRange_minn[col] + Math.random()
							* (ParRange_maxn[col] - ParRange_minn[col]);
				}
				if (y[row][col] > ParRange_maxn[col]) {
					y[row][col] = ParRange_minn[col] + Math.random()
							* (ParRange_maxn[col] - ParRange_minn[col]);
				}
			}
		}

		return y;
	}

}