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
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.boussinesq;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_KEYWORDS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_LABEL;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_LICENSE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_NAME;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSSHALSTAB_STATUS;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;

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

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;

@Description(OMSSHALSTAB_DESCRIPTION)
@Author(name = "Formetta Serafin Cordano")
@Keywords(OMSSHALSTAB_KEYWORDS)
@Label(OMSSHALSTAB_LABEL)
@Name(OMSSHALSTAB_NAME)
@Status(OMSSHALSTAB_STATUS)
@License(OMSSHALSTAB_LICENSE)
public class ComputeBeqInput extends JGTModel {

	@Description("Digital Elevation Model of the bedrock")
	@In
	public GridCoverage2D inDtm = null;

	@Description("Porosity Map")
	@In
	public GridCoverage2D inPorosity = null;

	@Description("Porosity Map")
	@In
	public double valueSourceV = -9999;

	@Description("Porosity Map")
	@In
	public double valueEtaInitialCondVtoAddtoDem = -9999;

	@Description("Porosity Map")
	@In
	public double valuePorosityV = -9999;

	@Description("Porosity Map")
	@In
	public double valueCV = -9999;
	@Description("Porosity Map")
	@In
	public double valueMV = -9999;
	@Description("Porosity Map")
	@In
	public double valueEtaDricheletV = -9999;

	@Description("If true do simmetric")
	@In
	public boolean doSimmetricMatrix = true;

	@Description("If true do simmetric")
	@In
	public String pPath = null;

	@Description("Numbered Polygon")
	@Out
	public GridCoverage2D outShalstab = null;

	@Description("Numbered Polygon in rcform")
	@Out
	public int[] Mp = null;

	@Description("Row pointer  in rcform ")
	@Out
	public int[] Mj = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] Ml = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vSource = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vEtaInitialCondV = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vPorosity = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vC = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vM = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vEtaDrichelet = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vBedrock = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vPlanarArea = null;

	@Description("Shared sides  in rcform")
	@Out
	public double[] vHydraulicCondictivity = null;

	public final double EPS = 0.01;

	@Execute
	public void process() throws Exception {
		// if (!concatOr(outShalstab == null, doReset)) {
		// return;
		// }
		checkNull(inDtm);

		RenderedImage slopeRI = inDtm.getRenderedImage();

		qcrit(slopeRI);
	}

	/**
	 * Calculates the trasmissivity in every pixel of the map.
	 * 
	 * @throws IOException
	 */
	private void qcrit(RenderedImage slope) throws IOException {

		HashMap<String, Double> regionMap = CoverageUtilities
				.getRegionParamsFromGridCoverage(inDtm);
		int cols = regionMap.get(CoverageUtilities.COLS).intValue();
		int rows = regionMap.get(CoverageUtilities.ROWS).intValue();
		WritableRaster classiWR = CoverageUtilities.createDoubleWritableRaster(
				cols, rows, null, null, null);
		WritableRandomIter classiIter = RandomIterFactory.createWritable(
				classiWR, null);
		RandomIter slopeRI = RandomIterFactory.create(slope, null);
		int contaPoligoni = 0;
		LinkedList<Integer> rowP = new LinkedList<Integer>();
		LinkedList<Integer> colP = new LinkedList<Integer>();
		LinkedList<Integer> valueP = new LinkedList<Integer>();

		LinkedList<Double> valuePlanArea = new LinkedList<Double>();
		LinkedList<Double> valueSource = new LinkedList<Double>();
		LinkedList<Double> valueEtaInitialCond = new LinkedList<Double>();
		LinkedList<Double> valuePorosity = new LinkedList<Double>();
		LinkedList<Double> valueC = new LinkedList<Double>();
		LinkedList<Double> valueM = new LinkedList<Double>();
		LinkedList<Double> valueEtaDrcihelet = new LinkedList<Double>();
		LinkedList<Double> valueBedrock = new LinkedList<Double>();

		for (int j = 0; j < rows; j++) {
			for (int i = 0; i < cols; i++) {

				double dtmValue = slopeRI.getSampleDouble(i, j, 0);
				double valueEtaInitialCondV = dtmValue
						+ valueEtaInitialCondVtoAddtoDem;

				if (!isNovalue(dtmValue)) {
					contaPoligoni += 1;
					rowP.add(j);
					colP.add(i);
					valueP.add(contaPoligoni);
					valueC.add(valueCV);
					valueM.add(valueMV);
					valueEtaDrcihelet.add(valueEtaDricheletV);
					valueEtaInitialCond.add(valueEtaInitialCondV);
					valuePlanArea.add(regionMap.get(CoverageUtilities.XRES)
							.doubleValue()
							* regionMap.get(CoverageUtilities.YRES)
									.doubleValue());
					valuePorosity.add(valuePorosityV);
					valueSource.add(valueSourceV);
					valueBedrock.add(dtmValue);

					classiIter.setSample(i, j, 0, contaPoligoni);

				} else {
					classiIter.setSample(i, j, 0, doubleNovalue);
				}
			}
		}

		vSource = new double[valueP.size()];

		vEtaInitialCondV = new double[valueP.size()];

		vPorosity = new double[valueP.size()];

		vC = new double[valueP.size()];

		vM = new double[valueP.size()];

		vEtaDrichelet = new double[valueP.size()];
		vBedrock = new double[valueP.size()];
		vPlanarArea = new double[valueP.size()];

		for (int kkk = 0; kkk < vPlanarArea.length; kkk++) {
			vSource[kkk] = valueSource.get(kkk);
			vEtaInitialCondV[kkk] = valueEtaInitialCond.get(kkk);
			vPorosity[kkk] = valuePorosity.get(kkk);
			vC[kkk] = valueC.get(kkk);
			vM[kkk] = valueM.get(kkk);
			vEtaDrichelet[kkk] = valueEtaDrcihelet.get(kkk);
			vBedrock[kkk] = valueBedrock.get(kkk);
			vPlanarArea[kkk] = valuePlanArea.get(kkk);

		}

		// pm.done();

		outShalstab = CoverageUtilities.buildCoverage("classi", classiWR,
				regionMap, inDtm.getCoordinateReferenceSystem());

		LinkedList<Integer> rowL = new LinkedList<Integer>();
		LinkedList<Integer> colL = new LinkedList<Integer>();
		LinkedList<Integer> valueL = new LinkedList<Integer>();

		RenderedImage p = outShalstab.getRenderedImage();
		RandomIter pRI = RandomIterFactory.create(p, null);

		for (int j = 0; j < rows; j++) {
			for (int i = 0; i < cols; i++) {
				System.out.println(pRI.getSampleDouble(i, j, 0));
				if (!isNovalue(pRI.getSampleDouble(i, j, 0))) {
					if ((i + 1) < cols) {

						if (pRI.getSampleDouble(i + 1, j, 0) != doubleNovalue) {

							rowL.add((int) pRI.getSampleDouble(i, j, 0));
							colL.add((int) pRI.getSampleDouble(i + 1, j, 0));
							valueL.add(2);

						}
					}
					if ((i - 1) > 0) {
						if (pRI.getSampleDouble(i - 1, j, 0) != doubleNovalue) {

							rowL.add((int) pRI.getSampleDouble(i, j, 0));
							colL.add((int) pRI.getSampleDouble(i - 1, j, 0));
							valueL.add(2);

						}
					}
					if ((j + 1) < rows) {
						if (pRI.getSampleDouble(i, j + 1, 0) != doubleNovalue) {

							rowL.add((int) pRI.getSampleDouble(i, j, 0));
							colL.add((int) pRI.getSampleDouble(i, j + 1, 0));

							valueL.add(2);

						}
					}
					if ((j - 1) > 0) {
						if (pRI.getSampleDouble(i, j - 1, 0) != doubleNovalue) {

							rowL.add((int) pRI.getSampleDouble(i, j, 0));
							colL.add((int) pRI.getSampleDouble(i, j - 1, 0));
							valueL.add(2);

						}
					}

				}
			}
		}

		if (doSimmetricMatrix) {
			for (int i = 0; i < colL.size(); i++) {
				if (rowL.get(i) > colL.get(i)) {
					System.out.println(rowL.get(i));
					System.out.println(colL.get(i));

					// colL.remove(i);
					// rowL.remove(i);
					// valueL.remove(i);
					valueL.set(i, -9999);
					// rowL.set(i, -9999);
					// colL.set(i, -9999);

				}
			}

		}
		int cnt = 0;
		for (int i = 0; i < valueL.size(); i++) {
			if (valueL.get(i) != -9999) {
				cnt++;
				valueL.set(i, cnt);
			}

		}

		LinkedList<Integer> rowLnew = new LinkedList<Integer>();
		LinkedList<Integer> colLnew = new LinkedList<Integer>();
		LinkedList<Integer> valueLnew = new LinkedList<Integer>();

		if (doSimmetricMatrix) {
			for (int i = 0; i < valueL.size(); i++) {
				if (valueL.get(i) != -9999) {
					rowLnew.add(rowL.get(i));
					colLnew.add(colL.get(i));
					valueLnew.add(valueL.get(i));

					rowLnew.add(colL.get(i));
					colLnew.add(rowL.get(i));
					valueLnew.add(valueL.get(i));

				}

			}

		} else {
			rowLnew = rowL;
			colLnew = colL;
			valueLnew = valueL;
		}

		for (int i = 1; i < contaPoligoni + 1; i++) {
			rowLnew.add(i);
			colLnew.add(i);
			valueLnew.add(-1);
		}

		int colLnewVect[] = new int[colLnew.size()];

		int rowLnewVect[] = new int[colLnew.size()];
		double[] valueLnewVect = new double[colLnew.size()];

		for (int i = 0; i < valueLnewVect.length; i++) {
			colLnewVect[i] = colLnew.get(i);
			rowLnewVect[i] = rowLnew.get(i);
			valueLnewVect[i] = valueLnew.get(i);

		}
		computeMp(colLnewVect, rowLnewVect, valueLnewVect, contaPoligoni);

		if (pPath != null && pPath.length() > 0) {
			String s = (pPath + "/vSource");
			FileWriter Rstatfile = new FileWriter(s);
			PrintWriter errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vSource.length; i++) {
				errestat.println(vSource[i]);
			}
			Rstatfile.close();

			s = (pPath + "/vEtaInitialCondV");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vEtaInitialCondV.length; i++) {
				errestat.println(vEtaInitialCondV[i]);
			}
			Rstatfile.close();

			s = (pPath + "/vPorosity");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vPorosity.length; i++) {
				errestat.println(vPorosity[i]);
			}
			Rstatfile.close();

			s = (pPath + "/vC");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vC.length; i++) {
				errestat.println(vC[i]);
			}
			Rstatfile.close();

			s = (pPath + "/vM");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vM.length; i++) {
				errestat.println(vM[i]);
			}
			Rstatfile.close();

			s = (pPath + "/vEtaDrichelet");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vEtaDrichelet.length; i++) {
				errestat.println(vEtaDrichelet[i]);
			}
			Rstatfile.close();

			s = (pPath + "/vBedrock");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vBedrock.length; i++) {
				errestat.println(vBedrock[i]);
			}
			Rstatfile.close();

			s = (pPath + "/vPlanarArea");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < vPlanarArea.length; i++) {
				errestat.println(vPlanarArea[i]);
			}
			Rstatfile.close();

			s = (pPath + "/Mp");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < Mp.length; i++) {
				errestat.println(Mp[i]);
			}
			Rstatfile.close();

			s = (pPath + "/Mj");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < Mj.length; i++) {
				errestat.println(Mj[i]);
			}
			Rstatfile.close();

			s = (pPath + "/Ml");
			Rstatfile = new FileWriter(s);
			errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < Ml.length; i++) {
				errestat.println(Ml[i]);
			}
			Rstatfile.close();

		}

	}

	/**
	 * Compute Mp.
	 * 
	 * @desc
	 * 
	 * @param Mj
	 *            the mj
	 * @param Np
	 *            the np
	 * @return the int[]
	 */
	public void computeMp(int[] pMi, int[] pMj, double[] pMl, int Np) {

		Mp = new int[Np + 1];
		int[] Mi = new int[pMi.length];
		Mj = new int[pMj.length];
		Ml = new double[pMl.length];

		int index = 0;
		int min = 0;
		boolean intoLoop = false;
		// Mp[index] = 0;
		// index++;

		ArrayList<Integer> alMi = new ArrayList<Integer>();
		ArrayList<Integer> alMj = new ArrayList<Integer>();
		ArrayList<Double> alMl = new ArrayList<Double>();

		for (int i = 0; i < Np; i++) {

			for (int j = 0; j < pMi.length; j++) {

				if (pMi[j] == (i + 1)) {

					alMi.add(pMi[j]);
					alMl.add(pMl[j]);
					alMj.add(pMj[j]);

					intoLoop = true;

				}
			}

			if (intoLoop) {

				Mp[i] = index;
				intoLoop = false;

				// System.out.println(index);

				while (alMi.size() > 0) {

					min = alMj.indexOf(Collections.min(alMj));
					Mi[index] = alMi.get(min) - 1;
					Mj[index] = alMj.get(min) - 1;
					Ml[index] = alMl.get(min);

					index++;

					alMi.remove(min);
					alMj.remove(min);
					alMl.remove(min);

				}

				// alMi.clear();
				// alMj.clear();
				// alMl.clear();

			}

		}

		Mp[Np] = pMi.length;

		// return Mp;

	}
}
