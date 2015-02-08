/*
 * JGrass - Free Open Source Java GIS http://www.jgrass.org 
 * 
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Library General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Library General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Library General Public License
 * along with this library; if not, write to the Free Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.jgrasstools.hortonmachine.modules.hydrogeomorphology.boussinesq;

import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_AUTHORCONTACTS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_KEYWORDS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_LABEL;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_LICENSE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_NAME;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_STATUS;

import java.io.IOException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

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

import org.jgrasstools.gears.libs.modules.JGTModel;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

@Description(OMSADIGE_DESCRIPTION)
@Author(name = "Serafin, Formetta, Cordano", contact = OMSADIGE_AUTHORCONTACTS)
@Keywords(OMSADIGE_KEYWORDS)
@Label(OMSADIGE_LABEL)
@Name(OMSADIGE_NAME)
@Status(OMSADIGE_STATUS)
@License(OMSADIGE_LICENSE)
public class OmsBoussinesq extends JGTModel {

	@Description("Vector the planar Area of each polygon")
	@In
	@Unit("m2")
	public double[] planArea;

	@Description("Integration timestep in second")
	@In
	@Unit("s")
	public int TIMESTEP;

	@Description("Vector of piezometric head in the Drichlet cells: all novalues but the Drichlet cells")
	@Unit("m")
	@In
	public double[] etaDrichelet;

	@Description("Vector of piezometric head in the cells")
	@Unit("m")
	@In
	public double[] eta;

	@Description("No value")
	@Unit("-")
	@In
	public double NOVALUE;

	@Description("Row pointer vector")
	@Unit("-")
	@In
	public int[] Mi = null;

	@Description("Index cells vector")
	@Unit("-")
	@In
	public int[] Mp = null;

	@Description("Shared sides")
	@Unit("-")
	@In
	public double[] Ml = null;

	@Description("Euclidean distances between cells")
	@Unit("m")
	@In
	public double[] euclideanDistance = null;

	@Description("Hydr Conductivity on the sides")
	@Unit("m/s")
	@In
	public double[] hydrConductivity = null;

	@Description("length of shared sides")
	@Unit("m")
	@In
	public double[] lengthSides = null;

	@Description("Bottom elevation")
	@Unit("m")
	@In
	public double[] bedrockElevation = null;

	@Description("Cell Porosity")
	@Unit("-")
	@In
	public double[] porosity = null;

	@Description("Rainfall-Source term")
	@Unit("m/s")
	@In
	public double[] source = null;

	@Description("OutFlow multiplicative factor vector")
	@Unit("m^(1-d)/s")
	@In
	public double[] c = null;

	@Description("OutFlow exponent vector")
	@Unit("-")
	@In
	public double[] d = null;

	@Description("OutFlow exponent vector")
	@Unit("-")
	@Out
	public double[] outflow = null;

	@Description("Aquifer Thickness")
	@Unit("m")
	@Out
	public double[] aquiferThickness = null;

	@SuppressWarnings("nls")
	public static String boundaryConditions;
	public int Np = planArea.length;

	public double[] volumeSource;

	protected int[] indexDiag;
	protected double tolerance;

	public double matT[];
	public double matTDirichlet[];

	double[] volume;

	double volumeOld = 0;
	double volumeNew = 0;

	double[] matTNoDirichlet;

	double[] arrb;

	public double[] arrbNoDrichelet;
	public RCConjugateGradient cg = new RCConjugateGradient(Np);

	@Execute
	public void process() throws Exception {

		aquiferThickness = new double[Np];
		boundaryConditions = "NoDirichlet";

		// search a Dirichlet cell into the array of eta of Dirichlet
		for (int i = 0; i < etaDrichelet.length; i++) {

			if (etaDrichelet[i] != NOVALUE) {

				boundaryConditions = "Dirichlet";
				break;// go out the loop at the first Dirichlet cell
			}

		}

		// the simulation type is shown by the video output

		// choose the type of simulation at run time
		if (boundaryConditions.equals("Dirichlet")) {

			computeBEqDrichelet(boundaryConditions);

		} else {

			computeBEqNoDrichelet(boundaryConditions);

		}

	}

	/**
	 * @brief Initial method, to compute the Water Table without Dirichlet's BC
	 * It calls all the methods useful to compute the Water Table in a domain
	 * without Dirichlet's BC
	 * @param
	 * @return
	 */
	public void computeBEqNoDrichelet(String boundaryCondition)
			throws IOException, IterativeSolverDoubleNotConvergedException {
		// allocate the memory for eta array

		firstThings();

		computeInitialVolume();

		// initialize eta array
		System.arraycopy(eta, 0, eta, 0, eta.length);// la eta che copiavo era
														// in computational
														// domain quindi diverso
														// dalla eta in cui sto
														// copiando

		outflow = new double[Np];
		computeBEqArraysNoDrichelet(eta);

		eta = solutionMethod(eta, matT, arrb);

		computeOutputFeatures(eta);

		computeVolumeConservation();

		// countLoop++;

	}

	/**
	 * @brief A method to compute the two elements of the linear system
	 * @param[in] eta The piezometric head at the previous timestep
	 */
	public void computeBEqArraysNoDrichelet(double[] eta) {

		matT = computeT(eta);
		arrb = computeBNoDrichelet(eta);

	}

	/**
	 * @brief Eq. (25) of [ Cordano & Rigon, 2012 ]
	 *
	 * @param[in] eta The piezometric head at the previous timestep
	 * @return arrB The array of the known terms
	 */
	public double[] computeBNoDrichelet(double[] eta) {

		double[] arrB = new double[Np];

		for (int i = 0; i < Np; i++) {

			// water volume stored in the cell at the previous timestep
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], bedrockElevation[i], porosity[i], planArea[i]);

			// TIMESTEP has to be less then 1/c
			arrB[i] = volume + TIMESTEP * planArea[i] * source[i] - TIMESTEP
					* planArea[i] * c[i] * Math.pow(volume / planArea[i], d[i]);
			outflow[i] = TIMESTEP * planArea[i] * c[i]
					* Math.pow(volume / planArea[i], d[i]);

			if (arrB[i] < 0) {

				System.out.println("WARNING!!!\nThe element " + i
						+ " of the array of known terms is NEGATIVE");

			}

		}

		return arrB;
	}

	public double[] computeT(double[] eta) {

		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */
		double rowSum = 0;

		/* to identify the diagonal entry of matrix T in row-compressed form */
		int index = 0;

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < Np; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
			for (int j = Mp[i]; j < Mp[i + 1]; j++) {

				if (Mi[j] != i) {
					// equation (21)
					arrayT[j] = -TIMESTEP
							* (1 / euclideanDistance[(int) Ml[j] - 1])
							* hydrConductivity[(int) Ml[j] - 1]
							* lengthSides[(int) Ml[j] - 1]
							* Math.max(
									Math.max(0, eta[Mi[j]]
											- bedrockElevation[Mi[j]]),
									Math.max(0, eta[i] - bedrockElevation[i]));

					rowSum += -arrayT[j];

				} else {
					index = j;
				}

			}
			// equation (20)
			if (rowSum == 0) {

				arrayT[index] = 1;

			} else {

				arrayT[index] = rowSum;
				rowSum = 0;

			}
		}

		return arrayT;
	}

	public double[] newtonIteration(double[] arrb, double[] arrT,
			int[] indexDiag, double[] eta, RCConjugateGradient cg,
			double tolerance) throws IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		SparseDoubleMatrix1D matrixr;

		double maxResidual = 10;

		do {

			// compute Jr
			double[] jr = computeJr(indexDiag, arrT, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixJr = new SparseRCDoubleMatrix2D(Np, Np, Mp, Mi, jr);

			// compute the residual function
			double[] r = computeR(arrT, arrb, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixr = new SparseDoubleMatrix1D(r);

			cg.solverCG(matrixr, matrixJr);

			// compute the new eta for every cell
			for (int i = 0; i < Np; i++) {
				eta[i] = eta[i] - cg.matSol.get(i);
			}

			// compute the max residual
			maxResidual = Math.max(Math.abs(cg.matSol.getMaxLocation()[0]),
					Math.abs(cg.matSol.getMinLocation()[0]));

			// System.out.println(maxResidual);

		} while (maxResidual > tolerance * 1000);

		return eta;

	}

	public double[] solutionMethod(double[] etaOld, double[] matT, double[] arrb)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newtonIteration(arrb, matT, indexDiag, etaOld, cg, tolerance);

		return eta;

	}

	public double[] solutionMethodDrichelet(double[] etaOld, double[] matT,
			double[] arrb) throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newtonIterationDrichelet(arrb, matT, indexDiag, etaOld, cg,
				tolerance);

		return eta;

	}

	public double[] newtonIterationDrichelet(double[] arrb, double[] arrT,
			int[] indexDiag, double[] eta, RCConjugateGradient cg,
			double tolerance) throws IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		SparseDoubleMatrix1D matrixr;

		double maxResidual = 10;

		do {

			// compute Jr
			double[] jr = computeJrDrichelet(indexDiag, arrT, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixJr = new SparseRCDoubleMatrix2D(Np, Np, Mp, Mi, jr);

			// compute the residual function
			double[] r = computeRDrichelet(arrT, arrb, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixr = new SparseDoubleMatrix1D(r);

			cg.solverCG(matrixr, matrixJr);

			// compute the new eta for every cell
			for (int i = 0; i < eta.length; i++) {
				eta[i] = eta[i] - cg.matSol.get(i);
			}

			// compute the max residual
			maxResidual = Math.max(Math.abs(cg.matSol.getMaxLocation()[0]),
					Math.abs(cg.matSol.getMinLocation()[0]));

			// System.out.println("Residual: " + maxResidual);

		} while (maxResidual > tolerance * 100);

		// System.out.println("Compute time: " +
		// (ComputeBEqDirichlet.timeCompute)/ 1000000000.0 + " [s]");
		// System.out.println("Solver time: " +
		// (ComputeBEqDirichlet.timeSolver)/ 1000000000.0 + " [s]");

		return eta;

	}

	public double[] computeJrDrichelet(int[] indexDiag, double[] arrT,
			double[] eta) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		double[] arrJr = new double[arrT.length];

		System.arraycopy(arrT, 0, arrJr, 0, arrT.length);

		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < indexDiag.length; i++) {

			if (isNoValue(etaDrichelet[i], NOVALUE)) {
				// non Dirichlet cells
				// equation (A6)
				arrJr[indexDiag[i]] = arrT[indexDiag[i]]
						+ PolygonGeometricalWetProperties.computeWetArea(
								eta[i], bedrockElevation[i], porosity[i],
								planArea[i]);

			} else {
				// Dirichlet cells
				arrJr[indexDiag[i]] = arrT[indexDiag[i]];
			}
		}

		return arrJr;
	}

	public double[] computeRDrichelet(double[] arrT, double[] arrb, double[] eta) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		double[] arrR = new double[Np];

		for (int i = 0; i < Np; i++) {
			if (isNoValue(etaDrichelet[i], NOVALUE)) {

				// non Dirichlet cells
				for (int j = Mp[i]; j < Mp[i + 1]; j++) {
					sum += arrT[j] * eta[Mi[j]];
				}

				double waterVolume = PolygonGeometricalWetProperties
						.computeWaterVolume(eta[i], bedrockElevation[i],
								porosity[i], planArea[i]);
				// equation (A3)
				arrR[i] = waterVolume + sum - arrb[i];

				sum = 0;
			} else {

				// Dirichlet cells
				arrR[i] = 0;
			}
		}

		return arrR;
	}

	public void firstThings() {
		RCIndexDiagonalElement rcIndexDiagonalElement = new RCIndexDiagonalElement();
		indexDiag = rcIndexDiagonalElement.computeIndexDiag(Np, Mp, Mi);
		MachineEpsilon cMEd = new MachineEpsilon();
		tolerance = cMEd.computeMachineEpsilonDouble();

	}

	public void computeBEqDrichelet(String boundaryCondition)
			throws IOException, IterativeSolverDoubleNotConvergedException {

		firstThings();

		computeInitialVolume();

		System.arraycopy(eta, 0, eta, 0, eta.length);
		outflow = new double[Np];

		eta = etaInitialization(eta);

		computeBEqArraysDrichelet(eta);

		eta = solutionMethodDrichelet(eta, matTNoDirichlet, arrb);

		computeOutputFeatures(eta);

		// }

		computeVolumeConservation();

	}

	public void computeBEqArraysDrichelet(double[] eta) {

		matT = computeT(eta);
		matTDirichlet = computeTDirichlet(matT);
		matTNoDirichlet = computeTNoDirichlet(matT, indexDiag);
		arrb = computeBDrichelet(eta, matTDirichlet);

	}

	public double[] etaInitialization(double[] eta) {

		// initialize eta array
		for (int i = 0; i < eta.length; i++) {
			if (isNoValue(etaDrichelet[i], NOVALUE)) {

				// not Dirichlet cells
				eta[i] = eta[i];
			} else {

				// Dirichlet cells
				eta[i] = etaDrichelet[i];
			}
		}

		return eta;

	}

	public double[] computeTDirichlet(double[] T) {

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < Np; i++) {

			if (!isNoValue(etaDrichelet[i], NOVALUE)) {

				// Dirichlet cells
				for (int j = Mp[i]; j < Mp[i + 1]; j++) {
					arrayT[j] = T[j];
				}
			} else {

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
				for (int j = Mp[i]; j < Mp[i + 1]; j++) {

					if (!isNoValue(etaDrichelet[Mi[j]], NOVALUE)) {

						// adjacent Dirichlet cell
						arrayT[j] = T[j];
					} else {

						// adjacent non Dirichlet cell
						arrayT[j] = 0.0;
					}
				}

			}

		}

		return arrayT;
	}

	public boolean isNoValue(double x, double noValue) {

		if (x <= noValue) {
			return true;
		} else {
			return false;
		}
	}

	public double[] computeTNoDirichlet(double[] T, int[] indexDiag) {

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < Np; i++) {

			if (isNoValue(etaDrichelet[i], NOVALUE)) {

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
				for (int j = Mp[i]; j < Mp[i + 1]; j++) {

					if (isNoValue(etaDrichelet[Mi[j]], NOVALUE)) {

						// adjacent non Dirichlet cells
						arrayT[j] = T[j];
					} else {

						// adjacent Dirichlet cells
						arrayT[j] = 0.0;
					}

					if (j == indexDiag[i] && arrayT[j] == 0.0) {

						arrayT[j] = 1;

					}

				}

			} else {

				for (int j = Mp[i]; j < Mp[i + 1]; j++) {

					if (j == indexDiag[i] && arrayT[j] == 0.0) {

						arrayT[j] = 1;

					}

				}

			}

		}

		return arrayT;
	}

	public double[] convertDoubles(LinkedList<Double> doubles) {
		double[] ret = new double[doubles.size()];
		Iterator<Double> iterator = doubles.iterator();
		int i = 0;
		while (iterator.hasNext()) {
			ret[i] = iterator.next().doubleValue();
			i++;
		}
		return ret;
	}

	// declaration of the array that holds the known terms of the linear
	// system
	public double[] computeB(double[] eta) {
		double[] arrB = new double[Np];

		for (int i = 0; i < Np; i++) {
			// compute the water volume stored in the cell
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], bedrockElevation[i], porosity[i], planArea[i]);

			// delta t deve essere minore di 1/c
			arrB[i] = volume + TIMESTEP * planArea[i] * source[i] - TIMESTEP
					* planArea[i] * c[i] * Math.pow(volume / planArea[i], d[i]);

			outflow[i] = TIMESTEP * planArea[i] * c[i]
					* Math.pow(volume / planArea[i], d[i]);

			// TextIO.putln(ComputationalDomain.outflow[i]);

			if (arrB[i] < 0) {

				System.out.println("WARNING!!!\nThe element " + i
						+ " of the array of known terms is NEGATIVE");

			}

		}

		return arrB;
	}

	public double[] computeBDrichelet(double[] eta, double[] Tdirichlet) {

		// declaration of the array that holds the known terms of the linear
		// system
		double[] arrB = new double[Np];

		for (int i = 0; i < Np; i++) {
			// compute the water volume stored in the cell
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], bedrockElevation[i], porosity[i], planArea[i]);
			// equation (19)
			double sum = 0;
			for (int j = Mp[i]; j < Mp[i + 1]; j++) {
				sum += Tdirichlet[j] * etaDrichelet[Mi[j]];
			}

			// delta t deve essere minore di 1/c
			arrB[i] = volume + TIMESTEP * planArea[i] * source[i] - sum
					- TIMESTEP * planArea[i] * c[i]
					* Math.pow(volume / planArea[i], d[i]);
			//
			// if (arrB[i] < 0){
			//
			// TextIO.putln("WARNING!!!\nThe element " + i +
			// " of the array of known terms is NEGATIVE");
			//
			// }

		}

		return arrB;
	}

	public double computeVolume(int index, double eta) {

		double volume;

		volume = PolygonGeometricalWetProperties.computeWaterVolume(eta,
				bedrockElevation[index], porosity[index], planArea[index]);

		volume = volume - TIMESTEP * planArea[index] * source[index] + TIMESTEP
				* planArea[index] * c[index]
				* Math.pow(volume / planArea[index], d[index]);

		return volume;

	}

	public void computeOutputFeatures(double[] eta) {

		for (int j = 0; j < eta.length; j++) {

			volumeSource[j] = source[j] * planArea[j];
			aquiferThickness[j] = Math.max(0, eta[j] - bedrockElevation[j]);
			volume[j] = computeVolume(j, eta[j]);
			volumeNew = volumeNew + volume[j];

		}

	}

	public void computeInitialVolume() {

		for (int i = 0; i < Np; i++) {

			volume[i] = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], bedrockElevation[i], porosity[i], planArea[i]);
			volumeOld = volumeOld + volume[i];

		}

		System.out.println("Initial volume: " + volumeOld);

	}

	public void computeVolumeConservation() {

		if (Math.abs(volumeNew - volumeOld) > Math.pow(10, -6)) {

			System.out.println("WARNING!!! The system is losing mass");
			System.out
					.println("The difference between initial volume and compute volume is: "
							+ Math.abs(volumeNew - volumeOld));

		}

		volumeNew = 0;

	}

	/**
	 * Compute Jr.
	 * 
	 * @desc this method computes the Jacobian matrix of the water volume stored
	 *       into every cell. In this case the array Jr, in Row Compressed Form,
	 *       is evaluated like sum between array T and the wet area, according
	 *       the equation (A6) and (A7) of [Cordano & Rigon, 2012]. The array Jr
	 *       is a copy of T where only diagonal entries are summed to P, because
	 *       P is a diagonal matrix in Row Compressed Form too. These operations
	 *       are made only in case the Jacobian is computed only in a non
	 *       Dirichlet cell. Otherwise the volume of water stored is constant
	 *       with eta and P is equal to zero.
	 * 
	 * @param indexDiag
	 *            the array of the indices of the diagonal entries
	 * @param arrT
	 *            the array of T in Row Compressed Form
	 * @param eta
	 *            the piezometric head
	 * @param zetaBedrock
	 *            the zeta bedrock
	 * @param porosity
	 *            the porosity
	 * @param planimetricArea
	 *            the planimetric area of a cell
	 * @param etaDirichlet
	 *            the eta of Dirichlet cells
	 * @param NOVALUE
	 *            the novalue
	 * 
	 * @return the Jacobian array of water volume stored in Row Compressed Form
	 */
	public double[] computeJr(int[] indexDiag, double[] arrT, double[] eta) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		double[] arrJr = new double[arrT.length];

		System.arraycopy(arrT, 0, arrJr, 0, arrT.length);

		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < indexDiag.length; i++) {

			// equation (A6)
			arrJr[indexDiag[i]] = arrT[indexDiag[i]]
					+ PolygonGeometricalWetProperties.computeWetArea(eta[i],
							bedrockElevation[i], porosity[i], planArea[i]);

		}

		return arrJr;
	}

	public double[] computeR(double[] arrT, double[] arrb, double[] eta) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		double[] arrR = new double[Np];

		for (int i = 0; i < Np; i++) {

			for (int j = Mp[i]; j < Mp[i + 1]; j++) {
				sum += arrT[j] * eta[Mi[j]];
			}

			double waterVolume = PolygonGeometricalWetProperties
					.computeWaterVolume(eta[i], bedrockElevation[i],
							porosity[i], planArea[i]);
			// equation (A3)
			arrR[i] = waterVolume + sum - arrb[i];

			sum = 0;

		}

		return arrR;
	}

}
