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

import cern.colt.Arrays;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleDiagonal;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleICC;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleILU;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoublePreconditioner;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

/**
 * The Class RCConjugateGradient.
 */
public class RCConjugateGradient {

	/** The matrix_ a. */
	SparseRCDoubleMatrix2D matrix_A;

	/** The matrix_x. */
	DenseDoubleMatrix1D matrix_x;

	/** The mat sol. */
	public DoubleMatrix1D matSol;

	DoubleMatrix1D prova;
	
	DoublePreconditioner dd;
	
	/**
	 * Instantiates a new rC conjugate gradient.
	 * 
	 * @param SIZE
	 *            the size
	 * @param Mp
	 *            the mp
	 * @param Mi
	 *            the mi
	 * @param Ml
	 *            the ml
	 */
	public RCConjugateGradient(int SIZE) {

		matrix_x = new DenseDoubleMatrix1D(SIZE);
		matSol = new DenseDoubleMatrix1D(SIZE);
		dd = new DoubleICC(SIZE);

//		matrix_x.set(0, 0);
//		matrix_x.set(0, 1);
		
	}

	/**
	 * Solver cg.
	 * 
	 * @param matrix_b
	 *            the matrix_b
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public void solverCG(DoubleMatrix1D matrix_b,
			SparseRCDoubleMatrix2D matrix_A)
			throws IterativeSolverDoubleNotConvergedException {

		dd.setMatrix(matrix_A);
		
		DoubleCG conjugateGradient = new DoubleCG(matrix_x);
		
		
		
		conjugateGradient.setPreconditioner(dd);
		
//		DoubleCG conjugateGradient = new DoubleCG(matrix_b);
		matSol = conjugateGradient.solve(matrix_A, matrix_b, matrix_x);
		
	}

	/**
	 * The main method.
	 * 
	 * @param args
	 *            the arguments
	 * @throws IterativeSolverDoubleNotConvergedException 
	 */
	public static void main(String[] args) throws IterativeSolverDoubleNotConvergedException {

		double[] b = { 1, 2 };
		int[] Mp = { 0, 2, 4 };
		int[] Mi = { 0, 1, 0, 1 };
		double[] Ml = { 4, 1, 1, 3 };

		DenseDoubleMatrix1D matrix_b = new DenseDoubleMatrix1D(b);
		SparseRCDoubleMatrix2D matrix_A = new SparseRCDoubleMatrix2D(b.length,
				b.length, Mp, Mi, Ml);
		// Mesh grid1 = new Mesh("Song");
		RCConjugateGradient cg = new RCConjugateGradient(b.length);

		cg.solverCG(matrix_b, matrix_A);
		
		System.out.println(cg.matSol.get(0));

	}

}
