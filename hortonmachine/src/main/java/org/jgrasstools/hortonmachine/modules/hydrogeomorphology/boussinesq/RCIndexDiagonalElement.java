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

public class RCIndexDiagonalElement {

	/**
	 * Compute index of diagonal.
	 * 
	 * @desc this method computes the indices of the diagonal of the adjacency
	 *       matrix; this matrix and all the sparse matrices are stored in Row
	 *       Compressed Form. More information at the web site
	 *       https://en.wikipedia.org/wiki/Sparse_matrix
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * 
	 * @return the array that holds the indices of the diagonal entries of the
	 *         sparse adjacency matrix in Row Compressed Form
	 */
	public int[] computeIndexDiag(int size, int[] Mp, int[] Mi) {

		// declaration of the array that holds the indices of diagonal entries
		int[] indexDiag = new int[size];

		/* for-loop to analyze the matrix cell by cell */
		for (int i = 0; i < size; i++) {
			/*
			 * nested for-loop to analyze diagonal entries, which are identified
			 * by a negative number
			 */
			for (int j = Mp[i]; j < Mp[i + 1]; j++) {

				if (Mi[j] == i) {
					indexDiag[i] = j;
				}

			}
		}

		return indexDiag;
	}
	
}
