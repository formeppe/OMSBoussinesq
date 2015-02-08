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

public class PolygonGeometricalWetProperties {

	/**
	 * Compute wet area.
	 * 
	 * @desc this method computes the wet area of the cell; if the piezometric
	 *       head is less than the bedrock elevation, the wet area is null,
	 *       otherwise is equal to the porosity for planimetric area of the
	 *       cell.
	 * 
	 * @param eta
	 *            the piezometric head
	 * @param zetaBedrock
	 *            the bedrock elevation
	 * @param porosity
	 *            the porosity
	 * @param planimetricArea
	 *            the planimetric area of the cell
	 * 
	 * @return the wet area of the cell
	 */
	public static double computeWetArea(double eta, double zetaBedrock,
			double porosity, double planimetricArea) {

		// the wet area is the variable that the method returns
		double wetArea = 0;

		if (eta > zetaBedrock) {
			wetArea = porosity * planimetricArea;
		} else {
			wetArea = 0;
		}

		return wetArea;
	}

	/**
	 * Compute water volume.
	 * 
	 * @desc this method computes the volume of stored water in the cell like
	 *       multiplication between the wet area and the thickness of the water
	 *       table
	 * 
	 * @param eta
	 *            the piezometric head
	 * @param zetaBedrock
	 *            the bedrock elevation
	 * @param planimetricArea
	 *            the planimetric area of the cell
	 * 
	 * @return the volume of the stored water in the cell
	 */
	public static double computeWaterVolume(double eta, double zetaBedrock,
			double porosity, double planimetricArea) {

		double volume = computeWetArea(eta, zetaBedrock, porosity, planimetricArea) * (eta - zetaBedrock);

		return volume;

	}

}
