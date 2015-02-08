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

public class MachineEpsilon {

	/**
	 * Calculate machine epsilon double.
	 * 
	 * @desc this method compute the tolerance of the machine. For more info go
	 *       to
	 *       https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java
	 *       . In c/c++ section there's write that:
	 * 
	 *       In such languages as C or C++ when you do something like while( 1.0
	 *       + eps > 1.0 ) the expression is calculated not with 64 bits
	 *       (double) but with the processor precision (80 bits or more depends
	 *       on the processor and compile options). Below program calculates
	 *       exactly on 32 bits (float) and 64 bits (double)
	 * 
	 * 
	 * @return the tolerance of the machine for double
	 */
	public double computeMachineEpsilonDouble() {

		// machine tolerance
		double machEps = 1.0d;

		do
			machEps /= 2.0d;
		while ((double) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

	/**
	 * Calculate machine epsilon float.
	 * 
	 * @desc this method compute the tolerance of the machine. For more info go
	 *       to
	 *       https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java
	 *       . In c/c++ section there's write that:
	 * 
	 *       In such languages as C or C++ when you do something like while( 1.0
	 *       + eps > 1.0 ) the expression is calculated not with 64 bits
	 *       (double) but with the processor precision (80 bits or more depends
	 *       on the processor and compile options). Below program calculates
	 *       exactly on 32 bits (float) and 64 bits (double)
	 * 
	 * 
	 * @return the tolerance of the machine for float
	 */
	public float computeMachineEpsilonFloat() {
		float machEps = 1.0f;

		do
			machEps /= 2.0f;
		while ((float) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

	public static void main(String[] args) {

		double doubleTolerance;
		float floatTolerance;

		MachineEpsilon cMEd = new MachineEpsilon();
		MachineEpsilon cMEf = new MachineEpsilon();
		
		doubleTolerance = cMEd.computeMachineEpsilonDouble();
		floatTolerance = cMEf.computeMachineEpsilonFloat();

		System.out.println("The machine precision for double is : "
				+ doubleTolerance);
		System.out.println("The machine precision for float is : "
				+ floatTolerance);

	}

}
