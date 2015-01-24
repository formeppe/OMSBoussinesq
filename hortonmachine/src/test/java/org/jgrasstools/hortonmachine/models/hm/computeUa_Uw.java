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
package org.jgrasstools.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.geomorphology.ab.OmsAb;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test ab.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class computeUa_Uw extends HMTestCase {
	public final static double GAMMAW_KNM3 = 9810;

	public static void main(String[] Args) {

		double nVG = 5.1270;
		double alphaVG = 0.14;
		double ks = 1.0e-6;
		double ua_uw[] = new double[500];
		double tt[] = new double[500];
		double zz[] = new double[500];
		double q = -1.5e-7;
		double thetaR=0.2;
		double thetaS=0.6;
		for (int i = 0; i < 500; i++) {
			double z = i/100.0 ;
			zz[i]=z;
			ua_uw[i] = compUa_Uw(alphaVG, q, ks, GAMMAW_KNM3, z);
			tt[i]=computeTheta( ua_uw[i],  alphaVG,  nVG,thetaR,  thetaS);
			System.out.println(z+" ,"+ua_uw[i]);
		}
		System.out.println("ciao");

	}

	private static double compUa_Uw(double alpha, double q, double ks,
			double GAMMAW_KNM3, double z) {
		double ua_uw = 0;
		
		
		ua_uw = -(1 / alpha) * Math.log((1.0 + q / ks)
				* Math.exp(-(GAMMAW_KNM3/1000) * alpha * z)-q / ks);
		return ua_uw;
	}

	private static double computeTheta(double ua_uw, double alphaVG, double nVG,
			double thetaR, double thetaS) {
		double theta = 0;
		double ttt = 0;
		double a = Math.pow(alphaVG * ua_uw, nVG);
		ttt = Math.pow(1 / (1 + a), 1 - 1 / nVG);
		theta = thetaR + (thetaS - thetaR) * ttt;

		return ttt;
	}

}
