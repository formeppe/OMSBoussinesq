package org.jgrasstools.hortonmachine.models.hm;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.sql.rowset.Joinable;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.jgrasstools.hortonmachine.modules.statistics.cb.Cb;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

/* Test for the Cb module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class aaaaGRAZIA extends HMTestCase {

	
	
	public static boolean verifica(int M[][]){
		for(int i=0;i<M.length;i++){
			boolean primo=eprimo(M[i][i]);
			if(primo==false ||M[i][i]!=M[M.length-1-i][i]){
				return false;
			}
		}
		return true;
	}
	public static boolean eprimo(int n) {
		if (n == 2)
			return true;
		if (n % 2 == 0)
			return false;
		for (int i = 3; i * i <= n; i += 2) {
			if (n % i == 0)
				return false;
		}
		return true;
	}

	public static void main(String args[]) {
		int n = 1;
		System.out.println((eprimo(n)));
		int M[][]=new int [][]{{5,2,3,4,5,17},{1,2,3,4,13,6},{1,2,3,11,5,6},{1,2,3,11,5,6},{1,2,3,4,13,6},{5,2,3,4,13,17}};
		for(int i=0;i<M.length;i++){
			for(int j=0;j<M[0].length;j++){
				System.out.print(M[i][j]+" ");
			}
			System.out.println(" ");
		}
		System.out.println((verifica(M)));

		
	}

}
