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
public class Cb_for_me_LU_PLANAR extends HMTestCase {

	public void testCb() throws Exception {
		String[] nomimappe = new String[] { "thetaliqL0010N0003.asc",
				"thetaliqL0010N0006.asc", "thetaliqL0010N0012.asc",
				"thetaliqL0010N0024.asc", "thetaliqL0010N0047.asc",
				"thetaliqL0020N0003.asc", "thetaliqL0020N0006.asc",
				"thetaliqL0020N0012.asc", "thetaliqL0020N0024.asc",
				"thetaliqL0020N0047.asc", "thetaliqL0030N0003.asc",
				"thetaliqL0030N0006.asc", "thetaliqL0030N0012.asc",
				"thetaliqL0030N0024.asc", "thetaliqL0030N0047.asc",
				"pressureLIQL0010N0003.asc", "pressureLIQL0010N0006.asc",
				"pressureLIQL0010N0012.asc", "pressureLIQL0010N0024.asc",
				"pressureLIQL0010N0047.asc", "pressureLIQL0020N0003.asc",
				"pressureLIQL0020N0006.asc", "pressureLIQL0020N0012.asc",
				"pressureLIQL0020N0024.asc", "pressureLIQL0020N0047.asc",
				"pressureLIQL0030N0003.asc", "pressureLIQL0030N0006.asc",
				"pressureLIQL0030N0012.asc", "pressureLIQL0030N0024.asc",
				"pressureLIQL0030N0047.asc", "FS_NL10N0003.asc",
				"FS_NL10N0006.asc", "FS_NL10N0012.asc", "FS_NL10N0024.asc",
				"FS_NL10N0048.asc", "FS_NL20N0003.asc", "FS_NL20N0006.asc",
				"FS_NL20N0012.asc", "FS_NL20N0024.asc", "FS_NL20N0048.asc",
				"FS_NL30N0003.asc", "FS_NL30N0006.asc", "FS_NL30N0012.asc",
				"FS_NL30N0024.asc", "FS_NL30N0048.asc", "sigmaS_NL10N0003.asc",
				"sigmaS_NL10N0006.asc", "sigmaS_NL10N0012.asc",
				"sigmaS_NL10N0024.asc", "sigmaS_NL10N0048.asc",
				"sigmaS_NL20N0003.asc", "sigmaS_NL20N0006.asc",
				"sigmaS_NL20N0012.asc", "sigmaS_NL20N0024.asc",
				"sigmaS_NL20N0048.asc", "sigmaS_NL30N0003.asc",
				"sigmaS_NL30N0006.asc", "sigmaS_NL30N0012.asc",
				"sigmaS_NL30N0024.asc", "sigmaS_NL30N0048.asc" };

		String[] pathToDir = new String[] {
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I1/DRY/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I2/DRY/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I3/DRY/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I1/MEDIUM/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I2/MEDIUM/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I3/MEDIUM/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I1/WET/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I2/WET/planar/output_maps/",
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/I3/WET/planar/output_maps/" };

		int c = 0;
		String[] pathToDirOK = new String[pathToDir.length * nomimappe.length];
		for (int i = 0; i < pathToDir.length; i++) {
			for (int j = 0; j < nomimappe.length; j++) {

				System.out.println(c + "  " + pathToDir[i] + nomimappe[j]);
				pathToDirOK[c] = pathToDir[i] + nomimappe[j];
				c = c + 1;
			}

		}
		String path = new File(
				"/Users/giuseppeformetta/Desktop/MINES/PAPER_SIMONI/loc_4806/planar/cell/slopePlanarSites")
				.getAbsolutePath();
		OmsRasterReader reader = new OmsRasterReader();
		reader.file = path;
		reader.fileNovalue = -9999.0;
		reader.geodataNovalue = Double.NaN;
		reader.process();
		GridCoverage2D slope = reader.outRaster;

		String a = "/Users/giuseppeformetta/Desktop/XX_CB/";

		for (int ii = 0; ii < pathToDirOK.length; ii++) {

			try {
				String filename = pathToDirOK[ii];
				FileWriter fw = new FileWriter(filename, true); // the true will
				// append the
				// new data
				fw.write("\n");// appends the string to the file
				fw.close();
			} catch (IOException ioe) {
				System.err.println("IOException: " + ioe.getMessage());
			}

			String path2 = pathToDirOK[ii];
			OmsRasterReader reader2 = new OmsRasterReader();
			reader2.file = path2;
			reader2.fileNovalue = -9999.0;
			reader2.geodataNovalue = Double.NaN;
			reader2.process();

			GridCoverage2D h2cd = reader2.outRaster;

			PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(
					System.out, System.out);

			Cb cb = new Cb();
			cb.pBins = 100;
			cb.pFirst = 1;
			cb.pLast = 4;
			cb.inMap1 = slope;
			cb.inMap2 = h2cd;
			cb.pm = pm;
			cb.process();

			String[] parts = pathToDirOK[ii].split("/");
			String p = parts[parts.length - 1].split(".asc")[0];

			double[][] moments = cb.outCb;

			// System.out.println(a.concat(listOfFiles[ii].getName()));
			String ppp = a.concat(parts[6]).concat("/").concat(parts[7])
					.concat("/").concat("planar/").concat("CB_").concat(p);
			//
			FileWriter Rstatfile = new FileWriter(ppp);
			PrintWriter errestat = new PrintWriter(Rstatfile);
			for (int i = 0; i < (moments.length); i++) {
				for (int j = 0; j < (moments[0].length); j++) {

					errestat.print(moments[i][j] + " ");
					System.out.print(moments[i][j] + " ");

				}
				errestat.println();
				System.out.println();
			}

			Rstatfile.close();
		}

		System.out.print("ciao");

	}
		

		

}
