package eu.hydrologis.jgrass.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import eu.hydrologis.jgrass.hortonmachine.libs.monitor.PrintStreamProgressMonitor;
import eu.hydrologis.jgrass.hortonmachine.modules.geomorphology.draindir.DrainDir;
import eu.hydrologis.jgrass.hortonmachine.utils.HMTestCase;
import eu.hydrologis.jgrass.hortonmachine.utils.HMTestMaps;
import eu.hydrologis.jgrass.hortonmachine.utils.coverage.CoverageUtilities;

/**
 * Test the {@link DrainDir} module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestDrain extends HMTestCase {

    public void testDrain() throws Exception {
        HashMap<String, Double> envelopeParams = HMTestMaps.envelopeParams;
        CoordinateReferenceSystem crs = HMTestMaps.crs;

        double[][] pitfillerData = HMTestMaps.pitData;
        GridCoverage2D pitfillerCoverage = CoverageUtilities.buildCoverage("pitfiller",
                pitfillerData, envelopeParams, crs);
        double[][] flowData = HMTestMaps.flowData;
        GridCoverage2D flowCoverage = CoverageUtilities.buildCoverage("flow", flowData,
                envelopeParams, crs);

        PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(System.out, System.out);

        DrainDir drainDir = new DrainDir();
        // drainDir.doLad = false;
        drainDir.pLambda = 1;
        drainDir.inPit = pitfillerCoverage;
        drainDir.inFlow = flowCoverage;
        drainDir.pm = pm;

        drainDir.process();

        GridCoverage2D draindirCoverage = drainDir.outFlow;
        GridCoverage2D tcaCoverage = drainDir.outTca;

        checkMatrixEqual(draindirCoverage.getRenderedImage(), HMTestMaps.drainData1);
        checkMatrixEqual(tcaCoverage.getRenderedImage(), HMTestMaps.mtcaData);
    }

}