package org.jgrasstools.hortonmachine.models.hm;

import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.geomorphology.aspect.Aspect;
import org.jgrasstools.hortonmachine.utils.HMTestCase;
import org.jgrasstools.hortonmachine.utils.HMTestMaps;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 * Test the {@link Aspect} module.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public class TestAspect extends HMTestCase {
    public void testAspectDegrees() throws Exception {
        double[][] pitData = HMTestMaps.pitData;
        HashMap<String, Double> envelopeParams = HMTestMaps.envelopeParams;
        CoordinateReferenceSystem crs = HMTestMaps.crs;
        GridCoverage2D pitCoverage = CoverageUtilities.buildCoverage("pit", pitData, envelopeParams, crs, true);

        PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(System.out, System.out);

        Aspect aspect = new Aspect();
        aspect.inDem = pitCoverage;
        aspect.doRound = true;
        aspect.pm = pm;

        aspect.process();

        GridCoverage2D aspectCoverage = aspect.outAspect;

        checkMatrixEqual(aspectCoverage.getRenderedImage(), HMTestMaps.aspectDataDegrees, 0.01);
    }

    public void testAspectRadiants() throws Exception {
        double[][] pitData = HMTestMaps.pitData;
        HashMap<String, Double> envelopeParams = HMTestMaps.envelopeParams;
        CoordinateReferenceSystem crs = HMTestMaps.crs;
        GridCoverage2D pitCoverage = CoverageUtilities.buildCoverage("pit", pitData, envelopeParams, crs, true);
        
        PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(System.out, System.out);
        
        Aspect aspect = new Aspect();
        aspect.inDem = pitCoverage;
        aspect.doRadiants = true;
        aspect.pm = pm;
        
        aspect.process();
        
        GridCoverage2D aspectCoverage = aspect.outAspect;
        
        checkMatrixEqual(aspectCoverage.getRenderedImage(), HMTestMaps.aspectDataRadiants, 0.01);
    }

}