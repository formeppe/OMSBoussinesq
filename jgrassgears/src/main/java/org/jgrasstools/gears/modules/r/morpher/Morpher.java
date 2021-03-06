/*
 * This file is part of JGrasstools (http://www.jgrasstools.org)
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
package org.jgrasstools.gears.modules.r.morpher;

import static org.jgrasstools.gears.libs.modules.JGTConstants.doubleNovalue;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;
import static org.jgrasstools.gears.libs.modules.Variables.CLOSE;
import static org.jgrasstools.gears.libs.modules.Variables.DILATE;
import static org.jgrasstools.gears.libs.modules.Variables.ERODE;
import static org.jgrasstools.gears.libs.modules.Variables.OPEN;
import static org.jgrasstools.gears.libs.modules.Variables.SKELETONIZE;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;

import javax.media.jai.iterator.RandomIter;
import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

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
import oms3.annotations.UI;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.libs.exceptions.ModelsIllegalargumentException;
import org.jgrasstools.gears.libs.modules.GridNode;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.monitor.DummyProgressMonitor;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.modules.utils.BinaryFast;
import org.jgrasstools.gears.utils.RegionMap;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;

@Description("Mophologic binary operations")
@Author(name = "Simon Horne, Andrea Antonello", contact = "http://homepages.inf.ed.ac.uk/rbf/HIPR2/, www.hydrologis.com")
@Keywords("Dilation, Erosion, Skeletonize, Open, Close, Raster")
@Label(JGTConstants.RASTERPROCESSING)
@Name("morpher")
@Status(Status.DRAFT)
@License("http://www.gnu.org/licenses/gpl-3.0.html")
public class Morpher extends JGTModel {

    @Description("The map to morph.")
    @In
    public GridCoverage2D inMap = null;

    @Description("The progress monitor.")
    @In
    public IJGTProgressMonitor pm = new DummyProgressMonitor();

    @Description("A kernel to use instead of the default.")
    @In
    public int[] pKernel = null;

    @Description("Process in binary mode (novalue = 0, valid value = 1.")
    @In
    public boolean doBinary = true;

    @Description("The operation type to perform (dilate, erode, skeletonize, open, close)")
    @UI("combo:" + DILATE + "," + ERODE + "," + SKELETONIZE + "," + OPEN + "," + CLOSE)
    @In
    public String pMode = DILATE;

    @Description("The resulting map.")
    @Out
    public GridCoverage2D outMap = null;

    @Execute
    public void process() throws Exception {
        if (!concatOr(outMap == null, doReset)) {
            return;
        }

        if (pKernel == null) {
            pKernel = MorpherHelp.DEFAULT3X3KERNEL;
        }

        if (pMode.equals(SKELETONIZE)) {
            skeletonize();
        } else {
            RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(inMap);
            WritableRaster inWR = CoverageUtilities.renderedImage2WritableRaster(inMap.getRenderedImage(), false);
            WritableRaster[] wrHolder = new WritableRaster[1];
            outMap = CoverageUtilities.createCoverageFromTemplate(inMap, doubleNovalue, wrHolder);
            WritableRaster outWR = wrHolder[0];
            if (pMode.equals(DILATE)) {
                dilate(inWR, regionMap, outWR, pKernel, doBinary);
            } else if (pMode.equals(ERODE)) {
                erode(inWR, regionMap, outWR, pKernel);
            } else if (pMode.equals(OPEN)) {
                open(inWR, regionMap, outWR, pKernel, doBinary);
            } else if (pMode.equals(CLOSE)) {
                close(inWR, regionMap, outWR, pKernel, doBinary);
            } else {
                throw new ModelsIllegalargumentException("Could not recognize mode.", this);
            }
        }

    }

    /**
     * Morphologically dilates an input raster by a given kernel. 
     * 
     * @param inWR the input raster.
     * @param regionMap the {@link RegionMap}.
     * @param outWR the raster to modify.
     * @param kernelArray the kernel to use.
     * @param binary if <code>true</code>, binary mode is used.
     */
    public static void dilate( WritableRaster inWR, RegionMap regionMap, WritableRaster outWR, int[] kernelArray, boolean binary ) {
        int cols = regionMap.getCols();
        int rows = regionMap.getRows();
        double xres = regionMap.getXres();
        double yres = regionMap.getYres();

        WritableRandomIter inIter = RandomIterFactory.createWritable(inWR, null);
        WritableRandomIter outIter = RandomIterFactory.createWritable(outWR, null);

        int[][] kernel = MorpherHelp.getSquareKernelMatrix(kernelArray);

        for( int c = 0; c < cols; c++ ) {
            for( int r = 0; r < rows; r++ ) {
                GridNode node = new GridNode(inIter, cols, rows, xres, yres, c, r);
                if (!node.isValid()) {
                    double[][] nodeNeighbours = node.getWindow(kernel.length, false);
                    // apply the kernel
                    boolean set = false;
                    boolean doBreak = false;
                    double max = Double.NEGATIVE_INFINITY;
                    for( int kr = 0; kr < kernel.length; kr++ ) {
                        for( int kc = 0; kc < kernel[0].length; kc++ ) {
                            if (kernel[kr][kc] == 1) {
                                // valid kernel value
                                if (!isNovalue(nodeNeighbours[kr][kc]) && nodeNeighbours[kr][kc] > max) {
                                    max = nodeNeighbours[kr][kc];
                                    set = true;
                                    if (binary) {
                                        break;
                                    }
                                }
                            }
                        }
                        if (doBreak) {
                            break;
                        }
                    }
                    if (set) {
                        node.setValueInMap(outIter, max);
                    } else {
                        node.setValueInMap(outIter, JGTConstants.doubleNovalue);
                    }
                    continue;
                } else {
                    node.setValueInMap(outIter, node.elevation);
                }
            }
        }
    }

    /**
     * Morphologically erodes an input raster by a given kernel. 
     * 
     * @param inWR the input raster.
     * @param regionMap the {@link RegionMap}.
     * @param outWR the raster to modify.
     * @param kernelArray the kernel to use.
     */
    public static void erode( WritableRaster inWR, RegionMap regionMap, WritableRaster outWR, int[] kernelArray ) {
        int cols = regionMap.getCols();
        int rows = regionMap.getRows();
        double xres = regionMap.getXres();
        double yres = regionMap.getYres();

        WritableRandomIter inIter = RandomIterFactory.createWritable(inWR, null);
        WritableRandomIter outIter = RandomIterFactory.createWritable(outWR, null);

        int[][] kernel = MorpherHelp.getSquareKernelMatrix(kernelArray);

        for( int c = 0; c < cols; c++ ) {
            for( int r = 0; r < rows; r++ ) {
                GridNode node = new GridNode(inIter, cols, rows, xres, yres, c, r);
                if (node.isValid()) {
                    double[][] nodeNeighbours = node.getWindow(kernel.length, false);
                    // apply the kernel
                    boolean set = false;
                    double min = Double.POSITIVE_INFINITY;
                    for( int kr = 0; kr < kernel.length; kr++ ) {
                        for( int kc = 0; kc < kernel[0].length; kc++ ) {
                            if (kernel[kr][kc] == 1) {
                                // valid kernel value
                                if (isNovalue(nodeNeighbours[kr][kc])) {
                                    // novalue is the absolute min
                                    min = doubleNovalue;
                                    set = true;
                                    break;
                                } else if (nodeNeighbours[kr][kc] < min) {
                                    min = nodeNeighbours[kr][kc];
                                    set = true;
                                }
                            }
                        }
                        if (isNovalue(min)) {
                            break;
                        }
                    }
                    if (set) {
                        node.setValueInMap(outIter, min);
                    } else {
                        node.setValueInMap(outIter, node.elevation);
                    }
                    continue;
                }
            }
        }
    }

    /**
     * Morphologically opens an input raster by a given kernel. 
     * 
     * @param inWR the input raster.
     * @param regionMap the {@link RegionMap}.
     * @param outWR the raster to modify.
     * @param kernelArray the kernel to use.
     * @param binary if <code>true</code>, binary mode is used.
     */
    public static void open( WritableRaster inWR, RegionMap regionMap, WritableRaster outWR, int[] kernelArray, boolean binary ) {
        erode(inWR, regionMap, outWR, kernelArray);
        inWR.setDataElements(0, 0, outWR);
        clearRaster(regionMap, outWR);
        dilate(inWR, regionMap, outWR, kernelArray, binary);
    }

    /**
     * Morphologically closes an input raster by a given kernel. 
     * 
     * @param inWR the input raster.
     * @param regionMap the {@link RegionMap}.
     * @param outWR the raster to modify.
     * @param kernelArray the kernel to use.
     * @param binary if <code>true</code>, binary mode is used.
     */
    public static void close( WritableRaster inWR, RegionMap regionMap, WritableRaster outWR, int[] kernelArray, boolean binary ) {
        dilate(inWR, regionMap, outWR, kernelArray, binary);
        inWR.setDataElements(0, 0, outWR);
        clearRaster(regionMap, outWR);
        erode(inWR, regionMap, outWR, kernelArray);
    }

    private void skeletonize() {
        final RenderedImage renderedImage = inMap.getRenderedImage();
        RandomIter iter = RandomIterFactory.create(renderedImage, null);
        int width = renderedImage.getWidth();
        int height = renderedImage.getHeight();
        int[][] data = new int[width][height];
        for( int c = 0; c < width; c++ ) {
            for( int r = 0; r < height; r++ ) {
                double value = iter.getSampleDouble(c, r, 0);
                data[c][r] = BinaryFast.BACKGROUND;
                if (isNovalue(value)) {
                    continue;
                } else {
                    data[c][r] = BinaryFast.FOREGROUND;
                }
            }
        }

        BinaryFast binaryData = new BinaryFast(data);
        if (pKernel != null) {
            pm.errorMessage("Thinning doesn't permit for custom kernels. Using standard kernel.");
        }
        new Thin().process(binaryData);
        int[] values = binaryData.getValues();
        WritableRaster dataWR = CoverageUtilities.createWritableRasterFromArray(width, height, values);
        RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(inMap);
        outMap = CoverageUtilities.buildCoverage("morphed", dataWR, regionMap, inMap.getCoordinateReferenceSystem()); //$NON-NLS-1$
    }

    private static void clearRaster( RegionMap regionMap, WritableRaster outWR ) {
        // clear raster
        for( int c = 0; c < regionMap.getCols(); c++ ) {
            for( int r = 0; r < regionMap.getRows(); r++ ) {
                outWR.setSample(c, r, 0, doubleNovalue);
            }
        }
    }

}
