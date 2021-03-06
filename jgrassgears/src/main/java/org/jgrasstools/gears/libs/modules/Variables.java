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
package org.jgrasstools.gears.libs.modules;


/**
 * Variable names that also need translations and used in the modules.
 * 
 * @author Andrea Antonello (www.hydrologis.com)
 */
public interface Variables {
    public String DEFAULT = "default";

    public String PROGRESS_MONITOR_EN = "The progress monitor.";

    public String TPS = "TPS";
    public String IDW = "IDW";

    public String NEAREST_NEIGHTBOUR = "nearest neightbour";
    public String BILINEAR = "bilinear";
    public String BICUBIC = "bicubic";

    public String INTERSECTION = "intersection";
    public String UNION = "union";
    public String DIFFERENCE = "difference";
    public String SYMDIFFERENCE = "symdifference";

    public String TCA = "tca";
    public String TCA_SLOPE = "tca and slope";
    public String TCA_CONVERGENT = "tca in convergent sites";

    public String FIXED_NETWORK = "with fixed network";

    public String CAP_ROUND = "round";
    public String CAP_FLAT = "flat";
    public String CAP_SQUARE = "square";
    public String JOIN_ROUND = "round";
    public String JOIN_MITRE = "mitre";
    public String JOIN_BEVEL = "bevel";

    public String DILATE = "dilate";
    public String ERODE = "erode";
    public String SKELETONIZE = "skeletonize";
    public String OPEN = "open";
    public String CLOSE = "close";

    public String CUSTOM = "custom";
    public String DECIDUOUS = "deciduous";
    public String CONIFER = "conifer";
    public String MIXED_PINES_AND_DECIDUOUS = "mixed pines and deciduous trees";

    public String BINARY = "binary";
    public String COSINE = "cosine";
    public String DISTANCE = "distance";
    public String EPANECHNIKOV = "epanechnikov";
    public String GAUSSIAN = "gaussian";
    public String INVERSE_DISTANCE = "inverse_distance";
    public String QUARTIC = "quartic";
    public String TRIANGULAR = "triangular";
    public String TRIWEIGHT = "triweight";

}
