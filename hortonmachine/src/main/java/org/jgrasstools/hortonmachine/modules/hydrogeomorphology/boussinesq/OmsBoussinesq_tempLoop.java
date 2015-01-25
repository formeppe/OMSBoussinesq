/*
 * JGrass - Free Open Source Java GIS http://www.jgrass.org 
 * (C) HydroloGIS - www.hydrologis.com 
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

import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_AUTHORCONTACTS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_DESCRIPTION;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_KEYWORDS;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_LABEL;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_LICENSE;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_NAME;
import static org.jgrasstools.hortonmachine.i18n.HortonMessages.OMSADIGE_STATUS;

import java.io.IOException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

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
import oms3.annotations.Unit;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.joda.time.DateTime;
import org.joda.time.DateTimeZone;
import org.joda.time.format.DateTimeFormat;
import org.joda.time.format.DateTimeFormatter;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

@Description(OMSADIGE_DESCRIPTION)
@Author(name = "Serafin, Formetta, Cordano", contact = OMSADIGE_AUTHORCONTACTS)
@Keywords(OMSADIGE_KEYWORDS)
@Label(OMSADIGE_LABEL)
@Name(OMSADIGE_NAME)
@Status(OMSADIGE_STATUS)
@License(OMSADIGE_LICENSE)
public class OmsBoussinesq_tempLoop extends JGTModel {

	@Description("Start date")
	@In
	@Unit("")
	public String StartDate;

	@Description("End date")
	@In
	@Unit("")
	public String EndDate;

	@Description("TimeStep")
	@In
	@Unit("s")
	public String TimeStep;

	@Execute
	public void process() throws Exception {
		DateTimeFormatter formatter = DateTimeFormat.forPattern(
				"yyyy-MM-dd HH:mm").withZone(DateTimeZone.UTC);
		DateTime startcurrentDatetime = formatter.parseDateTime(StartDate);

		DateTime endcurrentDatetime = formatter.parseDateTime(EndDate);
		DateTime timeStep = formatter.parseDateTime(TimeStep);

		long diff = 0;
		diff = (endcurrentDatetime.getMillis() - startcurrentDatetime
				.getMillis()) / ((timeStep.getMillis()));

		for (int i = 0; i < diff; i++) {
			String pathToSlope = "/Users/giuseppeformetta/Desktop/testMesh/pitLW100.asc";

			OmsRasterReader reader = new OmsRasterReader();
			reader.file = pathToSlope;
			reader.fileNovalue = -9999.0;
			reader.geodataNovalue = Double.NaN;
			reader.process();
			GridCoverage2D slopeCoverage = reader.outRaster;
			OmsBoussinesq b = new OmsBoussinesq();

			ComputeBeqInput c = new ComputeBeqInput();
			c.inDtm = slopeCoverage;
			c.pPath = "/Users/giuseppeformetta/Desktop/BEQ";

			// in computebeqinput cambiare valueEtaInitialCondVtoAddtoDem cn una
			// mappa
			// if(i==0){
			// dagli la mappa di condizioni iniziali (t=0)
			// c.valueEtaInitialCondVtoAddtoDem = 0.1;
			// }else{
			// dagli la eta in output di b
			// c.valueEtaInitialCondVtoAddtoDem = b.e;
			// }
			c.valueCV = 10.0;
			c.valueMV = 0.3;
			c.valueEtaDricheletV = -9999;
			c.valuePorosityV = 0.4;
			c.valueSourceV = 0.1;
			c.process();

			b.Mi = c.Mj;
			b.Ml = c.Ml;
			b.Mp = c.Mp;
			b.bedrockElevation = c.vBedrock;
			b.c = c.vC;
			b.d = c.vM;
			b.eta = c.vEtaInitialCondV;
			b.etaDrichelet = c.vEtaDrichelet;
			// DA CREARE b.euclideanDistance=c.
			b.hydrConductivity = c.vHydraulicCondictivity;
			// DA CREARE b.lengthSides=c.
			b.planArea = c.vPlanarArea;
			// DA CREARE b.NOVALUE=
			b.porosity = c.vPorosity;
			b.source = c.vSource;
			b.process();

			// passare dalle eta in rcf alla eta sulla mappa, usando come
			// mappatura il file dei poligoni numerati

		}

	}

}