package eu.hydrologis.edc.annotatedclasses.timeseries.netshortwaveradiation;

import javax.persistence.Entity;
import javax.persistence.Table;
import static eu.hydrologis.edc.utils.Constants.*;

import eu.hydrologis.edc.annotatedclasses.timeseries.SeriesMonitoringPointsTable;

@Entity
@Table(name = "series_netshortwaveradiation_2013", schema = "edcseries")
@org.hibernate.annotations.Table(appliesTo = "series_netshortwaveradiation_2013", 
        indexes = @org.hibernate.annotations.Index(
                name = "IDX_TIMESTAMP_MONPOINT_series_netshortwaveradiation_2013",
                columnNames = {TIMESTAMPUTC, MONITORINGPOINTS_ID}
))
public class SeriesNetshortwaveradiation2013 extends SeriesMonitoringPointsTable {
}