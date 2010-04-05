package eu.hydrologis.edc.annotatedclasses.timeseries.shortwaveradiationdirect;

import javax.persistence.Entity;
import javax.persistence.Table;
import static eu.hydrologis.edc.utils.Constants.*;

import eu.hydrologis.edc.annotatedclasses.timeseries.SeriesMonitoringPointsTable;

@Entity
@Table(name = "series_shortwaveradiationdirect_2012", schema = "edcseries")
@org.hibernate.annotations.Table(appliesTo = "series_shortwaveradiationdirect_2012", 
        indexes = @org.hibernate.annotations.Index(
                name = "IDX_TIMESTAMP_MONPOINT_series_shortwaveradiationdirect_2012",
                columnNames = {TIMESTAMPUTC, MONITORINGPOINTS_ID}
))
public class SeriesShortwaveradiationdirect2012 extends SeriesMonitoringPointsTable {
}