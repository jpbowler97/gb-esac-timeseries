package gb.esac.timeseries;

import gb.esac.binner.BinningException;
import java.io.IOException;


public interface ITimeSeriesFileReader {

    TimeSeries readTimeSeriesFile(String filename) throws TimeSeriesFileException, TimeSeriesException, BinningException, IOException ;

}
