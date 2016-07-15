package gb.esac.timeseries;

import nom.tam.fits.Fits;



public class TestTimeSeries {

    public static void main(String[] args) throws Exception  {


 	String filename = "flare.fits";
    	filename="empty.dat";
      	filename="/Users/gbelanger/Documents/integ/idx/date2rev.sh";
      	filename="/Users/gbelanger/Documents/integ/idx/point.lis";
     	filename="/Users/gbelanger/Documents/integ/idx/GNRL-SCWG-GRP-IDX.fits.gz";
     	filename="/Users/gbelanger/Documents/pubs/apj/belanger10/qdp/gx_10s_lc.fits";
//      	filename="/Users/gbelanger/Documents/pubs/apj/belanger10/qdp/gx_100s_lc.fits";
//     	filename="/Users/gbelanger/Documents/pubs/apj/belanger10/qdp/windowing_evlist.fits";
//  	filename="/Users/gbelanger/Documents/pubs/apj/belanger10/qdp/GX_1+4_deltat10.iilight.lc.fits";
//  	filename="simEvlist.fits";
//  	filename="binCentresAndCounts.qdp";

// 	filename = "/Users/gbelanger/Documents/astroData/rxte/saxj1808_lc1.fits";
// 	TimeSeries lc1 = TimeSeriesMaker.makeTimeSeries(filename);
// 	filename = "/Users/gbelanger/Documents/astroData/rxte/saxj1808_lc3.fits";
// 	TimeSeries lc2 = TimeSeriesMaker.makeTimeSeries(filename);
// 	filename = "/Users/gbelanger/Documents/astroData/rxte/saxj1808_lc2.fits";
// 	TimeSeries lc3 = TimeSeriesMaker.makeTimeSeries(filename);
// 	TimeSeries lc = TimeSeriesOperations.combine(new TimeSeries[] {lc1, lc2, lc3});


 	filename = "/Users/gbelanger/Documents/astroData/rxte/saxj1808_lightcurve_16msec.lc";
	TimeSeries lc = TimeSeriesMaker.makeTimeSeries(filename);

	double tStart = 0;
	double tStop = 10;
	TimeSeries seg = TimeSeriesOperations.getSegment(lc, tStart, tStop);
	seg.writeCountsAsQDP("seg.qdp");


	//  There is a problem with the resampling !!!

	//TimeSeries reb = TimeSeriesResampler.resampleToClosestPowerOfTwo(lc);
	//reb.writeRatesAsQDP("rates.qdp");

	

    }
}