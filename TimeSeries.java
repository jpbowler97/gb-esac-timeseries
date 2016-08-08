package gb.esac.timeseries;


import java.io.IOException;
import java.util.Arrays;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import org.apache.log4j.Logger;

import gb.esac.tools.BasicStats;
import gb.esac.tools.DataUtils;


/**
 *

The class <code>TimeSeries</code> is an object used to represent binned data ordered in time. 
It should be viewed as a histogram whose x-axis is in units of time (seconds).

The time axis is discretized into bins defined by two edges per bin; we consider that the bin width 
is constant if the variance of the bin widths is smaller than 1e-6. Gaps between bins are also defined 
by two gap edges; we consider that there is a gap in the time series if there is at least one gap that 
is longer than 1e-6 s.

The Y-axis or the height of each bin shows the total number of counts per bin, and there are therefore 
no errors on the bin heights. The rates are defined as the bin height divided by the bin width.
Upon constructing the time series, if the errors on the rates are set, then the time series is 
considered as non-Poissonian, and the errors are used to derive weights in all statistical operations.

All TimeSeries instances are immutable. The constructors are package-private and used by the public classes 
TimeSeriesMaker and TimeSeriesFileReader. All the setters are private and are therefore used internally to 
define all the fields and properties of a TimeSeries instance. The getters are public and return copies of the
internal objects like binCentres, binWidths, and binHeights for example. There are no static class variables
other than the logger, and so all are instance variables.

In June 2016
- Added new "instnace variables" and the corresponding methods to access these values
String instrument (which is null by default, and the only variable with a public setter method
double errorOnMeanRate 
double[] weightsOnRates
double weightedMeanRate
double errorOnWeightedMean


 *
 * @author <a href="mailto:guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (last modfied: August 2015, ESAC)
 */
public class TimeSeries {

    private static Logger logger  = Logger.getLogger(TimeSeries.class);
    private static String instrument = null;
    // about bins
    private double[] binEdges;
    private double[] leftBinEdges;
    private double[] rightBinEdges;
    private int nBins;
    private double tStart;
    private double tStop;
    private double tMid;
    private double duration;
    private double[] binCentres;
    private double[] binWidths;
    private double[] halfBinWidths;
    private boolean binWidthIsConstant = false;
    private double minBinWidth;
    private double maxBinWidth;
    private double avgBinWidth;
    private double sumOfBinWidths;
    private double binCentreAtMinBinHeight;
    private double binCentreAtMaxBinHeight;
    // gaps
    private double[] gapEdges;
    private double[] gapLengths;
    private int nGaps;
    private int nNonNaNs;
    private double minGap;
    private double maxGap;
    private double meanGap;
    private double sumOfGaps;
    private boolean thereAreGaps = false; 
    // sampling function
    private int nSamplingFunctionBins;
    private double[] samplingFunctionValues;
    private double[] samplingFunctionBinEdges;
    // about bin heights
    private double[] binHeights;
    private double minBinHeight;
    private double maxBinHeight;
    private double meanBinHeight;
    private double sumOfBinHeights;
    private double sumOfSquaredBinHeights;
    private double varianceInBinHeights;
    private double meanDeviationInBinHeights;
    private double skewnessInBinHeights;
    private double kurtosisInBinHeights;
    // about rates
    private double[] rates;
    private double[] errorsOnRates;
    private double[] weightsOnRates;
    private boolean errorsAreSet = false;
    private double minRate;
    private double maxRate;
    private double meanRate;
    private double weightedMeanRate;
    private double errorOnMeanRate;
    private double errorOnWeightedMeanRate;
    private double sumOfRates;
    private double sumOfSquaredRates;
    private double sumOfWeightsOnRates;
    private double varianceInRates;
    private double meanDeviationInRates;
    private double skewnessInRates;
    private double kurtosisInRates;
    private double skewnessStandardError;
    private double kurtosisStandardError;

    //  Package-private constructors
    TimeSeries(TimeSeries ts) {
	if ( ts.errorsAreSet() ) {
	    setBinEdges(ts.tStart(), ts.getBinEdges());
	    setRates(ts.getRates());
	    setErrorsOnRates(ts.getErrorsOnRates());
	}
	else {
	    setBinEdges(ts.tStart(), ts.getBinEdges());
	    setCounts(ts.getBinHeights());
	}
	printRateInfo();
    }

    TimeSeries(double tStart, double[] binEdges, double[] counts) {
	setBinEdges(tStart, binEdges);
	setCounts(counts);
	printRateInfo();
    }
    
    TimeSeries(double tStart, double[] binEdges, double[] rates, double[] errorsOnRates) {
	setBinEdges(tStart, binEdges);
	setRates(rates);
	setErrorsOnRates(errorsOnRates);
	printRateInfo();
    }

    // private info-printing 
    private void printRateInfo() {
	logger.info("Intensities are defined");
	logger.info("  Sum of binHeights = "+this.sumOfBinHeights);
	logger.info("  Mean rate = "+this.meanRate);
	logger.info("  Min rate = "+this.minRate);
	logger.info("  Max rate = "+this.maxRate);
	logger.info("  Variance in rates = "+this.varianceInRates);
    }

    //  Private setters
    private void setBinEdges(double tStart, double[] binEdges) {
	this.binEdges = new double[binEdges.length];
	this.leftBinEdges = new double[binEdges.length/2];
	this.rightBinEdges = new double[binEdges.length/2];
	this.nBins = binEdges.length/2;
	for ( int i=0; i < this.nBins; i++ ) {
	    this.binEdges[2*i] = binEdges[2*i];
	    this.binEdges[2*i+1] = binEdges[2*i+1];
	    this.leftBinEdges[i] = binEdges[2*i];
	    this.rightBinEdges[i] = binEdges[2*i+1];
	}
	this.tStart = tStart;
	this.duration = this.binEdges[this.binEdges.length-1] - this.binEdges[0];
	this.tStop = this.tStart + this.duration;
	this.tMid = (this.tStart + this.tStop)/2;
	logger.info("TimeSeries has "+this.nBins+" bins");
	logger.info("  TStart = "+this.tStart);
    	logger.info("  TStop = "+this.tStop);
	//logger.info("  TMid = "+this.tMid);
	logger.info("  Duration = "+this.duration);
	this.binCentres = new double[this.nBins];
	this.binWidths = new double[this.nBins];
	this.halfBinWidths = new double[this.nBins];
	double min = Double.MAX_VALUE;
	double max = -Double.MAX_VALUE;
	double avg = 0;
	for ( int i=0; i < this.nBins; i++ ) {
	    this.binCentres[i] = (this.binEdges[2*i] + this.binEdges[2*i+1])/2;
	    this.binWidths[i] = this.binEdges[2*i+1] - this.binEdges[2*i];
	    this.halfBinWidths[i] = this.binWidths[i]/2.0;
	    min = Math.min(min, binWidths[i]);
	    max = Math.max(max, binWidths[i]);
	    avg += binWidths[i];
	}
	this.minBinWidth = min;
	this.maxBinWidth = max;
	this.sumOfBinWidths = avg;
	this.avgBinWidth = avg/this.nBins;
	//  Check if bin width is constant, excluding the last bin
	double[] widths = getBinWidths();
	double[] w = new double[widths.length-1];
	for ( int i=0; i < widths.length-1; i++ ) {
	    w[i] = widths[i];
	}
	double var = BasicStats.getVariance(w);
	if ( var < 1e-10 || Double.isNaN(var) ) {
	    this.binWidthIsConstant = true;
	    logger.info("  Bin width is constant = "+this.binWidths[0]);
	}
	else {
	    this.binWidthIsConstant = false;
	    logger.warn("  Bin width is not constant");
	    logger.info("  Min bin width = "+this.minBinWidth);
	    logger.info("  Max bin width = "+this.maxBinWidth);
	    logger.info("  Average bin width = "+this.avgBinWidth);
	}
	// Define gapEdges and sampling function
	this.gapEdges = new double[2*(this.nBins-1)];
	this.gapLengths = new double[this.nBins-1];
	DoubleArrayList samplingFuncValuesList = new DoubleArrayList();
	DoubleArrayList samplingFuncEdgesList = new DoubleArrayList();
	samplingFuncEdgesList.add(this.binEdges[0]);
	samplingFuncEdgesList.add(this.binEdges[1]);
	samplingFuncValuesList.add(1); // time series never starts with a gap, because if there is one, we take it out
	double minGap = Double.MAX_VALUE;
	double maxGap = -Double.MAX_VALUE;
	int nGaps = 0;
	double sumOfGaps = 0;
	for ( int i=1; i < this.nBins; i++ ) {
 	    double gap = this.binEdges[2*i] - this.binEdges[2*i-1];
	    if ( gap > Math.ulp(2*this.binEdges[2*i]) ) {
		nGaps++;
		sumOfGaps += gap;
		samplingFuncEdgesList.add(this.binEdges[2*i-1]);
		samplingFuncEdgesList.add(this.binEdges[2*i]);
		samplingFuncValuesList.add(0);
	    }
	    samplingFuncEdgesList.add(this.binEdges[2*i]);
	    samplingFuncEdgesList.add(this.binEdges[2*i+1]);
	    samplingFuncValuesList.add(1);
	    minGap = Math.min(minGap, gap);
	    maxGap = Math.max(maxGap, gap);
 	    this.gapLengths[i-1] = gap;
	    this.gapEdges[2*(i-1)] = this.binEdges[2*i-1];
	    this.gapEdges[2*(i-1)+1] = this.binEdges[2*i];
	}
	if ( maxGap > Math.ulp(2*this.binEdges[binEdges.length-1]) ) {
	    this.thereAreGaps = true;
	    this.nGaps = nGaps;
	    this.sumOfGaps = sumOfGaps;
	    this.meanGap = sumOfGaps/nGaps;
	    this.maxGap = maxGap;
	    this.minGap = minGap;
	    logger.warn("There are "+nGaps+" gaps in timeline");
	    logger.info("  Total gap time = "+sumOfGaps);
	    logger.info("  Gap fraction wrt duration = "+(sumOfGaps/this.duration));
	    logger.info("  Mean gap = "+meanGap);
	    logger.info("  Max gap = "+maxGap);
	}
	else {
	    this.thereAreGaps = false;
	    this.nGaps = 0;
	    this.sumOfGaps = 0;
	    this.meanGap = 0;
	    this.maxGap = 0;
	    this.minGap = 0;
	    logger.info("No gaps in timeline");
	}
	samplingFuncValuesList.trimToSize();
	samplingFuncEdgesList.trimToSize();
	this.samplingFunctionValues = samplingFuncValuesList.elements();
	this.samplingFunctionBinEdges = samplingFuncEdgesList.elements();
	this.nSamplingFunctionBins = (this.samplingFunctionValues).length;
	logger.info("Sampling function is defined");
	logger.info("  nZeros = "+this.nGaps);
	logger.info("  nOnes = "+this.nBins);
    }


    private void setCounts(double[] counts) {
	this.binHeights = new double[this.nBins];
	this.rates = new double[this.nBins];
	double minBinHeight = Double.MAX_VALUE; 
	double maxBinHeight = -Double.MAX_VALUE; 
	double sumOfBinHeights = 0;
	double sumOfSquaredBinHeights = 0;
	double sumOfRates = 0;
	double sumOfSquaredRates = 0;
	double minRate = Double.MAX_VALUE;
	double maxRate = -Double.MAX_VALUE; 
	int nNaNs = 0;
	int nNonNaNs = 0;
	for ( int i=0; i < this.nBins; i++ ) {
	    this.binHeights[i] = counts[i];
	    this.rates[i] = this.binHeights[i]/this.binWidths[i];
	    if ( Double.isNaN(this.binHeights[i]) ) {
		//logger.warn("NaN encountered in binHeights: index "+i+". Excluding from calculations.");
		nNaNs++;
		this.thereAreGaps = true;
	    }
	    else {
		minBinHeight = Math.min(minBinHeight, this.binHeights[i]);
		maxBinHeight = Math.max(maxBinHeight, this.binHeights[i]);
		sumOfBinHeights += this.binHeights[i];
		sumOfSquaredBinHeights += this.binHeights[i]*this.binHeights[i];
		minRate = Math.min(minRate, this.rates[i]);
		maxRate = Math.max(maxRate, this.rates[i]);
		sumOfRates += this.rates[i];
		sumOfSquaredRates += this.rates[i]*this.rates[i];
		nNonNaNs++;
	    }
	}
	this.nNonNaNs = nNonNaNs;
	this.minRate = minRate;
	this.maxRate = maxRate;
	this.minBinHeight = minBinHeight;
	this.maxBinHeight = maxBinHeight;
	this.sumOfBinHeights = sumOfBinHeights;
	this.sumOfSquaredBinHeights = sumOfSquaredBinHeights;
	this.sumOfSquaredRates = sumOfSquaredRates;
	setStatsOnIntensities();
    }

    public void setRates(double[] r) {
	this.rates = new double[this.nBins];
	this.binHeights = new double[this.nBins];
	double minRate = Double.MAX_VALUE;
	double maxRate = -Double.MAX_VALUE;
	double sumOfRates = 0;
	double sumOfSquaredRates = 0;
	double minBinHeight = Double.MAX_VALUE; 
	double maxBinHeight = -Double.MAX_VALUE; 
	double sumOfBinHeights = 0;
	double sumOfSquaredBinHeights = 0;
	int nNaNs = 0;
	int nNonNaNs = 0;
	for ( int i=0; i < this.nBins; i++ ) {
	    //  Rate
	    double rate = r[i];
	    this.rates[i] = rate;
	    double counts = this.rates[i]*this.binWidths[i];
	    this.binHeights[i] = counts;
	    if ( Double.isNaN(rate) ) {
		nNaNs++;
		thereAreGaps = true;
	    }
	    else {
		minRate = Math.min(minRate, rate);
		maxRate = Math.max(maxRate, rate);
		sumOfRates += rate;
		sumOfSquaredRates += rate*rate;
		minBinHeight = Math.min(minBinHeight, this.binHeights[i]);
		maxBinHeight = Math.max(maxBinHeight, this.binHeights[i]);
		sumOfBinHeights += counts;
		sumOfSquaredBinHeights += counts*counts;
		nNonNaNs++;
	    }
	}
	if ( thereAreGaps ) {
	    logger.warn("There are "+nNaNs+" NaN values in the RATE column");
	}
	this.nNonNaNs = nNonNaNs;
	this.minRate = minRate;
	this.maxRate = maxRate;
	this.minBinHeight = minBinHeight;
	this.maxBinHeight = maxBinHeight;
	this.sumOfRates = sumOfRates;
	this.sumOfBinHeights = sumOfBinHeights;
	this.sumOfSquaredBinHeights = sumOfSquaredBinHeights;
	this.sumOfSquaredRates = sumOfSquaredRates;

    }

    public void setErrorsOnRates(double[] errors) {
	this.errorsOnRates = new double[this.nBins];
	this.weightsOnRates = new double[this.nBins];
	double sum = 0;
	for ( int i=0; i < this.nBins; i++ ) {
	    if ( Double.isNaN(errors[i]) ) {
		if ( ! Double.isNaN(this.rates[i]) ) {
		    logger.warn("There is a NaN value in errors whose corresponding rate is not NaN. Setting error from mean counts per bin.");
		    double uncertainty = Math.sqrt(this.meanBinHeight);
		    this.errorsOnRates[i] = uncertainty/this.binWidths[i];
		}
	    }
	    else {
		this.errorsOnRates[i] = errors[i];
	    }
	    this.weightsOnRates[i] = 1./Math.pow(this.errorsOnRates[i], 2);
	    sum += this.weightsOnRates[i];
	}
	this.sumOfWeightsOnRates = sum;
	this.weightedMeanRate = Descriptive.weightedMean(new DoubleArrayList(this.rates), new DoubleArrayList(this.weightsOnRates));
	this.errorsAreSet = true;
	setStatsOnIntensities();
    }

    private void setStatsOnIntensities() {
	this.binCentreAtMinBinHeight = this.binCentres[DataUtils.getIndex(this.minBinHeight, this.binHeights)];
	this.binCentreAtMaxBinHeight = this.binCentres[DataUtils.getIndex(this.maxBinHeight, this.binHeights)];
	this.meanBinHeight = this.sumOfBinHeights/this.nNonNaNs;
	this.meanRate = this.sumOfBinHeights/this.sumOfBinWidths;
 	this.varianceInBinHeights = Descriptive.sampleVariance(this.nNonNaNs, this.sumOfBinHeights, this.sumOfSquaredBinHeights);
 	this.varianceInRates = Descriptive.sampleVariance(this.nNonNaNs, this.sumOfRates, this.sumOfSquaredRates);
	this.errorOnMeanRate = Math.sqrt(this.varianceInRates/this.nNonNaNs);
	this.errorOnWeightedMeanRate = 1./Math.sqrt(this.sumOfWeightsOnRates);
	this.meanDeviationInBinHeights = Descriptive.meanDeviation(new DoubleArrayList(this.binHeights), this.meanBinHeight);
	this.meanDeviationInRates = Descriptive.meanDeviation(new DoubleArrayList(this.rates), this.meanRate);
	this.skewnessInBinHeights = Descriptive.sampleSkew(new DoubleArrayList(this.binHeights), this.meanBinHeight, this.varianceInBinHeights);
	this.skewnessInRates = Descriptive.sampleSkew(new DoubleArrayList(this.rates), this.meanBinHeight, this.varianceInBinHeights);
	this.skewnessStandardError = Descriptive.sampleSkewStandardError(this.nBins);
	this.kurtosisInBinHeights = Descriptive.sampleKurtosis(new DoubleArrayList(this.binHeights), this.meanBinHeight, this.varianceInBinHeights);
	this.kurtosisInRates = Descriptive.sampleKurtosis(new DoubleArrayList(this.rates), this.meanBinHeight, this.varianceInBinHeights);
	this.kurtosisStandardError =  Descriptive.sampleKurtosisStandardError(this.nBins);
    }

    // The only public setter methods
    public void setInstrument(String instrument) {
	this.instrument = instrument;
    }
    public String instrument() {
	return new String(this.instrument);
    }
    
    //  About Bins
    public int nBins() { return this.nBins; }
    public double tStart() { return this.tStart; }
    public double tStop() { return this.tStop; }
    public double tMid() { return this.tMid; }
    public double duration() { return this.duration; }
    public double[] getBinCentres() { return Arrays.copyOf(this.binCentres, this.binCentres.length); }
    public double[] getBinWidths() { return Arrays.copyOf(this.binWidths, this.binWidths.length); }
    public double[] getHalfBinWidths() { return Arrays.copyOf(this.halfBinWidths, this.halfBinWidths.length); }
    public double[] getBinEdges() { return Arrays.copyOf(this.binEdges, this.binEdges.length); }
    public double[] getLeftBinEdges() { return Arrays.copyOf(this.leftBinEdges, this.leftBinEdges.length); }
    public double[] getRightBinEdges() { return Arrays.copyOf(this.rightBinEdges, this.rightBinEdges.length); }
    public double binCentreAtMinBinHeight() { return this.binCentreAtMinBinHeight; }
    public double binCentreAtMaxBinHeight() { return this.binCentreAtMaxBinHeight; }
    public double minBinWidth() { return this.minBinWidth; }
    public double maxBinWidth() { return this.maxBinWidth; }
    public double avgBinWidth() { return this.avgBinWidth; }
    public double sumOfBinWidths() { return this.sumOfBinWidths; }
    public double binWidth()  throws TimeSeriesException {
	if ( this.binWidthIsConstant )
	    return binWidths[0];
	else {
	    throw new TimeSeriesException("BinWidth is not constant. Use getBinWidths()");
	}
    }

    //  About Gaps
    public int nGaps() { return this.nGaps; }
    public double[] getGapEdges() { return Arrays.copyOf(this.gapEdges, this.gapEdges.length); }
    public double[] getGapLengths() { return Arrays.copyOf(this.gapLengths, this.gapLengths.length); }
    public double meanGap() { return this.meanGap; }
    public double minGap() { return this.minGap; }
    public double maxGap() { return this.maxGap; }
    public double sumOfGaps() { return this.sumOfGaps; }
    public int nSamplingFunctionBins() { return this.nSamplingFunctionBins; }
    public double[] getSamplingFunctionValues() { return Arrays.copyOf(this.samplingFunctionValues, this.samplingFunctionValues.length); }
    public double[] getSamplingFunctionBinEdges() { return Arrays.copyOf(this.samplingFunctionBinEdges, this.samplingFunctionBinEdges.length); }

    //  About Intensities
    public double[] getBinHeights() { return Arrays.copyOf(this.binHeights, this.binHeights.length); }
    public double sumOfBinHeights() { return this.sumOfBinHeights; }
    public double meanBinHeight() { return this.meanBinHeight; }
    public double minBinHeight() { return this.minBinHeight; }
    public double maxBinHeight() { return this.maxBinHeight; }
    public double varianceInBinHeights() { return this.varianceInBinHeights; }
    public double meanDeviationInBinHeights() { return this.meanDeviationInBinHeights; }
    public double kurtosisInBinHeights() { return this.kurtosisInBinHeights; }
    public double kurtosisStandardError() { return this.kurtosisStandardError; }
    public double skewnessInBinHeights() { return this.skewnessInBinHeights; }
    public double skewnessStandardError() { return this.skewnessStandardError; }
    public double[] getRates() { return Arrays.copyOf(this.rates, this.rates.length); }
    public double meanRate() { return this.meanRate; }
    public double minRate() { return this.minRate; }
    public double maxRate() { return this.maxRate; }
    public double errorOnMeanRate() { return this.errorOnMeanRate; }
    public double weightedMeanRate() { return this.weightedMeanRate; }
    public double errorOnWeightedMeanRate() { return this.errorOnWeightedMeanRate; }
    public double varianceInRates() { return this.varianceInRates; }
    public double meanDeviationInRates() { return this.meanDeviationInRates; }
    public double kurtosisInRates() { return this.kurtosisInRates; }
    public double skewnessInRates() { return this.skewnessInRates; }
    public double[] getErrorsOnRates() {
	if ( errorsAreSet ) 
	    return Arrays.copyOf(this.errorsOnRates, this.errorsOnRates.length);
	else {
	    double[] errorsOnRates = new double[this.nBins];
	    double uncertainty = Math.sqrt(this.meanBinHeight);
	    for ( int i=0; i < this.nBins; i++ ) {
		//double uncertainty = Math.sqrt(this.binHeights[i]);
		errorsOnRates[i] = uncertainty/this.binWidths[i];
	    }
	    return errorsOnRates;
	}
    }
    public double[] getMeanSubtractedRates() { 
	double[] meanSubRates = new double[this.nBins];
	for ( int i=0; i < this.nBins; i++ ) {
	    meanSubRates[i] = this.rates[i] - this.meanRate;
	}
	return meanSubRates;
    }
    public double[] getMeanSubtractedBinHeights() { 
	double[] meanSubBinHeights = new double[this.nBins];
	for ( int i=0; i < this.nBins; i++ ) {
	    meanSubBinHeights[i] = this.binHeights[i] - this.meanBinHeight;
	}
	return meanSubBinHeights;
    }

    //  Boolean checkers
    public boolean binWidthIsConstant() { return this.binWidthIsConstant; }
    public boolean thereAreGaps() { return this.thereAreGaps; }
    public boolean errorsAreSet() { return this.errorsAreSet; }

    //  Write to ascii file
    public void writeCountsAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeCountsAsQDP(this, filename);
    }

    public void writeCountsAsQDP(double[] function, String filename) throws IOException {
	TimeSeriesWriter.writeCountsAsQDP(this, function, filename);
    }

    public void writeRatesAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAsQDP(this, filename);
    }

    public void writeRatesAsQDP(double[] function, String filename) throws IOException {
	TimeSeriesWriter.writeRatesAsQDP(this, function, filename);
    }

    public void writeCountsAndSamplingAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeCountsAndSamplingAsQDP(this, filename);
    }

    public void writeRatesAndSamplingAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAndSamplingAsQDP(this, filename);
    }

}
