package gb.esac.timeseries;


import gb.esac.binner.BinningUtils;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.tools.BasicStats;
import gb.esac.tools.MinMax;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import org.apache.log4j.Logger;
import java.util.Date;


final class TimeSeriesWriter {

    private static Logger logger  = Logger.getLogger(TimeSeriesWriter.class);
    private static String classname = (TimeSeriesWriter.class).getCanonicalName();
    private static int bufferSize = 256000;
    private static DecimalFormat decimals = new DecimalFormat("0.00");
    private static DecimalFormat stats = new DecimalFormat("0.00E00");
    private static DecimalFormat number = new DecimalFormat("0.000");
    private static DecimalFormat noDigits = new DecimalFormat("0");
    private static DecimalFormat oneDigit = new DecimalFormat("0.0");
    private static DecimalFormat twoDigits = new DecimalFormat("0.00");
    private static DecimalFormat threeDigits = new DecimalFormat("0.000");
    private static DecimalFormat timeFormat = new DecimalFormat("0.000E0");

    // METHODS

    ////  Counts
    public static void writeCountsAsQDP(TimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	String[] header = makeQDPHeader(ts, "counts", filename);
	AsciiDataFileWriter writer = new AsciiDataFileWriter(filename);
	writer.writeData(header, ts.getBinCentres(), ts.getHalfBinWidths(), ts.getBinHeights());
	logger.info(tsClassName+" in counts written to "+filename);
    }

    public static void writeCountsAsQDP(TimeSeries ts, double[] function, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	String[] header = makeQDPHeader(ts, "counts", filename);
	AsciiDataFileWriter writer = new AsciiDataFileWriter(filename);
	writer.writeData(header, ts.getBinCentres(), ts.getHalfBinWidths(), ts.getBinHeights(), function);
	logger.info(tsClassName+" in counts written to "+filename);
    }

    public static void writeCountsAndSamplingAsQDP(TimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filename), bufferSize));
	String[] header = makeQDPHeader(ts, "counts", filename);
	//  Write the header 
	for ( int i=0; i < header.length; i++ )  pw.println(header[i]);
	//  Write the TimeSeries 
	double[] c1 = ts.getBinCentres();
	double[] c2 = ts.getHalfBinWidths();
	double[] c3 = ts.getBinHeights();
	int[] lengths = new int[] {c1.length, c2.length, c3.length};
	double var = BasicStats.getVariance(lengths);
	if ( var != 0 ) {
	    logger.warn("input column data of different lengths. Using min.");
	}
	int nbins = (new Double(MinMax.getMin(lengths))).intValue();
	for ( int i=0; i < nbins; i++ ) {
	    pw.println(c1[i] +"\t"+ c2[i] +"\t"+ c3[i] +"\t");
	}
	pw.println("NO NO NO");
	//  Write the sampling function
	double[] edges = ts.getSamplingFunctionBinEdges();
	double[] centres = BinningUtils.getBinCentresFromBinEdges(edges);
	double[] halfWidths = BinningUtils.getHalfBinWidthsFromBinEdges(edges);
	double[] func = ts.getSamplingFunctionValues();
	for ( int i=0; i < func.length; i++ ) {
	    pw.println(centres[i] +"\t"+ halfWidths[i] +"\t"+ func[i] +"\t" );
	}
	pw.close();
	logger.info(tsClassName+" in counts and sampling function written to "+filename);
    }

    ////  Rates
    public static void writeRatesAsQDP(TimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	String[] header = makeQDPHeader(ts, "rates", filename);
	AsciiDataFileWriter writer = new AsciiDataFileWriter(filename);
	writer.writeData(header, ts.getBinCentres(), ts.getHalfBinWidths(), ts.getRates(), ts.getErrorsOnRates());
	logger.info(tsClassName+" in rates written to "+filename);
    }
    
    public static void writeRatesAsQDP(TimeSeries ts, double[] function, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	String[] header = makeQDPHeader(ts, "rates", filename);
	AsciiDataFileWriter writer = new AsciiDataFileWriter(filename);
	writer.writeData(header, ts.getBinCentres(), ts.getHalfBinWidths(), ts.getRates(), ts.getErrorsOnRates(), function);
	logger.info(tsClassName+" in rates written to "+filename);
    }

    public static void writeRatesAndSamplingAsQDP(TimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filename), bufferSize));
	String[] header = makeQDPHeader(ts, "rates", filename);
	//  Write the header 
	for ( int i=0; i < header.length; i++ )  pw.println(header[i]);
	//  Write the TimeSeries 
	double[] c1 = ts.getBinCentres();
	double[] c2 = ts.getHalfBinWidths();
	double[] c3 = ts.getRates();
	double[] c4 = ts.getErrorsOnRates();
	int[] lengths = new int[] {c1.length, c2.length, c3.length, c4.length};
	double var = BasicStats.getVariance(lengths);
	if ( var != 0 ) {
	    logger.warn("input column data of different lengths. Using min.");
	}
	int nbins = (new Double(MinMax.getMin(lengths))).intValue();
	for ( int i=0; i < nbins; i++ ) {
	    pw.println(c1[i] +"\t"+ c2[i] +"\t"+ c3[i] +"\t"+ c4[i] +"\t" );
	}
	pw.println("NO NO NO");
	//  Write the sampling function
	double[] edges = ts.getSamplingFunctionBinEdges();
	double[] centres = BinningUtils.getBinCentresFromBinEdges(edges);
	double[] halfWidths = BinningUtils.getHalfBinWidthsFromBinEdges(edges);
	double[] func = ts.getSamplingFunctionValues();
	for ( int i=0; i < func.length; i++ ) {
	    pw.println(centres[i] +"\t"+ halfWidths[i] +"\t"+ func[i] +"\t");
	}
	pw.close();
	logger.info(tsClassName+" in rates and sampling function written to "+filename);
    }


    //  Methods acting on CodedMaskTimeSeries (only rates)

    public static void writeRatesAsQDP(CodedMaskTimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	String[] header = makeCodedMaskQDPHeader(ts, filename);
	AsciiDataFileWriter writer = new AsciiDataFileWriter(filename);
	writer.writeData(header, ts.getBinCentres(), ts.getHalfBinWidths(), ts.getRates(), ts.getErrorsOnRates());
	logger.info(tsClassName+" in rates written to "+filename);
    }

    public static void writeRatesAndDistToAxisAsQDP(CodedMaskTimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	String[] header = makeCodedMaskQDPHeader(ts, filename);
	AsciiDataFileWriter writer = new AsciiDataFileWriter(filename);
	writer.writeData(header, ts.getBinCentres(), ts.getHalfBinWidths(), ts.getRates(), ts.getErrorsOnRates(), ts.getDistToPointingAxis());
	logger.info(tsClassName+" in rates written to "+filename);
    }

    public static void writeRatesAsPLT(CodedMaskTimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filename), bufferSize));
	String[] header = makeCodedMaskPLTHeader(ts, filename);
	//  Write the header 
	for ( int i=0; i < header.length; i++ ) pw.println(header[i]);
	double[] t = ts.getBinCentres();
	double[] rates = ts.getRates();
	double[] rateErrs = ts.getErrorsOnRates();
	for ( int j=0; j < rates.length; j++ ) {
	    if ( !Double.isNaN(rates[j]) && rates[j] != 0.0 ) {
		double low = rates[j] - rateErrs[j];
		double hi = rates[j] + rateErrs[j];
		pw.println(t[j] +", "+ rates[j] +", "+ low +", "+ hi);
	    }
	}
	pw.flush();
	pw.close();
	logger.info(tsClassName+" in rates written to "+filename);
    }

    public static void writeRatesAsXML(CodedMaskTimeSeries ts, String filename) throws IOException {
	String tsClassName = ts.getClass().getCanonicalName();
	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filename), bufferSize));
	String[] header = makeCodedMaskXMLHeader(ts, filename);
	//  Write the header 
	for ( int i=0; i < header.length; i++ ) pw.println(header[i]);
	double[] t = ts.getBinCentres();
	double[] rates = ts.getRates();
	double[] rateErrs = ts.getErrorsOnRates();
	for ( int j=0; j < rates.length; j++ ) {
	    if ( !Double.isNaN(rates[j]) && rates[j] != 0.0 ) {
		double low = rates[j] - rateErrs[j];
		double hi = rates[j] + rateErrs[j];
		pw.println( "<p x='"+ t[j] +"' y='"+ rates[j] +"' "+ "lowErrorBar='"+low+ "' " + "highErrorBar='"+hi+ "'/>");
	    }
	}
	pw.println("</dataset>");
	pw.println("</plot>");
	pw.flush();
	pw.close();
	logger.info(tsClassName+" in rates written to "+filename);	
    }

    //  Methods to make headers
    
    private static String[] makeQDPHeader(TimeSeries ts, String type, String filename) throws IOException {
	// type can be "counts" or "rates"
	double[] binEdges = ts.getBinEdges();
	double xmin =  binEdges[0] - 0.05*ts.duration();
	double xmax =  binEdges[binEdges.length-1] + 0.05*ts.duration();
	//  Default is type.equals("counts")
	String yLabel = "Intensity (cts)";
	String line = "STEP ON 2";
	String serr = "SERR 1";
	double min = ts.minBinHeight();
	double max = ts.maxBinHeight();
	double binWidth = 0;
	try { binWidth = ts.binWidth(); }
	catch ( TimeSeriesException e ) { }	    
	if ( binWidth != 0 ) {
	    yLabel = "Intensity (cts per "+decimals.format(binWidth)+" s)";
	}
	if ( type.equals("rates") ) {
	    line = "OFF 2";
	    serr = "SERR 1 2";
	    double mean = ts.meanRate();
	    double sigma = Math.sqrt(ts.varianceInRates());
	    min = mean - 5*sigma;
	    max = mean + 5*sigma;
	    yLabel="Intensity (s\\u-1\\d)";
	}
	double yRange = max - min;
	double ymin = min; // -0.1*yRange;
	double ymax = max; // +0.1*yRange;
	String[] header = new String[] {
	    "! Filename: " + filename,
	    "! Produced by: "+ classname,
	    "! Date: "+new Date(),
	    "! Author: G. Belanger - ESA/ESAC",
	    "!",
	    "DEV /XS",
	    "READ "+serr,
	    "LAB T", "LAB F",
	    "TIME OFF",
	    "LINE "+line, 
	    "LINE STEP ON 3",
	    "LW 4", "CS 1.5",
	    "LAB X Time (s)  since "+ts.tStart(),
	    "LAB Y "+yLabel,
	    "VIEW 0.1 0.2 0.9 0.8",
	    "R X "+twoDigits.format(xmin)+" "+twoDigits.format(xmax),
	    "R Y "+twoDigits.format(ymin)+" "+twoDigits.format(ymax),
	    "!"
	};
	return header;
    }

    private static String[] makeCodedMaskQDPHeader(CodedMaskTimeSeries ts, String filename) throws IOException {
	String serr = "SERR 1 2";
	String instrument = ts.instrument();
	double emin = ts.energyRangeMin();
	double emax = ts.energyRangeMax();
	String eminStr = null;
	String emaxStr = null;
	if ( Math.round(emin) == emin && Math.round(emax) == emax ) {
	    eminStr = noDigits.format(emin);
	    emaxStr = noDigits.format(emax);
	}
	else {
	    eminStr = String.valueOf(emin);
	    emaxStr = String.valueOf(emax);
	}
	double ra = ts.targetRA();
	double dec = ts.targetDec();
	String oTitle = instrument+" Time Series ("+eminStr+"-"+emaxStr+" keV)";
	String title = "(RA="+ra+", Dec="+dec+")";
	String y2Label="Intensity (s\\u-1\\d)";
	String y3Label="Off-axis (deg)";
	double fracFC = ts.fullyCodedFraction();
	double fracPC = 1. - fracFC;
	int nGood = ts.nNonNaN_exposuresOnTarget();
	int duration = (int) ts.duration();
	int ontime = (int) ts.sumOfEffectivePointingDurations();
	int effectiveExposure = (int) ts.sumOfExposuresOnTarget();
	double weightedMean = ts.weightedMeanRate();
	double errorOnWMean = ts.errorOnWeightedMeanRate();
	double signif = weightedMean/errorOnWMean;
	double[] binEdges = ts.getBinEdges();
	double xmin =  binEdges[0] - 0.05*ts.duration();
	double xmax =  binEdges[binEdges.length-1] + 0.05*ts.duration();
	double xRange = xmax - xmin;
	double mean = ts.weightedMeanRate();
	double sigma = Math.sqrt(ts.varianceInRates());
	double min = mean - 5*sigma;
	double max = mean + 5*sigma;
	double yRange = max - min;
	double ymin = min; // -0.05*yRange;
	double ymax = max; // +0.05*yRange;
	String[] header = new String[] {
	    "! Filename: " + filename,
	    "! Produced by: "+ classname,
	    "! Date: "+new Date(),
	    "! Author: G. Belanger - ESA/ESAC",
	    "!",
	    "DEV /XS",
	    "READ "+serr,
	    "PLOT VERT",
	    "LAB F",
	    "TIME OFF",
	    "LW 3", 
	    "CS 1.0",
	    "LAB OT "+oTitle, 
	    "LAB T "+title,
	    "MA 1 ON",
	    "CO 2 ON 3",
	    "LAB X Time (s) since MJD "+(ts.tStart()/86400),
	    "LAB Y2 "+y2Label,
	    "LAB Y3 "+y3Label,
	    "R Y2 "+twoDigits.format(ymin)+" "+twoDigits.format(ymax), 
	    "R Y3 -2 17",
	    "R X "+twoDigits.format(xmin)+" "+twoDigits.format(xmax),
	    "VIEW 0.1 0.2 0.9 0.8",
	    "WIN 3",
	    "LOC 0 0.1 1 0.4",
	    "LAB 100 POS "+twoDigits.format(xmin + 0.01*xRange)+" 4.15 LINE 0 0.98 \"",
	    "LAB 100 LS 4 JUST CEN",
	    "LAB 101 VPOS 0.12 0.23 \""+(int)Math.round(fracPC*100)+"%\" CS 0.55 JUST CEN",
	    "LAB 102 VPOS 0.12 0.207 \""+(int)Math.round(fracFC*100)+"%\" CS 0.55 JUST CEN",
	    //"LAB 101 POS "+twoDigits.format(xmin + 0.025*xRange)+" 5.35 \""+(int)Math.round(fracPC*100)+"%\" CS 0.55 JUST CEN",
	    //"LAB 102 POS "+twoDigits.format(xmin + 0.025*xRange)+" 2.95 \""+(int)Math.round(fracFC*100)+"%\" CS 0.55 JUST CEN",
	    "WIN 2",
	    "LOC 0 0.22 1 0.92",
	    "LAB 201 VPOS 0.12 0.75 \"Weighted mean rate = "+threeDigits.format(weightedMean)+" +/- "+threeDigits.format(errorOnWMean)+" s\\u-1\\d\" CS 0.55 JUST LEFT",
	    "LAB 202 VPOS 0.12 0.73 \"Detection significance = "+threeDigits.format(signif)+"\" CS 0.55 JUST LEFT",
	    "LAB 203 VPOS 0.88 0.75 \"Number of observations = "+nGood+"\" CS 0.55 JUST RIGHT",
	    "LAB 204 VPOS 0.88 0.73 \"Time series duration = "+duration+" s\" CS 0.55 JUST RIGHT",
	    "LAB 205 VPOS 0.88 0.71 \"Sum of pointings = "+ontime+" s\" CS 0.55 JUST RIGHT",
	    "LAB 206 VPOS 0.88 0.69 \"Exposure on target = "+effectiveExposure+" s\" CS 0.55 JUST RIGHT",
	    "!"
	};
	return header;
    }

    private static String[] makeCodedMaskPLTHeader(CodedMaskTimeSeries ts, String filename) throws IOException {
	int nGood = ts.nNonNaN_exposuresOnTarget();
	int duration = (int) ts.duration();
	int effectiveExposure = (int) ts.sumOfExposuresOnTarget();
	double weightedMean = ts.weightedMeanRate();
	double errorOnWMean = ts.errorOnWeightedMeanRate();
	double signif = weightedMean/errorOnWMean;
	String instrument = ts.instrument();
	double emin = ts.energyRangeMin();
	double emax = ts.energyRangeMax();
	String eminStr = null;
	String emaxStr = null;
	if ( Math.round(emin) == emin && Math.round(emax) == emax ) {
	    eminStr = noDigits.format(emin);
	    emaxStr = noDigits.format(emax);
	}
	else {
	    eminStr = String.valueOf(emin);
	    emaxStr = String.valueOf(emax);
	}
	double ra = ts.targetRA();
	double dec = ts.targetDec();
	String titleText = instrument+" Time Series ("+eminStr+"-"+emaxStr+" keV) at RA="+ra+", Dec="+dec;
	String[] header = new String[]{
	    "# Data file for Ptplot 5.3",
	    "# Id: "+ filename +", "+ (new Date()) +" ptII Exp",
	    "#",
	    "# Filename: "+ filename,
	    "# Produced by: "+ classname,
	    "# Author: G. Belanger - ESA/ESAC ",
	    "# Summary:",
	    "#   Selected Observations = "+nGood,
	    "#   Time series duration = "+ duration +" s",
	    "#   Effective exposure = "+ effectiveExposure +" s",
	    "#   Weighted mean count rate = "+ weightedMean +" +/- "+ errorOnWMean,
	    "#   Detection significance = "+ signif,
	    "#",
	    "TitleText: "+ titleText,
	    "Marks: dots",
	    "Lines: off",
	    "XLabel: Time (s) since MJD "+(ts.tStart()/86400),
	    "YLabel: Count rate (cts/s)",
	    "DataSet:"
	};
	return header;
    }

    private static String[] makeCodedMaskXMLHeader(CodedMaskTimeSeries ts, String filename) throws IOException {
	int nGood = ts.nNonNaN_exposuresOnTarget();
	int duration = (int) ts.duration();
	int effectiveExposure = (int) ts.sumOfExposuresOnTarget();
	double weightedMean = ts.weightedMeanRate();
	double errorOnWMean = ts.errorOnWeightedMeanRate();
	double signif = weightedMean/errorOnWMean;
	String instrument = ts.instrument();
	double emin = ts.energyRangeMin();
	double emax = ts.energyRangeMax();
	String eminStr = null;
	String emaxStr = null;
	if ( Math.round(emin) == emin && Math.round(emax) == emax ) {
	    eminStr = noDigits.format(emin);
	    emaxStr = noDigits.format(emax);
	}
	else {
	    eminStr = String.valueOf(emin);
	    emaxStr = String.valueOf(emax);
	}
	double ra = ts.targetRA();
	double dec = ts.targetDec();
	String titleText = instrument+" Time Series ("+eminStr+"-"+emaxStr+" keV) at RA="+ra+", Dec="+dec;
	String[] header = new String[] {
	    "<?xml version ='1.0' standalone='no'?>",
	    "<!DOCTYPE model PUBLIC '-//UC Berkeley//DTD PlotML 1//EN'",
	    "'http://ptolemy.eecs.berkeley.edu/xml/dtd/PlotML_1.dtd'>",
	    "<plot>",
	    "<!-- Ptolemy plot, version 5.3, PlotML format. -->",
	    "<!-- Filename: "+ filename +" -->",
	    "<!-- Produced by: "+ classname +"-->",
	    "<!-- Date: "+ new Date() +"-->",
	    "<!-- Author: G. Belanger - ESA/ESAC -->",
	    "<!-- Summary: -->",
	    "<!--   Selected Observations = "+ nGood +" -->",
	    "<!--   Time series duration = "+ duration +" s -->",
	    "<!--   Effective exposure = "+ effectiveExposure +" s -->",
	    "<!--   Weighted mean count rate = "+ weightedMean +" +/- "+ errorOnWMean +" -->",
	    "<!--   Detection significance = "+ signif +" -->",
	    "<!-- -->",
	    "<title>"+ titleText +"</title>",
	    "<xLabel>Time (s)  since MJD "+(ts.tStart()/86400)+"</title>",
	    "<yLabel>Count rate (cts/s)</yLabel>",
	    "<dataset marks='dots' connected='no' stems='no'>"
	};
	return header;
    }
}
