package gb.esac.timeseries;

import java.io.IOException;
import java.util.Arrays;

import org.apache.log4j.Logger;
import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import java.awt.geom.Point2D;
import jsky.coords.WorldCoords;

/**
 *

The class <code>CodedMaskTimeSeries</code> extends <code>TimeSeries</code> to store the information
that is relevant and important for a time series made from a coded mask instrument. This includes
coordinates of the target, the energy range used to make the  time series, the angular distance
from the pointing axis to the target, and the maximum distance for full coding.

See documentation in <code>TimeSeries</code> for details on the parent class.

 *
 * @author <a href="mailto:guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (last modfied: June 2016, ESAC)

 */


public class CodedMaskTimeSeries extends TimeSeries {

    //  Class variables
    private static Logger logger  = Logger.getLogger(CodedMaskTimeSeries.class);
    private double targetRA;
    private double targetDec;
    private double energyRangeMin;
    private double energyRangeMax;
    private double maxDistForFullCoding;
    private int nFullyCoded;
    private double fullyCodedFraction;
    
    // about Ra, Dec of pointing
    private Point2D.Double[] raDecsOfPointings;
    private double[] rasOfPointings;
    private double[] decsOfPointings;
    private double minRaOfPointings;
    private double maxRaOfPointings;
    private double minDecOfPointings;
    private double maxDecOfPointings;
    // private double meanRaDecOfPointings;
    // private double meanRaOfPointings;
    // private double meanDecOfPointings;
    // private double varianceInRasOfPointings;
    // private double varianceInDecsOfPointings;
    // private double meanDeviationInRasOfPointings;
    // private double meanDeviationInDecsOfPointings;
    private int nNonNaN_raDecs;
    
    // about pointing durations
    private double[] effectivePointingDurations;
    private double[] deadTimeDurations;
    private double[] liveTimeFractions;
    private double[] deadTimeFractions;
    private double minEffectivePointingDuration;
    private double maxEffectivePointingDuration;
    private double meanEffectivePointingDuration;
    private double varianceInEffectivePointingDuration;
    private double meanDeviationInEffectivePointingDuration;    
    private double sumOfEffectivePointingDurations;
    private double sumOfSquaredEffectivePointingDurations;
    private double minDeadTimeDuration;
    private double maxDeadTimeDuration;
    private double meanDeadTimeDuration;
    private double sumOfDeadTimeDurations;
    private double minLiveTimeFraction;
    private double maxLiveTimeFraction;
    private double meanLiveTimeFraction;
    private double sumOfLiveTimeFractions;
    private double minDeadTimeFraction;
    private double maxDeadTimeFraction;
    private double meanDeadTimeFraction;
    private double sumOfDeadTimeFractions;
    private int nNonNaN_pointingDurations;
    
    // about distances to pointing axis
    private double[] distToPointingAxis;
    private double minDistToPointingAxis;
    private double maxDistToPointingAxis;
    private double sumOfDistToPointingAxis;
    private double sumOfSquaredDistToPointingAxis;
    private double meanDistToPointingAxis;
    private double varianceInDistToPointingAxis;
    private double meanDeviationInDistToPointingAxis;
    private int nNonNaN_distToPointingAxis;
    
    // about exposure on target
    private double[] exposuresOnTarget;
    private double minExposureOnTarget;
    private double maxExposureOnTarget;
    private double sumOfExposuresOnTarget;
    private double sumOfSquaredExposuresOnTarget;
    private double meanExposureOnTarget;
    private double varianceInExposureOnTarget;
    private double meanDeviationInExposureOnTarget;
    private int nNonNaN_exposuresOnTarget;
    
    // Package-private constructor specific to CodedMaskTimeSeries
    //   with additional attributes of raDecsOfPointings, distToPointingAxis, and exposuresOnTarget
    CodedMaskTimeSeries(Point2D.Double targetRaDec, Point2D.Double energyRangeMinMax, String instrument, double maxDistForFullCoding, double tStart, double[] binEdges, double[] effectivePointingDurations, double[] rates, double[] errors, Point2D.Double[] raDecsOfPointings, double[] exposuresOnTarget) throws IllegalArgumentException {
	super(tStart, binEdges, rates, errors);
	if ( Double.isNaN(targetRaDec.getX()) || Double.isNaN(targetRaDec.getY()) ) {
	    throw new IllegalArgumentException("CodedMaskTimeSeries requires attributes: targetRA, targetDec");
	}
	if ( Double.isNaN(energyRangeMinMax.getX()) || Double.isNaN(energyRangeMinMax.getY()) ) {
	    throw new IllegalArgumentException("CodedMaskTimeSeries requires attributes: energyRangeMin, energyRangeMax");
	}
	if ( Double.isNaN(maxDistForFullCoding) ) {
	    throw new IllegalArgumentException("CodedMaskTimeSeries requires attribute maxDistForFullCoding");
	}
	if ( effectivePointingDurations == null ) {
	    throw new IllegalArgumentException("CodedMaskTimeSeries requires effective durations of pointings");
	}
	if ( raDecsOfPointings == null ) {
	    throw new IllegalArgumentException("CodedMaskTimeSeries requires RA and Dec of pointings");
	}
	if ( exposuresOnTarget == null ) {
	    throw new IllegalArgumentException("CodedMaskTimeSeries requires effective exposures on target");
	}
	this.setInstrument(instrument);
	setAttributes(targetRaDec, energyRangeMinMax, maxDistForFullCoding);
	setRaDecsOfPointings(raDecsOfPointings);
	setPointingDurations(effectivePointingDurations);
	setExposures(exposuresOnTarget);
	printCodedMaskInfo();
    }

    // private info-printing 
    private void printCodedMaskInfo() {
	logger.info("CodedMaskTimeSeries extends TimeSeries");
	logger.info("RA of target = "+this.targetRA);
	logger.info("Dec of target = "+this.targetDec);
	logger.info("Energy range = "+this.energyRangeMin+" - "+this.energyRangeMax);
	logger.info("Max distance for full coding = "+this.maxDistForFullCoding);
	logger.info("Instrument = "+this.instrument());
	logger.info("Number of non-NaN pointing directions (RA, Dec) = "+this.nNonNaN_raDecs);
	logger.info("Number of non-NaN pointing durations = "+this.nNonNaN_pointingDurations);
	logger.info("  Total = "+this.sumOfEffectivePointingDurations);
	logger.info("  Mean effective duration = "+this.meanEffectivePointingDuration);
	logger.info("  Min = "+this.minEffectivePointingDuration);
	logger.info("  Max = "+this.maxEffectivePointingDuration);
	logger.info("  Mean deviation = "+this.meanDeviationInEffectivePointingDuration);
	logger.info("Dead time duration:");
	logger.info("  Total = "+this.sumOfDeadTimeDurations);
	logger.info("  Mean deadtime duration = "+this.meanDeadTimeDuration);
	logger.info("  Min = "+this.minDeadTimeDuration);
	logger.info("  Max = "+this.maxDeadTimeDuration);
	logger.info("Live time fraction:");
	logger.info("  Mean livetime fraction = "+this.meanLiveTimeFraction);
	logger.info("  Min = "+this.minLiveTimeFraction);
	logger.info("  Max = "+this.maxLiveTimeFraction);
	logger.info("Dead time fraction:");
	logger.info("  Mean deadtime fraction = "+this.meanDeadTimeFraction);
	logger.info("  Min = "+this.minDeadTimeFraction);
	logger.info("  Max = "+this.maxDeadTimeFraction);
	logger.info("Number of non-NaN distances to pointing axis = "+this.nNonNaN_distToPointingAxis);	
	logger.info("  Fully coded points = "+this.nFullyCoded);
	logger.info("  Fully coded fraction = "+this.fullyCodedFraction);
	logger.info("  Mean distance = "+this.meanDistToPointingAxis);
	logger.info("  Min = "+this.minDistToPointingAxis);
	logger.info("  Max = "+this.maxDistToPointingAxis);
	logger.info("  Mean deviation = "+this.meanDeviationInDistToPointingAxis);
	logger.info("Number of non-NaN exposures on target = "+this.nNonNaN_exposuresOnTarget);
	logger.info("  Total = "+this.sumOfExposuresOnTarget);
	logger.info("  Mean exposure = "+this.meanExposureOnTarget);
	logger.info("  Min = "+this.minExposureOnTarget);
	logger.info("  Max = "+this.maxExposureOnTarget);
	logger.info("  Mean deviation = "+this.meanDeviationInExposureOnTarget);
    }

    private void setAttributes(Point2D.Double targetRaDec, Point2D.Double energyRangeMinMax, double maxDistForFullCoding) {
	this.targetRA = targetRaDec.getX();
	this.targetDec = targetRaDec.getY();
	this.energyRangeMin = energyRangeMinMax.getX();
	this.energyRangeMax = energyRangeMinMax.getY();
	this.maxDistForFullCoding = maxDistForFullCoding;
    }
    
    private void setRaDecsOfPointings(Point2D.Double[] raDecsOfPointings) {
	this.raDecsOfPointings = new Point2D.Double[raDecsOfPointings.length];
	this.rasOfPointings = new double[raDecsOfPointings.length];
	this.decsOfPointings = new double[raDecsOfPointings.length];
	this.distToPointingAxis = new double[raDecsOfPointings.length];
	WorldCoords targetCoords = new WorldCoords(this.targetRA, this.targetDec);
	int nNonNaN_raDecs = 0;
	int nNonNaN_distToPointingAxis = 0;
	int nFullyCoded = 0;
	double minRaOfPointings = Double.MAX_VALUE;
	double maxRaOfPointings = -Double.MAX_VALUE;
	double minDecOfPointings = Double.MAX_VALUE;
	double maxDecOfPointings = -Double.MAX_VALUE;
	double minDistToPointingAxis = Double.MAX_VALUE;
	double maxDistToPointingAxis = -Double.MAX_VALUE;
	double[] rates = this.getRates();
	double sumOfDist = 0;
	double sum2OfDist = 0;
	for ( int i=0; i < raDecsOfPointings.length; i++ ) {
	    this.raDecsOfPointings[i] = raDecsOfPointings[i];
	    this.rasOfPointings[i] = raDecsOfPointings[i].getX();
	    this.decsOfPointings[i] = raDecsOfPointings[i].getY();
	    if ( Double.isNaN(this.rasOfPointings[i]) || Double.isNaN(this.decsOfPointings[i]) ) {
		if ( ! Double.isNaN(rates[i]) ) {
		    logger.warn("There is a NaN value in RA or Dec whose corresponding rate is not NaN.");
		}
	    }
	    else {
		nNonNaN_raDecs++;
		minRaOfPointings = Math.min(minRaOfPointings, this.rasOfPointings[i]);
		maxRaOfPointings = Math.max(maxRaOfPointings, this.rasOfPointings[i]);
		minDecOfPointings = Math.min(minDecOfPointings, this.decsOfPointings[i]);
		maxDecOfPointings = Math.max(maxDecOfPointings, this.decsOfPointings[i]);
	    }
	    WorldCoords pointingCoords = new WorldCoords(this.rasOfPointings[i], this.decsOfPointings[i]);
	    double dist = targetCoords.dist(pointingCoords); // returns arc minutes
	    double distInDeg = dist/60.;
	    this.distToPointingAxis[i] = distInDeg;
	    if ( Double.isNaN(this.distToPointingAxis[i]) ) {
		if ( ! Double.isNaN(rates[i]) ) {
		    logger.warn("There is a NaN value in distance from target to pointing axis whose corresponding rate is not NaN.");
		}
	    }
	    else {
		nNonNaN_distToPointingAxis++;
		minDistToPointingAxis = Math.min(minDistToPointingAxis, distInDeg);
		maxDistToPointingAxis = Math.max(maxDistToPointingAxis, distInDeg);
		sumOfDist += distInDeg;
		sum2OfDist += distInDeg*distInDeg;
		if ( distInDeg <= this.maxDistForFullCoding ) {
		    nFullyCoded++;
		}
	    }
	}
	this.nNonNaN_raDecs = nNonNaN_raDecs;
	this.nNonNaN_distToPointingAxis = nNonNaN_distToPointingAxis;
	this.nFullyCoded = nFullyCoded;
	this.fullyCodedFraction = (double)this.nFullyCoded/(double)this.nNonNaN_raDecs;
	this.minRaOfPointings = minRaOfPointings;
	this.maxRaOfPointings = maxRaOfPointings;
	this.minDecOfPointings = minDecOfPointings;
	this.maxDecOfPointings = maxDecOfPointings;
	this.minDistToPointingAxis = minDistToPointingAxis;
	this.maxDistToPointingAxis = maxDistToPointingAxis;
	this.sumOfDistToPointingAxis = sumOfDist;
	this.sumOfSquaredDistToPointingAxis = sum2OfDist;
	//setStatsOnRaDecsOfPointings();
	setStatsOnDistToPointingAxis();
    }

    // private void setStatsOnRaDecsOfPointings() {
    // 	//  Calculate means and variances RA and Dec for collection of pointings
    // double meanRaOfPointingsDecOfPointings = ;
    // double meanDeviationInRasOfPointings = ;
    // double meanDeviationInDecsOfPointings = ;    
    // double varianceInRasOfPointings = ;
    // double varianceInDecsOfPointings = ;
    // }

    private void setStatsOnDistToPointingAxis() {
	this.meanDistToPointingAxis = this.sumOfDistToPointingAxis/this.nNonNaN_distToPointingAxis;
 	this.varianceInDistToPointingAxis = Descriptive.sampleVariance(this.nNonNaN_distToPointingAxis, this.sumOfDistToPointingAxis, this.sumOfSquaredDistToPointingAxis);
	this.meanDeviationInDistToPointingAxis = Descriptive.meanDeviation(new DoubleArrayList(this.distToPointingAxis), this.meanDistToPointingAxis);
    }

    private void setPointingDurations(double[] effectivePointingDurations) {
	this.effectivePointingDurations = new double[effectivePointingDurations.length];
	this.deadTimeDurations = new double[effectivePointingDurations.length];
	this.liveTimeFractions = new double[effectivePointingDurations.length];
	this.deadTimeFractions = new double[effectivePointingDurations.length];
	double[] rates = this.getRates();
	double[] binWidths = this.getBinWidths();
	double min = Double.MAX_VALUE;
	double max = -Double.MAX_VALUE;
	double sum = 0;
	double sum2 = 0;
	double minDeadTime = Double.MAX_VALUE;
	double maxDeadTime = -Double.MAX_VALUE;
	double sumDeadTime = 0;
	double minLiveFraction = Double.MAX_VALUE;
	double maxLiveFraction = -Double.MAX_VALUE;
	double sumLiveFraction = 0;
	double minDeadFraction = Double.MAX_VALUE;
	double maxDeadFraction = -Double.MAX_VALUE;
	double sumDeadFraction = 0;
	int n = 0;
	for ( int i=0; i < this.nBins(); i++ ) {
	    this.effectivePointingDurations[i] = effectivePointingDurations[i];
	    if ( Double.isNaN(effectivePointingDurations[i]) ) {
		if ( ! Double.isNaN(rates[i]) ) {
		    logger.warn("There is a NaN value in RA or Dec whose corresponding rate is not NaN.");
		}
		this.deadTimeDurations[i] = Double.NaN;
		this.liveTimeFractions[i] = Double.NaN;
		this.deadTimeFractions[i] = Double.NaN;
	    }
	    else {
		n++;
		// effective duration
		sum += effectivePointingDurations[i];
		sum2 += effectivePointingDurations[i]*effectivePointingDurations[i];		
		min = Math.min(min, effectivePointingDurations[i]);
		max = Math.max(max, effectivePointingDurations[i]);
		// dead time
		double deadTime = binWidths[i] - effectivePointingDurations[i];
		this.deadTimeDurations[i] = deadTime;
		minDeadTime = Math.min(minDeadTime, deadTime);
		maxDeadTime = Math.max(maxDeadTime, deadTime);
		sumDeadTime += deadTime;
		// live time fraction
		double liveTimeFraction = effectivePointingDurations[i]/binWidths[i];
		this.liveTimeFractions[i] = liveTimeFraction;
		minLiveFraction = Math.min(minLiveFraction, liveTimeFraction);
		maxLiveFraction = Math.max(maxLiveFraction, liveTimeFraction);
		sumLiveFraction += liveTimeFraction;
		// dead time fraction
		double deadTimeFraction = 1. - liveTimeFraction; 
		this.deadTimeFractions[i] = deadTimeFraction;
		minDeadFraction = Math.min(minDeadFraction, deadTimeFraction);
		maxDeadFraction = Math.max(maxDeadFraction, deadTimeFraction);
		sumDeadFraction += deadTimeFraction;
	    }
	}
	this.nNonNaN_pointingDurations = n;
	this.sumOfEffectivePointingDurations = sum;
	this.sumOfSquaredEffectivePointingDurations = sum2;
	this.minEffectivePointingDuration = min;
	this.maxEffectivePointingDuration = max;

	this.sumOfDeadTimeDurations = sumDeadTime;
	this.minDeadTimeDuration = minDeadTime;
	this.maxDeadTimeDuration = maxDeadTime;

	this.sumOfLiveTimeFractions = sumLiveFraction;
	this.minLiveTimeFraction = minLiveFraction;
	this.maxLiveTimeFraction = maxLiveFraction;
	
	this.sumOfDeadTimeFractions = sumDeadFraction;
	this.minDeadTimeFraction = minDeadFraction;
	this.maxDeadTimeFraction = maxDeadFraction;
	setStatsOnPointingDurations();
    }

    private void setStatsOnPointingDurations() {
	this.meanEffectivePointingDuration = this.sumOfEffectivePointingDurations/this.nNonNaN_pointingDurations;
	this.varianceInEffectivePointingDuration = Descriptive.sampleVariance(this.nNonNaN_pointingDurations, this.sumOfEffectivePointingDurations, this.sumOfSquaredEffectivePointingDurations);
	this.meanDeviationInEffectivePointingDuration = Descriptive.meanDeviation(new DoubleArrayList(this.effectivePointingDurations), this.meanEffectivePointingDuration);	
	this.meanDeadTimeDuration = this.sumOfDeadTimeDurations/this.nNonNaN_pointingDurations;
	this.meanLiveTimeFraction = this.sumOfLiveTimeFractions/this.nNonNaN_pointingDurations;
	this.meanDeadTimeFraction = this.sumOfDeadTimeFractions/this.nNonNaN_pointingDurations;
    }
    
    private void setExposures(double[] exposuresOnTarget) {
	this.exposuresOnTarget = new double[exposuresOnTarget.length];
	double min = Double.MAX_VALUE;
	double max = -Double.MAX_VALUE;
	double sum = 0;
	double sum2 = 0;
	int n = 0;
	double[] rates = this.getRates();
	for ( int i=0; i < exposuresOnTarget.length; i++ ) {
	    this.exposuresOnTarget[i] = exposuresOnTarget[i];
	    if ( Double.isNaN(exposuresOnTarget[i]) ) {
		if ( ! Double.isNaN(rates[i]) ) {
		    logger.warn("There is a NaN value in effective exposures whose corresponding rate is not NaN.");
		}
	    }
	    else {
		min = Math.min(min, exposuresOnTarget[i]);
		max = Math.max(max, exposuresOnTarget[i]);
		sum += exposuresOnTarget[i];
		sum2 += exposuresOnTarget[i]*exposuresOnTarget[i];
		n++;
	    }
	}
	this.nNonNaN_exposuresOnTarget = n;
	this.minExposureOnTarget = min;
	this.maxExposureOnTarget = max;
	this.sumOfExposuresOnTarget = sum;
	this.sumOfSquaredExposuresOnTarget = sum2;
	setStatsOnExposures();
    }
    
    private void setStatsOnExposures() {
	this.meanExposureOnTarget = this.sumOfExposuresOnTarget/this.nNonNaN_exposuresOnTarget;
	this.varianceInExposureOnTarget = Descriptive.sampleVariance(this.nNonNaN_exposuresOnTarget, this.sumOfExposuresOnTarget, this.sumOfSquaredExposuresOnTarget);
	this.meanDeviationInExposureOnTarget = Descriptive.meanDeviation(new DoubleArrayList(this.exposuresOnTarget), this.meanExposureOnTarget);	
    }

    //  About attributes
    public double targetRA() { return this.targetRA; }
    public double targetDec() { return this.targetDec; }
    public double energyRangeMin() { return this.energyRangeMin; }
    public double energyRangeMax() { return this.energyRangeMax; }
    public double maxDistForFullCoding() { return this.maxDistForFullCoding; }
    public int nFullyCoded() { return this.nFullyCoded; }
    public double fullyCodedFraction() { return this.fullyCodedFraction; }

    
    // About RA and Dec
    public Point2D.Double[] getRaDecsOfPointings() { return Arrays.copyOf(this.raDecsOfPointings, this.raDecsOfPointings.length); }
    public int nNonNaN_raDecs() { return this.nNonNaN_raDecs; }
    public double[] getRasOfPointings() { return Arrays.copyOf(this.rasOfPointings, this.rasOfPointings.length); }
    public double[] getDecsOfPointings() { return Arrays.copyOf(this.decsOfPointings, this.decsOfPointings.length); }
    public double minRaOfPointings() { return this.minRaOfPointings; };
    public double maxRaOfPointings() { return this.maxRaOfPointings; };
    public double minDecOfPointings() { return this.minDecOfPointings; };
    public double maxDecOfPointings() { return this.maxDecOfPointings; };
    // public double meanRaDecOfPointings() { return this.meanRaDecOfPointings; }
    // public double meanDeviationInRasOfPointings(); { return this.meanDeviationInRasOfPointings(); }
    // public double meanDeviationInDecsOfPointings(); { return this.meanDeviationInDecsOfPointings(); }
    // public double varianceInRasOfPointings() { return this.varianceInRasOfPointings(); }
    // public double varianceInDecsOfPointings() { retutn this.varianceInDecsOfPointings();}

    //  About pointing durations
    public double[] getEffectivePointingDurations() { return Arrays.copyOf(this.effectivePointingDurations, this.effectivePointingDurations.length); }
    public double minEffectivePointingDuration() { return this.minEffectivePointingDuration; }
    public double maxEffectivePointingDuration() { return this.maxEffectivePointingDuration; }
    public double meanEffectivePointingDuration() { return this.meanEffectivePointingDuration; }
    public double sumOfEffectivePointingDurations() { return this.sumOfEffectivePointingDurations; }
    public double varianceInEffectivePointingDuration() { return this.varianceInEffectivePointingDuration; }
    public double meanDeviationInEffectivePointingDuration() { return this.meanDeviationInEffectivePointingDuration; }    
    public double[] getDeadTimeDurations() { return Arrays.copyOf(this.deadTimeDurations, this.deadTimeDurations.length); }
    public double minDeadTimeDuration() { return this.minDeadTimeDuration; }
    public double maxDeadTimeDuration() { return this.maxDeadTimeDuration; }
    public double meanDeadTimeDuration() { return this.meanDeadTimeDuration; }
    public double sumOfDeadTimeDurations() { return this.sumOfDeadTimeDurations; }
    public double[] getLiveTimeFractions() { return Arrays.copyOf(this.liveTimeFractions, this.liveTimeFractions.length); }
    public double minLiveTimeFraction() { return this.minLiveTimeFraction; }
    public double maxLiveTimeFraction() { return this.maxLiveTimeFraction; }
    public double meanLiveTimeFraction() { return this.meanLiveTimeFraction; }
    public double sumOfLiveTimeFractions() { return this.sumOfLiveTimeFractions; }
    public double[] getDeadTimeFractions() { return Arrays.copyOf(this.deadTimeFractions, this.deadTimeFractions.length); }    
    public double minDeadTimeFraction() { return this.minDeadTimeFraction; }
    public double maxDeadTimeFraction() { return this.maxDeadTimeFraction; }
    public double meanDeadTimeFraction() { return this.meanDeadTimeFraction; }
    public double sumOfDeadTimeFractions() { return this.sumOfDeadTimeFractions; }
    
    // About distance from pointing axis
    public double[] getDistToPointingAxis() { return Arrays.copyOf(this.distToPointingAxis, this.distToPointingAxis.length); }
    public int nNonNaN_distToPointingAxis() { return this.nNonNaN_distToPointingAxis; }
    public double minDistToPointingAxis() { return this.minDistToPointingAxis; }
    public double maxDistToPointingAxis() { return this.maxDistToPointingAxis; }
    public double meanDistToPointingAxis() { return this.meanDistToPointingAxis; }
    public double varianceInDistToPointingAxis() { return this.varianceInDistToPointingAxis; }
    public double meanDeviationInDistToPointingAxis() { return this.meanDeviationInDistToPointingAxis; }
    public double sumOfDistToPointingAxis() { return this.sumOfDistToPointingAxis; }
    public double sumOfSquaredDistToPointingAxis() { return this.sumOfSquaredDistToPointingAxis; }

    // About effective exposures
    public double[] getExposuresOnTarget() { return Arrays.copyOf(this.exposuresOnTarget, this.exposuresOnTarget.length); }
    public int nNonNaN_exposuresOnTarget() { return this.nNonNaN_exposuresOnTarget; }
    public double minExposureOnTarget() { return this.minExposureOnTarget; }
    public double maxExposureOnTarget() { return this.maxExposureOnTarget; }
    public double meanExposureOnTarget() { return this.meanExposureOnTarget; }
    public double varianceInExposureOnTarget() { return this.varianceInExposureOnTarget; }
    public double meanDeviationInExposureOnTarget() { return this.meanDeviationInExposureOnTarget; }
    public double sumOfExposuresOnTarget() { return this.sumOfExposuresOnTarget; }
    public double sumOfSquaredExposuresOnTarget() { return this.sumOfSquaredExposuresOnTarget; }

    //  Write to ascii file
    public void writeRatesAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAsQDP(this, filename);
    }

    public void writeRatesAndDistToAxisAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAndDistToAxisAsQDP(this, filename);
    }

    public void writeRatesAsPLT(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAsQDP(this, filename);
    }
    
    public void writeRatesAsXML(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAsQDP(this, filename);
    }

}
