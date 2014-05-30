/*
 * 1.3
 *
 * @(#)ABIChromatogram.java       1.0     2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.chrom;

import java.awt.*;
import java.lang.*;
import java.io.*;
import com.paracel.tt.io.*;
import com.paracel.tt.util.*;

/**
 * An implementation of an <code>AbstractChromatogram</code>.  This is used
 * for chromatograms generated from ABI files.
 */
public class ABIChromatogram extends AbstractChromatogram {

    /** The raw data traces array */
    private Trace[] rawDataTraces;

    /** The total number of raw data traces */
    private int numOfRawDataTraces;

    private int rawDataMaxIntensity;

    private int rawDataTraceLength;

    /** The called peak locations array */
    private short[] calledPeakLocations;

    /** The edited base call peak locations array */
    private short[] editedPeakLocations;

    /** The called bases */
    private String calledBases;

    /** The edited bases */
    private String editedBases;

    /** The fifth dye data trace */
    private Trace fifthDyeTrace;

    /** The flag to indicate whether fifth dye data exists */
    private boolean fifthDyeTraceExist;

    /** Creates an ABIChromatogram object */
    protected ABIChromatogram() {}

    /**
     * Creates an ABIChromatogram object from the specified ABI file name
     * @param   abiFileName	The name of the ABI file.
     * @exception FileNotFoundException
     * @exception IOException
     * @see     java.io.FileNotFoundException
     * @see     java.io.IOException
     */
    public ABIChromatogram(String abiFileName)
		throws FileNotFoundException, IOException {
	ABIFile abiFile = new ABIFile(abiFileName);
	initialize(abiFile);
    }

    /** Creates an ABIChromatogram object from the specified ABIFile object */
    public ABIChromatogram(ABIFile abiFile) {
	initialize(abiFile);
    }

    /** Creates an ABIChromatogram object from the specified ABIFile object */
    protected void initialize(ABIFile abiFile) {
	numOfAnalyzedTraces = 4; //ACGT
	calledBases = abiFile.getCalledBases().toUpperCase();
	editedBases = abiFile.getEditedBases().toUpperCase();
	calledPeakLocations = abiFile.getCalledPeakLocations();
	editedPeakLocations = abiFile.getEditedPeakLocations();

	// Get fifth dye data 
	int[] fifthDyeData;

	fifthDyeData = abiFile.getFifthDyeData();
	if (fifthDyeData.length > 0) {
	    fifthDyeTraceExist = true;
	    fifthDyeTrace = new Trace("Fifth Dye",
				      Trace.findTraceColor("Fifth Dye"),
				      fifthDyeData);
	    numOfRawDataTraces = 5; //ACGT5
	} else {
	    fifthDyeTraceExist = false;
	    numOfRawDataTraces = 4; //ACGT
	}

	// Get analyzed data
	analyzedTraces = new Trace[numOfAnalyzedTraces];
	String name;
	for (int i = 0; i < numOfAnalyzedTraces; i++) {
	    name = abiFile.getBaseFromDyeIndex(i+1);
	    analyzedTraces[i] = new Trace(name,
					  Trace.findTraceColor(name),
					  abiFile.getAnalyzedData(i+1));
	} 

	// Get raw data
	rawDataTraces = new Trace[numOfRawDataTraces];
	rawDataTraceLength = 0;
	for (int i = 0; i < numOfRawDataTraces; i++) {
	    if (i == 4) {
		rawDataTraces[i] = fifthDyeTrace;
	    } else {
	        name = abiFile.getBaseFromDyeIndex(i+1);
	        rawDataTraces[i] = new Trace(name,
					 Trace.findTraceColor(name),
					 abiFile.getRawData(i+1));
	    }
	    int length = rawDataTraces[i].getLength();
	    if (length > rawDataTraceLength) {
		rawDataTraceLength = length;
	    }

	    int intensity = rawDataTraces[i].getMaxIntensity();
	    if (intensity > rawDataMaxIntensity) {
		rawDataMaxIntensity = intensity;
	    }
	} 

    }

    /**
     * Identifies whether or not fifth dye data exists
     * @return  true if fifth dye data exists; false otherwise.
     */
    public boolean hasFifthDyeData() {
	return fifthDyeTraceExist;
    }

    /**
     * Gets the called peak locations
     * @return	An short int array that contains all called peak locations.
     */
    public short[] getCalledPeakLocations() {
	return calledPeakLocations;
    }

    /**
     * Gets the edited peak locations
     * @return	An short int array that contains all edited peak locations.
     */
    public short[] getEditedPeakLocations() {
	return editedPeakLocations;
    }

    /**
     * Gets the edited bases
     * @return	A string that represents all edited bases.
     */
    public String getEditedBases() {
	return editedBases;
    }

    /**
     * Gets the fifth dye trace
     * @return	The fifth dye trace.
     */
    public Trace getFifthDyeTrace() {
	return fifthDyeTrace;
    }

    /**
     * Gets the called bases
     * @return	A string that represents all called bases.
     */
    public String getCalledBases() {
	return calledBases;
    }

    /**
     * Gets the raw data traces
     * @return	A <code>Trace</code> array that contains all raw data traces.
     */
    public Trace[] getRawDataTraces() {
	return rawDataTraces;
    }

    /**
     * Gets the raw data trace at the specified index
     * @return	The raw data trace at the specified index.
     * @param	The specified index.
     * @exception  IndexOutOfBoundsException if index < 0 or >= total number
     *		   of raw data traces.
     * @see	java.lang.IndexOutOfBoundsException
     */
    public Trace getRawDataTraceAt(int index) {
        if ( (index >= length) || (index < 0) ) {
            throw new IndexOutOfBoundsException(
				"Wrong raw data trace array index");
	} else {
	    return rawDataTraces[index];
	}
    }

    /**
     * Gets the length of the raw data trace.
     * @return	the length of the raw data trace.
     */
    public int getRawDataTraceLength() {
	return rawDataTraceLength;
    }

    /**
     * Gets the maximum intensity of the raw data trace.
     * @return	the maximum intensity of the raw data trace.
     */
    public int getRawDataMaxIntensity() {
	return rawDataMaxIntensity;
    }
}
