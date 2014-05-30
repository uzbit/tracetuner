/*
 * 1.2
 *
 * @(#)AbstractChromatogram.java       1.0     2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.chrom;

import java.awt.*;
import java.io.*;
import java.lang.*;
import com.paracel.tt.io.*;
import com.paracel.tt.util.*;

/**
 * Defines common behaviors for chromatograms.
 */
public class AbstractChromatogram implements Chromatogram {

    public static final int ABI = 1;
    public static final int SCF = 2;

    /** The total number of traces */
    protected int numOfAnalyzedTraces;

    /** The maximum intensity among all the traces */
    protected int maxIntensity;

    /** The length of the chromatogram */
    protected int length;

    /** The array of all analyzed traces */
    protected Trace[] analyzedTraces;

    /** Returns the file type of the specified file. */
    public static int getChromType(String chromFileName)
		throws FileNotFoundException, IOException {
	ABIFile file = new ABIFile(chromFileName);
	return (file.getType());
    }

    /** Returns the file type of the specified ABIFile object. */
    public static int getChromType(ABIFile file) {
	return (file.getType());
    }

    /**
     * Gets the length of the chromatogram, which is the max length among
     * the analyzed traces
     * @return an int specifying the length of the chromatogram
     */
    public int getLength() {
	// iterate through all the traces in this chromatogram, find the
	// maximum trace length.
	int max = 0;
	int length;

	for (int i = 0; i < numOfAnalyzedTraces; i++) {
	    length = analyzedTraces[i].getLength();
	    if (length > max) {
		max = length;
	    }
	}

	this.length = max;
	return max;
    }

    /**
     * Gets total number of analyzed traces
     * @return an int specifying the total number of analyzed traces
     */
    public int getNumOfAnalyzedTraces() { 
	return numOfAnalyzedTraces; 
    }

    /**
     * Gets the maximum intensity among all the analyzed traces
     * @return an int specifying the maximum intensity among all the 
     *         analyzed traces
     */
    public int getMaxIntensity() {
        // iterate through all the analyzed traces in this chromatogram,
	// find the maximum intensity.
        maxIntensity = 0;
	int intensity;

        for (int i = 0; i < numOfAnalyzedTraces; i++) {

	    intensity = analyzedTraces[i].getMaxIntensity();
            if (intensity > maxIntensity) {
                maxIntensity = intensity;
            }
        }

        return maxIntensity;
    }

    /**
     * Gets the analyzed traces
     * @return  A <code>Trace</code> array that contains all analyzed traces.
     */
    public Trace[] getAnalyzedTraces() { 
	return analyzedTraces;
    }

    /**
     * Gets the trace at the specified index
     * @return a <code>Trace</code> object at the specified index
     * @exception IndexOutOfBoundsException if index value < 0
     *            or >= the total number of traces
     * @see       java.lang.IndexOutOfBoundsException
     */
    public Trace getAnalyzedTraceAt(int index)
		throws IndexOutOfBoundsException {
	if (index >= numOfAnalyzedTraces) {
	    throw new IndexOutOfBoundsException(
				"Wrong analyzed trace array index");
	} else {
	    return analyzedTraces[index];
	}
    }
}
