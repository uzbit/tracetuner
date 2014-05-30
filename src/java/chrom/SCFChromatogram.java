/*
 * 1.3
 *
 * @(#)SCFChromatogram.java       1.0     2000/10/25
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
 * for chromatograms generated from SCF files.
 *
 * <p>Comparing to ABI chromatograms, SCF chromatograms do not have raw data.
 *    Also SCF files do not contain data for edited base calls or their peak
 *    locations.
 * </p>
 */
public class SCFChromatogram extends AbstractChromatogram {

    /** The called peak locations array */
    private short[] calledPeakLocations;

    /** The called bases */
    private String calledBases;

    /** Creates an SCFChromatogram object */
    protected SCFChromatogram() {}

    /**
     * Creates an SCFChromatogram object from the specified SCF file name
     * @param   scfFileName	The name of the SCF file.
     * @exception FileNotFoundException
     * @exception IOException
     * @see     java.io.FileNotFoundException
     * @see     java.io.IOException
     */
    public SCFChromatogram(String scfFileName)
		throws FileNotFoundException, IOException {
	ABIFile scfFile = new ABIFile(scfFileName);
	initialize(scfFile);
    }

    /** Creates an ABIChromatogram object from the specified ABIFile object */
    public SCFChromatogram(ABIFile scfFile) {
	initialize(scfFile);
    }

    /** Creates an ABIChromatogram object from the specified ABIFile object */
    protected void initialize(ABIFile scfFile) {
	numOfAnalyzedTraces = 4; //ACGT
	analyzedTraces = new Trace[numOfAnalyzedTraces];
	calledBases = scfFile.getCalledBases().toUpperCase();
	calledPeakLocations = scfFile.getCalledPeakLocations();

	// Get analyzed data
	String name;
	for (int i = 0; i < numOfAnalyzedTraces; i++) {
	    name = scfFile.getBaseFromDyeIndex(i+1);
	    analyzedTraces[i] = new Trace(name,
					  Trace.findTraceColor(name),
					  scfFile.getAnalyzedData(i+1));
	} 
    }

    /**
     * Gets the called peak locations
     * @return	An short int array that contains all called peak locations.
     */
    public short[] getCalledPeakLocations() {
	return calledPeakLocations;
    }

    /**
     * Gets the called bases
     * @return	A string that represents all called bases.
     */
    public String getCalledBases() {
	return calledBases;
    }
}
