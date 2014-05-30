/*
 * @(#)Chromatogram.java	1.0	2000/10/24
 *
 * Copyright 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc. 
 * Use is subject to license terms.
 *
 */

package com.paracel.tt.chrom;

/**
 * The chromatogram interface for various types of chromatograms, such as
 * ABIChromatogram, SCFChromatogram, etc.
 */
public interface Chromatogram {
    /**
     * Gets the length of the chromatogram, which is the max length among
     * the traces
     * @return an int specifying the length of the chromatogram
     */
    public int getLength();

    /**
     * Gets total number of analyzed traces
     * @return an int specifying the total number of traces
     */
    public int getNumOfAnalyzedTraces();

    /**
     * Gets the maximum intensity among all the traces
     * @return an int specifying the maximum intensity among all the traces
     */
    public int getMaxIntensity();

    /**
     * Gets the analyzed traces
     * @return  A <code>Trace</code> array that contains all analyzed traces.
     */
    public Trace[] getAnalyzedTraces();

    /**
     * Gets the analyzed trace at the specified index
     * @return a <code>Trace</code> object at the specified index
     * @exception IndexOutOfBoundsException if index value < 0 
     *            or >= the total number of traces
     * @see       java.lang.IndexOutOfBoundsException
     */
    public Trace getAnalyzedTraceAt(int index) throws IndexOutOfBoundsException;
}
