/*
 * 1.2
 *
 * @(#)Trace.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.chrom;

import java.awt.*;
import java.lang.*;
import com.paracel.tt.util.*;

/**
 * This class creates a Trace object, which contains the name, color, data
 * points of the trace.
 */
public class Trace implements Cloneable {

    /** The name of the trace.  */
    private String name;
   
    /** The color of the trace.  */
    private Color color;
   
    /** The intensity value array of the trace. */
    private int[] intensity;

    /** The length of the trace. */
    private int length;

    /**
     * Creates a Trace.
     * @param     name The trace name
     * @param     color The trace color
     * @param     intensity The trace intensity value array
     */
    public Trace (String name, Color color, int[] intensity) {
	this.name = name;
	this.color = color;

	this.intensity = intensity;

	this.length = intensity.length;
    }

    /**
     * Creates a clone of the specified trace.
     * @return    a clone of the specified trace
     * @param     trace The specified trace to be cloned
     */
    public Object clone() {
	return new Trace(this.name, this.color, this.intensity);
    }

    /**
     * Releases the memory resources the specified Trace object holds.
     * @param	trace	The specified trace to be disposed
     */
    public void dispose() {
	this.name = null;
	this.color = null;
	this.intensity = null;
    }

    /**
     * Finds the trace color representing the trace with the specified name.
     * @return	  the color corresponding to the specified trace name
     */
   static public Color findTraceColor(String name) {
	if (name == null) {
	    return null;
	} else if (name.equalsIgnoreCase("A")) {
	    return Constants.green;
	} else if (name.equalsIgnoreCase("C")) {
	    return Color.blue;
	} else if (name.equalsIgnoreCase("G")) {
	    return Color.black;
	} else if (name.equalsIgnoreCase("T")) {
	    return Color.red;
	} else if (name.equalsIgnoreCase("Fifth Dye")) {
	    return Color.lightGray;
	} else {
            System.err.println("No color defined for trace under the name: "
                                + name);
            return null;
	}
   }

    /**
     * Gets the trace color.
     * @return    the <code>Color</code> value of the trace color
     */
    public Color getColor() {
	return this.color;
    }

    /**
     * Gets the trace name.
     * @return    the <code>String</code> value of the trace name
     */
    public String getName() {
	return this.name;
    }

    /**
     * Gets the intensity array.
     * @return    the array of the intensity values.
     */    
    public int[] getIntensityArray() {
	return intensity;
    }

    /**
     * Gets the intensity value at the specified index.
     * @param     index  The index position
     * @return    the int value of the intensity
     * @exception IndexOutOfBoundsException if index value < 0 or
     *  	  >= the trace length
     * @see	  java.lang.IndexOutOfBoundsException
     */
    public int getIntensityAt(int index) throws IndexOutOfBoundsException {
	if ( (index >= length) || (index < 0) ) {
            throw new IndexOutOfBoundsException("Wrong trace index: " + index);
        } else {
	    return this.intensity[index];
	}
    }

    /**
     * Gets the length of the trace.
     * @return    the int value of the trace length
     */
    public int getLength() {
	return this.length;
    }

    /**
     * Gets the maximum intensity value.
     * @return    the int value of the maximum intensity
     */
    public int getMaxIntensity() {
	int max = 0;
	for (int i = 0; i < length; i++) {
	    if (intensity[i] > max) {
		max = intensity[i];
	    }
	}
	return max;
    }

    /**
     * Sets the intensity value at the specified index to the specified value.
     */
    public void setIntensityAt(int index, int value) {
	intensity[index] = value;
    }
}
