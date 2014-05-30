/*
 * 1.1
 *
 * @(#)DisplaySignalValueEvent.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.event;

import java.awt.*;
import java.awt.event.*;

/**
 * This class creates a DisplaySignalValueEvent object. Such an event is
 * fired when the user clicks on the chromatogram to view the signal values.
 */
public class DisplaySignalValueEvent extends AWTEvent {

    /** The int array that contains the signal value related data. */
    private int[] signalValues;

    /** The flag that indicates whether or not fifth dye data exists. */
    private boolean fifthDyeExists;

    /** 
     * Creates a DisplaySignalValueEvent object from the specified data.
     * @param	o	the object where the event originated
     * @param	location	the x location of this signal value set.
     * @param	AValue		the signal value of base 'A'.
     * @param	CValue		the signal value of base 'C'.
     * @param	GValue		the signal value of base 'G'.
     * @param	TValue		the signal value of base 'T'.
     */
    public DisplaySignalValueEvent(Object o, int location, int AValue, 
				   int CValue, int GValue, int TValue) {
	super(o, AWTEvent.RESERVED_ID_MAX + 4);
	signalValues = new int[6];
	signalValues[0] = location;
	signalValues[1] = AValue;
	signalValues[2] = CValue;
	signalValues[3] = GValue;
	signalValues[4] = TValue;
	this.fifthDyeExists = false;
    }

    /** 
     * Creates a DisplaySignalValueEvent object from the specified data.
     * @param	o	the object where the event originated
     * @param	location	the x location of this signal value set.
     * @param	AValue		the signal value of base 'A'.
     * @param	CValue		the signal value of base 'C'.
     * @param	GValue		the signal value of base 'G'.
     * @param	TValue		the signal value of base 'T'.
     * @param	fithDyeValue	the signal value of the fifth dye trace.
     */
    public DisplaySignalValueEvent(Object o, int location, int AValue, 
				   int CValue, int GValue, 
				   int TValue, int fifthDyeValue) {
	this(o, location, AValue, CValue, GValue, TValue);
	signalValues[5] = fifthDyeValue;
	this.fifthDyeExists = true;
    }

    /**
     * Gets the signal value array.
     * @return the signal value array.
     */
    public int[] getSignalValues() { return signalValues; }

    /**
     * Indicates whether or not the fifth dye data exists.
     * @return	true if the fifth dye data exists, false otherwise.
     */
    public boolean fifthDyeDataExists() { return fifthDyeExists; }
}
