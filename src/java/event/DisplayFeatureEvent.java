/*
 * 1.3
 *
 * @(#)DisplayFeatureEvent.java	1.0	2000/10/24
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
 * This class creates a DisplayFeatureEvent object. A DisplayFeatureEvent is
 * fired when the user changes the display feature in the TraceTuner
 * viewer. For an example: the user clicks on "Display Original Bases" 
 * to view the original base calls.
 */
public class DisplayFeatureEvent extends AWTEvent {

    /** Whether or not to display raw data */
    private boolean displayRaw;

    /** Whether or not to display original base calls */
    private boolean displayOrig;

    /** Whether or not to display alternative base calls */
    private boolean displayAbc;

    /** Whether or not to display alignment with consensus */
    private boolean displayAln;

    /** Whether or not to display intrinsic peaks */
    private boolean displayTip;

    /** Whether or not to display total intrinsic signal */
    private boolean displayTotalIS;

    /** Whether or not to mark trimmed regions */
    private boolean displayTrim;

    /** Whether or not this display data event was triggered 
	by file index change. */
    private boolean fileIndexChanged;

    /**
     * Creates a DisplayFeatureEvent object.
     * @param	o	the object where the event originated
     */
    public DisplayFeatureEvent(Object o, 
			    boolean displayRaw,
			    boolean displayOrig,
			    boolean displayAbc,
			    boolean displayAln,
			    boolean displayTip,
			    boolean displayTotalIS,
			    boolean displayTrim,
			    boolean fileIndexChanged) {
	super(o, AWTEvent.RESERVED_ID_MAX + 2);
	this.displayRaw = displayRaw;
	this.displayOrig = displayOrig;
	this.displayAbc = displayAbc;
	this.displayAln = displayAln;
	this.displayTip = displayTip;
	this.displayTotalIS = displayTotalIS;
	this.displayTrim = displayTrim;
	this.fileIndexChanged = fileIndexChanged;
    }

    /** Returns true if raw data is being displayed; false, otherwise. */
    public boolean displayRaw() { return displayRaw; }

    /** Returns true if original base calls are being displayed; 
        false, otherwise. */
    public boolean displayOrig() { return displayOrig; }

    /** Returns true if alternative base calls are being displayed; 
        false, otherwise. */
    public boolean displayAbc() { return displayAbc; }

    /** Returns true if alignment with consensus is being displayed; 
        false, otherwise. */
    public boolean displayAln() { return displayAln; }

    /** Returns true if intrinsic peaks are being displayed; 
        false, otherwise. */
    public boolean displayTip() { return displayTip; }

    /** Returns true if total intrinsic signal is being displayed. */
    public boolean displayTotalIS() { return displayTotalIS; }

    /** Returns true if trimming effects are being displayed; 
        false, otherwise. */
    public boolean displayTrim() { return displayTrim; }

    /** Returns true if the event originated from file index change. */
    public boolean fileIndexChanged() { return fileIndexChanged; }
}
