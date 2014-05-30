/*
 * 1.2
 *
 * @(#)ZoomEvent.java	1.0	2000/10/24
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
 * This class creates a ZoomEvent object. Such an event is
 * fired when the user changes the zoom factor on either x or y axis.
 */
public class ZoomEvent extends AWTEvent {

    /** The axis. (X or Y) */
    private char axis;

    /** The old zoom factor value. */
    private float oldValue;

    /** The new zoom factor value. */
    private float newValue;

    /**
     * Creates a ZoomEvent object.
     * @param	o	the object where the event originated
     * @param	axis	the axis.
     * @param	oldValue	the old zoom factor.
     * @param	newValue	the new zoom factor.
     */
    public ZoomEvent(Object o, char axis, float oldValue, float newValue) {
	super(o, AWTEvent.RESERVED_ID_MAX + 1);
	this.axis = axis;
	this.oldValue = oldValue;
	this.newValue = newValue;
    }

    /**
     * Gets the axis of which the zoom factor was changed. 
     * @return	axis	the axis.
     */
    public char getAxis() { return axis; }

    /**
     * Gets the old zoom factor.
     * @return	oldValue	the old zoom factor.
     */
    public float getOldValue() { return oldValue; }

    /**
     * Gets the new zoom factor.
     * @return	newValue	the new zoom factor.
     */
    public float getNewValue() { return newValue; }
}
