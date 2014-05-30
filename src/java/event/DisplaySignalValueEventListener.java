/*
 * 1.1
 *
 * @(#)DisplaySignalValueEventListener.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.event;

import java.awt.*;
import java.awt.event.*;
import java.util.EventListener;

/**
 * The listener interface for receiving notification of 
 * DisplaySignalValueEvent.
 */
public interface DisplaySignalValueEventListener extends EventListener {
   
    /** Invoked when the user single-click on the chromatogram. */
    public void signalValueSet(DisplaySignalValueEvent event);
}
