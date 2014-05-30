/*
 * 1.1
 *
 * @(#)DisplayFeatureEventListener.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.event;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
 * The listener interface for receiving notification of DisplayFeatureEvent.
 */
public interface DisplayFeatureEventListener extends EventListener {

    /** Invoked when the display features are changed. */
    public void displayFeatureChanged(DisplayFeatureEvent event);
}
