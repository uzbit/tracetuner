/*
 * 1.1
 *
 * @(#)SearchEventListener.java	1.0	2001/03/21
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.event;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
 * The listener interface for receiving notification of SearchEvent.
 */
public interface SearchEventListener extends EventListener {

    /** Invoked when the search event fired. */
    public void search(SearchEvent event);
}
