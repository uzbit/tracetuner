/*
 * 1.1
 *
 * @(#)SearchEvent.java	1.0	2001/03/21
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.event;

import java.awt.*;
import java.awt.event.*;
import com.paracel.tt.util.SearchConstants;

/**
 * This class creates a SearchEvent object. Such an event is fired when
 * the user searches for a certain pattern in the TT base calls, ABI base
 * calls, or reference sequence.
 */
public class SearchEvent extends AWTEvent implements SearchConstants {

    /** The search pattern string (search for what). */
    private String pattern;

    /** The search subject (where to search). */
    private int subject;

    /** The searching direction. */
    private int direction, startIndex;

    /** Creates a SearchEvent object. */
    public SearchEvent(Object o) {
	super(o, AWTEvent.RESERVED_ID_MAX + 5);
    }

    /** Creates a SearchEvent object. */
    public SearchEvent(Object o, int subject, String pattern,
		       int direction) {
	this(o);
	this.subject = subject;
	this.pattern = pattern;
	this.direction = direction;
	this.startIndex = -1;
    }

    /** Creates a SearchEvent object. */
    public SearchEvent(Object o, int subject, String pattern,
		       int direction, int startIndex) {
	this(o, subject, pattern, direction);
	this.startIndex = startIndex;
    }

    /** Returns the subject of the search (where to search). */
    public int getSubject() { return subject; }

    /** Returns the search pattern. */
    public String getPattern() { return pattern; }

    /** Returns the search direction. */
    public int getDirection() { return direction; }

    /** Returns the search start index. */
    public int getStartIndex() { return startIndex; }
}
