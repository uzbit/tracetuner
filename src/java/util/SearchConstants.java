/*
 * 1.1 
 *
 * @(#)SearchConstants.java	2001/03/21
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */
package com.paracel.tt.util;


/** 
 * The interface holds constant data used for searching, such as search
 * subject type, search direction, etc.
 */
public interface SearchConstants {
    /** Search subject type constants. */
    public static final int NONE = -1;
    public static final int TT = 0;
    public static final int REF = 1;
    public static final int ABI = 2;

    /** Search direction constants. */
    public static final int FORWARD = 0;
    public static final int BACKWARD = 1;
}
