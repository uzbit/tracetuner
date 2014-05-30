/*
 * 1.3
 *
 * @(#)BaseColorIdentifier.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.util;

import java.awt.*;
import java.lang.*;

/** This class is used to assign a Color to each base code.  */
public class BaseColorIdentifier {

    /**
     * Gets the corresponding color of the specified base code.
     * @return  the correspoding color.
     * @param	base	The nucleic acid code.
     */
    public static Color getColor(char base) {
	switch (Character.toUpperCase(base)) {
	case 'A':
	    return Constants.green;
	case 'C':
	    return Color.blue;
	case 'G':
	    return Color.black;
	case 'T':
	    return Color.red;
	case 'N':
	case 'Y':
	case 'R':
	case 'K':
	case 'M':
	case 'S':
	case 'W':
	default:
	    return Constants.gray148;
	}
    }

}
