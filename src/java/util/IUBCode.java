/*
 * 1.1
 *
 * @(#)IUBCode.java	1.0	2001/05/15
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.util;

import java.lang.*;
import java.util.*;

/** This class is used for recognition of IUB codes.  */
public class IUBCode {

    /**
     * Returns true if the specified base is a valid IUB code; false otherwise.
     */
    public static boolean isValid(char base) {
	switch (Character.toUpperCase(base)) {
	case 'A': case 'B': case 'C': case 'D': case 'G':
	case 'H': case 'K': case 'M': case 'N': case 'R':
	case 'S': case 'T': case 'U': case 'V': case 'W': case 'Y':
	    return true;
	default:
	    return false;
	}
    }

    /**
     * Returns true if all bases in the specified sequence are valid
     * IUB codes; false otherwise.
     */
    public static boolean isValid(String sequence) {
	char c;
	for (int i = 0; i < sequence.length(); i++) {
	    c = sequence.charAt(i);
	    if (Character.isWhitespace(c)) {
		continue;
	    } else if (!isValid(sequence.charAt(i))) {
		return false;
	    }
	}
	return true;
    }
}
