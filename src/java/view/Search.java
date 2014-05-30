/*
 * 1.1
 *
 * @(#)Search.java       1.0     2001/03/21
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.lang.*;
import java.util.*;
import com.paracel.tt.util.*;

/** This class searches for a substring in a subject string. */
public class Search implements SearchConstants {

    /** The search result. */
    private SearchResult searchResult;

    /** All patterns to search for. */
    private Vector patterns;

    /** The subject type: (TT, ABI, or REF). */
    private int subject;

    /** The subject string. */
    private String sbjString;

    /** The starting index and the direction of the search. */
    private int startIndex, direction;

    /** Returns the search result of the search performed on the specified
	sbjString, for the specified pattern string (pat), starting from
	the specified startIndex, in the specified direction. */
    public SearchResult search(int subject, String sbjString, String pat,
			     int startIndex, int direction) {
	this.startIndex = startIndex;
	this.direction = direction;
	this.subject = subject;
	this.sbjString = sbjString.toUpperCase();
	patterns = parsePatternString(pat.toUpperCase());

	searchResult = new SearchResult();
	searchResult.patternFound = false;

	SearchResult temp;
	for (int i = 0; i < patterns.size(); i++) {
	    temp = searchForPattern((String)patterns.elementAt(i));
	    if (!temp.patternFound) { continue; }
	    if (!searchResult.patternFound) {
		searchResult = temp;
	    } else {
		if (direction == FORWARD && temp.begin < searchResult.begin) {
		    searchResult = temp;
		} else if (direction == BACKWARD 
			   && temp.begin > searchResult.begin) {
		    searchResult = temp;
		}
	    }
	}
	return searchResult;
    }

    /** Parses the specified pat string.  If the pat string contains delimiters
	such as whitespaces, ','s, and ';'s, the pat string is considered to be
	multiple patterns that the user looks for.  eg. if the user speicifies
	a pattern string of "A TTT", that means the user is looking for "A" or
	"TTT".  This function returns a vector containing all the patterns
	parsed from the specified pat string.*/
    protected Vector parsePatternString(String pat) {
	Vector temp = new Vector();
	StringTokenizer st1 = new StringTokenizer(pat);
	StringTokenizer st2, st3;
	while (st1.hasMoreTokens()) {
	    st2 = new StringTokenizer(st1.nextToken(), ";");
	    while (st2.hasMoreTokens()) {
	    	st3 = new StringTokenizer(st2.nextToken(), ",");
	    	while (st3.hasMoreTokens()) {
	    	    temp.add(st3.nextToken());
		}
	    }
	}
	return temp;
    }

    /** Returns the search result for the specified pattern (parsed individual
	subpattern). */
    protected SearchResult searchForPattern(String pattern) {
    	SearchResult temp = new SearchResult();

	int index;
	if (direction == FORWARD) {
	    index = sbjString.indexOf(pattern, startIndex);
	} else {
	    index = sbjString.lastIndexOf(pattern, startIndex);
	}
	if (index < 0) {
	    temp.patternFound = false; 
	    temp.nextStart = index;
	}  else {
	    temp.patternFound = true;
	    temp.subject = subject;
	    temp.pattern = pattern;
	    temp.begin = index;
	    temp.end = temp.begin + pattern.length() - 1;
	    temp.nextStart = (direction == FORWARD)
				? temp.begin + 1 
				: Math.max(temp.begin - 1, 0);
	    temp.direction = direction;
	}
	return temp;
    }

    /** An inner class fo <code>Search</code> for search results. */
    public class SearchResult {

	/** Constructor. */
	public SearchResult() {}

	/** The search subject type: TT, ABI, or REF. */
	public int subject;

	/** The begin index, end index of the found search result. */
	public int begin, end;

	/** The starting index of the next search. */
	public int nextStart;
	
	/** The search direction: FORWARD or BACKWARD. */
	public int direction;
	
	/** The searched pattern of this search result. */
	public String pattern;

	/** Whether or not the pattern was found. */
	public boolean patternFound;
    }
}
