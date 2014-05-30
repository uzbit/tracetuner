/*
 * 1.1
 *
 * @(#)CommandBuilder.java	1.0	2001/03/13
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.util;

import java.util.*;

/** This class builds a command array, which can be passed to
    java.lang.Runtime.exec().  The reason of writing this class
    is to build a command array considering the platform issue
    when the command option parameters contain spaces.  This class
    glues the arguments containing white spaces.  For an
    example, a file or directory name with a space in it. 
    <p>To use this class, one calls append() method to pass the
    arguments one by one.  When it's done, getCommandArray() method
    can be used to return the resulted command array.
    <p><pre>For an example:
	If one wants to run the command: (someApp "1 2" 3)
	Just do the following:
		CommandBuilder cb = new CommandBuilder();
		cb.append("someApp");
		cb.append("1 2");
		cb.append("3");
		String[] cmdArray = cb.getCommandArray();
		Process process = Runtime.getRuntime().exec(cmdArray);
    </pre></p>
 */
public class CommandBuilder {

    /** The command array. */
    private ArrayList cmdArray;

    /** Whether the application is running on windows system. */ 
    private boolean isWindows;

    /** Constructor. */
    public CommandBuilder() {
	isWindows = (Constants.OS_NAME.toLowerCase().trim()
	    		.indexOf("windows") >= 0) ? true : false;
	cmdArray = new ArrayList(20);
    }

    /** Append ONE argument as an entity into the command array.
        NOTE: Passing two or more arguments together will cause the command
	builder to MIS-FUNCTION. */
    public void append(String s) {
	if (s == null || s.trim().length() <= 0) {
	    return;
	}
	s = s.trim();

	/* If the argument ends with a backward quote, it messes things up
	   under Windows OSes when the argument contains empty spaces.  */
	while (s.endsWith("\\")) {
	    s = s.substring(0, s.length() - 1).trim();
	}
	if (s.length()<= 0) {
	    return;
	}
	/* If the argument contains empty space(s) and the application runs
	   on Windows, wrap the argument with double quotes. */
	if (s.indexOf(" ") >= 0 && isWindows) {
	    s = "\"" + s + "\"";
	}
	cmdArray.add(s);
    }

    /** Returns true if the built command array contains the specified
	argument (case sensitive). */
    public boolean contains(String s) {
	if (s == null) {
	    return false;
	}
	s = s.trim();
	if (s.length() <= 0) {
	    return false;
	}
	String temp;
	for (int i = 0; i < cmdArray.size(); i++) {
	    temp = (String) cmdArray.get(i);
	    if (s.equals(temp)) {
		return true;
	    }
	}
	return false;
    }
    
    /** Returns the built command array. */
    public String[] getCommandArray() {
	String[] result = new String[cmdArray.size()];
	for (int i = 0; i < result.length; i++) {
	    result[i] = (String) cmdArray.get(i);
	}
	return result;
    }

    /** Returns the command string which wraps the space-containing arguments
	with double quotes.  If this string is passed to Runtime.exec() and it
	contains double quotes, it will work on Windows OS, but fail on 
	Solaris at least for jdk1.3.0. */
    public String getCommandString() {
	StringBuffer buf = new StringBuffer();
	String temp;
	for (int i = 0; i < cmdArray.size(); i++) {
	    temp = (String) cmdArray.get(i);
	    if (temp.indexOf(" ") >= 0 && (!isWindows)) {
		temp = "\"" + temp + "\"";
	    }
	    buf.append(temp + " ");
	}
	return buf.toString();
    }

    /** Returns true if the command array contains nothing. */
    public boolean isEmpty() {
	return (cmdArray.size() == 0);
    }
}
