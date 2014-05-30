/**
 * 1.1
 *
 * @(#)HelpAction.java	
 * 
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 *
 */

package com.paracel.tt.util;

import java.io.*;
//import com.apple.mrj.MRJFileUtils;
//import com.paracel.shared.util.*;

/** This class lets you open a URL in your default browser for whatever platform you're on.
 *  Its purpose is to encapsulate all machine-specific browser knowledge.
 *
 * <p>
 * <b>Design Pattern:</b>  Adapter
 *
 * @author Michael Gibson
 * @version 1.1
 */
public class DefaultBrowser {

    public static String OS_NAME = System.getProperty("os.name");

    /** Whether the program is run on a windows OS. */
    public static boolean isWindows() {
        return OS_NAME.trim().toLowerCase().startsWith("windows");
    }

    /** Whether the program is run on a Linux OS. */
    public static boolean isLinux() {
        return OS_NAME.trim().toLowerCase().startsWith("linux");
    }

    /** Whether the program is run on a Sun OS. */
    public static boolean isSunOS() {
        return OS_NAME.trim().toLowerCase().startsWith("sunos");
    }

    /** Whether the program is run on Mac OS. */
    public static boolean isMacOS() {
        // Apple-recommended method for detecting Mac OS in Java.
        return System.getProperty("mrj.version") != null;
    }


    /** This is a temporary kludge, because 1.3 doesn't contain regular expressions.
     *  when we switch to 1.4, this can all be replaced by:
     *  return Pattern.matches("^.*://.*$", string);
     */
    private static boolean containsProtocol(String string) {
	int colonIndex = string.indexOf(':');
	if (colonIndex == -1) {
	    return false;
	}

	return string.regionMatches(false, colonIndex, new String("://"), 0, 3);
    }

    /// Open default web browser to show the content of the url
    public static void openURL(String urlString) throws java.io.IOException {
        //if (Utilities.isMacOS()) {
            //MRJFileUtils.openURL(urlString);
	//} else if (Utilities.isWindows()) {
	if (isWindows()) {
	    windowsOpenURL(urlString);
	} else if ( ! isMacOS() ) {
	    String[] commandArray = new String[] {
					"netscape", 
					urlString };
	    Process process = Runtime.getRuntime().exec(commandArray);
	}
    }

    private static void windowsOpenURL(String urlString)
	    throws java.io.IOException {

	// Sometimes, you can open a url in Windows by using
	//     rundll32 url.dll,FileProtocolHandler <urlString>
	// However, this doesn't work for
	//	urlStrings with spaces
	//	urlStrings that end in .html or .htm (go figure!)
	// Anyhow, the code below works for those cases
	File shortcut = File.createTempFile("InternetShortCut", ".url");
	shortcut = shortcut.getCanonicalFile();
	shortcut.deleteOnExit();
	PrintWriter output = new PrintWriter(new FileWriter(shortcut));
	output.println("[InternetShortcut]");
	output.println("URL=" + urlString);
	output.close();

	String cmd = "rundll32 url.dll,FileProtocolHandler " +
			shortcut.getCanonicalPath();
	Process process = Runtime.getRuntime().exec(cmd);

	// FUTURE: you can also open the URL by just executing the command
	// 	name_of_url_file
	// On my machine, running IE, that works fine.
	// On William's machine, running Netscape, both the rundll version
	// and the file only version give error messages, but open the
	// URL anyway.  So I'm not going to change the code at this time.
	// Some day, it might be good to figure out the Netscape problem,
	// which may involve getting rid of rundll.
    }

    /// Check the URL, then open default web browser to show the content of the URL
    public static void checkAndOpenURL(String urlString) throws java.io.IOException {
	// Check whether urlString starts with http:// or something like that
	if (!containsProtocol(urlString)) {
	    // If not, add http://
	    urlString = "http://" + urlString;
	}

	openURL(urlString);
    }
}
