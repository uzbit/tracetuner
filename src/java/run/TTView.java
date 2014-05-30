/*
 * 1.4
 *
 * @(#)TTView.java       1.0     2000/10/25
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.run;

import javax.swing.*;
import com.paracel.tt.view.*;
import com.paracel.tt.util.*;

/**
 * The class is provided for running stand-alone Viewer.  It is just
 * a wrapper class of TTViewFrame.
 */
public class TTView {

    public TTView(String args[]) {
	String chromatFileName = null;
	String phdFileName = null;
	String secondPhdFileName = null;
	int numOfArgs = args.length;
	int realArgStart = 0;

	if (numOfArgs == 0) {
	    System.out.println("Usage: java TTView "
                                + "sample-file [phd-file] [2nd-phd-file]");
	    System.exit(1);
	} else  {
	    if (args[0].equalsIgnoreCase("-dev")) {
		Constants.runAsDev = true;
		realArgStart = 1;
		numOfArgs--;
	    }
	    if (numOfArgs <= 0) {
	    	System.out.println("Usage: java TTView "
                                + "sample-file [phd-file] [2nd-phd-file]");
	    	System.exit(1);
	    } else if (numOfArgs == 1) {
	    	chromatFileName = args[realArgStart];
	    	phdFileName = chromatFileName + ".phd.1";
	    	secondPhdFileName = "";
            } else if (numOfArgs == 2) {
	    	chromatFileName = args[realArgStart];
	    	phdFileName = args[realArgStart + 1];
	    	secondPhdFileName = "";
	    } else if (numOfArgs == 3) {
	    	chromatFileName = args[realArgStart];
	    	phdFileName = args[realArgStart + 1];
	    	secondPhdFileName = args[realArgStart + 2];
	    } else {
	    	System.out.println("Usage: java TTView "
                                + "sample-file [phd-file] [2nd-phd-file]");
	    	System.exit(1);
	    }
	}

	String[] cNames = new String[1];
	String[] pNames = new String[1];

	cNames[0] = chromatFileName;
	pNames[0] = phdFileName;

        new TTViewFrame(cNames, pNames, secondPhdFileName, true);
    }

    public static void main(String args[]) {

	// Redirect System.err to a file.
	ErrorLog.turnOn(Constants.ERROR_LOG_FILE);
	ErrorLog.printDate();

	// Plug in the Windows Look&Feel if running on windows
	try {
	    if (Constants.OS_NAME.toLowerCase().trim().indexOf("windows")>=0) {
		UIManager.setLookAndFeel(
			"com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
	    }
	} catch (Exception e) {
	    System.err.println("Couldn't load Window look and feel: " + e);
	}

	new TTView(args);
    }
}
