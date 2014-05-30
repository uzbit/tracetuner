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

package com.paracel.tt.view;

import java.io.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import com.paracel.tt.util.DefaultBrowser;

public class HelpAction extends AbstractAction {

    /** Help Docs directory.  This is passed into JVM using
        "java -Dhelpdocs.dir=foo -jar ..."
     */
    public static String HELP_DOCS_DIR = System.getProperty("helpdocs.dir");

    /** Help Docs master html file name */
    public static String HELP_DOCS_MASTER_HTML_NAME = "wwhelp.htm";

    /** File separator */
    public static String FILE_SEPARATOR = System.getProperty("file.separator");

    /** Help Docs master html file path */
    public static String HELP_DOCS_MASTER_HTML
        = HELP_DOCS_DIR + FILE_SEPARATOR + HELP_DOCS_MASTER_HTML_NAME;


    public HelpAction(String title) {
	super(title);
    }

    public void actionPerformed(ActionEvent event) {

	// 1. checks whether the HTML master file exists
	String masterHTML = HELP_DOCS_MASTER_HTML;
	try {
	    File temp = new File(masterHTML);
	    if (! temp.exists()) {
		errorMessage("Cannot view Help Contents.\n"
				    + "File \"" + masterHTML
				    + "\" does not exist.");
		return;
	    }
	    masterHTML = temp.getAbsolutePath();
	} catch (Exception e) {
	    errorMessage("Cannot view Help Contents.\n"
				+ "File \"" + masterHTML
				+ "\" might have been corrupted.");
	    return;
	}

	// 2. brings up the HTML help in system default browser
	try {
	    // No Mac OS support for now.
	    //if ( isMacOS() )
		//DefaultBrowser.openURL("file://localhost" + masterHTML);
	    //else
		DefaultBrowser.openURL(masterHTML);
	} catch (Exception e) {
	    errorMessage("Cannot view Help Contents.\n\n"
				//+ e.getMessage()
				+ "Please make sure a web browser is"
				+ " installed on your system.\n");
	}
    }

    private void errorMessage(String message) {
        JOptionPane.showMessageDialog(null, message, "Error",
                                      JOptionPane.ERROR_MESSAGE);
    }
}

