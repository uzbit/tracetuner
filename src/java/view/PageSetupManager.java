/*
 * 1.2
 *
 * @(#)PageSetupManager.java	2001/03/06
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */
package com.paracel.tt.view;

import java.awt.print.*;
import com.paracel.tt.util.*;

/**
 * This class keeps track of the application-wise page setup paramters, 
 * including Paper Size, Orientation, Margins, and two parameters for Print
 * All: Approximate Number of Bases per Row, and Number of Rows per Page.
 * The default value of these parameters are: US Letter, Portrait, 1.0 inch,
 * 80, and 6, respectively.
 */
public class PageSetupManager {

    /** The current manager.  A class variable. */
    private static PageSetupManager manager;

    /** The current page format. */
    private static PageFormat pageFormat = PageSetupDialog.defaultPageFormat();

    /** The current Print All base num per row. */
    private static int paBaseNum = Constants.PA_DEFAULT_BASE_NUM_PER_ROW;

    /** The current Print All row num per page. */
    private static int paRowNum = Constants.PA_DEFAULT_ROW_NUM_PER_PAGE;

    /** Creates the manager with default settings. */
    protected PageSetupManager() {}

    /** Returns the current manager. */
    public static PageSetupManager currentManager() {
	if (manager == null) {
	    manager = new PageSetupManager();
	}
	return manager;
    }

    /** Invokes the page setup dailog, updates the page setup parameters
	application-wise when the user alters them.  It calls the 
        PageSetupDialog class method pageDialog(...) to do the actual work. */
    public static void pageDialog() {
	Object[] o = PageSetupDialog.pageDialog(pageFormat, 
						paBaseNum, paRowNum); 
	pageFormat = (PageFormat) o[0];
	paBaseNum = ((Integer) o[1]).intValue();
	paRowNum = ((Integer) o[2]).intValue();
    }

    /** Returns the default page format.  A wrapper function for 
        java.awt.print.TTPrintJob.defaultPage(). */
    public static PageFormat defaultPageFormat() { 
	return PageSetupDialog.defaultPageFormat();
    }

    /** Returns the current page format. */
    public static PageFormat getPageFormat() { return pageFormat; }

    /** Returns the current Print All row number per page. */
    public static int getPARowNum() { return paRowNum; }

    /** Returns the current Print All base number per row. */
    public static int getPABaseNum() { return paBaseNum; }
}
