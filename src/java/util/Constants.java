/*
 * 1.12
 *
 * @(#)Constants.java	2001/01/12
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */
package com.paracel.tt.util;

import java.awt.*;

/** 
 * The class holds constant data used accross the Viewer and Launcher.
 * such as log file name, ttuner command name, ini file name, version
 * string, phd file suffix, alignment file suffix, etc.
 */
public class Constants {
    public static final String OS_NAME = System.getProperty("os.name");

    public static final String CLASS_PATH = System.getProperty("java.class.path");
    public static final String USER_DIR = System.getProperty("user.dir");
    public static final String FILE_SEPARATOR = System.getProperty("file.separator");
    public static final String NEW_LINE = System.getProperty("line.separator");

    public static final String COMMAND_NAME = "ttuner";
    public static final String INIT_FILE = "ttuner.ini";
    public static final String FILES_NAME_FILE = "ttfiles.txt";
    public static final String QV_REPORT_FILE = "qvreport.txt";
    public static final String ERROR_LOG_FILE = "ttstderr.log";
    public static final String LOG_FILE = "ttlog.txt";

    public static final String LOOKUP_TABLE_SUFFIX = ".tbl";
    public static final String DEFAULT_CALIBRATION = "Automatic selection";
    public static final String ABI3730_POP7_LOOKUP_TABLE = "3730_Pop7_BigDye";
    public static final String ABI3700_POP5_LOOKUP_TABLE = "3700_Pop5_BigDye";
    public static final String ABI3700_POP6_LOOKUP_TABLE = "3700_Pop6_BigDye";
    public static final String ABI3100_POP6_LOOKUP_TABLE = "3100_Pop6_BigDye";
    public static final String MegaBACE_LOOKUP_TABLE = "MegaBACE";

    public static final float DEFAULT_MIN_RATIO = 0.1f; 
    public static final int DEFAULT_TRIM_THRESHOLD = 20; 
    public static final int DEFAULT_TRIM_WINDOW_SIZE = 10; 

    public static final String PHD_FILE_SUFFIX = ".phd.1";
    public static final String TAL_FILE_SUFFIX = ".tal";
    public static final String TAB_FILE_SUFFIX = ".tab";
    public static final String TIP_FILE_SUFFIX = ".tip";

    public static final String VERSION = "3.0";

    public static final Color gray208 = new Color(208, 208, 208);
    public static final Color gray148 = new Color(148, 148, 148);
    public static final Color gray108 = new Color(108, 108, 108);
    public static final Color green = new Color(0, 218, 0);

    /** Default row num per page for Print All print-outs. */
    public static final int PA_DEFAULT_ROW_NUM_PER_PAGE = 6;

    /** Maximum row num per page for Print All print-outs. */
    public static final int PA_MAX_ROW_NUM_PER_PAGE = 8;

    /** Default base num per row for Print All print-outs. */
    public static final int PA_DEFAULT_BASE_NUM_PER_ROW = 80;

    /** Average scans for each base call. */
    public static final int AVERAGE_SCANS_PER_BASE = 12;

    /** Whether or not TT is brought up for internal development purpose. */
    public static boolean runAsDev = false;
    
    /** whether or not the "Edited bases" Advanced Option is true. */
    public static boolean fromEdited = false;
}
