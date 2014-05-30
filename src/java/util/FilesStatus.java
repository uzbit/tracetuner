/*
 * 1.5
 *
 * @(#)FilesStatus.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.util;

import java.lang.*;
import java.awt.*;
import com.paracel.tt.io.*;

/**
 * This class creates a FilesStatus object, which contains the index of the
 * sample file in the sample files list, the total number of sample files
 * in the sample files list, the name and status of the sample files,
 * the phd file, the tab file, the tal file, and the tip file.
 */
public class FilesStatus {
    /** The valid file status constant. */
    public static final int VALID = 0;

    /** The does-not-exist file status constant. */
    public static final int DOES_NOT_EXIST = -1;

    /** The can-not-read file status constant. */
    public static final int CAN_NOT_READ = -2;

    /** For tip file status and tal file status:
       This constant applies when the file is not needed. */ 
    public static final int EXEMPT = -3;

    /** The status for an tal file which has no good alignments. */
    public static final int NO_GOOD_ALIGN = -4;

    /** The status for an tal file which has no alignments
	due to bad processing. */
    public static final int BAD_PROCESSING = -5;

    /** The status for an tal file which has no alignments
	due to possible repeats. */
    public static final int POSSIBLE_REPEATS = -6;

    /** The status for an tal file which has no alignments
	due to invalid alignment data. */
    public static final int INVALID_ALIGN_DATA = -7;

    /** The status for an tal file which has no alignments
	due to unknown reason. */
    public static final int UNKNOWN = -8;

    /** The status constant applies when the file contains invalid data. */
    public static final int INVALID_DATA = -9;

    /** The status for an abc file which has no alternative base calls. */
    public static final int NO_ABC = -10;

    /** The status for an abc file which has invalid alternative base 
	call data. */
    public static final int INVALID_ABC_DATA = -11;

    /** The valid-but-stale file status constant for valid .tab, .tal, 
	and .tip files, whose last modification time is more than 1-minute
	older than that of the corresponding .phd.1 file. */
    public static final int VALID_BUT_STALE = -12;

    /** The index of the sample file in the file list. */
    private int index;

    /** The total number of sample files in the sample files batch. */
    private int numOfFiles;

    /** The sample file. */
    private ABIFile sampleFile;

    /** The phd file. */
    private PhdFile phdFile;

    /** The tab file. */
    private TabFile tabFile;

    /** The tal file. */
    private TalFile talFile;

    /** The tip file name. */
    private String sampleFileName, phdFileName, tabFileName, 
		   talFileName, tipFileName;
   
    /** The sample file status. */
    private int   sampleFileStatus;
   
    /** The phd file status. */
    private int   phdFileStatus;
   
    /** The abc file status. */
    private int   tabFileStatus;
   
    /** The tal file status. */
    private int   talFileStatus;
   
    /** The tip file status. */
    private int   tipFileStatus;

    /**
     * Creates a FilesStatus object from the specified parameters.
     * @param	i	The index of the sample file in the file list.
     * @param	num	The number of the sample files in the file list.
     * @param	s	The sample file name.
     * @param	p	The phd file name.
     * @param	a	The tal file name.
     * @param	t	The tip file name.
     * @param	ss	The sample file status.
     * @param	ps	The phd file status.
     * @param	as	The tal file status.
     * @param	ts	The tip file status.
     */
    public FilesStatus(int i, int num,
		       ABIFile sf, PhdFile pf, TabFile bf, TalFile af, 
		       String s, String p, String b, String a,
		       String t, int ss, int ps, int bs, int as, int ts) {
	index = i;
	numOfFiles = num;
	sampleFile = sf;
	phdFile = pf;
	tabFile = bf;
	talFile = af;
	sampleFileName = s;
	phdFileName = p;
	tabFileName = b;
	talFileName = a;
	tipFileName = t;
	sampleFileStatus = ss;
	phdFileStatus = ps;
	tabFileStatus = bs;
	talFileStatus = as;
	tipFileStatus = ts;
    }

    /** Returns the sample file index. */
    public int getIndex() { return index; }

    /** Returns the total number of sample files. */
    public int getNumOfFiles() { return numOfFiles; }

    /** Returns the sample file name. */
    public String getSampleFileName() { return sampleFileName; }

    /** Returns the sample file ABIFile object. */
    public ABIFile getSampleFile() { return sampleFile; }

    /** Returns the phd file name. */
    public String getPhdFileName() { return phdFileName; }

    /** Returns the phd file PhdFile object. */
    public PhdFile getPhdFile() { return phdFile; }

    /** Returns the .tab file name. */
    public String getTabFileName() { return tabFileName; }

    /** Returns the .tab file TabFile object. */
    public TabFile getTabFile() { return tabFile; }

    /** Returns the tal file name. */
    public String getTalFileName() { return talFileName; }

    /** Returns the tal file TalFile object. */
    public TalFile getTalFile() { return talFile; }

    /** Returns the tip file name. */
    public String getTipFileName() { return tipFileName; }

    /** Returns the sample file status. */
    public int getSampleFileStatus() { return sampleFileStatus; }

    /** Returns the phd file status. */
    public int getPhdFileStatus() { return phdFileStatus; }

    /** Returns the .tab file status. */
    public int getTabFileStatus() { return tabFileStatus; }

    /** Returns the .tal file status. */
    public int getTalFileStatus() { return talFileStatus; }

    /** Returns the .tip file status. */
    public int getTipFileStatus() { return tipFileStatus; }

    /** Returns an approapriate text string for the specified status code. */
    public static String getFileStatusString(int status) {
	switch (status) {
	case VALID:
	    return "OK";
	case VALID_BUT_STALE:
	    return "OK (stale)";
	case DOES_NOT_EXIST:
	    return "Does not exist";
	case CAN_NOT_READ:
	    return "Can not read";
	case INVALID_DATA:
	    return "Invalid data";
	case NO_GOOD_ALIGN:
	    return "No good alignments found";
	case POSSIBLE_REPEATS:
	    return "Possible repeats";
	case BAD_PROCESSING:
	    return "Bad processing";
	case INVALID_ALIGN_DATA:
	    return "Invalid alignment data";
	case NO_ABC:
	    return "No alternative base calls";
	case INVALID_ABC_DATA:
	    return "Invalid alternative base call data";
	case UNKNOWN:
	default:
	    return "Invalid for unknown reason";
	}
    }

    /** Returns an appropriate Color for the speicified status code. */
    public static Color getFileStatusStringColor(int status) {
	if (status == VALID || status == VALID_BUT_STALE) {
	    return (Color.green).darker();
	} else {
	    return (Color.red).darker();
	}	
    }

    public void releaseResources() {
	if (sampleFile != null) {
	    sampleFile.releaseResources();
	    sampleFile = null;
	}
	if (phdFile != null) {
	    phdFile.releaseResources();
	    phdFile = null;
	}
	if (tabFile != null) {
	    tabFile.releaseResources();
	    tabFile = null;
	}
	if (talFile != null) {
	    talFile.releaseResources();
	    talFile = null;
	}
    }
}
