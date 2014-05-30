/*
 * 1.12
 *
 * @(#)Options.java	1.0	2001/01/12
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.io.*;
import com.paracel.tt.util.*;

/**
 * This class contains the options data presented in the Options dialog
 * window and the Advanced Options dialog window, including the output format,
 * ttuner run options, output directory, etc.
 */
public class Options {

    /** The name of the .ini file used by the TraceTuner Launcher to
	store the selected options. */
    public static final String INI_FILE = Constants.USER_DIR
    				+ Constants.FILE_SEPARATOR
				+ Constants.INIT_FILE;

    /** The status of the ini file: exisits, is readable and writable,
	and has valid data. */
    public static int OK = 0;

    /** The status of the ini file: does not exist. */
    public static int DOES_NOT_EXIST = 1;

    /** The status of the ini file: not readable. */
    public static int CAN_NOT_READ = 2;

    /** The status of the ini file: exception thrown during i/o. */
    public static int EXCEPTION_THROWN = 3;

    /** The status of the ini file: has in valid data. */
    public static int INVALID_DATA = 4;

    /** The status of the ini file: can't be created. */
    public static int CAN_NOT_CREATE = 5;

    /** The status of the ini file: not writable. */
    public static int CAN_NOT_WRITE = 6;

    /** The flags for output formats. */
    private boolean phd, qual, qa, scf, seq, sa, tip, tal, tab;

    /** The names of the -qa file, and the -sa file. */
    private String qaFileName, saFileName;

    /** The names of the specified lookuptable and selected lookuptable. */
    private String otherLookupTable, lookupTable;

    /** The ttuner run conditions. */
    private boolean noCall, recallN, het, mix, edBase, shift;

    /** The trim parameters. */
    private String trimWindowSize, trimThreshold;

    /** The minimum ratio parameter. */
    private String minRatio;

    /** The consensus file name. */
    private String consensus;

    /** The context table file name. */
    private String contextTable;

    /** The output directory. */
    private String outputDir;

    /** The status of the ini file. */
    private int status;

    /** The error message. */
    private String errMesg;

    /**
     * The constructor.  Set the lookuptable to default lookup table,
     * and Phd Files as the default output format.
     */
    public Options() {
	this(null);
    }

    public Options(Options o) {
	super();
	phd = true;
	lookupTable = Constants.DEFAULT_CALIBRATION;

	if (o == null) {
	    return;
	}

        phd       = o.phd; 
        qual      = o.qual; 
        qa        = o.qa;   
        scf       = o.scf; 
        seq       = o.seq;
	sa        = o.sa;   
        tip       = o.tip;   
        tal       = o.tal; 
        tab       = o.tab;
	noCall    = o.noCall;
	recallN   = o.recallN;
	het       = o.het;
        mix       = o.mix;
	edBase    = o.edBase;
	shift     = o.shift;
	status    = o.status;

	// make deep copy of String fields
	trimWindowSize = (o.trimWindowSize == null) ? null
					: new String(o.trimWindowSize);
	trimThreshold = (o.trimThreshold == null) ? null
					: new String(o.trimThreshold);
	minRatio = (o.minRatio == null) ? null
					: new String(o.minRatio);
	consensus = (o.consensus == null) ? null
					: new String(o.consensus);
	contextTable = (o.contextTable == null) ? null
					: new String(o.contextTable);
	outputDir = (o.outputDir == null) ? null
					: new String(o.outputDir);
	errMesg = (o.errMesg == null) ? null
					: new String(o.errMesg);
	qaFileName = (o.qaFileName == null) ? null
					: new String(o.qaFileName);
	saFileName = (o.saFileName == null) ? null
					: new String(o.saFileName);
	otherLookupTable = (o.otherLookupTable == null) ? null
					: new String(o.otherLookupTable);
	lookupTable = (o.lookupTable == null) ? null
					: new String(o.lookupTable);
    }

    /** Loads the options from the ini file. */
    public void load() {
	File f = new File(INI_FILE);
	if (!f.exists()) {
	    status = DOES_NOT_EXIST;
	    return; 
	}
	if (!f.canRead()) {
	    errMesg = "Read acces of : " + INI_FILE + " was denied.\n"
		      + "Default setting will be used.";
	    status = CAN_NOT_READ;
	    return; 
	}
	try {
	    BufferedReader reader = new BufferedReader(new FileReader(f));
	    String s;
	    String line = reader.readLine();
	    while (line != null && !line.equals("END_TTUNER_INI")) {
		if ((s = parse(line, "OUTPUT_FORMAT_PHD:")) != null) {
		    phd = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "OUTPUT_FORMAT_QUAL:")) != null) {
		    qual = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "OUTPUT_FORMAT_QA:")) != null) {
		    qa = (s.startsWith("true")) ? true : false;
		    qaFileName = qa ? s.substring(5) : null;
		} else if ((s = parse(line, "OUTPUT_FORMAT_SCF:")) != null) {
		    scf = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "OUTPUT_FORMAT_SEQ:")) != null) {
		    seq = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "OUTPUT_FORMAT_SA:")) != null) {
		    sa = (s.startsWith("true")) ? true : false;
		    saFileName = sa ? s.substring(5) : null;
		} else if ((s = parse(line, "OUTPUT_FORMAT_TIP:")) != null) {
		    tip = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "OUTPUT_FORMAT_TAL:")) != null) {
		    tal = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "OUTPUT_FORMAT_TAB:")) != null) {
		    tab = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "DYE_LOOKUP_DIR:")) != null) {
		    otherLookupTable = new String(s);
		} else if ((s = parse(line, "DYE_LOOKUP_LIST:")) != null) {
		    lookupTable = new String(s);
		} else if ((s = parse(line, "PROCESS_OPTION_NOCALL:")) 
			   != null) {
		    noCall = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "PROCESS_OPTION_RECALL_N:"))
			   != null) {
		    recallN = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "PROCESS_OPTION_HET:"))!=null) {
		    het = (s.equals("true")) ? true : false;
                } else if ((s = parse(line, "PROCESS_OPTION_MIX:"))!=null) {
                    mix = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "EDITED_BASES:"))!=null) {
		    edBase = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "MOBILITY_SHIFT:"))!=null) {
		    shift = (s.equals("true")) ? true : false;
		} else if ((s = parse(line, "CHG_TRIM_WINDOW_SIZE:"))!=null) {
		    trimWindowSize = (s.startsWith("true")) ? s.substring(5) 
							    : null;
		} else if ((s = parse(line, "CHG_TRIM_THRESHOLD:"))!=null) {
		    trimThreshold = (s.startsWith("true")) ? s.substring(5) 
							   : null;
		} else if ((s = parse(line, "CHG_MIN_RATIO:"))!=null) {
		    minRatio = (s.startsWith("true")) ? s.substring(5) : null;
		} else if ((s = parse(line, "CONSENSUS:"))!=null) {
		    consensus = (s.startsWith("true")) ? s.substring(5) : null;
		} else if ((s = parse(line, "CONTEXT_TABLE:"))!=null) {
		    contextTable = (s.startsWith("true")) ? s.substring(5)
							  : null;
		} else if ((s = parse(line, "CHG_OUTPUT:"))!=null) {
		    outputDir = (s.startsWith("true")) ? s.substring(5) : null;
		}

		line = reader.readLine();
	    }
	    if (readValidate()) {
		status = OK;
	    } else {
		status = INVALID_DATA;
	    }
	    reader.close();
	} catch (Exception e) {
	    errMesg = "Error occured while reading " + INI_FILE 
		      + "\nDefault setting will be used.";
	    System.err.println("Options--load:" + e.toString());
	    e.printStackTrace(System.err);
	    status = EXCEPTION_THROWN;
	}
    }

    /** Resets the options to default values. */
    public void reset() {
	phd = true;
	qual = scf = seq = tip = tal = tab = qa = sa = false;
	noCall = recallN = het = mix = edBase = shift = false;
	lookupTable = Constants.DEFAULT_CALIBRATION;
	otherLookupTable = trimWindowSize = trimThreshold = minRatio
			 = outputDir = consensus = contextTable = null;
    }

    /** Saves the options to the ini file. */
    public void save() {
	File f = new File(INI_FILE);
	try {
	    if (!f.exists()) {
	    	f.createNewFile();
	    }
	} catch (Exception e) {
	    errMesg = "Can not create the file: " + INI_FILE;
	    status = CAN_NOT_CREATE;
	    return;
	}
	if (!f.canWrite()) {
	    errMesg = "Can not write to the file: " + INI_FILE;
	    status = CAN_NOT_WRITE;
	    return; 
	}
	try {
	    BufferedWriter writer = new BufferedWriter(new FileWriter(f));
	    writer.write("BEGIN_TTUNER_INI" + Constants.NEW_LINE);
	    writer.write("OUTPUT_FORMAT_PHD: " + phd + Constants.NEW_LINE);
	    writer.write("OUTPUT_FORMAT_QUAL: " + qual + Constants.NEW_LINE);
	    writer.write("OUTPUT_FORMAT_SCF: " + scf + Constants.NEW_LINE);
	    writer.write("OUTPUT_FORMAT_SEQ: " + seq + Constants.NEW_LINE);
	    writer.write("OUTPUT_FORMAT_TIP: " + tip + Constants.NEW_LINE);
	    writer.write("OUTPUT_FORMAT_TAL: " + tal + Constants.NEW_LINE);
	    writer.write("OUTPUT_FORMAT_TAB: " + tab + Constants.NEW_LINE);
	    if (qaFileName != null && qa) {
	    	writer.write("OUTPUT_FORMAT_QA: true " + qaFileName 
			     + Constants.NEW_LINE);
	    } else {
		writer.write("OUTPUT_FORMAT_QA: false " + Constants.NEW_LINE);
	    }
	    if (saFileName != null && sa) {
	    	writer.write("OUTPUT_FORMAT_SA: true " + saFileName 
			     + Constants.NEW_LINE);
	    } else {
		writer.write("OUTPUT_FORMAT_SA: false " + Constants.NEW_LINE);
	    }

	    if (lookupTable != null) {
	    	writer.write("DYE_LOOKUP_LIST: " + lookupTable 
			     + Constants.NEW_LINE);
	    } else {
	    	writer.write("DYE_LOOKUP_DIR: " + otherLookupTable
			     + Constants.NEW_LINE);
	    }
	    writer.write("PROCESS_OPTION_NOCALL: " + noCall 
			 + Constants.NEW_LINE);
	    writer.write("PROCESS_OPTION_RECALL_N: " + recallN 
			 + Constants.NEW_LINE);
	    writer.write("PROCESS_OPTION_HET: " + het 
			 + Constants.NEW_LINE);
            writer.write("PROCESS_OPTION_MIX: " + mix
                         + Constants.NEW_LINE);
	    writer.write("EDITED_BASES: " + edBase
			 + Constants.NEW_LINE);
	    writer.write("MOBILITY_SHIFT: " + shift
			 + Constants.NEW_LINE);
	    if (trimWindowSize != null) {
	    	writer.write("CHG_TRIM_WINDOW_SIZE: true " + trimWindowSize 
			     + Constants.NEW_LINE);
	    } else {
	    	writer.write("CHG_TRIM_WINDOW_SIZE: false " 
			     + Constants.NEW_LINE);
	    }
	    if (trimThreshold != null) {
	    	writer.write("CHG_TRIM_THRESHOLD: true " + trimThreshold 
			     + Constants.NEW_LINE);
	    } else {
	    	writer.write("CHG_TRIM_THRESHOLD: false " 
			     + Constants.NEW_LINE);
	    }
	    if (minRatio != null) {
	    	writer.write("CHG_MIN_RATIO: true " + minRatio 
			     + Constants.NEW_LINE);
	    } else {
	    	writer.write("CHG_MIN_RATIO: false " 
			     + Constants.NEW_LINE);
	    }
	    if (consensus != null) {
	    	writer.write("CONSENSUS: true " + consensus 
			     + Constants.NEW_LINE);
	    } else {
		writer.write("CONSENSUS: false " + Constants.NEW_LINE);
	    }
	    if (contextTable != null) {
	    	writer.write("CONTEXT_TABLE: true " + contextTable 
			     + Constants.NEW_LINE);
	    } else {
		writer.write("CONTEXT_TABLE: false " + Constants.NEW_LINE);
	    }
	    if (outputDir != null) {
	    	writer.write("CHG_OUTPUT: true " + outputDir 
			     + Constants.NEW_LINE);
	    } else {
	    	writer.write("CHG_OUTPUT: false " + Constants.NEW_LINE);
	    }

	    writer.write("END_TTUNER_INI" + Constants.NEW_LINE);
	    status = OK;
	    writer.close();
	} catch (Exception e) {
	    errMesg = "Error saving " + INI_FILE + ":\n"+ e.getMessage();
	    status = EXCEPTION_THROWN;
	}
    }

    /** 
     * Validates the options read from the ini file.
     * @return true if valid; false otherwise.
     */
    public boolean readValidate() {
	if (!(phd || qual || qa || scf || seq || sa || tip || tal || tab)) {
	    errMesg = "At least one output format is required.";
	    return false;
	}
	if (tab && (!phd)) {
	    errMesg = "Phd Files option is required when Tab Files option"
			+ " is selected.";
	    return false;
	}
	if (tal && (!phd)) {
	    errMesg = "Phd Files option is required when Tal Files option"
			+ " is selected.";
	    return false;
	}
	if (tip && (!phd)) {
	    errMesg = "Phd Files option is required when Tip Files option"
			+ " is selected.";
	    return false;
	}
	if ((noCall && recallN) || (recallN && het) || (noCall && het) ||
            (noCall && mix    ) || (recallN && mix) || (mix    && het)) 
        {
	    errMesg = "Only one option allowed among noCall, "
			+ "recallN, het and mix.";
	    return false;
	}
	if (lookupTable == null && otherLookupTable == null) {
	    errMesg = "At least one lookup table is required.";
	    return false;
	}
	return true;
    }

    /** 
     * Parses the specified string according to the specified pattern. 
     * @return the rest of the string after the specified pattern;
     *	       null, if the specified string does not start with the
     *	       specified pattern.
     */
    protected String parse(String s, String pat) {
	if (!s.startsWith(pat)) {
	    return null;
	}
	String ss = s.substring(pat.length() + 1, s.length()).trim();
	return ss;
    }

    public boolean outputPhd() { return phd; }
    public boolean outputQual() { return qual; }
    public boolean outputQa() { return qa; }
    public boolean outputScf() { return scf; }
    public boolean outputSeq() { return seq; }
    public boolean outputSa() { return sa; }
    public boolean outputTip() { return tip; }
    public boolean outputTal() { return tal; }
    public boolean outputTab() { return tab; }
    public String getOtherLookupTable() { return otherLookupTable; }
    public String getLookupTable() { return lookupTable; }
    public String getTrimThreshold() { return trimThreshold; }
    public String getTrimWindowSize() { return trimWindowSize; }
    public String getMinRatio() { return minRatio; }
    public String getConsensus() { return consensus; }
    public String getContextTable() { return contextTable; }
    public String getOutputDir() { return outputDir; }
    public String getQAFileName() { return qaFileName; }
    public String getSAFileName() { return saFileName; }
    public boolean noCall() { return noCall; }
    public boolean recallN() { return recallN; }
    public boolean het() { return het; }
    public boolean mix() { return mix; }
    public boolean editedBase() { return edBase; }
    public boolean shift() { return shift; }
    public int getStatus() { return status; }
    public String getErrorMessage() { return errMesg; }

    public void setOutputPhd(boolean b) { phd = b; }
    public void setOutputQual(boolean b) { qual = b; }
    public void setOutputQa(boolean b) { qa = b; }
    public void setOutputScf(boolean b) { scf = b; }
    public void setOutputSeq(boolean b) { seq = b; }
    public void setOutputSa(boolean b) { sa = b; }
    public void setOutputTip(boolean b) { tip = b; }
    public void setOutputTal(boolean b) { tal = b; }
    public void setOutputTab(boolean b) { tab = b; }
    public void setOtherLookupTable(String s) { otherLookupTable = s; }
    public void setLookupTable(String s) { lookupTable = s; }
    public void setTrimThreshold(String s) { trimThreshold = s; }
    public void setTrimWindowSize(String s) { trimWindowSize = s; }
    public void setMinRatio(String s) { minRatio = s; }
    public void setConsensus(String s) { consensus = s; }
    public void setContextTable(String s) { contextTable = s; }
    public void setOutputDir(String s) { outputDir = s; }
    public void setQAFileName(String s) { qaFileName = s; }
    public void setSAFileName(String s) { saFileName = s; }
    public void setNoCall(boolean b) { noCall = b; }
    public void setRecallN(boolean b) { recallN = b; }
    public void setHet(boolean b) { het = b; }
    public void setMix(boolean b) { mix = b; }
    public void setEdBase(boolean b) { edBase = b; }
    public void setShift(boolean b) { shift = b; }

    /**
     * Builds the command array using the specified string as
     * the output directory.
     * @return the command builder.
     */
    public CommandBuilder buildCommand(String outDir) {

	CommandBuilder cb = new CommandBuilder();
	outDir = outDir.trim();
	boolean hasSpace = (outDir.indexOf(" ") >= 0) ? true : false;
	
	cb.append(Constants.COMMAND_NAME);
	String qvFile = outDir + Constants.FILE_SEPARATOR
			       + Constants.QV_REPORT_FILE;
	String qaFile = outDir + Constants.FILE_SEPARATOR + qaFileName;
	String saFile = outDir + Constants.FILE_SEPARATOR + saFileName;
	cb.append("-qr");
	cb.append(qvFile);

	if (noCall) {
	    cb.append("-nocall");
	}
	if (recallN) {
	    cb.append("-recalln");
	}
	if (het) {
	    cb.append("-het");
	}
        if (mix) {
            cb.append("-mix");
        }
	if ((het || mix) && minRatio != null) {
	    cb.append("-min_ratio");
	    cb.append(minRatio+"");
	}
	if (contextTable != null && Constants.runAsDev) {
	    cb.append("-ct");
	    cb.append(contextTable);
	}
	if (edBase) {
	    cb.append("-edited_bases");
	}
	if (shift) {
	    cb.append("-shift");
	}
	if (otherLookupTable != null) {
	    cb.append("-t");
	    cb.append(otherLookupTable);
	} else if (lookupTable != null) {
        if (lookupTable.equals(Constants.ABI3700_POP5_LOOKUP_TABLE)) {
                        cb.append("-3730");
        } else if (lookupTable.equals(Constants.ABI3700_POP5_LOOKUP_TABLE)) {
	    	cb.append("-3700pop5");
	    } else if (lookupTable.equals(
				Constants.ABI3700_POP6_LOOKUP_TABLE)) {
	    	cb.append("-3700pop6");
	    } else if (lookupTable.equals(
				Constants.ABI3100_POP6_LOOKUP_TABLE)) {
	    	cb.append("-3100");
	    } else if (lookupTable.equals(Constants.MegaBACE_LOOKUP_TABLE)) {
	    	cb.append("-mbace");
	    } else if (!lookupTable.equals(Constants.DEFAULT_CALIBRATION)) {
		cb.append("-t");
		cb.append(lookupTable);
	    }
	}

	if (trimWindowSize != null) {
	    cb.append("-trim_window");
	    cb.append(trimWindowSize + "");
	}
	if (trimThreshold != null) {
	    cb.append("-trim_threshold");
	    cb.append(trimThreshold + "");
	}
	if (phd) {
	    cb.append("-pd");
	    cb.append(outDir);
	}
	if (qual) {
	    cb.append("-qd");
	    cb.append(outDir);
	}
	if (qa) {
	    cb.append("-qa");
	    cb.append(qaFile);
	}
	if (seq) {
	    cb.append("-sd");
	    cb.append(outDir);
	}
	if (sa) {
	    cb.append("-sa");
	    cb.append(saFile);
	}
	if (scf) {
	    cb.append("-cd");
	    cb.append(outDir);
	}
	if (tip) {
	    cb.append("-tipd");
	    cb.append(outDir);
	}
	if (tal) {
	    if (consensus != null) {
		cb.append("-C");
		cb.append(consensus);
	    }
	    cb.append("-tald");
	    cb.append(outDir);
	}
	if (tab) {
	    cb.append("-tabd");
	    cb.append(outDir);
	}

	return cb;
    }

}
