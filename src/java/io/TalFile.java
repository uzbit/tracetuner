/* 1.7
 * 
 * @(#) TalFile.java	2000/12/05
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 *
 */

package com.paracel.tt.io;

import java.io.*;
import java.util.*;

/**
 *   This class reads and holds in memory the data from a .tal file.  All
 *   get operations are done from memory so they should be fairly quick
 *   once the file is read in.  This class will not notice if the file is
 *   deleted once the class is instantiated.
 */
public class TalFile extends RandomAccessFile {
    /** The constant which means the tal file has valid alignment data. */
    public final static int OK = 0;

    /** The status constant which means the tal file has no good alignments. */
    public final static int NO_GOOD_ALIGN = 1;

    /** The status constant which means the tal file has no good data
	due to possible repeats. */
    public final static int POSSIBLE_REPEATS = 2;

    /** The status constant which means the tal file has no good data
	due to bad processing. */
    public final static int BAD_PROCESSING = 3;

    /** The status constant which means the tal file has no good data
	due to an unknown reason. */
    public final static int UNKNOWN = -1;

    /** Whether the file contains invalid data: such as ALIGN_LENGTH is
	missing or negative, ALIGN_LENGTH does not agree with the actual
	alignment length present in the file, etc. */
    private boolean invalidFileData = false;
	
    /** The alignment statuts. */
    private int alignStatus;

    /**
     * The flag shows whether Trace Tuner was used to generate the 
     * base calls of the alignment in the .tal file.
     */ 
    private boolean isTT = false;

    /**
     * The version string of the software used to generate the base calls.
     */
    private String version = "unknown";

    /** The chromat file name. */
    private String chromatFileName = "unknown";

    /** The consensus file name. */
    private String consensusFileName = "unknown";

    /** The length of the alignment, if any. */
    private int alignLength = -1;

    /** The left trim, read from the .tal file. */
    private int leftTrim    = 0;

    /** The right trim, read from the .tal file. */
    private int rightTrim    = 0; 

    /** 
     * The array holds the data of the fragment base and the corresponding
     * consensus pairs in the alignment. 
     */
    private AlignData[]  alignData;

    /** Whether or not there exists invalid data in the alignment, if any. */
    private boolean invalidAlignData;


    /**
     * Creates the TalFile object from the specified .tal file.
     * @param	name 	The .tal file name.
     * @exception FileNotFoundException
     * @exception IOException
     * @see     java.io.FileNotFoundException
     * @see     java.io.IOException
     */
    public TalFile(String name) throws FileNotFoundException, IOException {

	super(name,"r");

	/*
         *  Read the tal file pre-amble and pick up as many of the fields
	 *  we're interested in as are present.
	 */
	String s;
	String line = readLine();
	if (line != null) {
	    line = line.trim();
	}

	alignStatus = UNKNOWN;
	invalidAlignData = false;

	while(line != null && line.charAt(0) == '#') {

	    if (line.length() == 1) {
		// Skip the lines only containing "#"
	    	line = readLine();
	    	if (line != null) {
	    	    line = line.trim();
	        }
		continue;
	    }

	    // trim '#' from the beginning of the string
	    line = line.substring(1, line.length()).trim();

	    if((s = parse(line, "CHROMAT_FILE:")) != null) {
		chromatFileName = s;
	    }
	    
	    if((s = parse(line, "CONSENSUS_FILE:")) != null) {
		consensusFileName = s;
	    } else {
		/* the sample file contains the "aSEQ" tag, which
		   means the sample file self-contains the consensus
		   or rather, the correct sequence. */
		consensusFileName = chromatFileName;
	    }

	    if((s = parse(line, "SOFTWARE_VERSION:")) != null) {
		version = s;
		if(version.substring(0,2).equals("TT")) {
		    isTT = true;
		} else {
		    isTT = false;
		}
	    }

	    if((s = parse(line, "ALIGN_STATUS:")) != null) {
		s = s.toLowerCase();
		if (s.indexOf("ok") >= 0) {
		    alignStatus = OK;
		} else if (s.indexOf("no good alignment") >= 0) {
		    alignStatus = NO_GOOD_ALIGN;
		} else if (s.indexOf("possible repeats") >= 0) {
		    alignStatus = POSSIBLE_REPEATS;
		} else if (s.indexOf("bad processing") >= 0) {
		    alignStatus = BAD_PROCESSING;
		} else {
		    alignStatus = UNKNOWN;
		}
	    }

	    if((s = parse(line, "ALIGN_LENGTH:")) != null) {
		try {
		    alignLength = Integer.parseInt(s);
		} catch (NumberFormatException e) {
		    invalidFileData = true;
		    close();
		    return;
		}
	    }

	    line = readLine();
	    if (line != null) {
	    	line = line.trim();
	    }
	}

	if (alignStatus != OK) {
	    close();
	    return; 
	}

	if(line == null || alignLength < 0) {
	    close();
	    invalidFileData = true;
	    return;
	}

	alignData = new AlignData[alignLength];

	// Read align data
	int i = 0;
	while(line != null) {
	    line = line.trim();
	    if (i >= alignLength) {
		invalidFileData = true; 
		close();
		return;
	    }
	    alignData[i] = parseAlignData(line);
	    if (alignData[i] == null) {
		invalidAlignData = true;
		close();
		return;
	    }
	    i++;
	    line = readLine();
	}

	if (i != alignLength) {
	    invalidFileData = true;
	}

	close();
    }


    /**
     * Parses each comment line (line starting with "#") read from 
     * the .tal file.
     * @param	s	The line string.
     * @param	pat	The pattern string.
     * @return		null, if s does not start with pat or s == pat;
     *			the rest of s after pat, otherwise.
     */
    private String parse(String s, String pat) {

	if(!s.startsWith(pat)) {
	    return null;
	}
	if (s.length() == pat.length()) {
	    return null;
	}
	s = s.substring(pat.length()+1,s.length());
	return s.trim();
    }

    /**
     * Parses each alignment data line read from the .tal file.
     * @param	s	The data line string.
     * @return		null, if this line contains invalid data;
     *			an AlignData object holding the data, otherwise.
     */
    private AlignData parseAlignData(String s) {
	StringTokenizer st = new StringTokenizer(s, "\t");
	//  If number of tokens is not 5, we have a bad line.
	if (st.countTokens() < 5) {
	    // bad line
	    return null;
	}
	AlignData temp = new AlignData();
	try {
	    temp.fragPosition   = Integer.parseInt(st.nextToken());
	    temp.fragBase = st.nextToken().charAt(0);
	    temp.consPosition   = Integer.parseInt(st.nextToken());
	    temp.consBase = st.nextToken().charAt(0);
	    int i = Integer.parseInt(st.nextToken());
	    if (i == 1) {
		temp.isMatch = true;
	    } else if (i == 0) {
		temp.isMatch = false;
	    } else {
		// Invalid line
		return null;
	    }
	} catch (Exception e) {
	    return null;
	}
	if (temp.fragPosition < 0 || temp.consPosition < 0) {
	    return null;
	}
	return temp;
    }

    /** Attempts to release the resources. */
    public void releaseResources() {
	if (alignData != null) {
	    for (int i = 0; i < alignData.length; i++) {
	    	alignData[i] = null;
	    }
	    alignData = null;
	}
    }


    /**
     * Gets the length of the alignment.
     * @return	the length of the alignment.
     */
    public int      getAlignLength()     { return alignLength;   }

    /**
     * Identifies whether or not the bases were called by TraceTuner.
     * @return	true, if the bases are called by TraceTuner;
     *		false, otherwise.
     */
    public boolean  isTT()            { return isTT;        }

    /**
     * Gets the version of the software used to call bases.
     * @return	the version string.
     */
    public String   getVersion()      { return version;     }

    /**
     * Gets the chromat file name.
     * @return the chromatogram file name.
     */
    public String   getChromatFileName()  { return chromatFileName;}

    /**
     * Gets the consensus file name.
     * @return the consensus file name.
     */
    public String   getConsensusFileName()  { return consensusFileName;}

    /**
     * Gets the alignment status.
     * @return the alignment status.
     */
    public int      getAlignStatus() { return alignStatus; }

    /**
     * Identifies whether or not the alignment contains invalid data.
     * @return true, if there are invalid data; false, otherwise.
     */
    public boolean  hasInvalidAlignData() { return invalidAlignData; }

    /**
     * Identifies whether or not the file contains invalid data, such as:
     * ALIGN_LENGTH is negative or missing, ALIGN_LENGTH does not agree
     * with the actual alignment length present in the file, etc.
     * @return true, if there are invalid data; false, otherwise.
     */
    public boolean  hasInvalidFileData() { return invalidFileData; }

    /**
     * Gets the AlignData array that holds the data of each aligned pair.
     * @return 	An AlignData array of all aligned base pair data.
     */
    public AlignData[] getAlignData() { return alignData; }

    /**
     * main() function for testing only.  Normally this class wouldn't 
     * be stand-alone.
     */
    public static void main(String[] args) {
	TalFile p;
	int  i, alignLength;
	AlignData[] alignData;

	if(args.length != 1) {
	    System.out.println("Usage: java TalFile <filename>");
	    System.exit(0);
	}
	try {
	    p = new TalFile(args[0]);
	}
	catch (Exception e) {
	    System.out.println("Error: "+e.toString());
	    e.printStackTrace(System.out);
	    return;
	}

	System.out.println("Input file      is " + args[0]);
	System.out.println("length of alignment is " + p.getAlignLength());
	System.out.println("version         is " + p.getVersion() + 
			   (p.isTT() ? "  (Tracetuner)" : "  (Phred)"));

	System.out.println("chromat_file    is " + p.getChromatFileName());
	System.out.println("consensus_file    is " + p.getConsensusFileName());

	alignData = p.getAlignData();

	System.out.println("The first 30 aligns are: ");
	for(i=0; i < 30; i++) {
	    System.out.print(alignData[i].getFragPosition()+"\t");
	    System.out.print(alignData[i].getFragBase()+"\t");
	    System.out.print(alignData[i].getConsPosition()+"\t");
	    System.out.print(alignData[i].getConsBase()+"\t");
	    if (alignData[i].isMatch()) {
	        System.out.print("1\n");
	    } else {
	        System.out.print("0\n");
	    }
	}

	System.out.println("The last 20 aligns are: ");
	for(i=p.getAlignLength() - 20; i < p.getAlignLength(); i++) {
	    System.out.print(alignData[i].getFragPosition()+"\t");
	    System.out.print(alignData[i].getFragBase()+"\t");
	    System.out.print(alignData[i].getConsPosition()+"\t");
	    System.out.print(alignData[i].getConsBase()+"\t");
	    if (alignData[i].isMatch()) {
	        System.out.print("1\n");
	    } else {
	        System.out.print("0\n");
	    }
	}
    }


    /**
     * A public inner class of TalFile.
     * This class creates an AlignData object, which contains information
     * about one align data in the alignment, including: the fragment base,
     * the fragment position, the consensus base, the consensus position,
     * and whether this fragment base matchs this consensus base.
     */
    public class AlignData implements Cloneable {

	/** The position of the fragment base in the fragment. */
	private int fragPosition;

	/** The fragment base code. */
	private char fragBase;

	/** The position of the consensus base in the consensus. */
	private int  consPosition;

	/** The fragment base code. */
	private char consBase;

	/** 
	 * The flag indicates whether this fragment base matchs this
	 * consensus base.
	 */
	private boolean isMatch;

	/** Constructor */
	public AlignData() {}

	/** Constructor */
	public AlignData(int fragPosition, char fragBase,
			 int consPosition, char consBase, boolean isMatch) {
	    this.fragPosition = fragPosition;
	    this.fragBase = fragBase;
	    this.consPosition = consPosition;
	    this.consBase = consBase;
	    this.isMatch = isMatch;
	}

	/** Constructor */
	public AlignData(AlignData a) {
	    this(a.fragPosition, a.fragBase, 
		 a.consPosition, a.consBase, a.isMatch);
	}

	/** Creates a clone. */
	public Object clone() {
	    try {
		return super.clone();
	    } catch (CloneNotSupportedException e) {
		throw new InternalError(e.toString());
	    }
	}

	/**
	 * Gets the position of this fragment base in the fragment sequence.
	 * @return the position.
	 */
	public int getFragPosition() { return fragPosition; }

	/**
	 * Gets the fragment base code.
	 * @return the base code.
	 */
	public char getFragBase() { return fragBase; }

	/**
	 * Gets the position of this consensus base in the consensus sequence.
	 * @return the position.
	 */
	public int getConsPosition() { return consPosition; }

	/**
	 * Gets the consensus base code.
	 * @return the base code.
	 */
	public char getConsBase() { return consBase; }

	/**
	 * Identifies whether this fragment base matches this consensus base.
	 * @return	true, if it is a match; false, otherwise.
	 */
	public boolean isMatch() { return isMatch; }
    }
}
