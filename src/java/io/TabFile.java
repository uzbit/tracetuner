/* 1.5
 * 
 * @(#) TabFile.java	2001/01/04
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */
package com.paracel.tt.io;

import java.io.*;
import java.lang.*;
import java.util.*;

/**
 *   This class reads and holds in memory the data from a .tab file.  All
 *   get operations are done from memory so they should be fairly quick
 *   once the file is read in.  This class will not notice if the file is
 *   deleted once the class is instantiated.
 */
public class TabFile extends RandomAccessFile {
    /** The type constant for substitution alternative base calls */
    public static int SUBSTITUTION = 1;

    /** The type constant for deletion abcs. */
    public static int DELETION = 2;

    /** Total number of alternative base calls. */
    private int numOfAbcs, numOfHighestQVAbcs;

    /** Total number of substitution alternative base calls. */
    private int numOfSubAbcs;

    /** Total number of deletion alternative base calls. */
    private int numOfDelAbcs;

    /**
     * The flag shows whether Trace Tuner was used to generate the 
     * base calls of the alignment in the .tab file.
     */ 
    private boolean isTT = false;

    /**
     * The version string of the software used to generate the base calls.
     */
    private String version = "unknown";

    /** The chromat file name. */
    private String chromatFileName = "unknown";

    /** The array holds the data of all highestQV alternative base calls.
        For each TT base call that has alternative base calls (abc), there
        may be more than 1 alternative base calls, highestQV abc means
        the one with the highest QV. */
    private AbcData[] highestQVAbcData;

    /** Whether or not there exists invalid data, such as: NUM_ABC does not
	agree with the acutal number of alternatice base calls presented in
	the file, NUM_SUBSTITUTIONS does not agree with the actual number of
	substitution abcs present in the file, etc. */
    private boolean invalidFileData;

    /** Whether or not there exists invalid abc data, if any.  For example:
     *  Any of the abc data rows misses a column, or has wrong format, etc.
     */
    private boolean invalidAbcData;

    /**
     * Creates the TabFile object from the specified .abc file.
     * @param	name 	The .abc file name.
     * @exception FileNotFoundException
     * @exception IOException
     * @see     java.io.FileNotFoundException
     * @see     java.io.IOException
     */
    public TabFile(String name) throws FileNotFoundException, IOException {
	super(name,"r");

        /* Read the abc file pre-amble and pick up as many of the fields
	   we're interested in as are present. */
	String s;
	String line = readLine();
	if (line != null) {
	    line = line.trim();
	}

	numOfAbcs = -1;
	numOfSubAbcs = -1;
	numOfDelAbcs = -1;
	invalidAbcData = false;

	while(line != null && line.charAt(0) == '#') {

	    if (line.length() == 1) {
		// Skip the lines only containing "#"
	    	line = readLine();
	    	if (line != null) {
	    	    line = line.trim();
	        }
		continue;
	    }

	    // trim "# " from the beginning of the string
	    line = line.substring(1, line.length()).trim();

	    if((s = parse(line, "CHROMAT_FILE:")) != null) {
		chromatFileName = s;
	    }

	    if((s = parse(line, "SOFTWARE_VERSION:")) != null) {
		version = s;
		if(version.substring(0,2).equals("TT")) {
		    isTT = true;
		} else {
		    isTT = false;
		}
	    }

	    if ((s = parse(line, "ABC_STATUS:")) != null) {
		if (s.equals("No alternative base calls")) {
		    invalidFileData = false;
		    numOfAbcs = 0;
		    numOfSubAbcs = 0;
		    numOfDelAbcs = 0;
		    close();
		    return;
		}
	    }

	    if((s = parse(line, "NUM_ABC:")) != null) {
		try {
		    numOfAbcs = Integer.parseInt(s);
		} catch (NumberFormatException e) {
		    invalidFileData = true;
		    close();
		    return;
		}
	    }

	    if((s = parse(line, "NUM_SUBSTITUTIONS:")) != null) {
		try {
		    numOfSubAbcs = Integer.parseInt(s);
		} catch (NumberFormatException e) {
		    invalidFileData = true;
		    close();
		    return;
		}
	    }

	    if((s = parse(line, "NUM_DELETIONS:")) != null) {
		try {
		    numOfDelAbcs = Integer.parseInt(s);
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

	if (line == null || numOfAbcs < 0
	        || numOfSubAbcs < 0 || numOfDelAbcs < 0
	    	|| (numOfSubAbcs + numOfDelAbcs) != numOfAbcs) {
	    invalidFileData = true;
	    close();
	    return;
	}

	AbcData[] abcData = new AbcData[numOfAbcs];

	// Read align data
	int iHighestQVAbcs = 0;
	int iSubAbcs = 0;
	int iDelAbcs = 0;
	AbcData currHighestQVAbc = null, currAbc;
	while(line != null) {
	    line = line.trim();
	    if (iHighestQVAbcs >= numOfAbcs
		    || iSubAbcs > numOfSubAbcs
		    || iDelAbcs > numOfDelAbcs) {
		invalidFileData = true;
		close();
		return;
	    }
	    currAbc = parseAbcData(line);
	    if (currAbc == null) {
		invalidAbcData = true;
		close();
		return;
	    } 
	    if (currAbc.type == SUBSTITUTION) {
		iSubAbcs++;
	    } else {
		//presumably it is DELETION
		iDelAbcs++; 
	    }
	    // store the highestQV abcs
	    if (currHighestQVAbc == null) {
		currHighestQVAbc = currAbc;
	    } else {
		if (currAbc.index != currHighestQVAbc.index) {
	    	    abcData[iHighestQVAbcs++] = currHighestQVAbc;
		    currHighestQVAbc = currAbc;
	    	} else {
		    // compare the QVs for the abcs of the same phd base call,
		    // set the one with the highest QV to be the highestQV abc.
		    if (currAbc.qv > currHighestQVAbc.qv) {
		    	currHighestQVAbc = currAbc;
		    }
	    	}
	    }
	    line = readLine();
	    if (line == null && currHighestQVAbc != null) {
	    	abcData[iHighestQVAbcs++] = currHighestQVAbc;
	    }
	}

	if (iSubAbcs != numOfSubAbcs || iDelAbcs != numOfDelAbcs) {
	    invalidFileData = true;
	}

	if (iHighestQVAbcs == numOfAbcs) {
	    highestQVAbcData = abcData;
	} else {
	    highestQVAbcData = new AbcData[iHighestQVAbcs];
	    System.arraycopy(abcData, 0, highestQVAbcData, 0, iHighestQVAbcs);
	}
	numOfHighestQVAbcs = iHighestQVAbcs;

	close();
    }


    /**
     * Parses each comment line (line starting with "#") read from 
     * the .tab file.
     * @param	s	The line string.
     * @param	pat	The pattern string.
     * @return		null, if s does not start with pat or s == pat;
     *			the rest of s after pat, otherwise.
     */
    private String parse(String s, String pat) {

	if (!s.startsWith(pat)) {
	    return null;
	}
	if (s.length() == pat.length()) {
	    return null;
	}
	s = s.substring(pat.length()+1,s.length());
	return s.trim();
    }

    /**
     * Parses each alignment data line read from the .tab file.
     * @param	s	The data line string.
     * @return		null, if this line contains invalid data;
     *			an AbcData object holding the data, otherwise.
     */
    private AbcData parseAbcData(String s) {
	
	StringTokenizer st = new StringTokenizer(s);
	//  If number of tokens is not 4, we have a bad line.
	if (st.countTokens() < 4) {
	    // bad line
	    return null;
	}
	AbcData temp = new AbcData();
	try {
	    temp.altBase = Character.toUpperCase(st.nextToken().charAt(0));
	    temp.qv = Integer.parseInt(st.nextToken());
	    temp.position   = Integer.parseInt(st.nextToken());
	    temp.index = Integer.parseInt(st.nextToken());
	} catch (Exception e) {
	    return null;
	}
	if (temp.qv < 0 || temp.position < 0) {
	    return null;
	}

	if (temp.index < 0) {
	    temp.type = DELETION;
	} else {
	    temp.type = SUBSTITUTION;
	}
	return temp;
    }

    /** Attempts to release the resources. */
    public void releaseResources() {
	if (highestQVAbcData != null) {
	    for (int i = 0; i < highestQVAbcData.length; i++) {
	    	highestQVAbcData[i] = null;
	    }
	    highestQVAbcData = null;
	}
    }

    /** Returns the total number of alternative base calls. */
    public int getNumOfAbcs() { return numOfAbcs;   }

    /** Returns the number of highestQV alternative base calls. */
    public int getNumOfHighestQVAbcs() { return numOfHighestQVAbcs;   }

    /** Returns the total number of substitution alternative base calls. */
    public int getNumOfSubAbcs() { return numOfSubAbcs;   }

    /** Returns the total number of deletion alternative base calls. */
    public int getNumOfDelAbcs() { return numOfDelAbcs;   }

    /**
     * Identifies whether or not the bases were called by TraceTuner.
     * @return	true, if the bases are called by TraceTuner;
     *		false, otherwise.
     */
    public boolean  isTT()            { return isTT;        }

    /** Returns the version of the software used to call bases. */
    public String   getVersion()      { return version;     }

    /** Returns the chromat file name. */
    public String   getChromatFileName()  { return chromatFileName;}

    /** Returns true if there exists invalid data, such as: NUM_ABC does not
	agree with the acutal number of alternatice base calls presented in
	the file, NUM_SUBSTITUTIONS does not agree with the actual number of
	substitution abcs present in the file, etc;  false, otherwise. */
    public boolean  hasInvalidFileData() { return invalidFileData; }

    /** Returns true if any row of the abc data is invalid; false, otherwise.
     */
    public boolean  hasInvalidAbcData() { return invalidAbcData; }

    /** Returns the AbcData array that holds all highestQV abc data. */
    public AbcData[] getHighestQVAbcData() { return highestQVAbcData; }

    /**
     * main() function for testing only.  Normally this class wouldn't 
     * be stand-alone.
     */
    public static void main(String[] args) {
	TabFile p;
	int  i, numOfAbcs, numOfHighestQVAbcs;
	AbcData[] abcData;

	if(args.length != 1) {
	    System.out.println("Usage: java TabFile <filename>");
	    System.exit(0);
	}
	try {
	    p = new TabFile(args[0]);
	    numOfAbcs = p.getNumOfAbcs();
	    numOfHighestQVAbcs = p.getNumOfHighestQVAbcs();
	}
	catch (Exception e) {
	    System.out.println("Error: "+e.toString());
	    e.printStackTrace(System.out);
	    return;
	}

	System.out.println("Input file      is " + args[0]);
	System.out.println("num of abcs is " + numOfAbcs);
	System.out.println("num of highestQV abcs is " + numOfHighestQVAbcs);
	System.out.println("num of sub abcs is " + p.getNumOfSubAbcs());
	System.out.println("num of del abcs is " + p.getNumOfDelAbcs());
	System.out.println("version         is " + p.getVersion() + 
			   (p.isTT() ? "  (Tracetuner)" : "  (Phred)"));

	System.out.println("chromat_file    is " + p.getChromatFileName());

	abcData = p.getHighestQVAbcData();
	if (p.hasInvalidFileData()) {
	    System.out.println("has invalid file data");
	}
	if (p.hasInvalidAbcData()) {
	    System.out.println("has invalid abc data");
	}

	System.out.println("The first 30 highestQV abcs are: ");
	for (i=0; i < 30; i++) {
	    if (i >= numOfHighestQVAbcs) {
		break;
	    }
	    System.out.print(abcData[i].getAltBase()+"\t");
	    System.out.print(abcData[i].getQv()+"\t");
	    System.out.print(abcData[i].getPosition()+"\t");
	    System.out.print(abcData[i].getIndex()+"\t");
	    if (abcData[i].getType() == SUBSTITUTION) {
	        System.out.print("sub\n");
	    } else {
	        System.out.print("del\n");
	    }
	}

	System.out.println("The last 20 highestQV abcs are: ");
	for (i=numOfHighestQVAbcs - 20; i < numOfHighestQVAbcs; i++) {
	    if (i < 0) {
		break;
	    }
	    System.out.print(abcData[i].getAltBase()+"\t");
	    System.out.print(abcData[i].getQv()+"\t");
	    System.out.print(abcData[i].getPosition()+"\t");
	    System.out.print(abcData[i].getIndex()+"\t");
	    if (abcData[i].getType() == TabFile.SUBSTITUTION) {
	        System.out.print("sub\n");
	    } else {
	        System.out.print("del\n");
	    }
	}
    }


    /**
     * A public inner class of TabFile.
     * This class creates an AbcData object, which contains information
     * about one abc data in the alignment, including: the alt base,
     * the quality value, the position in the fragment, the index.
     */
    public class AbcData implements Cloneable {

	/** The alternative base code. */
	private char altBase;

	/** The quality value of this alternative base call. */
	private int qv;

	/** The position of the base. */
	private int position;

	/** The index of the base this alternative base was called on. */
	private int index;

	/** The type of this alternative base call: SUBSTITUTION or DELETION */
	private int type;

	/** Default constructor. */
	public AbcData() {}

	/** Constructor. */
	public AbcData(char altBase, int qv, int pos, int index, int type) {
	    this.altBase = altBase;
	    this.qv = qv;
	    this.position = pos;
	    this.index = index;
	    this.type = type;
	}

	/** Constructor. */
	public AbcData(AbcData a) {
	    this(a.altBase, a.qv, a.position, a.index, a.type);
	}

	/** Creates a clone. */
	public Object clone() {
	    return new AbcData(this);
	}

	/**
	 * Gets the alternative base code.
	 * @return the base code.
	 */
	public char getAltBase() { return altBase; }

	/** 
	 * Gets the position of this alternative base in the fragment sequence.
	 * @return the base position.
	 */
	public int getPosition() { return position; }

	/**
	 * Gets the quality value of this alternative base call.
	 * @return the quality value.
	 */
	public int getQv() { return qv; }

	/**
	 * Gets the index of the base on which this alternative base 
	 * was called.
	 * @return the index.
	 */
	public int getIndex() { return index; }

	/**
	 * Gets the type of this alternative base call: SUBSITUTION or
	 * DELETION.
	 * @return the type.
	 */
	public int getType() { return type; }
    }
}
