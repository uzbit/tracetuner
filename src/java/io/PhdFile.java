/* 1.5
 * PhdFile.java - Reader for .phd files as from Phred or TraceTuner.
 *              - MFM Oct-99
 *
 *   This class reads and holds in memory the data from a phd file.  All
 *   get operations are done from memory so they should be fairly quick
 *   once the file is read in.  This class will not notice if the file is
 *   deleted once the class is instantiated.
 */

package com.paracel.tt.io;

import java.io.*;
import java.util.*;

public class PhdFile
{
    private static final String MIXED_BASES = "RYKMSW";

    private boolean isTT         = false;
    private String  version      = "unknown";
    private String  call_method  = "unknown";
    private String  chromat_file = "unknown";
    private int     trace_min    = 0;
    private int     trace_max    = 0; 
    private int     num_bases    = 0;
    private char[]  bases        = null;
    private int[]   qv_array     = null;
    private int[]   peak_array   = null;
    private int     numSNPs    = 0;
    private char[]  hetBases        = null;
    private int[]   hetIndices     = null;
    private int[]   hetPositions   = null;
    private int[]   hetQVs   = null;

    /** Left trim point.  Qualified region starts from this base index. 
        (the base indices starts from 0) */
    private int leftTrimPoint = -2;

    /** Right trim point.  Qualified region ends at this base index. 
        (the base indices starts from 0) */
    private int rightTrimPoint = -2;

    /** Trim quality value threshold. */
    private int trimThreshold = -1;

    /** Whether the trim Data were found. */
    private boolean trimDataFound;

    public PhdFile(String name)
	throws FileNotFoundException, IOException
    {
	RandomAccessFile f = new RandomAccessFile(name,"r");

	/*
         *  Read the phd file pre-amble and pick up as many of the fields
	 *  we're interested in as are present.
	 */
	String s;
	String line = f.readLine().trim();

	while(line != null && !line.equals("BEGIN_DNA")) {

	    if((s = parse(line, "CHROMAT_FILE:")) != null)
		chromat_file = s;
	    if((s = parse(line, "PHRED_VERSION:")) != null) {
		version = s;
		if(version.substring(0,2).equals("TT"))
		    isTT = true;
		else
		    isTT = false;
	    }
	    if((s = parse(line, "TT_VERSION:")) != null) {
		version = "TT_"+s;
		isTT = true;
	    }
	    /* QV_VERSION was used prior to TT release 1.0 */
	    if((s = parse(line, "QV_VERSION:")) != null) {
		version = s;
		isTT = true;
	    }
	    if((s = parse(line, "CALL_METHOD:")) != null)
		call_method = s;
	    if((s = parse(line, "TRACE_ARRAY_MIN_INDEX:")) != null)
		trace_min = Integer.parseInt(s);
	    if((s = parse(line, "TRACE_ARRAY_MAX_INDEX:")) != null)
		trace_max = Integer.parseInt(s);

	    if((s = parse(line, "TRIM:")) != null) {
		trimDataFound = true;
	    	StringTokenizer st = new StringTokenizer(s);
		int num = st.countTokens();
		if (num >= 3) {
		    try {
		    	leftTrimPoint = Integer.parseInt(st.nextToken());
		    	rightTrimPoint = Integer.parseInt(st.nextToken());
		    	double d = Double.parseDouble(st.nextToken());
		    	trimThreshold = (int) -Math.round(
						(Math.log(d)/Math.log(10))*10);
		    } catch (NumberFormatException e) {
		    }
		}
	    }

	    line = f.readLine().trim();
	}
	if(line == null) 
	    throw new IOException(name + 
			" does not seem to be a valid .phd file");

	/*
         *  Read the phd file base calls, quality values, and peak locations.
	 *  Make one pass through to see how many there are, reset the file
	 *  location, and make a second pass through to fetch them.
	 */
	long dna_section = f.getFilePointer();  // save this file location
	line = f.readLine().trim();
	char c;
	while(line != null && !line.equals("END_DNA")) {
	    num_bases++;
	    if (line.length() > 0) {
	    	c = Character.toUpperCase(line.charAt(0));
		if (MIXED_BASES.indexOf(c) >= 0) {
		    numSNPs++; 
		}
	    }
	    line = f.readLine().trim();
	}

	// Allocate memory for now known number of bases
	bases      = new char[num_bases];
	qv_array   = new int[num_bases];
	peak_array = new int[num_bases];

	if (numSNPs > 0) {
	    hetBases = new char[numSNPs];
	    hetIndices = new int[numSNPs];
	    hetPositions = new int[numSNPs];
	    hetQVs = new int[numSNPs];
	}

	// Restore file location
	f.seek(dna_section);

	// Read bases
	int i = 0;
	int j = 0;
	line = f.readLine();
	if (line != null) {
	    line.trim();
	}
	boolean isSNP;
	while(line != null && !line.equals("END_DNA") && i < num_bases) {
	    isSNP = false;
	    StringTokenizer st = new StringTokenizer(line);
	    //  If number of tokens in not 3, we have a bad line.
	    //  Fill in as much as we can
	    if(st.hasMoreTokens()) {
		bases[i] = Character.toUpperCase(st.nextToken().charAt(0));
		if (MIXED_BASES.indexOf(bases[i]) >= 0) {
		    isSNP = true;
		    hetBases[j] = bases[i];
		    hetIndices[j] = i + 1;
		}
	    } else {
		bases[i]      = 'N';
	    }
	    try {
		if(st.hasMoreTokens()) {
		    qv_array[i]   = Integer.parseInt(st.nextToken());
		    if (isSNP) {
			hetQVs[j] = qv_array[i];
		    }
		}
		if(st.hasMoreTokens()) {
		    peak_array[i] = Integer.parseInt(st.nextToken());
		    if (isSNP) {
			hetPositions[j] = peak_array[i];
		    }
		}
	    }
	    catch (Exception e) { } // Just continue if we hit something weird
	    line = f.readLine();
	    if (line != null) {
		line.trim();
	    }
	    i++;
	    if (isSNP) {
	    	j++;
	    }
	}

	f.close();  // That's all folks!

	if (j < numSNPs) {
	    // marks the end of the valid het data in the het arrays.
	    hetBases[j] = '-';
	    numSNPs = j;
	}
    }

    private String parse(String s, String pat) {

	if(!s.startsWith(pat))
	    return null;
	s = s.substring(pat.length()+1,s.length());
	return s.trim();
    }


    public int      getNumBases()     { return num_bases;   }
    public boolean  isTT()            { return isTT;        }
    public String   getVersion()      { return version;     }
    public String   getCallMethod()   { return call_method; }
    public String   getChromatFile()  { return chromat_file;}
    public int      getTraceMin()     { return trace_min;   }
    public int      getTraceMax()     { return trace_max;   }
    /*
     *  Note that the following three array objects are returned directly.
     *  The calling code can write into the arrays and change the in-memory
     *  copy.  Subsequent calls to getBases() will return an altered result
     *  if this is done.  This is either a feature or a bug depending on
     *  your viewpoint.
     */
    public char[] getBases()          { return bases;       }
    public int[]  getQvValues()       { return qv_array;    }
    public int[]  getPeakLocations()  { return peak_array;  }

    public int getNumOfSNPs() { return numSNPs; }
    public char[] getSnpBases() { return hetBases; }
    public int[] getSnpIndices() { return hetIndices; }
    public int[] getSnpPositions() { return hetPositions; }
    public int[] getSnpQVs() { return hetQVs; }

    public int[] getTrimData() {
	if ((!trimDataFound) || leftTrimPoint <= -2 
	        || rightTrimPoint <= -2 || trimThreshold <= -1 
	        || leftTrimPoint >= num_bases || rightTrimPoint >= num_bases) {
	    // invalid trim data detected.
	    return null;
	} else {
	    int[] res = { leftTrimPoint, rightTrimPoint, trimThreshold };
	    return res;
	}
    }

    public void releaseResources() {
	bases = null;
	qv_array = null;
	peak_array = null;
	hetBases = null;
	hetIndices = null;
	hetPositions = null;
	hetQVs = null;
    }

    /*
     *  The following main routine is for testing only.  Normally
     *  this class wouldn't be stand-alone.
     */
    public static void main(String[] args)  // for testing only
    {
	PhdFile p;
	int  i, numbases;
	char bases[];
	int  qv_values[];
	int  peak_locations[];

	if(args.length != 1) {
	    System.out.println("Usage: java PhdFile <filename>");
	    System.exit(0);
	}
	try {
	    p = new PhdFile(args[0]);
	}
	catch (Exception e) {
	    System.err.println("Error: "+e.getMessage());
	    return;
	}

	System.out.println("Input file      is " + args[0]);
	System.out.println("number of bases is " + p.getNumBases());
	System.out.println("version         is " + p.getVersion() + 
			   (p.isTT() ? "  (Tracetuner)" : "  (Phred)"));

	System.out.println("chromat_file    is " + p.getChromatFile());
	System.out.println("call method     is " + p.getCallMethod());
	System.out.println("trace min       is " + p.getTraceMin());
	System.out.println("trace max       is " + p.getTraceMax());
	System.out.println("left trim       is " + p.getTrimData()[0]);
	System.out.println("right trim      is " + p.getTrimData()[1]);
	System.out.println("trim threshold  is " + p.getTrimData()[2]);

	bases = p.getBases();
	qv_values = p.getQvValues();
	peak_locations = p.getPeakLocations();
	numbases = p.getNumBases();

	System.out.println("Bases are: " + bases);
	System.out.println("First 20 Qv values are: ");
	for(i=0; i < 20; i++)
	    System.out.print(""+qv_values[i]+" ");
	System.out.println();
	System.out.println("First 20 Peak locations are: ");
	for(i=0; i < 20; i++)
	    System.out.print(""+peak_locations[i]+" ");
	System.out.println();
    }
}
