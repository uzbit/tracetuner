/* 1.8
 * TipFile.java - Reader for .tip files from TraceTuner.
 *              - GD Mar-2000
 *
 *   This class reads and holds in memory the data from a tip file. 
 *   The file contains information about all the intrinsic peaks       
 *   that TraceTuner computes from chematogram of a given sample file
 */
package com.paracel.tt.io;

import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;

public class TipFile
{
    private String file_name = "unknown";
    public  int num_peaks;
    public  IntrinsicPeak peak_list[];
    private long bytesRead;
    private long fileLength;
    private JProgressBar progressBar;
    boolean sleepNeeded;

    /** Total intrinsic signal peak list. */
    private IntrinsicPeak totalSignalPeaks[];
    private int maxX = 0;

    /* shiftX is the absolute value of the beginning x location of the
       first peak in the peak_list, if it is negative. */
    int shiftX;

    // Inner class 
    public class IntrinsicPeak implements Cloneable
    {
        public char base;
        public int type;
        public int length;
        public int x[]; 
        public int y[]; 
       
        // Constructor of inner class 

	public IntrinsicPeak() {}
        public IntrinsicPeak(char b, int typ, int beg, int sig[]) {
            base = b;  
            type   = typ;
	    length = sig.length;
            x      = new int[length];
            y      = new int[length];

	    for (int i = 0; i < length; i++) {
		x[i] = i + beg;
	    }
	    System.arraycopy(sig, 0, y, 0, length);
        }
	public IntrinsicPeak(IntrinsicPeak target) {
	    this(target.base, target.type, target.x[0], target.y);
	}

	public Object clone() {
	    return new IntrinsicPeak(this);
	}

	public void dispose() {
	    x = null;
	    y = null;
	}

        // Methods of the inner class
        public int getBase()     { return base;   }
        public int getType()      { return type;    }
        public int getLength()    { return length;  }
        public int[] getXValues()    { return x;  }
        public int[] getYValues()    { return y;  }
    }

    // Constructor of TipFile class 
    public TipFile(String name, JFrame frame)
	throws FileNotFoundException, IOException {
	JDialog progressWindow = null;
	if (frame != null) {
	/* create the progress bar to inform the user the progress of
	   tip file loading.  The reason is that tip files usually take
	   long time to load. */
	    progressWindow = new JDialog(frame, "Loading tip file");
	    progressWindow.setSize(300, 50);

	    UIManager.put("ProgressBar.selectionBackground", Color.black);
	    UIManager.put("ProgressBar.selectionForeground", Color.white);
	    UIManager.put("ProgressBar.foreground", new Color(8, 32, 128));

	    progressBar = new JProgressBar();
	    progressBar.setMinimum(0);
	    progressBar.setStringPainted(false);
	    progressBar.setSize(300, 20);

	    progressWindow.getContentPane().add(progressBar, 
						BorderLayout.CENTER);

	    // locate the progress window at the center of the frame
	    int x1 = (int) frame.getLocation().getX();
	    int y1 = (int) frame.getLocation().getY();
	    int w1 = frame.getWidth();
	    int h1 = frame.getHeight();
	    int w2 = progressWindow.getWidth();
	    int h2 = progressWindow.getHeight();
	    int x2 = x1 + (w1 - w2) / 2;
	    int y2 = y1 + (h1 - h2) / 2;
	    progressWindow.setLocation(x2, y2);
	    progressWindow.setResizable(false);
	    progressWindow.setVisible(true);

	    /* The sleeping is the work-around added for a JDK JProgressBar
	       bug:  If the progressing is too fast, the progress bar fails
	       to show the progress. */
	    String osName = System.getProperty("os.name");
	    if (osName != null) {
	    	osName = osName.toLowerCase().trim();
	    	if (osName.indexOf("windows") >= 0) {
		    sleepNeeded = true;
	    	} else {
		    sleepNeeded = false;
	    	}
	    }
	}

	RandomAccessFile f = new RandomAccessFile(name,"r");
	ArrayList peaks = new ArrayList(3000);

        int num_lines = 0;
        char b = 'N';
        int  typ=0, beg=0, len=0, ipos=0, num= 0, sig[];

	int counter = 0;
	bytesRead = 0;
	fileLength = f.length();
	if (progressBar != null) {
	    progressBar.setMaximum((int) fileLength);
	}

	int percentageComplete = 0;

	//  Read the tip file header 
	String s;
	String line = f.readLine();

	if (line == null || !line.startsWith(">")) {
	    throw new IOException(name + 
			" does not seem to be a valid .tip file");
        } else {
	    bytesRead += line.length();
            line = line.trim();
            num_lines++;
        }

	// Read intrinsic peaks
 	/* For each intrinsic peak, two lines of output will be produced. The 
     * first line contains the base character corresponding to the peak,
	 * intrinsic position, the number of integer values in the second
	 * line, apparent peak boundaries (beginning and end), and the peak
	 * type. The second line contains the y-values for the
	 * right half of intrinsic peak, rounded to integer, starting from
	 * the intrinsic peak position and ending with the first position
	 * where the rounded value is zero.
	 */
	int i = 0;
	line = f.readLine();
	bytesRead += line.length();
	if (line != null) {
	    line = line.trim();
	}
	while (line != null) {
            StringTokenizer st = new StringTokenizer(line);
	    // the data of each peak composes of two lines.
	    if (i % 2 == 0) {    
		// 1st line of intrinsic peak

                //  The first token is the peak base         
                if (st.hasMoreTokens()) {
                    b = st.nextToken().charAt(0);
                }
                try {
                    if (st.hasMoreTokens()) {
			//intrinsic peak position
                        ipos = Integer.parseInt(st.nextToken());
                    }
                    if (st.hasMoreTokens()) {
			// num of datapoints on the second line
                        num = Integer.parseInt(st.nextToken());
			len = Math.max(0, 2 * num - 1);
                    }
		    beg = ipos - num + 1;
		    // skip the next three fields
                    if (st.hasMoreTokens()) {
			// the apparent beginning of the intrinsic peak
			st.nextToken();
		    }
                    if (st.hasMoreTokens()) {
			// the apparent end of the intrinsic peak
			st.nextToken();
		    }
                    if (st.hasMoreTokens()) {
            // peak resolution
            st.nextToken();
            }
                    if (st.hasMoreTokens()) {
			// the type of the intrinsic peak
                        typ = Integer.parseInt(st.nextToken());
                    }
                }
                catch (Exception e) {
		    // Just continue if we hit something weird
		}
            } else {
                sig = new int[len];  
                for (int j = num - 1; j >= 0; j--) { 
                    try {
                        if (st.hasMoreTokens()) {
                            sig[j] = Integer.parseInt(st.nextToken());
			    sig[len - 1 - j] = sig[j];
                        }
                    }
                    catch (Exception e) { 
			// Just continue if we hit something weird
		    }
                }
                peaks.add(new IntrinsicPeak(b, typ, beg, sig));

		// find the maximum x location value
		if (maxX < (beg + sig.length - 1)) {
		    maxX = beg + sig.length - 1;
		}
            }
            i++;
	    line = f.readLine();
            if (line != null) {
		if (progressBar != null) {
		    bytesRead += line.length();
		    if (bytesRead * 100 / fileLength >= counter) {
		    	setProgressValue(bytesRead);
		    	counter += 10;
		    }
		}
                line = line.trim();
            }
        }
	f.close();
	
	// Convert the ArrayList to an IntrinsicPeak array
	num_peaks = peaks.size();
	peak_list = new IntrinsicPeak[num_peaks];
	for (i = 0; i < num_peaks; i++) {
	    peak_list[i] = (IntrinsicPeak) (peaks.get(i));
	}
	peaks = null;

	/* The rest of the constructor is to compute the total intrinsic
	   signal peaks. */

	/* Sometimes the first intrinsic peak begins with a negative x
	   location.  In this case, we expand the length of the int
	   arrays in order to hold all intrinsic signal data. */
	IntrinsicPeak firstPeak = peak_list[0];
	shiftX = (firstPeak.x[0] < 0) ? -firstPeak.x[0] : 0;
	int maxLen = maxX + 1 + shiftX;
	int[] ASignals = new int[maxLen];
	int[] CSignals = new int[maxLen];
	int[] GSignals = new int[maxLen];
	int[] TSignals = new int[maxLen];

	/* Calculates the total intrinsic signal value of each base type 
	   at each x point, and constructs one int array for each base type
	   to hold the total intrinsic signal value at all x points for
	   this base type. */
	IntrinsicPeak peak;
	int x;
	for (i = 0; i < peak_list.length; i++) {
	    peak = peak_list[i];
	    for (int j = 0; j < peak.x.length; j++) {
		x = peak.x[j] + shiftX;
		if (x < 0)  {
		    continue;
		}
		if (x >= maxLen)  {
		    break;
		}
	    	switch (peak.base) {
	    	case 'A':
		    ASignals[x] += peak.y[j]; break;
	    	case 'C':
		    CSignals[x] += peak.y[j]; break;
	    	case 'G':
		    GSignals[x] += peak.y[j]; break;
	    	case 'T':
		    TSignals[x] += peak.y[j]; break;
	    	}
	    }
	}

	/* Then from the four total intrinsic signal value arrays, we
	   compute the total intrinsic signal peak array.  The resulted
	   array is not ordered.
	   Since the number of total intrinsic signal peaks is no greater
	   than that of the individual intrinsic peaks.  We first use
	   an IntrinsicPeak array of the num of individual intrinsic
	   peaks to hold the computed total intrinsic signal peaks. */
	IntrinsicPeak[] totalSPeaks = new IntrinsicPeak[peak_list.length];
	IntrinsicPeak temp;
	int index = 0, numOfPeaksFound = 0;
	int AFromIndex = 0, CFromIndex = 0, GFromIndex = 0, TFromIndex = 0;
	while (true) {
	    numOfPeaksFound = 0;
	    temp = findPeak('A', ASignals, AFromIndex);
	    if (temp != null) {
		totalSPeaks[index++] = temp;
		AFromIndex = temp.x[temp.x.length - 1] + shiftX;
	    	numOfPeaksFound++;
	    } else {
		// No more peaks of A can be found.  By setting AFromIndex
		// to the length of ASignals, we guarantee the next findPeak
		// action will return right away with null value.
		AFromIndex = ASignals.length;
	    }
	    temp = findPeak('C', CSignals, CFromIndex);
	    if (temp != null) {
		totalSPeaks[index++] = temp;
		CFromIndex = temp.x[temp.x.length - 1] + shiftX;
	    	numOfPeaksFound++;
	    } else {
		CFromIndex = CSignals.length;
	    }
	    temp = findPeak('G', GSignals, GFromIndex);
	    if (temp != null) {
		totalSPeaks[index++] = temp;
		GFromIndex = temp.x[temp.x.length - 1] + shiftX;
	    	numOfPeaksFound++;
	    } else {
		GFromIndex = GSignals.length;
	    }
	    temp = findPeak('T', TSignals, TFromIndex);
	    if (temp != null) {
		totalSPeaks[index++] = temp;
		TFromIndex = temp.x[temp.x.length - 1] + shiftX;
	    	numOfPeaksFound++;
	    } else {
		TFromIndex = TSignals.length;
	    }
	    if (numOfPeaksFound == 0) {
		// No more peaks of any base type can be found.
		break; 
	    }
	    if (index >= totalSPeaks.length) {
		/* Normally the index should not exceed the length
		   of totalSPeaks.  This condition is added to ensure
		   that index will not exceed this limit, so to avoid
		   ArrayIndexOutOfBoundsException, if it does happen. */
		break; 
	    }
	}

	totalSignalPeaks = new IntrinsicPeak[index];
	System.arraycopy(totalSPeaks, 0, totalSignalPeaks, 0, index);
	totalSPeaks = null;

	if (frame != null) {
	    setProgressValue(fileLength);
	    progressWindow.setVisible(false);
	    progressWindow.dispose();
	    progressWindow = null;
	    System.gc();
	}
    }

    /* Returns the first total intrinsic signal peak found in the specified
       total intrinsic signal value array starting from the specified array
       index; returns null if no peak found. */
    protected IntrinsicPeak findPeak(char base, int[] signals, int fromIndex) {
	if (fromIndex < 0 || signals == null 
	        || fromIndex >= signals.length-1) {
	    return null;
	}
	int start = -1, end = -1;
	if (signals[fromIndex] > 0) {
	    start = fromIndex; 
	} else {
	    for (int i = fromIndex; i < signals.length - 1; i++) {
	   	if (signals[i] == 0 && signals[i+1] > 0) {
	            start = i;
	       	    break;
	   	}
	    }
	}
	if (start == -1) {
	    return null; 
	}
	for (int i = start + 1; i < signals.length - 1; i++) {
	   if (signals[i] > 0 && signals[i+1] == 0) {
	       end = i + 1;
	       break;
	   } else {
	       /* Enforce the size of the total signal peak to be less than
	          100, because too wide peaks slow down drawing */
	       if (i - start + 1 >= 100) {
		   end = i;
		   break;
	       }
	   }
	}
	if (end == -1 && signals[signals.length-1] > 0) {
	    end = signals.length - 1;
	}
	if (end < start) {
	    return null; 
	}
	int[] sig = new int[end - start + 1];
	System.arraycopy(signals, start, sig, 0, sig.length);
	return (new IntrinsicPeak(base, 0, start - shiftX, sig));
    }

    protected void sleep(long ms) {
	try {
	    Thread.currentThread().sleep(ms);
	} catch (Exception e) {
	    System.err.println("TipFile--Constructor" + e.toString());
	    e.printStackTrace(System.err);
	}
    }

    protected void setProgressValue(final long value) {
	if (progressBar == null) {
	    return;
	}
	if (sleepNeeded) {
	    sleep(1);
	}
	progressBar.setValue((int)value);
    }
    
    /** Attempts to release the resources. */
    public void releaseResources() {
	for (int i = 0; i < peak_list.length; i++) {
	     if (peak_list[i] != null) {
		 peak_list[i].dispose();
		 peak_list[i] = null;
	     }
	}
	peak_list = null;

	for (int i = 0; i < totalSignalPeaks.length; i++) {
	     if (totalSignalPeaks[i] != null) {
		 totalSignalPeaks[i].dispose();
		 totalSignalPeaks[i] = null;
	     }
	}
	totalSignalPeaks = null;
    }

    public String   getName()         { return file_name;   }
    public int      getNumPeaks()     { return num_peaks;   }
    public IntrinsicPeak[]      getPeaks()     { return peak_list;   }
    public IntrinsicPeak[] getTotalSignalPeaks() { return totalSignalPeaks; }
    public int getNumTotalSignalPeaks() { return totalSignalPeaks.length; }

    /*
     *  The following main routine is for testing only.  Normally
     *  this class wouldn't be stand-alone.
     */
    public static void main(String[] args)  // for testing only
    {
	TipFile t;
	int  i, num_peaks; 
        char b;
	int  typ, beg, len;
	int  sig[];              

	if(args.length != 1) {
	    System.out.println("Usage: java TipFile <filename>");
	    System.exit(0);
	}
	try {
	    t = new TipFile(args[0], null);
	}
	catch (Exception e) {
	    System.err.println("Error: "+e.getMessage());
	    e.printStackTrace(System.err);
	    return;
	}

	System.out.println("Input file is: " + args[0]);
	System.out.println("Number of peaks is " + t.getNumPeaks());
	System.out.println("Number of total intrinsic signal peaks is "
			   + t.getNumTotalSignalPeaks());

    }
}
