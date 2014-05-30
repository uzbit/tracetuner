/*
 * 1.12
 *
 * @(#)ChromSnapShot.java       1.0     2000/10/25
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.lang.*;
import java.awt.*;
import com.paracel.tt.chrom.*;
import com.paracel.tt.util.*;
import com.paracel.tt.io.*;

/**
 * This class creates a snapshot of the data involved in the panels of the
 * TTView window, such as the traces, quality values, peak locations, signal
 * values, file names, etc.  The snapshots are created for printing purpose.
 * The data are all made deep copies to ensure the printing will not be
 * affected while the user operates on the viewer.
 */
public class ChromSnapShot {
    
    /** The type of the print job for which this snapshot was created. */
    private int type;

    /** The product name. */
    private String productName;

    /** The sample file name and the phd file names. */
    private String chromFileName, phdFileName, phd2FileName;

    /** The trace lines and the analyzed trace lines. */
    private Trace[] traceLines, analyzedTraceLines;

    /** The phd base peak locations for 2 phd files. */
    private int[] phdPeakLocations, phdPeakLocations2;

    /** The quality values for 2 phd files. */
    private int[] qualityValues, qualityValues2;

    /** The phd base calls for 2 phd files. */
    private char[] phdBases, phdBases2;

    /** Whether the base calls are from TraceTuner for 2 phd files. */
    private boolean isTT, isTT2;

    /** The original base calls. */
    private String origBases;

    /** The original base call locations. */
    private short[] origBasePeakLocations;

    /** The trace lines length and the analyzed trace lines length */
    private int traceLength, analyzedTraceLength;

    /** The alternative base calls data. */
    private TabFile.AbcData[] abcData;

    /** The alignment data. */
    private TalFile.AlignData[] alignData;

    /** The intrinsic peak lines and the total intrinsic signal peak lines. */
    private TipFile.IntrinsicPeak[] intrinsicPeakLines, totalISPeakLines;

    /** The print view chrom layout and print all chrom layout. */
    private ChromLayout pvChromLayout, paChromLayout;

    /** The viewer window width and height. */
    private int viewWidth, viewHeight;

    /** The view position of the viewer window. */
    private Point viewPos;

    /** The closest peak mark. */
    private int[] markClosestPeak;

    /** The flags which determine display feature, therefore determine
	what to print. */
    private boolean displayRaw, displayOrigBases, displayAbc, displayAln, 
		    displayTip, displayTotalIS, displayTrim, displayPhd2,
		    signalValuesShown, fifthDyeDataShown;

    /** The signal values. */
    private int[] signalValues;

    /** The zoom value on x axis. */
    private float zoomX;

    /** The trim points and threshold. */
    private int leftTrimPoint, rightTrimPoint, trimThreshold;

    /** Whether sample file or phd file is invalid. */
    private boolean sampleOrPhdInvalid;

    /** The index of the sample file of this ChromSnapShot, and the total
	number of input files selected by the Viewer. */
    private int currFileIndex, totalNumOfFiles;

    /**
     * Creates a snapshot at the upper and bottom panels of the viewer
     * window according to the print job type.
     * @param	type	the print job type (TTPrintJob.PRINT_LOCAL,
     *			TTPrintJob.PRINT_SINGLE_GLOBAL or
     *			TTPrintJob.PRINT_ALL_GLOBAL).
     * @param	cp	the chrom view panel (upper panel) in the viewer.
     * @param	mp	the message panel (bottom panel) in the viewer.
     */
    public ChromSnapShot(int type, ChromViewPanel cp, MessagePanel mp) {
	this.type = type;

	int length = cp.getPhdPeakLocations().length;
	phdPeakLocations = new int[length];
	System.arraycopy(cp.getPhdPeakLocations(), 0, 
			 phdPeakLocations, 0, length);

	length = cp.getQualityValues().length;
	qualityValues = new int[length];
	System.arraycopy(cp.getQualityValues(), 0, 
			 qualityValues, 0, length);

	length = cp.getPhdBases().length;
	phdBases = new char[length];
	System.arraycopy(cp.getPhdBases(), 0, 
			 phdBases, 0, length);

	isTT = cp.isTT();

	String temp = cp.getCurrChromFileName();
	chromFileName = temp.substring(
			    temp.lastIndexOf(Constants.FILE_SEPARATOR) + 1, 
			    temp.length());
	temp = cp.getCurrPhdFileName();
	phdFileName = temp.substring(
			    temp.lastIndexOf(Constants.FILE_SEPARATOR) + 1, 
			    temp.length());
	temp = cp.getCurrPhd2FileName();
	if (temp == null) {
	    phd2FileName = null;
	} else {
	    phd2FileName = temp.substring(
			    temp.lastIndexOf(Constants.FILE_SEPARATOR) + 1, 
			    temp.length());
	}
	productName = new String(cp.getProductName());
	displayPhd2 = cp.displayingPhd2();
	
	if (displayPhd2) {
	    length = cp.getPhdPeakLocations2().length;
	    phdPeakLocations2 = new int[length];
	    System.arraycopy(cp.getPhdPeakLocations2(), 0, 
			     phdPeakLocations2, 0, length);

	    length = cp.getQualityValues2().length;
	    qualityValues2 = new int[length];
	    System.arraycopy(cp.getQualityValues2(), 0, 
			     qualityValues2, 0, length);

	    length = cp.getPhdBases2().length;
	    phdBases2 = new char[length];
	    System.arraycopy(cp.getPhdBases2(), 0, 
			     phdBases2, 0, length);

	    isTT2 = cp.isTT2();
	}
	
	displayAbc = cp.displayingAbc();
	if (displayAbc) {
	    length = cp.getAbcData().length;
	    abcData = new TabFile.AbcData[length];
	    for (int i = 0; i < length; i++) {
		abcData[i] = (TabFile.AbcData) (cp.getAbcData()[i].clone());
	    }
	}

	displayTrim = cp.displayingTrim();
	if (displayTrim) {
	    leftTrimPoint = cp.getLeftTrimPoint();
	    rightTrimPoint = cp.getRightTrimPoint();
	    trimThreshold = cp.getTrimThreshold();
	}

	if (type == TTPrintJob.PRINT_SINGLE_GLOBAL) {
	    paChromLayout = ChromLayout.createPrintAllLayout(displayAbc,
						     displayPhd2,
						     TTPrintJob.ROW_HEIGHT);

	    Trace[] analyzedTraces = cp.getAnalyzedTraces();
	    int analyzedTraceMaxIntensity = cp.getAnalyzedTraceMaxIntensity();
	    analyzedTraceLength = cp.getAnalyzedTraceLength();

	    length = analyzedTraces.length;
	    analyzedTraceLines = new Trace[length];

	    // Convert intensity traces to y pixel traces.
	    float yToPixelRatio = ((float)analyzedTraceMaxIntensity) 
	    			/ ((float)paChromLayout.chromHeight);

	    int tmp;
	    int[] yPoints;
	    for (int i = 0; i < length; i++) {
	    	yPoints = new int[analyzedTraces[i].getLength()];
	    	for (int j = 0; j <yPoints.length; j++) {
		    tmp = (int) (((float)analyzedTraces[i]
				   .getIntensityAt(j)) / yToPixelRatio);
		    if (tmp > paChromLayout.chromHeight) {
		    	tmp = paChromLayout.chromHeight;
		    }
		    yPoints[j] = paChromLayout.traceLinesY - tmp;
	    	}
	    	analyzedTraceLines[i] = new Trace(analyzedTraces[i].getName(),
	    			      		  analyzedTraces[i].getColor(),
						  yPoints);
	    }
	    return;
	}

	// The rest data are only needed for PRINT_LOCAL.
	viewPos = new Point(cp.getCurrViewPos());

	viewWidth = cp.getCurrViewWidth();
	viewHeight = cp.getCurrViewHeight();
	traceLength = cp.getTraceLength();

	zoomX = cp.getZoomX();

	displayRaw = cp.displayingRaw();
	displayOrigBases = cp.displayingOrigBases();
	displayAln = cp.displayingAln();
	displayTip = cp.displayingTip();
	displayTotalIS = cp.displayingTotalIS();
	markClosestPeak = new int[2];
	markClosestPeak[0] = cp.getMarkClosestPeak()[0];
	markClosestPeak[1] = cp.getMarkClosestPeak()[1];

	pvChromLayout = (ChromLayout) cp.getChromLayout().clone();

	origBases = new String(cp.getOrigBases());

	length = cp.getTraceLines().length;
	traceLines = new Trace[length];
	for (int i = 0; i < length; i++) {
	    traceLines[i] = (Trace) cp.getTraceLines()[i].clone();
	}

	if (displayOrigBases) {
	    length = cp.getOrigBasePeakLocations().length;
	    origBasePeakLocations = new short[length];
	    System.arraycopy(cp.getOrigBasePeakLocations(), 0, 
			     origBasePeakLocations, 0, length);
	}

	if (displayAln) {
	    length = cp.getAlignData().length;
	    alignData = new TalFile.AlignData[length];
	    for (int i = 0; i < length; i++) {
		alignData[i] = (TalFile.AlignData)
		    		(cp.getAlignData()[i].clone());
	    }
	}
	
	if (displayTip) {
	    length = cp.getIntrinsicPeakLines().length;
	    intrinsicPeakLines = new TipFile.IntrinsicPeak[length];
	    for (int i = 0; i < length; i++) {
		intrinsicPeakLines[i] = (TipFile.IntrinsicPeak)
		    		(cp.getIntrinsicPeakLines()[i].clone());
	    }
	}

	if (displayTotalIS) {
	    length = cp.getTotalISPeakLines().length;
	    totalISPeakLines = new TipFile.IntrinsicPeak[length];
	    for (int i = 0; i < length; i++) {
		totalISPeakLines[i] = (TipFile.IntrinsicPeak)
		    		(cp.getTotalISPeakLines()[i].clone());
	    }
	}

	// get data from MessagePanel
	signalValuesShown = mp.signalValuesShown();
	fifthDyeDataShown = mp.fifthDyeDataShown();
	if (signalValuesShown) {
	    int[] tmp = mp.getSignalValues();
	    signalValues = new int[tmp.length];
	    System.arraycopy(tmp, 0, signalValues, 0, tmp.length);
        }
    }

    /**
     * Creates a snapshot from the specified TTPrintJob type and the specified
     * FilesStatus object.
     */
    public ChromSnapShot(int type, FilesStatus status) {
	this.type = type;
	currFileIndex = status.getIndex();
	totalNumOfFiles = status.getNumOfFiles();
	if (status.getSampleFileStatus() != FilesStatus.VALID
	        || status.getPhdFileStatus() != FilesStatus.VALID) {
	    sampleOrPhdInvalid = true;
	    return;
	}
	ABIFile chromFile = status.getSampleFile();
	PhdFile phdFile = status.getPhdFile();

	phdPeakLocations = phdFile.getPeakLocations();
	qualityValues = phdFile.getQvValues();
	phdBases = LangTools.toUpperCase(phdFile.getBases());

	String temp = status.getSampleFileName();
	chromFileName = temp.substring(
			    temp.lastIndexOf(Constants.FILE_SEPARATOR) + 1, 
			    temp.length());
	temp = status.getPhdFileName();
	phdFileName = temp.substring(
			    temp.lastIndexOf(Constants.FILE_SEPARATOR) + 1, 
			    temp.length());
	productName = "Paracel TraceTuner Version " + phdFile.getVersion();
	displayPhd2 = false;
	displayAbc = false;
	displayTrim = false;

	paChromLayout = ChromLayout.createPrintAllLayout(displayAbc,
						     displayPhd2,
						     TTPrintJob.ROW_HEIGHT);

	int chromType = AbstractChromatogram.getChromType(chromFile);

	Chromatogram chromatogram;
	if (chromType == AbstractChromatogram.ABI) {
	    chromatogram = new ABIChromatogram(chromFile);
	} else {
	    chromatogram = new SCFChromatogram(chromFile);
	}

	Trace[] analyzedTraces = chromatogram.getAnalyzedTraces();
	int analyzedTraceMaxIntensity = chromatogram.getMaxIntensity();
	analyzedTraceLength = chromatogram.getLength();

	int length = analyzedTraces.length;
	analyzedTraceLines = new Trace[length];

	// Convert intensity traces to y pixel traces.
	float yToPixelRatio = ((float)analyzedTraceMaxIntensity) 
	    			/ ((float)paChromLayout.chromHeight);

	int tmp;
	int[] yPoints;
	for (int i = 0; i < length; i++) {
	    yPoints = new int[analyzedTraces[i].getLength()];
	    for (int j = 0; j <yPoints.length; j++) {
		tmp = (int) (((float)analyzedTraces[i]
				   .getIntensityAt(j)) / yToPixelRatio);
		if (tmp > paChromLayout.chromHeight) {
		    	tmp = paChromLayout.chromHeight;
		}
		yPoints[j] = paChromLayout.traceLinesY - tmp;
	    }
	    analyzedTraceLines[i] = new Trace(analyzedTraces[i].getName(),
	    			      		  analyzedTraces[i].getColor(),
						  yPoints);
	}
    }

    /** Releases the memory resource used by this ChromSnapShot object. */
    public void dispose() {
	if (traceLines != null) {
	    for (int i = 0; i < traceLines.length; i++) {
		traceLines[i].dispose();
		traceLines[i] = null;
	    }
	    traceLines = null;
	}
	if (analyzedTraceLines != null) {
	    for (int i = 0; i < analyzedTraceLines.length; i++) {
		analyzedTraceLines[i].dispose();
		analyzedTraceLines[i] = null;
	    }
	    analyzedTraceLines = null;
	}
	phdPeakLocations = null;
	qualityValues = null;
	phdBases = null;
	phdPeakLocations2 = null;
	qualityValues2 = null;
	phdBases2 = null;
	origBases = null;
	origBasePeakLocations = null;
	if (abcData != null) {
	    for (int i = 0; i < abcData.length; i++) {
		abcData[i] = null;
	    }
	    abcData = null;
	}
	if (alignData != null) {
	    for (int i = 0; i < alignData.length; i++) {
		alignData[i] = null;
	    }
	    alignData = null;
	}
	if (intrinsicPeakLines != null) {
	    for (int i = 0; i < intrinsicPeakLines.length; i++) {
		intrinsicPeakLines[i] = null;
	    }
	    intrinsicPeakLines = null;
	}
	if (totalISPeakLines != null) {
	    for (int i = 0; i < totalISPeakLines.length; i++) {
		totalISPeakLines[i] = null;
	    }
	    totalISPeakLines = null;
	}
	System.gc();
    }

    /**
     * Draws the local view of the content of this snapshot onto the
     * specified Graphics2D object.
     * @param	g	the Graphics2D object to draw onto.
     */
    public void drawContentPrintView(Graphics2D g) {
	Rectangle rect = g.getClipBounds();
	int xPos = (int) rect.getX();
	xPos = Math.max(xPos, 0);

	ChromDataPainter.paint(g, xPos, viewWidth, 0, traceLength, 
			       zoomX, pvChromLayout, traceLines, 
			       phdPeakLocations, phdBases, qualityValues,
			       phdPeakLocations2, phdBases2, 
			       qualityValues2, abcData, alignData, 
			       origBasePeakLocations, origBases,
			       intrinsicPeakLines, totalISPeakLines,
			       displayRaw, displayAbc,
			       displayAln, displayOrigBases, displayTip,
			       displayTotalIS,
			       displayTrim, leftTrimPoint, rightTrimPoint, 
			       trimThreshold, displayPhd2, markClosestPeak,
			       ChromDataPainter.PRINT_LOCAL, null,
			       isTT, isTT2);
    }

    /**
     * Draws the global view of the content of this snapshot onto the
     * specified Graphics2D object.  The snapshot to be printed starts
     * from the specified position.
     * @param	g	the Graphics2D object to draw onto.
     * @param	xPos	the starting point of the content to be printed.
     * @param	width	the width of the content to be printed.
     */
    public void drawContentPrintAll(Graphics2D g, int xPos, int width) {
	xPos = Math.max(xPos, 0);
	int[] noMarkForClosestPeak = { -1, -1 };
	ChromDataPainter.paint(g, xPos, width, 0, 
			       analyzedTraceLength, 1.0f, 
			       paChromLayout, analyzedTraceLines, 
			       phdPeakLocations, phdBases, qualityValues,
			       phdPeakLocations2, phdBases2, 
			       qualityValues2, abcData, 
			       null, null, null, null, null,
			       false, displayAbc, false, false, false, false,
			       displayTrim, leftTrimPoint, rightTrimPoint, 
			       trimThreshold,
			       displayPhd2, noMarkForClosestPeak, 
			       ChromDataPainter.PRINT_GLOBAL, null,
			       isTT, isTT2);

    }

    /* This last section of the code is to provide public access to
     * some of the data in the ChromSnapShot data, especially for TTPrintJob
     * class.
     */
    public int[] getTraceSignalValues() { return signalValues; }
    public boolean signalValuesShown() { return signalValuesShown; }
    public boolean fifthDyeDataShown() { return fifthDyeDataShown; }
    public String getChromFileName() { return chromFileName; }
    public String getPhdFileName() { return phdFileName; }
    public String getPhd2FileName() { return phd2FileName; }
    public String getProductName() { return productName; }
    public Point getViewPos() { return viewPos; }
    public int getViewWidth() { return viewWidth; }
    public int getViewHeight() { return viewHeight; }
    public int getType() { return type; }
    public int getAnalyzedTraceLength() { return analyzedTraceLength; }
    //public boolean displayTrim() { return displayTrim; }
    /** Returns null if trim regions are not being displayed;
	trim data otherwise. */
    public int[] getTrimData() { 
	if (!displayTrim) {
	    return null; 
	}
	int[] trim = {leftTrimPoint, rightTrimPoint, trimThreshold};
	return trim;
    }
    public boolean hasValidSampleAndPhd() { return !sampleOrPhdInvalid; }
    public int getCurrentFileIndex() { return currFileIndex; }
    public int getTotalNumOfFiles() { return totalNumOfFiles; }
}
