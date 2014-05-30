/*
 * 1.25
 *
 * @(#)ChromViewPanel.java       1.0     2000/10/25
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import com.paracel.tt.chrom.*;
import com.paracel.tt.util.*;
import com.paracel.tt.event.*;
import com.paracel.tt.io.*;

/**
 * This class creates the upper panel of TTView window, which contains the
 * actual drawing of the traces, quality values, etc.
 */
public class ChromViewPanel extends JPanel
			    implements ZoomEventListener,
				       DisplayFeatureEventListener,
				       FileNavigateEventListener,
				       SearchEventListener,
				       SearchConstants {
    final static Font baseFont = new Font("Courier", Font.PLAIN, 12);
    final static Font rowHeaderFont = new Font("Monospaced", Font.PLAIN, 11);

    /** The chromat view scroll panel. */
    private JScrollPane scrollPanel;

    /** The rowheader of the scroll panel. */
    private JLabel rowHeader;

    /** The TTViewFrame this panel resides in. */
    private JFrame frame;

    /** The chromatogram file name. */
    private String chromFileName;

    /** The phd file name. */
    private String phdFileName;

    /** The tab file name. */
    private String tabFileName;

    /** The tal file name. */
    private String talFileName;

    /** The tip file name. */
    private String tipFileName;

    /** The other phd file name. */
    private String phd2FileName;

    /** The chromatogram. */
    private AbstractChromatogram chromatogram;

    /** The sample file data. */
    private ABIFile chromFile;

    /** The phd data. */
    private PhdFile phdFile;

    /** The second phd data. */
    private PhdFile phdFile2;

    /** The tab file data. */
    private TabFile tabFile;

    /** The tal file data. */
    private TalFile talFile;

    /** The tip file data. */
    private TipFile tipFile;

    /** The chromatogram type: ABI or SCF. */
    private int chromType;

    /** The flag indicating whether to show another phd file. */
    private boolean displayPhd2;

    /** The flag to determin whether display raw data. */
    private boolean displayRaw;

    /** The flag to determin whether display original base calls. */
    private boolean displayOrig;
    
    /** The flag to determin whether display alternative base calls. */
    private boolean displayAbc;
    
    /** The flag to determin whether display consensus. */
    private boolean displayAln;
    
    /** The flag to determin whether display intrinsic peaks. */
    private boolean displayTip;
    
    /** The flag to determin whether display total intrinsic signal. */
    private boolean displayTotalIS;
    
    /** The flag to determin whether mark trimmed regions. */
    private boolean displayTrim;
    
    /** The flag indicating whether both sample file and phd file are valid. */
    private boolean spFilesAreValid;

    /** The maximum intensity of the traces currently shown in the viewer. */
    private int maxIntensity; 

    /** The maximum length of the traces currently shown in the viewer. */
    private int traceLength; 

    /** The maximum length of the analyzed traces. */
    private int analyzedTraceLength; 

    /** The maximum intensity of the analyzed traces. */
    private int analyzedTraceMaxIntensity; 

    /** The currently-shown trace lines which has y converted from intensity 
	to pixels. */
    private Trace[] traceLines;

    /** The currently-shown traces data. */
    private Trace[] tracesToDisplay;

    /** The analyzed traces. */
    private Trace[] analyzedTraces;

    /** The current and previous zoom values. */
    private float zoomX, zoomY, lastZoomX, lastZoomY;

    /** The phd base peak location array of 2 phd files. */
    private int[] phdPeakLocations, phdPeakLocations2;

    /** The quality value array of 2 phd files. */
    private int[] qualityValues, qualityValues2;

    /** The called bases array of 2 phd files. */
    private char[] phdBases, phdBases2;

    /** Whether the base calls are from TraceTuner for 2 phd files. */
    private boolean isTT, isTT2;

    /** The original bases. */
    private String origBases;

    /** The peak location array of the original bases. */
    private short[] origBasePeakLocations;

    /** The status of the tab file. */
    private int tabFileStatus;

    /** The abc data array. */
    private TabFile.AbcData[] abcData;

    /** The status of the tal file. */
    private int talFileStatus;

    /** The align data array. */
    private TalFile.AlignData[] alignData;

    /** The layout of the elements to be drawn on the panel. */
    private ChromLayout chromLayout;

    /** A 2-element int array for the closest peak to be marked on the graph.
        The first element is 0 for an ABI (original) base call peak, 1 for
        a TT base call peak, 2 for an alternative base call. The second 
	element is the actual x location of the peak or the base call. */
    private int[] markClosestPeak;

    /** The vector of the registered DisplaySignalValueEventListeners. */
    private Vector displaySigValEventListeners;

    /** The intrinsic peak list, and the total intrinsic signal peak list
	read from tip file, and their data after intensity --> pixel 
	conversion on y axis. */
    private TipFile.IntrinsicPeak[] intrinsicPeaks, intrinsicPeakLines,
    				    totalISPeaks, totalISPeakLines;

    /** Indicates whether data prepared for printing needs to be updated. */
    private boolean dataChanged;

    /** Trim points and threshold. */
    private int leftTrimPoint, rightTrimPoint, trimThreshold;

    /** The current file index. */
    private int fileIndex;

    /** The base strings for search. */
    private String phdSearchString, abiSearchString, refSearchString;

    /** The real indexes array of refSearchString.  eg. the consensus base of
	refSearchString[i] is actually the consensus base of
	alignData[refSearchRealIndexes[i]].   The reason of having this 
        referencing index array is the gaps ('-'s) in the consensus sequence
        of the alignment. */
    private int[] refSearchRealIndexes;

    /** The index where the current search starts. */
    private int currSearchIndex;

    /** The subject of the last search event. (search-what?) */
    private int lastSearchSubject = NONE;
    
    /** The Search object. */
    private Search search = new Search();
    
    /** The search result. */
    private Search.SearchResult searchResult;

    private static int MIN_HEIGHT;

    /** Constructor */
    public ChromViewPanel(JFrame f, String nameOfPhd2File, int phd2Status) {
	super();
	frame = f;

	zoomX = 1.0f;
	zoomY = 1.0f;
	
	phd2FileName = nameOfPhd2File;

	if (phd2FileName == null || phd2Status != FilesStatus.VALID) {
	    displayPhd2 = false;
	    MIN_HEIGHT = 300;
	} else {
	    displayPhd2 = true;
	    MIN_HEIGHT = 400;
	}

	displaySigValEventListeners = new Vector(1);

	addMouseListener(new MouseAdapter() {
	    public void mouseClicked(MouseEvent e) {
		if (!spFilesAreValid) {
		    return; 
		}
		int clickCount = e.getClickCount();
		int xPos = (int) e.getX();
		int yPos = (int) e.getY();
		if (clickCount == 1) {
		    repaint();
		    if ((yPos < chromLayout.traceLinesY 
			    	- chromLayout.chromHeight)
		    	    || (yPos > chromLayout.traceLinesY)) {
		    	// if the click is not on the chromatogram, do nothing.
		    	return; 
		    }
		    showSignalValues(xPos);
		} else if (clickCount == 2) {
		    if (displayRaw) {
			return;
		    }
		    showClosestPeak(xPos, yPos);
		    repaint();
		}
	    }
	});

	chromLayout = new ChromLayout();
	markClosestPeak = new int[2];

	scrollPanel = new JScrollPane(this);
	scrollPanel.setPreferredSize(new Dimension(f.getSize().width, 
					f.getSize().height - 150));
	final JScrollBar hBar = scrollPanel.getHorizontalScrollBar();
	hBar.setUnitIncrement(10);
	scrollPanel.getVerticalScrollBar().setUnitIncrement(5);

	hBar.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent e) {
		if (!(hBar.getValueIsAdjusting())) {
                    repaint();
                }
            }
        });

	hBar.addComponentListener(new ComponentAdapter() {
	    public void componentResized(ComponentEvent e) {
		/* Set the block increment to be 30 pixels less than
		   the visible amount for better visual effect: when
		   the user clicks on the paging area of the scroll bar,
		   the view changes by a page (the size of the current
		   view window) but has 30-pixel image overlapping with
		   the previous view for continuity. */
		int i = hBar.getVisibleAmount() - 30;
		if (i <= 0) {
		    i = hBar.getVisibleAmount();
		}
		hBar.setBlockIncrement(i);
	    }
	});

	// set up the row header for printing the labels
	rowHeader = new JLabel() {
	    public void paintComponent(Graphics g) {
		if (viewPositionAdjusting) {
		    return;
		}
		super.paintComponent(g);
		//Clipping for speed
		if (!dataLoaded || !spFilesAreValid || displayRaw) {
		    return; 
		}
		Rectangle r = g.getClipBounds();
		g.setFont(rowHeaderFont);
		g.setColor(Color.black);
		if (phdFile != null) {
		    String s1 = isTT ? "tt" : "ph"; 
		    String s2 = isTT ? "TT" : "PH"; 
		    g.drawString(s1, r.x + 6, chromLayout.indexY);
		    g.drawString(s2, r.x + 6, chromLayout.phdBasesY);
		}

		if (displayPhd2) {
		    String s3 = isTT2 ? "TT" : "PH"; 
		    g.drawString(s3, r.x + 6, chromLayout.phdBases2Y);
		}
		if (displayAln) {
		    g.drawString("ref", r.x + 2, chromLayout.consIndexY);
		    g.drawString("REF", r.x + 2, chromLayout.consBasesY);
		}

		if (displayOrig) {
		    g.drawString("ABI", r.x + 2, chromLayout.origBasesY);
		}

		if (displayAbc) {
		    g.drawString("ALT", r.x + 2, chromLayout.altBasesY);
		}
	    }

	    public Dimension getPreferredSize() {
		return new Dimension(25, 
			scrollPanel.getViewport().getViewSize().height);
	    }
	};
	rowHeader.setBorder(BorderFactory.createRaisedBevelBorder());
	rowHeader.setBackground(Constants.gray208);
	rowHeader.setOpaque(true);
	scrollPanel.setRowHeaderView(rowHeader);

	scrollPanel.addComponentListener(new ComponentAdapter() {
	    public void componentResized(ComponentEvent e) {
		scrollPanelResized();
	    }
	});
    }

    /** Returns the scroll panel that wraps the ChromViewPanel. */
    public JComponent getWrapperComponent() {
	return scrollPanel;
    }
    
    /** Attempts to release the resouces. */
    public void releaseResources() {
	scrollPanel = null;
	rowHeader = null;
	chromatogram = null;
	if (chromFile != null) {
	    chromFile.releaseResources();
	    chromFile = null;
	}
	if (phdFile != null) {
	    phdFile.releaseResources();
	    phdFile = null;
	}
	if (phdFile2 != null) {
	    phdFile2.releaseResources();
	    phdFile2 = null;
	}
	if (tabFile != null) {
	    tabFile.releaseResources();
	    tabFile = null;
	}
	if (talFile != null) {
	    talFile.releaseResources();
	    talFile = null;
	}
	if (tipFile != null) {
	    tipFile.releaseResources();
	    tipFile = null;
	}
	if (traceLines != null) {
	    for (int i = 0; i < traceLines.length; i++) {
		if (traceLines[i] != null) {
		    traceLines[i].dispose();
		    traceLines[i] = null;
		}
	    }
	}
	if (tracesToDisplay != null) {
	    for (int i = 0; i < tracesToDisplay.length; i++) {
		if (tracesToDisplay[i] != null) {
		    tracesToDisplay[i].dispose();
		    tracesToDisplay[i] = null;
		}
	    }
	}
	if (analyzedTraces != null) {
	    for (int i = 0; i < analyzedTraces.length; i++) {
		if (analyzedTraces[i] != null) {
		    analyzedTraces[i].dispose();
		    analyzedTraces[i] = null;
		}
	    }
	}
	phdPeakLocations = null;
	phdPeakLocations2 = null;
	qualityValues = null;
	qualityValues2 = null;
	phdBases = null;
	phdBases2 = null;
	origBases = null;
	origBasePeakLocations = null;
	chromLayout = null;
	intrinsicPeaks = null;
	intrinsicPeakLines = null;
	totalISPeaks = null;
	totalISPeakLines = null;
	phdSearchString = null;
	abiSearchString = null;
	refSearchString = null;
	refSearchRealIndexes = null;
	search = null;
	displaySigValEventListeners = null;
    }

    /* The flag indicating whether tip file data is in loading process. */
    boolean tipFileLoading;

    /* The flag indicating whether tip file data is loaded. */
    boolean tipFileLoaded;

    /* The flag indicating whether data is loaded. The reason of using this
     * flag is that paintComponent() function is called for some other reason,
     * which happened before the actual data loading finishes.
     */
    boolean dataLoaded;

    /**
     * Draws the traces, bases calls, quality values, etc.
     * @param	g   The <code>Graphics</code> object.
     */
    public void paintComponent(Graphics g) {
	rowHeader.repaint();
	Dimension d = getSize();
	g.setColor(Color.white);
	g.fillRect(0, 0, d.width, d.height);

	/* If data loading is not done, or sample file and phd file 
	   are invalid, do not paint anything. */
	if (!dataLoaded || !spFilesAreValid) {
	    return;
	}

	// get current clip of the chrom view window
	Rectangle rect = (Rectangle) g.getClipBounds();
	int currViewXPos = Math.max((int) rect.getX(), 0);
	int currViewWidth = (int) rect.getWidth();
	
	if (dataChanged) {
	    dataChanged = false;
	    layoutElements();
	    getTraceLinesData();
	    if (tipFileLoaded && (displayTip||displayTotalIS) && zoomX > 0.5) {
    	    	updateIntrPeakLinesY();
	    }
	}

	// Draw the elements onto the panel.
	try {
	    ChromDataPainter.paint((Graphics2D) g, currViewXPos, currViewWidth,
		    5, traceLength, zoomX, chromLayout, traceLines,
		    phdPeakLocations, phdBases, qualityValues,
		    phdPeakLocations2, phdBases2, qualityValues2,
		    abcData, alignData, origBasePeakLocations, origBases,
		    intrinsicPeakLines, totalISPeakLines,
		    displayRaw, displayAbc, displayAln, displayOrig,
		    tipFileLoaded && displayTip && (zoomX > 0.5),
		    tipFileLoaded && displayTotalIS && (zoomX > 0.5),
		    displayTrim, leftTrimPoint, rightTrimPoint, trimThreshold,
		    displayPhd2, markClosestPeak, ChromDataPainter.VIEW_PANEL,
		    searchResult, isTT, isTT2);
	} catch (ArrayIndexOutOfBoundsException ignored) {
	    /* This exception is ignored because it is caused by
	       repainting during the view position is adjusting 
	       when the zoom X value changes. */
	} catch (Exception e) {
	    //e.printStackTrace(System.err);
	}
    }

    /** Lays out the elements to be drawn on the panel. */
    protected void layoutElements() {
	int panelHeight = getSize().height;
	int baseLetterHeight = getFontMetrics(baseFont).getHeight();
	int offset = 18;
	int qvHeight = 100 + baseLetterHeight;

	chromLayout.traceLinesY = panelHeight - 30;
	chromLayout.peakLocY = panelHeight - 5;

	if (displayRaw) {
	   chromLayout.chromHeight = panelHeight - 80;
	} else {
	    if (displayOrig && displayAln) {
	    	chromLayout.consIndexY = baseLetterHeight - 2;
	    	chromLayout.indexY = chromLayout.consIndexY+baseLetterHeight-4;
	    	chromLayout.indexLineY = chromLayout.consIndexY - 4;
		chromLayout.consBasesY = chromLayout.indexY + 14;
		chromLayout.origBasesY = chromLayout.consBasesY + offset;
	    	chromLayout.phdBasesY = chromLayout.origBasesY + offset;
	    } else if (displayOrig) {
	    	chromLayout.indexY = baseLetterHeight - 2;
	    	chromLayout.indexLineY = chromLayout.indexY - 4;
		chromLayout.origBasesY = chromLayout.indexY + offset + 4;
	    	chromLayout.phdBasesY = chromLayout.origBasesY + offset;
	    } else if (displayAln) {
	    	chromLayout.consIndexY = baseLetterHeight - 2;
	    	chromLayout.indexY = chromLayout.consIndexY+baseLetterHeight-4;
	    	chromLayout.indexLineY = chromLayout.consIndexY - 4;
		chromLayout.consBasesY = chromLayout.indexY + 14;
		chromLayout.phdBasesY = chromLayout.consBasesY + offset;
	    } else {
	    	chromLayout.indexY = baseLetterHeight - 2;
	    	chromLayout.indexLineY = chromLayout.indexY - 4;
		chromLayout.phdBasesY = chromLayout.indexY + offset + 4;
	    } 
	    if (displayAbc) {
		chromLayout.altBasesY = chromLayout.phdBasesY;
		chromLayout.phdBasesY = chromLayout.altBasesY + offset + 8;
	    }

	    chromLayout.indexLineLength = baseLetterHeight;
	    chromLayout.qvBarY = chromLayout.phdBasesY + 2;
	    chromLayout.chromHeight = chromLayout.traceLinesY -
		    			  chromLayout.phdBasesY - qvHeight;
	    if (displayPhd2) {
		chromLayout.phdBases2Y = chromLayout.phdBasesY + qvHeight;
	        chromLayout.qvBar2Y = chromLayout.phdBases2Y + 2;
	        chromLayout.chromHeight = chromLayout.traceLinesY -
		    			  chromLayout.phdBases2Y - qvHeight;
	    }
	}
    }

    /** Mark the peak closest to the specified x location and y location. */
    protected void showClosestPeak(int xPos, int yPos) {
	int actualXPos = (int) (xPos/zoomX);
	int temp, temp1;
	if (displayOrig && yPos >= (chromLayout.origBasesY - 12)
	    	&& yPos <= (chromLayout.origBasesY + 4)) {
	    /* if double clicked on the ABI original base calls,
	       show the location of the closest original base peak. */
	    markClosestPeak[0] = 0; // 0 indicates ABI(original) base call peak
	    for (int i = 0; i < origBasePeakLocations.length; i++) {
	    	temp = origBasePeakLocations[i];
		if (temp == actualXPos) {
		    markClosestPeak[1] = temp;
		    return;
	    	} else if (temp > actualXPos) {
		    if (i == 0) {
		    	markClosestPeak[1] = temp;
			return;
		    } else {
		    	// check whether the previous peak is closer
		    	temp1 = origBasePeakLocations[i-1];
		    	markClosestPeak[1] = 
				((actualXPos-temp1)<=(temp-actualXPos)) ? temp1
			    						: temp;
		    	return;
		    }
	    	}
	    }
	    return;
	} 

	if (displayAbc && yPos >= (chromLayout.altBasesY - 12)
		   && yPos <= (chromLayout.altBasesY + 4)) {
	    /* If double clicked on an alternative base call,
	       show its location. */
	    markClosestPeak[0] = 2; // 2 indicates alternative base call peak
	    for (int i=0; i < abcData.length; i++) {
		temp = abcData[i].getPosition();
		if (temp == actualXPos) {
		    markClosestPeak[1] = temp;
		    return;
	    	} else if (temp > actualXPos) {
		    if (i == 0) {
		    	markClosestPeak[1] = temp;
			return;
		    } else {
		    	// check whether the previous base call is closer
		    	temp1 = abcData[i-1].getPosition();
		    	markClosestPeak[1] = 
				((actualXPos-temp1)<=(temp-actualXPos)) ? temp1
			    						: temp;
		    	return;
		    }
	    	}
	    }
	    return;
	} 

	if ((yPos >= (chromLayout.phdBasesY - 12) 
		   && yPos <= (chromLayout.phdBasesY + 4))
	    	   || (yPos >= chromLayout.traceLinesY - chromLayout.chromHeight)) {
	    /* If double clicked on a TT base call or the chromatogram,
	       show the location of the TT base call. */
	    markClosestPeak[0] = 1; // 1 indicates TT base call peak
	    for (int i = 0; i < phdPeakLocations.length; i++) {
	    	temp = phdPeakLocations[i];
	    	if (temp == actualXPos) {
		    markClosestPeak[1] = temp;
		    return;
	    	} else if (temp > actualXPos) {
		    if (i == 0) {
		    	markClosestPeak[1] = temp;
		    	return;
		    } else {
		    	// check whether the previous peak is closer
		    	temp1 = phdPeakLocations[i-1];
		    	markClosestPeak[1] = 
				((actualXPos-temp1)<=(temp-actualXPos)) ? temp1
			    						: temp;
		    	return;
		    }
	    	}
	    }
	}

	if (displayPhd2 && yPos >= (chromLayout.phdBases2Y - 12) 
		   && yPos <= (chromLayout.phdBases2Y + 4)) {
	    /* If double clicked on a TT base call of the second phd file,
	       show its location. */
	    markClosestPeak[0] = 3; // 3 indicates 2nd phd TT base call peak
	    for (int i = 0; i < phdPeakLocations2.length; i++) {
	    	temp = phdPeakLocations2[i];
	    	if (temp == actualXPos) {
		    markClosestPeak[1] = temp;
		    return;
	    	} else if (temp > actualXPos) {
		    if (i == 0) {
		    	markClosestPeak[1] = temp;
		    	return;
		    } else {
		    	// check whether the previous peak is closer
		    	temp1 = phdPeakLocations2[i-1];
		    	markClosestPeak[1] = 
				((actualXPos-temp1)<=(temp-actualXPos)) ? temp1
			    						: temp;
		    	return;
		    }
	    	}
	    }
	}
    }

    /**
     * Used to eliminate the effect caused by unnecessary calls to
     * paintComponent() when view position is adjusting because of zooming.
     *
     * Its value is set in adjustViewPositionAfterZoom().
     * In paintComponent(), its value is checked.  If it is true,
     * paintComponent() does nothing when called.
     */
    boolean viewPositionAdjusting;

    /**
     * Adjust the view position of the view port ( the left upper corner
     * of the TTView trace window ) when zoom factor is changed to ensure
     * the traces in the current window start with the same base as before
     * zooming.
     */
    protected void adjustViewPositionAfterZoom() {
        viewPositionAdjusting = true;
        Point currViewPos = scrollPanel.getViewport().getViewPosition();
        int x, y;
        if (zoomX <= 0.1) {
            // if zoom to 10% or less, set the window to show the traces
            // from the beginning of the traces.
            x = 0;
        } else {
            // otherwise, compute the correct view position according to
            // to the zooming factors.
            x =  (int) (((float)currViewPos.getX()) / lastZoomX * zoomX);
	    int panelWidth = (int) (traceLength * zoomX) + 5;
	    int restTraceLen = panelWidth - x;
	    int frameWidth = frame.getSize().width;
	    if ( (panelWidth - x) < (frameWidth)) {
		x = panelWidth - frameWidth;
	    }
        }
        y = (int)currViewPos.getY();

        // Note: the scrollRectToVisible() is essential to cooperate with
        //       setViewPosition() to force the window scroll to the right
        //       position, because of a possible JDK bug in setViewPosition():
        //       If the position you are setting is too far away from
        //       current position, the window only gets half way scrolled.
        scrollPanel.getViewport().scrollRectToVisible(
                                new Rectangle((int)x, y, 5, 5));
        scrollPanel.getViewport().setViewPosition(new Point((int)x, y));

        viewPositionAdjusting = false;
    }

    /** Fires a DisplaySignalValueEvent for the specified x location */
    protected void showSignalValues(int xPos) {
	int AValue = 0;
	int CValue = 0;
	int GValue = 0;
	int TValue = 0;
	int fifthDyeValue = 0;
	String name = "";
	boolean fifthDyeExists = false;
	// Round the x location to the closest value available in the
	// sample file, instead of using an interpolated value.
	int actualXPos = Math.round(((float)xPos)/zoomX);
	if (actualXPos < 0) {
	    actualXPos = 0;
	}
	
	for (int i = 0; i < tracesToDisplay.length; i++) {
	    // each trace
	    if (actualXPos >= tracesToDisplay[i].getLength()) {
		return;
	    }
	    name = tracesToDisplay[i].getName();
	    if (name.equalsIgnoreCase("A")) {
		AValue = tracesToDisplay[i].getIntensityAt(actualXPos);
	    } else if (name.equalsIgnoreCase("C")) {
		CValue = tracesToDisplay[i].getIntensityAt(actualXPos);
	    } else if (name.equalsIgnoreCase("G")) {
		GValue = tracesToDisplay[i].getIntensityAt(actualXPos);
	    } else if (name.equalsIgnoreCase("T")) {
		TValue = tracesToDisplay[i].getIntensityAt(actualXPos);
	    } else if (name.equalsIgnoreCase("Fifth Dye")) {
		fifthDyeValue = tracesToDisplay[i].getIntensityAt(actualXPos);
		fifthDyeExists = true;
	    }
	}

	DisplaySignalValueEvent e;
	if (fifthDyeExists) {
	    e = new DisplaySignalValueEvent(this, actualXPos, AValue, CValue, 
					    GValue, TValue, fifthDyeValue);
	} else {
	    e = new DisplaySignalValueEvent(this, actualXPos, AValue, CValue, 
					    GValue, TValue);
	}

	for (int i = 0; i < displaySigValEventListeners.size(); i++) {
	    ((DisplaySignalValueEventListener) (displaySigValEventListeners.
		elementAt(i))).signalValueSet(e);
	}
    }

    /** Converts the trace signal from intensity values to pixel values
        that suit the panel size. */
    protected void getTraceLinesData() {

	if (traceLines == null) {
	    traceLines = new Trace[tracesToDisplay.length];
	    int len;
	    for (int i = 0; i < tracesToDisplay.length; i++) {
		len = tracesToDisplay[i].getLength();
	    	traceLines[i] = new Trace(tracesToDisplay[i].getName(),
	    			      tracesToDisplay[i].getColor(),
	    			      (new int[len]));
	    }
	}

	// Convert intensity traces to y pixel traces.
	float yToPixelRatio = ((float)maxIntensity) 
	    			/ ((float)chromLayout.chromHeight);
	int temp;
	for (int i = 0; i < tracesToDisplay.length; i++) {
	    for (int j = 0; j <tracesToDisplay[i].getLength(); j++) {
		temp = (int) (((float)(tracesToDisplay[i].getIntensityAt(j))) 
				       * zoomY / yToPixelRatio);
		if (temp > chromLayout.chromHeight) {
		    temp = chromLayout.chromHeight;
		}
		traceLines[i].setIntensityAt(j,chromLayout.traceLinesY - temp);
	    }
	}
    }

    /** Generates the intrinsic peak lines data from the intrinsic
	peak data. */
    protected void getIntrinsicPeakLinesData() {
	int length = intrinsicPeaks.length;
	intrinsicPeakLines = null;
	intrinsicPeakLines = new TipFile.IntrinsicPeak[length];
	for (int i = 0; i < length; i++) {
	    intrinsicPeakLines[i] = tipFile.new IntrinsicPeak(
						intrinsicPeaks[i]);
	}

	length = totalISPeaks.length;
	totalISPeakLines = null;
	totalISPeakLines = new TipFile.IntrinsicPeak[length];
	for (int i = 0; i < length; i++) {
	    totalISPeakLines[i] = tipFile.new IntrinsicPeak(totalISPeaks[i]);
	}

	updateIntrPeakLinesY();
    }

    /** Updates the y values of the intrinsic peak lines data. */
    protected void updateIntrPeakLinesY() {
	// Convert intensity traces to y pixel traces.
	float yToPixelRatio = ((float)maxIntensity) 
	    			/ ((float)chromLayout.chromHeight);

	int length = intrinsicPeakLines.length;
	int len, temp;
	for (int i = 0; i < length; i++) {
	    len = intrinsicPeakLines[i].length;
	    for (int j = 0; j < len; j++) {
		temp = (int) (((float)(intrinsicPeaks[i].y[j])) * zoomY
			      / yToPixelRatio);
		if (temp > chromLayout.chromHeight) {
		    temp = chromLayout.chromHeight;
		}
		intrinsicPeakLines[i].y[j] = chromLayout.traceLinesY - temp;
	    }
	}

	length = totalISPeakLines.length;
	for (int i = 0; i < length; i++) {
	    len = totalISPeakLines[i].length;
	    for (int j = 0; j < len; j++) {
		temp = (int) (((float)(totalISPeaks[i].y[j])) * zoomY
			      / yToPixelRatio);
		if (temp > chromLayout.chromHeight) {
		    temp = chromLayout.chromHeight;
		}
		totalISPeakLines[i].y[j] = chromLayout.traceLinesY - temp;
	    }
	}
    }

    /** Update the view to show the specified files. */
    public void updateFiles(ABIFile chromatFile, PhdFile phdFile,
			    TabFile tabFile, TalFile talFile,
			    String nameOfChromatFile, String nameOfPhdFile,
			    String nameOfTabFile,
			    String nameOfTalFile, String nameOfTipFile) {
	dataLoaded = false;
	if (this.chromFile != null) {
	    this.chromFile.releaseResources();
	    this.chromFile = null;
	}
	if (this.phdFile != null) {
	    this.phdFile.releaseResources();
	    this.phdFile = null;
	}
	if (this.tabFile != null) {
	    this.tabFile.releaseResources();
	    this.tabFile = null;
	}
	if (this.talFile != null) {
	    this.talFile.releaseResources();
	    this.talFile = null;
	}
	this.chromFile = chromatFile;
	this.phdFile = phdFile;
	this.tabFile = tabFile;
	this.talFile = talFile;
	chromFileName = nameOfChromatFile;
	phdFileName = nameOfPhdFile;
	tabFileName = nameOfTabFile;
	talFileName = nameOfTalFile;
	tipFileName = nameOfTipFile;
	markClosestPeak[0] = -1;
	markClosestPeak[1] = -1;
	initialize();
	dataLoaded = true;
	dataChanged = true;
	revalidate();
	repaint();
	rowHeader.repaint();
    }

    /** Initialize the data from the files (sample file, phd file(s),
        tip file, tab file, tal file). */
    public void initialize() {
	chromType = AbstractChromatogram.getChromType(chromFile);

	if (chromType == AbstractChromatogram.ABI) {
	    chromatogram = new ABIChromatogram(chromFile);
	    if (Constants.fromEdited) {
	    	origBases = ((ABIChromatogram)
			    	chromatogram).getEditedBases();
	    	origBasePeakLocations = ((ABIChromatogram)
			    	chromatogram).getEditedPeakLocations();
	    } else {
	    	origBases = ((ABIChromatogram)
			    	chromatogram).getCalledBases();
	    	origBasePeakLocations = ((ABIChromatogram)
			    	chromatogram).getCalledPeakLocations();
	    }
	} else {
	    chromatogram = new SCFChromatogram(chromFile);
	    origBases = ((SCFChromatogram)
			    chromatogram).getCalledBases();
	    origBasePeakLocations = ((SCFChromatogram)
			    chromatogram).getCalledPeakLocations();
	}
	abiSearchString = null;
	abiSearchString = new String(origBases);

	analyzedTraces = chromatogram.getAnalyzedTraces();
	analyzedTraceMaxIntensity = chromatogram.getMaxIntensity();
	analyzedTraceLength = chromatogram.getLength();

	isTT = phdFile.isTT();
	phdBases = null;
	phdBases = LangTools.toUpperCase(phdFile.getBases());
	phdSearchString = null;
	phdSearchString = new String(phdBases);
	qualityValues = null;
	qualityValues = phdFile.getQvValues();
	phdPeakLocations = null;
	phdPeakLocations = phdFile.getPeakLocations();
	int trim[] = phdFile.getTrimData();
	if (trim != null) {
	    leftTrimPoint = trim[0];
	    rightTrimPoint = trim[1];
	    trimThreshold = trim[2];
	}

	if (displayPhd2) {
	    try {
	    	phdFile2 = new PhdFile(phd2FileName);
		isTT2 = phdFile2.isTT();
	    	phdBases2 = LangTools.toUpperCase(phdFile2.getBases());
	    	qualityValues2 = phdFile2.getQvValues();
	    	phdPeakLocations2 = phdFile2.getPeakLocations();
	    } catch (Exception e) {
	    	System.err.println(e.toString());
	    	displayPhd2 = false;
	    }
	}

	if (tabFileStatus == FilesStatus.VALID
	    	|| tabFileStatus == FilesStatus.VALID_BUT_STALE) {
	    abcData = tabFile.getHighestQVAbcData();
	}

	if (talFileStatus == FilesStatus.VALID
	    	|| talFileStatus == FilesStatus.VALID_BUT_STALE) {
	    alignData = talFile.getAlignData();
	    /* Builds the aligned consensus sequence string for search 
	       purpose.  '-'s are removed. */
	    StringBuffer sb = new StringBuffer();
	    int[] indexes = new int[alignData.length];
	    char c;
	    int count = 0;
	    for (int i = 0; i < alignData.length; i++) {
		c = alignData[i].getConsBase();
		if (c != '-') {
		    sb.append(c);
		    indexes[count++] = i;
		}
	    }
	    refSearchString = sb.toString();
	    refSearchRealIndexes = null;
	    refSearchRealIndexes = new int[count];
	    System.arraycopy(indexes, 0, refSearchRealIndexes, 0, count);
	}

	getTracesToDisplay();
	layoutElements();
	getTraceLinesData();

	loadTipFile();
    }

    /** Loads the tip file. */
    protected void loadTipFile() {
	if (tipFileLoaded || tipFileLoading) {
            frame.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	    return; 
	}
	
	if (displayTip || displayTotalIS) {    
	    tipFileLoading = true;
	    Thread runner = new Thread() {
	    	public void run() {
	    	    try {
            		frame.setCursor(Cursor.getPredefinedCursor(
						Cursor.WAIT_CURSOR));
            		setCursor(Cursor.getPredefinedCursor(
						Cursor.WAIT_CURSOR));
			if (tipFile != null) {
			    tipFile.releaseResources();
			    tipFile = null;
		        }
	    	    	tipFile = new TipFile(tipFileName, frame);
		    	intrinsicPeaks = tipFile.getPeaks();
		    	totalISPeaks = tipFile.getTotalSignalPeaks();
		    	getIntrinsicPeakLinesData();
		    	tipFileLoaded = true;
			tipFileLoading = false;
		    	repaint();
	    	    } catch (Exception e) {
		    	System.err.println("ChromViewPanel--initialize:"
						+e.toString());
		    	e.printStackTrace(System.err);
		    	tipFileLoaded = false;
			tipFileLoading = false;
	    	    } finally {
            		frame.setCursor(Cursor.getPredefinedCursor(
						Cursor.DEFAULT_CURSOR));
            		setCursor(Cursor.getPredefinedCursor(
						Cursor.DEFAULT_CURSOR));
		    }
	    	}
	    };
	    runner.start();
	} else {
            frame.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	}
    }

    /** Gets the right traces to display (analyzed or raw depending on
	the menu selection. */
    public void getTracesToDisplay() {
	if (traceLines != null) {
	    for (int i = 0; i < traceLines.length; i++) {
		if (traceLines[i] != null) {
	    	    traceLines[i].dispose();
		    traceLines[i] = null;
		}
	    }
	    traceLines = null;
	}

	if (!displayRaw) {
	    tracesToDisplay = chromatogram.getAnalyzedTraces();
	    maxIntensity = chromatogram.getMaxIntensity();
	    traceLength = chromatogram.getLength();
	} else {
	    tracesToDisplay = ((ABIChromatogram)
					chromatogram).getRawDataTraces();
	    maxIntensity = ((ABIChromatogram)
					chromatogram).getRawDataMaxIntensity();
	    traceLength = ((ABIChromatogram)
					chromatogram).getRawDataTraceLength();
	}

	int minHeight = Math.max(scrollPanel.getSize().height - 50,MIN_HEIGHT);
	int width = Math.max(traceLength, 
			     phdPeakLocations[phdPeakLocations.length - 1]);
	setPreferredSize(new Dimension((int)(width * zoomX) + 25, 
				       minHeight));
    }

    /** Invoked when zoom factor on either X or Y axis is changed. */
    public void zoomValueChanged(ZoomEvent e) {
	if (e.getAxis() == 'X') {
	    lastZoomX = zoomX;
	    zoomX = e.getNewValue();
	    int minHeight = Math.max(scrollPanel.getSize().height - 50,
				     MIN_HEIGHT);
	    int width = Math.max(traceLength, 
			     phdPeakLocations[phdPeakLocations.length - 1]);
	    setPreferredSize(new Dimension((int) (width * zoomX) + 25,
					   minHeight));
	    revalidate();
	    if (lastZoomX != zoomX) {
		adjustViewPositionAfterZoom();
	    }
	} else {
	    lastZoomY = zoomY;
	    zoomY = e.getNewValue();
	    dataChanged = true;
	}

	repaint();
    }

    /** Invoked when the display features are changed. */
    public void displayFeatureChanged(DisplayFeatureEvent e) {
	boolean oldDisplayRaw = displayRaw;
	displayRaw = e.displayRaw();
	displayOrig = e.displayOrig();
	displayAbc = e.displayAbc();
	displayAln = e.displayAln();
	displayTip = e.displayTip();
	displayTotalIS = e.displayTotalIS();
	displayTrim = e.displayTrim();

	if (!e.fileIndexChanged()) {
	    /* if the display features are changed because of file index
	       change, getTracesToDisplay() and loadTipFile() will be 
	       called by navigateTo().*/
	    if (oldDisplayRaw != displayRaw) {
	    	getTracesToDisplay();
	    }
	    loadTipFile();
	    dataChanged = true;
	    revalidate();
	    repaint();
	    rowHeader.repaint();
	} else {
	    /* revalidate() and repaint() will be called by updateFiles() */
	    dataChanged = true;
	}
    }

    /** Invoked when navigating to another file in the file list. */
    public void navigateTo(FileNavigateEvent e) {
	searchResult = null;
	lastSearchSubject = NONE;
	int newIndex = e.getFilesStatus().getIndex();
	tipFileLoaded = false;
	tipFileLoading = false;
	int sampleFileStatus = e.getFilesStatus().getSampleFileStatus();
	int phdFileStatus = e.getFilesStatus().getPhdFileStatus();
	tabFileStatus = e.getFilesStatus().getTabFileStatus();
	talFileStatus = e.getFilesStatus().getTalFileStatus();

	if (traceLines != null) {
	    for (int i = 0; i < traceLines.length; i++) {
		if (traceLines[i] != null) {
	    	    traceLines[i].dispose();
		    traceLines[i] = null;
		}
	    }
	    traceLines = null;
	}
	if ((sampleFileStatus != FilesStatus.VALID)
		|| (phdFileStatus != FilesStatus.VALID)) {
	    spFilesAreValid = false;
	    repaint();
            frame.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	    return;
	}

	spFilesAreValid = true;
	updateFiles(e.getFilesStatus().getSampleFile(),
		    e.getFilesStatus().getPhdFile(),
		    e.getFilesStatus().getTabFile(),
		    e.getFilesStatus().getTalFile(),
		    e.getFilesStatus().getSampleFileName(),
		    e.getFilesStatus().getPhdFileName(),
		    e.getFilesStatus().getTabFileName(),
		    e.getFilesStatus().getTalFileName(),
		    e.getFilesStatus().getTipFileName());
	if (newIndex != fileIndex) {
	    fileIndex = newIndex;
	    scrollPanel.getHorizontalScrollBar().setValue(0);
	}
    }

    /** Invoked when a SearchEvent is fired: "Find" button in the search
	dialog is clicked, "Find Again" button/option is clicked, or one
	data row in SNPs table of the message panel is double-clicked. */
    public void search(SearchEvent se) {
	if (!dataLoaded || !spFilesAreValid) {
	    return;
	}
	int subject = se.getSubject();
	String pattern = se.getPattern();
	int direction = se.getDirection();
	Object source = se.getSource();
	String sbjString;

	if (source instanceof MessagePanel) {
	    /* If the SearchEvent comes from the Message Panel, it is not a
	       "real" search event.  The SNP table borrows the SearchEvent
	       facilities to draw a black solid arrow above the SNP that was
	       double-clicked by the user in the table.  So no real search
	       action needs to be taken.  The "search result" is already known,
	       which is the SNP double-clicked in the table. */
	    if (lastSearchSubject == NONE) {
		lastSearchSubject = subject;
	    }
	    if (subject == TT && displayRaw) {
		JOptionPane.showMessageDialog(frame,
				"TraceTuner base calls are not shown.");
		return;
	    }
	    searchResult = search.new SearchResult();
	    searchResult.begin = searchResult.end = se.getStartIndex();
	    searchResult.subject = subject;
	    searchResult.pattern = pattern;
	    searchResult.direction = direction;
	    currSearchIndex = searchResult.nextStart = searchResult.begin + 1;
	    searchResult.patternFound = true;
	    paintSearchResult();
	    return;
	}

	// Code below are for "real" search events.
	if (subject == REF && (!displayAln)) {
	    JOptionPane.showMessageDialog(frame,
				"Reference/consensus sequence is not shown.");
	    return;
	}
	if (subject == ABI && (!displayOrig)) {
	    JOptionPane.showMessageDialog(frame,
				"Original base calls are not shown.");
	    return;
	}
	if (subject == TT && displayRaw) {
	    JOptionPane.showMessageDialog(frame,
				"TraceTuner base calls are not shown.");
	    return;
	}

	if (subject == TT) {
	    sbjString = phdSearchString;
	} else if (subject == ABI) {
	    sbjString = abiSearchString;
	} else {
	    // subject == REF
	    sbjString = refSearchString;
	}
	if (lastSearchSubject != subject || lastSearchSubject == NONE) {
	    currSearchIndex = (direction == FORWARD)
				? 0 : (sbjString.length() - 1);
	    lastSearchSubject = subject;
	}
	searchResult = search.search(subject, sbjString,
				     pattern, currSearchIndex, direction);
	currSearchIndex = searchResult.nextStart;
	if (searchResult.patternFound) {
	    paintSearchResult();
	} else {
	    String comment = (direction == FORWARD)
		? "End of sequence reached; continue from beginning?"
		: "Beginning of sequence reached; continue from end?";
	    int res = JOptionPane.showConfirmDialog(frame, comment, "Question",
					  	JOptionPane.OK_CANCEL_OPTION);
	    if (res == JOptionPane.OK_OPTION) {
		currSearchIndex = (direction == FORWARD)
		    			? 0 : sbjString.length() - 1;
		searchResult = search.search(subject, sbjString, pattern,
					     currSearchIndex, direction);
		currSearchIndex = searchResult.nextStart;
		if (searchResult.patternFound) {
	   	    paintSearchResult();
		} else {
	    	    JOptionPane.showMessageDialog(frame,
						 "Search sequence not found.");
		}
	    }
	}
    }

    /** Paints the search result: Places a black solid arrow above each
	found base.  This function is called by search() funtion. */
    protected void paintSearchResult() {
	if (!dataLoaded || !spFilesAreValid || displayRaw) {
	    return;
	}
	if ((searchResult.subject == REF && (!displayAln))
	    	|| (searchResult.subject == TT && displayRaw)
	    	|| (searchResult.subject == ABI && (!displayOrig))) {
	    return;
	}

	//check whether the search result is in current window.
	Point p = scrollPanel.getViewport().getViewPosition();
	int x1 = p.x;
	int x2 = x1 + scrollPanel.getViewport().getWidth();
	int begin, end;
	if (searchResult.subject == TT) {
	    begin = (int) (phdPeakLocations[searchResult.begin] * zoomX);
	    end = (int) (phdPeakLocations[searchResult.end] * zoomX);
	} else if (searchResult.subject == ABI) {
	    begin = (int) (origBasePeakLocations[searchResult.begin] * zoomX);
	    end = (int) (origBasePeakLocations[searchResult.end] * zoomX);
	} else {
	    //searchResult.subject == REF
	    /* Because of the possible gaps, the index returned by search()
	       might not be the index of alignData array.  So the actual
	       index which was stored in refSearchRealIndexes array is used. */
	    searchResult.begin = refSearchRealIndexes[searchResult.begin];
	    searchResult.end = refSearchRealIndexes[searchResult.end];
	    begin = (int) (findConsBaseLocation(searchResult.begin) * zoomX);
	    end = (int) (findConsBaseLocation(searchResult.end) * zoomX);
	}

	boolean inRange = (begin >= x2 || end <= x1) ? false : true;
	if (inRange) {
            repaint();
	} else {
	    /* If the search result is not in the range of the current view,
	       scroll the view until it shows up. */
	    int x;
	    int viewportWidth = scrollPanel.getViewport().getWidth();
	    if (searchResult.direction == FORWARD) {
	    	x = Math.min(begin - 50, getWidth() - viewportWidth);
	    } else {
	    	x = Math.min(end + 50 - viewportWidth,
			     getWidth() - viewportWidth);
	    }
	    if (x < 0) {
		x = 0;
	    }
            scrollPanel.getViewport().scrollRectToVisible(
                                new Rectangle(x, p.y, 5, 5));
	    scrollPanel.getViewport().setViewPosition(new Point(x, p.y));
	}
    }

    /** Calculates the location of the consensus base with the specified
	index, which should be an actual index of alignData array.  This
        function is called by paintSearchResult(). */
    protected int findConsBaseLocation(int ind) {
	int fragPos, peakLoc, delsAfter, delsBefore, totalDels, i1, i2;
	char fragBase;
        fragPos = alignData[ind].getFragPosition();
        fragBase = alignData[ind].getFragBase();
        if (fragBase == '-') {
            // consecutive deletions after this one
            delsAfter = 0;
            for (i1=ind+1; i1<alignData.length; i1++) {
                if (alignData[i1].getFragBase()!='-') {
                    break;
                }
                delsAfter++;
            }
            // consecutive deletions before this one
            delsBefore = 0;
            for (i2=ind-1; i2 >= 0; i2--) {
                if (alignData[i2].getFragBase()!='-') {
                    break;
                }
                delsBefore++;
            }
            totalDels = delsBefore + 1 + delsAfter;
            peakLoc = phdPeakLocations[fragPos - 1]
                          + ((phdPeakLocations[fragPos]
                             - phdPeakLocations[fragPos -1])
                          	* (delsBefore + 1)) / (totalDels + 1);
        } else {
            peakLoc = phdPeakLocations[fragPos - 1];
        }

	return peakLoc;
    }

    /** Invoked when the scrollPanel is resized. */
    protected void scrollPanelResized() {
	if (!dataLoaded || !spFilesAreValid) {
	    return; 
	}
	dataChanged = true;
	int minHeight = Math.max(scrollPanel.getSize().height - 50,MIN_HEIGHT);
	int width = Math.max(traceLength, 
			     phdPeakLocations[phdPeakLocations.length - 1]);
	setPreferredSize(new Dimension((int)(width * zoomX) + 25, 
				       minHeight));
	revalidate();
	repaint();
	rowHeader.revalidate();
	rowHeader.repaint();
    }

    /** Add the specified DisplaySignalValueEventListener to 
	the listener vector. */
    public void addDisplaySignalValueEventListener(
				DisplaySignalValueEventListener l) {
	displaySigValEventListeners.add(l);
    }

    /* The last section of code is to provide public access to most of
     *  The data in the ChromViewPanel object, especially for ChromSnapShot.
     */
    public boolean displayingRaw() { return displayRaw; } 
    public boolean displayingOrigBases() { return displayOrig; }
    public boolean displayingAbc() { return displayAbc; }
    public boolean displayingAln() { return displayAln; }
    public boolean displayingTrim() { return displayTrim; }
    public boolean displayingPhd2() { return displayPhd2; }
    public boolean displayingTip() { 
	return tipFileLoaded && displayTip && (zoomX > 0.5);
    }
    public boolean displayingTotalIS() { 
	return tipFileLoaded && displayTotalIS && (zoomX > 0.5);
    }
    
    public String getProductName() { 
	return ("Paracel TraceTuner Version " + phdFile.getVersion());
    }
    public String getCurrChromFileName() { return chromFileName; }
    public String getCurrPhdFileName() { return phdFileName; }
    public String getCurrPhd2FileName() { return phd2FileName; }

    public float getZoomX() { return zoomX; }
    public int[] getMarkClosestPeak() { return markClosestPeak; }
    public ChromLayout getChromLayout() {return chromLayout; }

    public Trace[] getTraceLines() { return traceLines; }
    public int getTraceLength() { return traceLength; }
    public int getAnalyzedTraceLength() { return analyzedTraceLength; }
    public int getAnalyzedTraceMaxIntensity() { return analyzedTraceMaxIntensity; }
    public Trace[] getAnalyzedTraces() { return analyzedTraces; }
    public int[] getPhdPeakLocations() { return phdPeakLocations; }
    public int[] getQualityValues() { return qualityValues; }
    public boolean isTT() { return isTT; }
    public char[] getPhdBases() { return phdBases; }
    public int[] getPhdPeakLocations2() { return phdPeakLocations2; }
    public int[] getQualityValues2() { return qualityValues2; }
    public char[] getPhdBases2() { return phdBases2; }
    public boolean isTT2() { return isTT2; }
    public String getOrigBases() { return origBases; }
    public short[] getOrigBasePeakLocations() { return origBasePeakLocations; }
    public TabFile.AbcData[] getAbcData() { return abcData; }
    public TalFile.AlignData[] getAlignData() { return alignData; }
    public TipFile.IntrinsicPeak[] getIntrinsicPeakLines() {
	return intrinsicPeakLines;
    }
    public TipFile.IntrinsicPeak[] getTotalISPeakLines() {
	return totalISPeakLines;
    }
    public int getLeftTrimPoint() { return leftTrimPoint; }
    public int getRightTrimPoint() { return rightTrimPoint; }
    public int getTrimThreshold() { return trimThreshold; }

    public int getCurrViewWidth() { return frame.getSize().width - 5; }
    public int getCurrViewHeight() { return getSize().height; }
    public Point getCurrViewPos() { 
	return scrollPanel.getViewport().getViewPosition();
    }
}
