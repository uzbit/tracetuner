/*
 * 1.5
 *
 * @(#)ChromLayout.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.awt.*;
import com.paracel.tt.util.*;

/**
 * This class contains the lay-out data of the elements, such as indices,
 * phd base calles, original base calls, quality values, chromatograms, etc.
 * This determines how these elements are layed out along y-axis in the viewer.
 */
public class ChromLayout implements Cloneable {

    /** The Y position of the index string */
    public int indexY;

    /** The Y position of the consensus base index string */
    public int consIndexY;

    /** The Y position of the index marker line */
    public int indexLineY;

    /** The length of the index marker line */
    public int indexLineLength;

    /** The Y position of the original base characters */
    public int origBasesY;

    /** The Y position of the consensus base characters */
    public int consBasesY;

    /** The Y position of the phd base characters */
    public int phdBasesY;

    /** The Y position of the alternative base characters */
    public int altBasesY;

    /** The Y position of the top of the quality value bars. */
    public int qvBarY;

    /** The Y position of the bottom of the trace lines. */
    public int traceLinesY;

    /** The height of the area for drawing the traces. */
    public int chromHeight;

    /** The Y position of the peak location string. */
    public int peakLocY;

    /** The Y position of the 2nd phd base characters if 2 phd files are
	shown in the same viewer. */
    public int phdBases2Y;

    /** The Y position of the bottom of the quality value bars 
	of the 2nd phd bases. */
    public int qvBar2Y;

    /** Creates a clone of the current ChromLayout object. */
    public Object clone() {
	ChromLayout newObject = new ChromLayout();
	newObject.indexY = indexY;
	newObject.origBasesY = origBasesY;
	newObject.consBasesY = consBasesY;
	newObject.phdBasesY = phdBasesY;
	newObject.qvBarY = qvBarY;
	newObject.altBasesY = altBasesY;
	newObject.phdBases2Y = phdBases2Y;
	newObject.qvBar2Y = qvBar2Y;
	newObject.traceLinesY = traceLinesY;
	newObject.chromHeight = chromHeight;
	newObject.peakLocY = peakLocY;
	newObject.consIndexY = consIndexY;
	newObject.indexLineY = indexLineY;
	newObject.indexLineLength = indexLineLength;

	return newObject;
    }

    /**
     * Creates a ChromLayout for Print-All.
     * @return the Print-All chrom layout.
     * @param	displayPhd2	whether showing 2 phd files in one viewer.
     */
    public static ChromLayout createPrintAllLayout(boolean displayAbc,
					    	   boolean displayPhd2,
						   int height) {
	    ChromLayout printAllLayout = new ChromLayout();
	    Font baseFont = Fonts.BASE_FONT;
	    int baseLetterHeight = Fonts.getFontMetrics(baseFont).getHeight();
	    int qvHeight = 33 + baseLetterHeight + 5;

	    printAllLayout.traceLinesY = height - 10;
	    printAllLayout.indexY = baseLetterHeight;
	    printAllLayout.indexLineLength = 5;
	    printAllLayout.indexLineY = printAllLayout.indexY 
					+ printAllLayout.indexLineLength;
	    if (displayAbc) {
	        printAllLayout.altBasesY = printAllLayout.indexLineY
					    + baseLetterHeight;
	        printAllLayout.phdBasesY = printAllLayout.indexLineY
					    + 2 * baseLetterHeight + 8;
	    } else {
	        printAllLayout.phdBasesY = printAllLayout.indexLineY
					    + baseLetterHeight;
	    }
	    printAllLayout.qvBarY = printAllLayout.phdBasesY + 2;

	    if (displayPhd2) {
	    	printAllLayout.phdBases2Y = printAllLayout.qvBarY
		    			+ qvHeight;
	    	printAllLayout.qvBar2Y = printAllLayout.phdBases2Y + 2;
	    	printAllLayout.chromHeight = printAllLayout.traceLinesY
		    			     - printAllLayout.qvBar2Y
		    			     - qvHeight;
	    } else {
	    	printAllLayout.chromHeight = printAllLayout.traceLinesY
		    			     - printAllLayout.qvBarY
		    			     - qvHeight;
	    }

	return printAllLayout;
    }
}
