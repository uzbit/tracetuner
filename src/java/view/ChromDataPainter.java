/*
 * 1.12
 *
 * @(#)ChromDataPainter.java       1.0     2001/02/09
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.lang.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import com.paracel.tt.chrom.*;
import com.paracel.tt.io.*;
import com.paracel.tt.util.*;

/**
 * This class performs the actual drawing of the traces, bases, QVs, etc.
 */
public class ChromDataPainter implements SearchConstants {
    final static Font baseFont = new Font("Courier", Font.PLAIN, 12);
    final static Font misMatchedBaseFont = new Font("Monospaced",
						    Font.BOLD, 12);
    final static Font indexFont = new Font("Arial", Font.PLAIN, 10);
    final static Font qvFont = new Font("Arial", Font.PLAIN, 8);
    final static Font labelFont = new Font("Monospaced", Font.BOLD, 8);

    final static float dash1[] = {5.0f};
    final static BasicStroke dashed = new BasicStroke(1.0f,
        BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, dash1, 2.5f);
    final static BasicStroke stroke = new BasicStroke(1.0f);
    final static BasicStroke wideStroke = new BasicStroke(14.0f);

    /** ChromDataPainter is used to print to the chrom view panel. */
    public static final int VIEW_PANEL = 0;

    /** ChromDataPainter is used to print local view to print-outs. */
    public static final int PRINT_LOCAL = 1;

    /** ChromDataPainter is used to print global view to print-outs. */
    public static final int PRINT_GLOBAL = 2;

    /**
     * Draws a clip of the traces, bases, indexes, intrinsic peak lines
     * onto the specified Graphics2D object.
     * <p>The parameters specify the clip, and the necessary data.
     * This function can be used for both printing view and printing all,
     * given correct parameters.
     */
    public static void paint(Graphics2D g, int currXPos, int currWidth,
			     int extra, int traceLen, float zoomXVal,
			     ChromLayout cLayout, Trace[] traces,
			     int[] pPeakLocs, char[] pBases, int[] qvs, 
			     int[] pPeakLocs2, char[] pBases2, 
			     int[] qvs2, 
			     TabFile.AbcData[] abcData,
			     TalFile.AlignData[] alnData,
			     short[] oPeakLocs, String oBases,
			     TipFile.IntrinsicPeak[] iPeakLines,
			     TipFile.IntrinsicPeak[] totalISPeakLines,
			     boolean drawRaw, boolean drawAbc,
			     boolean drawAln, boolean drawOrig, 
			     boolean drawTip, boolean drawTotalIS,
			     boolean drawTrim, 
			     int leftTrim, int rightTrim, 
			     int trimThreshold,
			     boolean drawPhd2, int[] markPeak,
			     int type, Search.SearchResult sr,
			     boolean isTT, boolean isTT2) {
	int chromStart = (int) (currXPos / zoomXVal);
	int chromWidth = (int) (currWidth / zoomXVal);
	int actualTraceWidth = Math.max(traceLen - chromStart, 0);
	int width = Math.min(chromWidth + extra, actualTraceWidth);
	// If Print View, print labels at the right side of the actual
	// chrom data. So the actual graph is right shifted for 30 pixels.
	int graphOffset = (type == PRINT_LOCAL) ? 30 : 0;

	// draw traces
	int[] xPoints, yPoints;
	if (traces == null) {
	    return;
	}
	if (type == PRINT_GLOBAL) {
	    // print every other point of the trace to speed up printing
	    // and decrease the printed data load
	    xPoints = new int[(width+1)>>1];
	    yPoints = new int[(width+1)>>1];
	    for (int i= 0; i < width; i+=2) {
		xPoints[i>>1] = (int)(zoomXVal * (i + chromStart))+graphOffset;
	    }
	    int[] intensity;
	    for (int i = 0; i < traces.length; i++) {
		intensity = traces[i].getIntensityArray();
		for (int j = 0; j < width; j+=2) {
		    yPoints[j>>1] = intensity[chromStart+j];
		}
	    	g.setColor(traces[i].getColor());
	    	g.drawPolyline(xPoints, yPoints, xPoints.length);
	    }
	} else {
	    xPoints = new int[width];
	    yPoints = new int[width];
	    for (int i= 0; i < width; i++) {
	    	xPoints[i] = (int)(zoomXVal * (i + chromStart)) + graphOffset;
	    }
	    for (int i = 0; i < traces.length; i++) {
		System.arraycopy(traces[i].getIntensityArray(), chromStart,
			     	 yPoints, 0, width);
		g.setColor(traces[i].getColor());
		g.drawPolyline(xPoints, yPoints, width);
	    }
	}

	if (drawRaw) {
	    return;
	}

	if (pBases == null || pPeakLocs == null || qvs == null) {
	    return;
	}

	if (drawTrim) {
	    // draw yellow background for the trimmed regions
	    int trimLineY = cLayout.phdBasesY - 4;
	    int currEndX = currXPos + currWidth;
	    int leftTrimStart = -1;
	    int leftTrimEnd = -1;
	    int rightTrimStart = -1;
	    int rightTrimEnd = -1;
	    if (leftTrim == -1 && rightTrim == -1) {
		leftTrimStart = (int) (pPeakLocs[0] * zoomXVal);
		leftTrimEnd = (int) (pPeakLocs[pPeakLocs.length - 1] 
				       		* zoomXVal);
	    } else {
	    	if (leftTrim > 0) {
		    leftTrimStart = (int) (pPeakLocs[0] * zoomXVal) ;
		    leftTrimEnd = (int) (pPeakLocs[leftTrim - 1] * zoomXVal);
	    	}
		if (rightTrim < pPeakLocs.length - 1) {
		    rightTrimStart = (int) (pPeakLocs[rightTrim + 1] 
					    	* zoomXVal);
		    rightTrimEnd = (int) (pPeakLocs[pPeakLocs.length - 1] 
				      		* zoomXVal);
		}
	    }
	    drawLine(g, Color.yellow, wideStroke, leftTrimStart, leftTrimEnd,
		     currXPos + graphOffset, currEndX, trimLineY);
	    drawLine(g, Color.yellow, wideStroke, rightTrimStart, rightTrimEnd,
		     currXPos + graphOffset, currEndX, trimLineY);

	    g.setStroke(stroke);
	}

	//draw phd bases
	int phdPeakLocation, qvBarHeight, offset;
	
	int indexInterval = 10;
	for (int i = 0; i < pBases.length; i++) {
	    phdPeakLocation = (int) (pPeakLocs[i] * zoomXVal);
	    if (phdPeakLocation < currXPos) {
		continue;
	    } else if (phdPeakLocation > (currXPos + currWidth)){
		break;
	    } else {
		if (zoomXVal <= 0.2) {
		    indexInterval = 100;
		} 
		if (zoomXVal >= 5.0) {
		    indexInterval = 5;
		} 
		if (zoomXVal >= 10.0) {
		    indexInterval = 2;
		} 
		if (zoomXVal >= 30.0) {
		    indexInterval = 1;
		} 
		if ((i % indexInterval) == (indexInterval - 1)) {
		   g.setFont(indexFont);
		   g.setColor(Color.black);
		   g.drawString(String.valueOf(i+1), 
				phdPeakLocation+3+graphOffset, cLayout.indexY);
		   g.drawLine(phdPeakLocation + graphOffset, 
			      cLayout.indexLineY 
			      + cLayout.indexLineLength,
			      phdPeakLocation+graphOffset, cLayout.indexLineY);
		}

		// draw base
		g.setFont(baseFont);
		g.setColor(BaseColorIdentifier.getColor(pBases[i]));
		g.drawChars(pBases, i, 1, 
			    phdPeakLocation-3+graphOffset, cLayout.phdBasesY);

		// draw qv value and qv bar
		qvBarHeight = (type == PRINT_GLOBAL) ? qvs[i] / 3 : qvs[i];
		qvBarHeight = Math.max(1, qvBarHeight);
		if (zoomXVal > 0.2 || i%2 == 0) {
		    // if zoomXVal <= 0.2, print qv bar every other base
		    g.setColor(findQvBarColor(qvs[i]));
		    g.fillRect(phdPeakLocation-1+graphOffset, cLayout.qvBarY, 
			       3, qvBarHeight);
		}

		if (zoomXVal > 0.2 || i%20 == 0) {
		    // if zoomXVal <= 0.2, print qv value every 20 bases
		    g.setFont(qvFont);
		    g.setColor(Color.black);
		    offset = (qvs[i] > 9) ? 3 : 1;
		    g.drawString(String.valueOf(qvs[i]),
			         phdPeakLocation-offset+graphOffset, 
			         cLayout.qvBarY + qvBarHeight + 10);
		}

		// draw triangular search result mark
		if (sr != null && sr.patternFound && sr.subject == TT
		    	&& i >= sr.begin && i <= sr.end) {
		    markSearchResult(g, phdPeakLocation + graphOffset,
				 cLayout.phdBasesY - 10);
		}
	    }
	}

	if (drawPhd2 && pBases2 != null
	    	&& pPeakLocs2 != null && qvs2 != null) {
	    int phdPeakLocation2, qvBarHeight2;
	    for (int i = 0; i < pBases2.length; i++) {
	    	phdPeakLocation2 = (int) (pPeakLocs2[i] * zoomXVal);
	    	if (phdPeakLocation2 < currXPos) {
		    continue;
	    	} else if (phdPeakLocation2 > (currXPos + currWidth)){
		    break;
	    	} else {
		    // draw second phd base
		    g.setFont(baseFont);
		    g.setColor(BaseColorIdentifier.getColor(pBases2[i]));
		    g.drawChars(pBases2, i, 1, 
			    	phdPeakLocation2-3+graphOffset,
				cLayout.phdBases2Y);

		    // draw qv value and qv bar
		    qvBarHeight2 = (type == PRINT_GLOBAL) ? qvs2[i] / 3
							  : qvs2[i];
		    qvBarHeight2 = Math.max(1, qvBarHeight2);
		    if (zoomXVal > 0.2 || i%2 == 0) {
		    	// if zoomXVal <= 0.2, print qv bar every other base
			g.setColor(findQvBarColor(qvs2[i]));
		    	g.fillRect(phdPeakLocation2 - 1 + graphOffset,
				   cLayout.qvBar2Y, 3, qvBarHeight2);
		    }

		    if (zoomXVal > 0.2) {
		    	// if zoomXVal <= 0.2, do not print qv value
		    	g.setFont(qvFont);
		    	g.setColor(Color.black);
		        offset = (qvs2[i] > 9) ? 3 : 1;
		    	g.drawString(String.valueOf(qvs2[i]),
			             phdPeakLocation2 - offset + graphOffset, 
			             cLayout.qvBar2Y + qvBarHeight2 + 10);
		    }
		}
	    }
	}

	// draw the alternative bases
	if (drawAbc && abcData != null && pPeakLocs != null) {
	    char altBase;
	    int qv, position, index, pos, phdBasePos;
	    for (int i=0; i < abcData.length; i++) {
		altBase = abcData[i].getAltBase();
		qv = abcData[i].getQv();
		position = abcData[i].getPosition();
		index = abcData[i].getIndex();
		pos = (int) (position * zoomXVal);
		if (pos < currXPos) {
		    continue;
		} else if (pos > (currXPos + currWidth)) {
		    break;
		} else {
		    // draw alternative base call
		    g.setFont(baseFont);
		    g.setColor(BaseColorIdentifier.getColor(altBase));
		    g.drawString(altBase + "", pos - 3 + graphOffset, 
				 cLayout.altBasesY);

		    g.setFont(qvFont);
		    g.setColor(Color.black);
		    offset = (qv > 9) ? 3 : 1;
		    g.drawString(qv+"", pos - offset + graphOffset,
				 cLayout.altBasesY + 8);

		    if (abcData[i].getType() == TabFile.SUBSTITUTION) {
		    	// mark the substituted base
		    	g.setColor(Color.orange);
		    	phdBasePos = (int) (pPeakLocs[index] * zoomXVal);
		    	g.drawRect(phdBasePos - 4 + graphOffset, 
				   cLayout.phdBasesY - 10, 8, 12);
		    }
		}
	    }
	}
		
	// draw the consensus bases and index
	if (drawAln && alnData != null && pPeakLocs != null) {
            /* display the alignment read from .tal file
             * The consensus sequence will be shown at the top;
             * the original base calls will be in the middle, if any;
             * and the called bases will be at the bottom.
             */
	    int fragPos, consPos, peakLoc, index; 
	    int delsAfter, delsBefore, totalDels, i1, i2;
	    char fragBase, consBase;
            for (int i=0; i < alnData.length; i++) {
		if (alnData[i] == null) {
		    continue;
		}
                fragPos = alnData[i].getFragPosition();
                fragBase = alnData[i].getFragBase();
                consPos = alnData[i].getConsPosition();
                consBase = alnData[i].getConsBase();
                if (fragBase == '-') {
                    // consecutive deletions after this one
                    delsAfter = 0;
                    for (i1=i+1; i1<alnData.length; i1++) {
                        if (alnData[i1].getFragBase()!='-') {
                            break;
                        }
                        delsAfter++;
                    }
                    // consecutive deletions before this one
                    delsBefore = 0;
                    for (i2=i-1; i2 >= 0; i2--) {
                        if (alnData[i2].getFragBase()!='-') {
                            break;
                        }
                        delsBefore++;
                    }
                    totalDels = delsBefore + 1 + delsAfter;
                    peakLoc = pPeakLocs[fragPos - 1] +
                                  ((pPeakLocs[fragPos] -
                                    pPeakLocs[fragPos -1])
                                  * (delsBefore + 1)) / (totalDels + 1);
                } else {
                    peakLoc = pPeakLocs[fragPos - 1];
                }

		peakLoc = (int) (peakLoc * zoomXVal);
	        if (peakLoc < currXPos) {
		    continue;
	    	} else if (peakLoc > (currXPos + currWidth)){
		    break;
	   	} else {
		    if (consBase == '-') {
		    	g.setColor(Color.black);
		    	g.fillRect(peakLoc - 3 + graphOffset,
				   cLayout.consBasesY - 5, 6, 2);
		    } else {
		    	g.setColor(BaseColorIdentifier.getColor(consBase));
		        if (consBase == fragBase) {
		    	    g.setFont(baseFont);
			    g.drawString(".", peakLoc - 2 + graphOffset,
				         cLayout.consBasesY - 3);
		    	} else {
		    	    g.setFont(misMatchedBaseFont);
		    	    g.drawString(consBase + "",
					 peakLoc - 3 + graphOffset,
					 cLayout.consBasesY);
			}
		    }

		    if (fragBase == '-') {
		    	g.setColor(Color.black);
		    	g.fillRect(peakLoc - 3 + graphOffset,
				   cLayout.phdBasesY - 5, 6, 2);
		    }

		    /* draw the index of the consensus base in the
                     * consensus sequence on the top of the index
                     * of the fragment base in the fragment sequence.
                     */
                    index = fragPos - 1;
		    if ((index % indexInterval) == (indexInterval - 1)) {
		   	g.setFont(indexFont);
		   	g.setColor(Constants.gray108);
		   	g.drawString(String.valueOf(consPos), 
                                (int) (pPeakLocs[fragPos - 1] *
                                	zoomXVal + 3) + graphOffset,
				cLayout.consIndexY);
                    }
	    	}

		// draw triangular search result mark
		if (sr != null && sr.patternFound && sr.subject == REF
		    	&& i >= sr.begin && i <= sr.end) {
		    markSearchResult(g, peakLoc + graphOffset,
				 cLayout.consBasesY - 10);
		}
            }           
	}
	
	// draw original bases
	if (drawOrig && oPeakLocs != null && oBases != null) {
	    int origPeakLocation;
	    for (int i = 0; i < oBases.length(); i++) {
	    	origPeakLocation = (int) (oPeakLocs[i] * zoomXVal);
	    	if (origPeakLocation < currXPos) {
		    continue;
	    	} else if (origPeakLocation > (currXPos + currWidth)){
		    break;
	    	} else {
		    g.setFont(baseFont);
		    g.setColor(BaseColorIdentifier.getColor(
							oBases.charAt(i)));
		    g.drawString(oBases.substring(i, i+1),
				 origPeakLocation - 3 + graphOffset,
				 cLayout.origBasesY);
		}

		// draw triangular search result mark
		if (sr != null && sr.patternFound && sr.subject == ABI
		    	&& i >= sr.begin && i <= sr.end) {
		    markSearchResult(g, origPeakLocation + graphOffset,
				 cLayout.origBasesY - 10);
		}
	    }
	}

	//draw intrinsic peaks
	if (drawTip && iPeakLines != null) {
	    g.setStroke(dashed);

	    int peakLen, begin, end;
	    int[] peakXPoints;
	    for(int i = 0; i < iPeakLines.length; i++) {
		peakLen = iPeakLines[i].length;
		begin = (int)(iPeakLines[i].x[0] * zoomXVal);
		end = (int)(iPeakLines[i].x[peakLen-1] * zoomXVal);
		if (end <= currXPos) {
		    continue;
		} else if (begin >= currXPos + currWidth) {
		    break;
		} else {
		    peakXPoints = new int[peakLen];
	            for (int j = 0; j < peakLen; j++) {
	    	    	peakXPoints[j] = (int)(iPeakLines[i].x[j] * zoomXVal)
			    			+ graphOffset;
	            }
		    g.setColor(BaseColorIdentifier.getColor(
					iPeakLines[i].base));
		    g.drawPolyline(peakXPoints, iPeakLines[i].y, peakLen);
		    peakXPoints = null;
		}
	    }
	    
	    g.setStroke(stroke);
	}

	//draw total intrinsic signal
	if (drawTotalIS && totalISPeakLines != null) {
	    g.setStroke(dashed);

	    int peakLen, begin, end;
	    int[] peakXPoints;
	    for (int i = 0; i < totalISPeakLines.length; i++) {
		peakLen = totalISPeakLines[i].length;
		begin = (int)(totalISPeakLines[i].x[0] * zoomXVal);
		end = (int)(totalISPeakLines[i].x[peakLen-1] * zoomXVal);
		if ((end < currXPos) || (begin > currXPos + currWidth)) {
		    continue;
		} else {
		    peakXPoints = new int[peakLen];
	            for (int j = 0; j < peakLen; j++) {
	    	    	peakXPoints[j] = (int)(totalISPeakLines[i].x[j]
					       * zoomXVal) + graphOffset;
	            }
		    g.setColor(BaseColorIdentifier.getColor(
					totalISPeakLines[i].base));
		    g.drawPolyline(peakXPoints, 
				   totalISPeakLines[i].y, peakLen);
		    peakXPoints = null;
		}
	    }
	    
	    g.setStroke(stroke);
	}

	//draw the closest peak location
	if (markPeak != null && markPeak[0] >= 0 && markPeak[1] >= 0) {
	    int drawXPos = (int) (markPeak[1] * zoomXVal);
	    int y = 0;
	    boolean errFound = false;
	    switch (markPeak[0]) {
	    case 0:
		y = cLayout.origBasesY; break;
	    case 1:
		y = cLayout.phdBasesY; break;
	    case 2:
		y = cLayout.altBasesY; break;
	    case 3:
		y = cLayout.phdBases2Y; break;
	    default:
		errFound = true; break;
	    }

	    if (!errFound) {
		g.setColor(Color.black);
		g.drawLine(drawXPos + graphOffset, cLayout.peakLocY,
		           drawXPos + graphOffset, y);

		g.setFont(indexFont);
		g.drawString(String.valueOf(markPeak[1]),
			     drawXPos + 2 + graphOffset, cLayout.peakLocY);
	    }
	}

	if (type == PRINT_LOCAL) {
	    if (drawTip || drawTotalIS) {
	    	/* Erase the extra pieces of the intrinsic and total intrinsic
	    	   lines at the left end.  The reason that intrinsic signals
		   have extra piece drawn at front is the way they were drawn:
		   As long as part of the peak is in the view range, the
		   entire peak is drawn. */
	    	g.setColor(Color.white);
	    	g.fillRect(currXPos, cLayout.traceLinesY-cLayout.chromHeight-2,
		           graphOffset, cLayout.chromHeight + 12);
	    }

	    //draw the labels indicating the role of each row
	    g.setFont(labelFont);
	    g.setColor(Color.black);
	    if (!drawRaw) {
		if (isTT) {
		    g.drawString("[tt]", currXPos + 2, cLayout.indexY);
		    g.drawString("[TT]", currXPos + 2, cLayout.phdBasesY);
		} else {
		    g.drawString("[ph]", currXPos + 2, cLayout.indexY);
		    g.drawString("[PH]", currXPos + 2, cLayout.phdBasesY);
		}
	    }
	    if (drawAln) {
		g.drawString("[ref]", currXPos, cLayout.consIndexY);
		g.drawString("[REF]", currXPos, cLayout.consBasesY);
	    }
	    if (drawOrig) {
		g.drawString("[ABI]", currXPos, cLayout.origBasesY);
	    }
	    if (drawAbc) {
		g.drawString("[ALT]", currXPos, cLayout.altBasesY);
	    }
	    if (drawPhd2) {
		if (isTT2) {
		    g.drawString("[TT]", currXPos + 2, cLayout.phdBases2Y);
		} else {
		    g.drawString("[PH]", currXPos + 2, cLayout.phdBases2Y);
		}
	    }
	}
    }

    /** Draws a line using the specified color and stroke on the
	specfied Graphics2D object starting at start and ending at
	end.  Consider the start and end of the current view window
	(viewStart and viewEnd). */
    protected static void drawLine(Graphics2D g, Color c, BasicStroke s,
			      int start, int end, 
			      int viewStart, int viewEnd,
			      int y) {
	if (start != -1 && end != -1) {
	    if (end < viewStart || start > viewEnd) {
		return; 
	    }
	    g.setColor(c);
	    g.setStroke(s);
	    if (start >= viewStart && end <= viewEnd) {
		g.drawLine(start, y, end, y);
	    } else if (start >= viewStart) {
		g.drawLine(start, y, viewEnd, y);
	    } else {
		g.drawLine(viewStart, y, end, y);
	    }
	}
    }

    /** Draws a black solid arrow above the bases within the search result. */
    protected static void markSearchResult(Graphics2D g, int x, int y){
	g.setColor(Color.black);
	int[] xPoints = { x - 4, x, x + 4 };
	int[] yPoints = { y - 4, y, y - 4 };
	g.fillPolygon(xPoints, yPoints, 3);
    }

    /** Returns the Color object corresponding to the specified qv. */
    protected static Color findQvBarColor(int qv) {
	if (qv < 20) {
	    return Color.lightGray;
	} else if (qv < 30) {
	    return Color.gray;
	} else {
	    return Color.darkGray;
	}
    }
}
