/*
 * 1.3
 *
 * @(#)TTPrintHeader.java       1.0     2000/11/22
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.io.*;
import java.util.*;
import java.lang.*;
import java.text.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.print.*;
import javax.swing.*;
import com.paracel.tt.util.*;

/**
 * This class encapsulates the printing of the top portion of TTPrintJob.
 * It contains info such as: Logo, Sample filename, Phd file name, date/time, 
 * page number etc.
 *
 * @see    	TTPrintJob 
 */
public class TTPrintHeader {

    /** The title font (for product name) */
    final static Font titleFont = new Font("Arial", Font.BOLD, 8);

    /** The content font */
    final static Font contentFont = new Font("Arial", Font.PLAIN, 8);

    /** The spacing constant between lines (or rows). */
    static final int lineSpace = -2;

    /** The extra spacing constant between lines (or rows). */
    static final int extraLineSpace = 2;

    /** The spacing constant between columns. */
    static final int columnSpace = 20;

    /** The logo image data (byte array) read from the JAR file. */
    private static byte[] logoImageData;

    /** The width and the height of the header area. */ 
    private int       	areaWidth, areaHeight;

    /** The x and the y of the origin. */ 
    private int       	originX, originY;

    /** The product name, sample file name, and phd file names. */
    private String 	productName, sampleFileName, phdFileName, phd2FileName;

    /** The current page number and total page number. */
    private int    	pageNumber, totalPageNumber;

    /** The time stamp string. */
    private String 	timeStamp;

    /** Indicates whether logo file exists. */
    private boolean	logoExists;

    /** The orientation of the page: Portrait or Landscape. */
    private int		pageOrientation;

    /** Indicates whether to print trace signal values. */
    private boolean	printTraceSignal;

    /** Indicates whether to print fifth dye singal value. */
    private boolean	printFifthDye;

    /** trace signal values. */
    private int[]	traceSignalValues;

    /** trim data (left trim point, right trim point, trim threshold) */
    private int[] trimData;

    /** The font metrics. */
    FontMetrics titleFontMetrics, contentFontMetrics;

    // Layout elements.
    int firstLineY, secondLineY, thirdLineY, fourthLineY, fifthLineY; 
    int column1StartX , column2StartX;

    /** Creates a TTPrintHeader object */
    public TTPrintHeader() {
	doRowLayout();
	setTimeStamp();
    }

    /** Lays out the rows: calculates the y position of each row. */
    protected void doRowLayout() {
	Container c = new Container();
	titleFontMetrics = c.getFontMetrics(titleFont);
    	contentFontMetrics = c.getFontMetrics(contentFont);
	c = null;

    	firstLineY = titleFontMetrics.getHeight() - 4;
    	secondLineY = firstLineY + (lineSpace + extraLineSpace) 
			+ contentFontMetrics.getHeight();
    	thirdLineY = secondLineY + lineSpace + contentFontMetrics.getHeight();
    	fourthLineY = thirdLineY + lineSpace + contentFontMetrics.getHeight();
    	fifthLineY = fourthLineY + lineSpace + contentFontMetrics.getHeight();
    }

    /** Prints the header onto the specified Graphics2D object. */
    public void print(Graphics2D g2D) {
	if (pageOrientation == PageFormat.LANDSCAPE) {
	    column2StartX = originX + areaWidth - 200;
	}
	drawLogo(g2D);
	drawProductName(g2D);
	drawSampleAndPhdFileName(g2D);
	if (phd2FileName != null) {
	    drawPhd2FileName(g2D);
	}
	/* Draws page number and timestamp before drawing trace signals
	 * because the position of trace signals depends on the page
	 * number and timestamp positions when the page orientation is
	 * LANDSCAPE.
	 */
	drawPageNumber(g2D);
	drawTimeStamp(g2D);
	if (printTraceSignal) {
	    drawTraceSignalValues(g2D);
	}
	if (trimData != null) {
	    drawTrimData(g2D);
	}
    }

    /** Getst the logo image data from JAR. */
    private static void getLogoImageData() {
	logoImageData = null;
	try {
	    InputStream in = LangTools.getZipEntryInputStream(
				"ttuner_tools.jar", "images/plogo.gif");
	    logoImageData = new byte[in.available()];
	    //in.read(logoImageData);
	    //This does not read all the data.
	    for (int i = 0; i < logoImageData.length; i++) {
		logoImageData[i] = (byte)in.read();
	    }
	} catch (Exception e) {} 
    }

    /** Draws the logo at the leftmost, if it exists. */
    private void drawLogo(Graphics2D g2D) {
	/* Draw the logo if it exists. */
	getLogoImageData();
	if (logoImageData == null) {
	    logoExists = false;
	} else {
	    ImageIcon tempLogo = new ImageIcon(logoImageData);
	    ImageIcon logo = new ImageIcon();

	    logoExists = true;
	    if (tempLogo.getIconWidth() > 75) {
		if (tempLogo.getIconHeight() > 35) {
		    logo = new ImageIcon(
				tempLogo.getImage().getScaledInstance(
				75, 35, Image.SCALE_DEFAULT));
		} else {
		    logo = new ImageIcon(
			    tempLogo.getImage().getScaledInstance(
			    75, tempLogo.getIconHeight(), 
			    Image.SCALE_DEFAULT));
		}
	    } else {
		 logo = tempLogo;
	    }

	    int xPosition = originX;
	    int yPosition = originY + (thirdLineY - logo.getIconHeight())/2;
	    g2D.drawImage(logo.getImage(), xPosition, yPosition, null);
	}
    }

 
    /** Draws the product Name to the right-top side of the logo. */
    private void drawProductName(Graphics2D g2D) {
	int xPosition = originX;
	int yPosition = originY + firstLineY;
	if (logoExists)
	    xPosition += 80;

	column1StartX = xPosition;
	int end = column1StartX + titleFontMetrics.stringWidth(productName);

	g2D.setColor(Color.black);
	g2D.setFont(titleFont);
	g2D.drawString(productName, xPosition, yPosition);
    }


    /** Draws the sample filename and phd filename below the product Name. */
    private void drawSampleAndPhdFileName(Graphics2D g2D) {
        int xPosition = column1StartX;
	int yPosition1 = originY + secondLineY;
	int yPosition2 = originY + thirdLineY;

	g2D.setColor(Color.black);
	g2D.setFont(contentFont);
	g2D.drawString(sampleFileName, xPosition, yPosition1);
	g2D.drawString(phdFileName, xPosition, yPosition2);
    }


    /** Draws the second phd filename below the phd file Name. */
    private void drawPhd2FileName(Graphics2D g2D) {
        int xPosition = column1StartX;
	int yPosition = originY + fourthLineY;

	g2D.setColor(Color.black);
	g2D.setFont(contentFont);
	g2D.drawString(phd2FileName, xPosition, yPosition);
    }

    /** Draws the page number at the rightmost. */
    private void drawPageNumber(Graphics2D g2D) {
	String pageNumString = "Page " + pageNumber + " of " + totalPageNumber;
	g2D.setColor(Color.black);
	g2D.setFont(contentFont);
        /* Two spaces are added to the end of "pageNumString" while 
	   calculating the xPosition in order to shift the rightmost column 
	   a little left.  This ensures all characters be printed out.
	 */
	int xPosition = originX + areaWidth 
			- contentFontMetrics.stringWidth(pageNumString+"   ");
	int yPosition = originY + firstLineY;

	g2D.drawString(pageNumString, xPosition, yPosition);
    }


    /** Draws the time stamp below the page number. */
    private void drawTimeStamp(Graphics2D g2D) {
	g2D.setColor(Color.black);
	g2D.setFont(contentFont);
        /* Three spaces are added to the end of "timeStamp" while calculating
           the xPosition in order to shift the rightmost column a little left.
           This ensures all characters be printed out.
	 */
	int xPosition = originX + areaWidth 
			- contentFontMetrics.stringWidth(timeStamp + "   ");
	int yPosition = originY + secondLineY;

	g2D.drawString(timeStamp, xPosition, yPosition);
    }


    /**
     * If the page orientation is Landscape:
     * <pre>Draws the signal values in the second column.
     * 	-- 1st row: "Signal A:125 C: 127 G:85 T:112 5thDye:1000"
     * 	-- 2nd row: "X Loc:377"</pre>
     * <p>
     * If the page orientation is Portrait:
     *  <pre>Draws the signal values in the frist column, 4th row.
     *  -- "Signal A:125 C:127 G:85 T:112 5thDye:1000  X Loc:377"</pre>
     */
    private void drawTraceSignalValues(Graphics2D g2D) {
	/* calculate the max width of this column (consider its first row 
	   as the longest */

	String locationXValue = "X Loc:" + traceSignalValues[0];
	String signals = "Signal A:" + traceSignalValues[1]
				 + " C:" + traceSignalValues[2]
				 + " G:" + traceSignalValues[3]
				 + " T:" + traceSignalValues[4];
	if (printFifthDye) {
	    signals += " 5thDye:" + traceSignalValues[5] + " ";
	}

	int widthFor3rdColumn = contentFontMetrics.stringWidth(timeStamp);
	int columnWidth = contentFontMetrics.stringWidth(signals);

	int xPosition1, xPosition2, yPosition1, yPosition2;

	if (pageOrientation == PageFormat.LANDSCAPE) {
	    xPosition1 = xPosition2 = column2StartX;
	    yPosition1 = originY + firstLineY;
	    yPosition2 = originY + secondLineY;
	} else {
	    // page orientation is PORTRAIT
	    xPosition1 = column1StartX;
	    xPosition2 = xPosition1 + contentFontMetrics.stringWidth(
						signals + "  ");
	    if (phd2FileName == null) {
	        yPosition1 = yPosition2 = fourthLineY;
	    } else {
	        yPosition1 = yPosition2 = fifthLineY;
	    }
	}

	drawSignalValues(g2D, xPosition1, yPosition1);
	
	g2D.drawString(locationXValue, xPosition2, yPosition2);
    }

    /** Draws trim data below the time stamp. */
    private void drawTrimData(Graphics2D g2D) {
	double prob = Math.pow(10.0, - ((double)trimData[2])/10.0);
	String temp = String.valueOf(prob);

	//Round the trim threshold value and make it an 8-char string.
	String trimThreshold;
	if (temp.length() == 8) {
	    trimThreshold = temp;
	} else if (temp.length() < 8) {
	    temp += "00000000";
	    trimThreshold = temp.substring(0, 8);
	} else {
	    char c = temp.charAt(8);
	    if (c >= '5') {
		char c1 = (char) (temp.charAt(7) + 1);
		trimThreshold = temp.substring(0, 7)  + c1;
	    } else {
		trimThreshold = temp.substring(0, 8);
	    }
	}

	String trimString = "Trim: " + trimData[0] + " "
	    			+ trimData[1] + " " + trimThreshold;
	int x = originX + areaWidth 
			- contentFontMetrics.stringWidth(trimString + "   ");
	int y = originY + thirdLineY;
	g2D.setFont(contentFont);
	g2D.drawString(trimString, x, y);
    }

    /** Draws the signal string for the specified base at the 
	specified (x,y). */
    private int drawSignalString(Graphics2D g, char base, int value,
				 int x, int y) {
	int xPos = x;

	// draw the base character
	String tempString = String.valueOf(base);
	g.setColor(BaseColorIdentifier.getColor(base));
	g.drawString(tempString, xPos, y);

	// draw the signal value of this base
	xPos += contentFontMetrics.stringWidth(tempString);
	tempString = ":" + value + " ";
	g.setColor(Color.black);
	g.drawString(tempString, xPos, y);

	// calculate the xPosition of the next string to be drawn
	xPos += contentFontMetrics.stringWidth(tempString);

	return xPos;
    }

    /** Draws the signal values at the specified (x, y). */
    private void drawSignalValues(Graphics2D g2D, int x, int y) {
	g2D.setFont(contentFont);

	// draw "Signal" in black
	int xPosition = x;

	g2D.setColor(Color.black);
	g2D.drawString("Signal ", xPosition, y);
	xPosition += contentFontMetrics.stringWidth("Signal ");

	xPosition = drawSignalString(g2D, 'A', traceSignalValues[1],
				     xPosition, y);
	xPosition = drawSignalString(g2D, 'C', traceSignalValues[2],
				     xPosition, y);
	xPosition = drawSignalString(g2D, 'G', traceSignalValues[3],
				     xPosition, y);
	xPosition = drawSignalString(g2D, 'T', traceSignalValues[4],
				     xPosition, y);
	if (printFifthDye) {
	    draw5thDyeSignalValues(g2D, xPosition, y);
	}
    }


    /** Draws the 5th Dye signal value, at the specified (x, y). */
    private void draw5thDyeSignalValues(Graphics2D g2D, int x, int y) {
	g2D.setFont(contentFont);

	int xPosition = x;
	String tempString = "5thDye";

	g2D.setColor(Color.lightGray);
	g2D.drawString(tempString, xPosition, y);

	xPosition += contentFontMetrics.stringWidth(tempString);
	tempString = ":" + traceSignalValues[5];

	g2D.setColor(Color.black);
	g2D.drawString(tempString, xPosition, y);
    }


    /** Sets the origin to the specified (x, y). */
    public void setOrigin(int x, int y) {
	this.originX = originX;
	this.originY= originY;
    }


    /** Sets the current page number to the specified value. */
    public void setPageNumber(int pageNumber) 
	    throws Exception {
	if (pageNumber <= 0)
	    throw new Exception("pageNumber must be a positive number");
	this.pageNumber = pageNumber;
    }


    /** Sets the header area height and width to the speicified values. */
    public void setArea(int areaHeight, int areaWidth)
	    throws Exception {
	if (areaHeight <= 0)
	    throw new Exception("areaHeight must be a positive number");
	if (areaWidth <= 0)
	    throw new Exception("areaWidth must be a positive number");
	this.areaHeight = areaHeight;
	this.areaWidth = areaWidth;
    }


    /** Sets the page orientation to the specified value. */
    public void setPageOrientation(int pageOrientation) {
	this.pageOrientation = pageOrientation;
    }

    /** Sets the total page number to the specified value. */
    public void setTotalPageNumber(int totalPageNumber)
	    throws Exception {
	if (totalPageNumber <= 0)
	    throw new Exception("totalPageNumber must be a positive number");
	this.totalPageNumber = totalPageNumber;
    }


    /** Sets the product name to the specified value. */
    public void setProductName(String pName) 
		throws Exception {
	if (pName == null)
	    throw new Exception("pName can not be null.");
	productName = pName;
    }


    /** Sets the sample file name to the specified value. */
    public void setSampleFileName(String sfName) 
		throws Exception {
	if (sfName == null)
	    throw new Exception("sfName can not be null.");
	sampleFileName = sfName;
    }


    /** Sets the phd file name to the specified value. */
    public void setPhdFileName(String pfName) 
		throws Exception {
	if (pfName == null)
	    throw new Exception("pfName can not be null.");
	phdFileName = pfName;
    }


    /** Sets the second phd file name to the specified value. */
    public void setPhd2FileName(String p2fName) {
	phd2FileName = p2fName;
    }


    /** Sets the time stamp to the current time (in String format). */
    public void setTimeStamp() {
	Date today = new Date();
	DateFormat dateFormat = DateFormat.getDateTimeInstance(
				DateFormat.SHORT, DateFormat.SHORT);
	timeStamp = dateFormat.format(new Date());
    }


    /**
     * Sets the trace signal printing flag and values.
     */
    public void setTraceSignalValues(boolean printTraceSignal,
				     boolean printFifthDye,
				     int[] traceSignalValues) 
		throws Exception {
	this.printTraceSignal = printTraceSignal;
	this.printFifthDye = printFifthDye;
	if (traceSignalValues != null && traceSignalValues.length != 6) {
	    throw new Exception("Array traceSignalValues has invalid length.");
	}
	this.traceSignalValues = traceSignalValues;
    }

    /** Sets the trim data. */
    public void setTrimData(int[] trim) { trimData = trim; }
}
