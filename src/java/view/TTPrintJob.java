/*
 * 1.7
 *
 * @(#)TTPrintJob.java       1.0     2000/11/22
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.awt.print.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;
import com.paracel.tt.util.*;

/**
 * This class creates a print job for the printing in the viewer.
 * There are two kinds of print jobs supported for TraceTuner Viewer:
 * PRINT_LOCAL and PRINT_SINGLE_GLOBAL.
 *
 * @see    	ChromSnapShot
 * @see    	TTPrintHeader 
 */
public class TTPrintJob implements Printable, Pageable {
    
    /** The print-global-view-of-current-file TTPrintJob type constant. */
    public static final int PRINT_SINGLE_GLOBAL = 0;

    /** The print-global-view-of-all-files type constant. */
    public static final int PRINT_ALL_GLOBAL = 1;

    /** The print-local-view TTPrintJob type constant. */
    public static final int PRINT_LOCAL = 2;

    /** The height of each row for print-all jobs. */
    public static final int ROW_HEIGHT = 210;

    /** The height of the header portion (which has the file name, signal
	value, page number, product name, time stamp information). */
    private int headerHeight = 45;

    /** The height of the extra space betweeen the header and body. */
    private int extraSpaceHeight = 10;

    /** The source of the data for print.  Here is the viewer window. */
    private Component source;

    /** The print job. */
    private PrinterJob printerJob;

    /** The snapshot to be printed. */
    private ChromSnapShot chromSnapShot;

    /** The total number of pages to be printed. */
    private int totalNumOfPages;

    /** The TTPrintJob type.  (PRINT_LOCAL, PRINT_SINGLE_GLOBAL,
	or PRINT_ALL_GLOBAL) */
    private int type;

    /** The page format used for printing. */
    private PageFormat pageFormat;

    /** The object which takes care of the header printing. */
    private TTPrintHeader printHeader;

    /** The label to show the printing status. */
    private JLabel printStatusLabel;

    /** Indicates whether the printing as cancelled. */
    private boolean printCancelled;

    /** The scaling factors on x and y axis to fit the drawing to the paper. */
    private double scaleX, scaleY;

    /** index of a row on a page */
    private int rowIndexOnPage;

    /** total number of rows */
    private int totalRowCounts;

    /** The number of rows to print on one page for print-all jobs. */
    private int paRowPerPage;

    /** The number of base calls to print on one row for print-all jobs. */
    private int paRowWidth;

    private boolean isFirstFile = true;

    /**
     * Creates a TTPrintJob from the specified ChromSnapShot.
     * @param	c	the source.
     * @param	css	the snapshot to be printed.
     * @param	pf	the page format to be used.
     */
    public TTPrintJob(Component c, ChromSnapShot css) {
	source = c;
	chromSnapShot = css;
	type = css.getType();

	PageSetupManager manager = PageSetupManager.currentManager();
	this.pageFormat = (PageFormat) manager.getPageFormat().clone();
	if (type == PRINT_SINGLE_GLOBAL || type == PRINT_ALL_GLOBAL) {
	    this.paRowPerPage = manager.getPARowNum();
	    this.paRowWidth = Constants.AVERAGE_SCANS_PER_BASE 
	    			* manager.getPABaseNum();
	}

	calcNumOfPages();
	printHeader = new TTPrintHeader();
	// Sets data that won't change from page to page
	try {
	    printHeader.setTotalPageNumber(totalNumOfPages);
	    printHeader.setProductName(chromSnapShot.getProductName());
	    printHeader.setSampleFileName(chromSnapShot.getChromFileName());
	    printHeader.setPhdFileName(chromSnapShot.getPhdFileName());
	    printHeader.setTrimData(chromSnapShot.getTrimData());
	    printHeader.setPhd2FileName(chromSnapShot.getPhd2FileName());
	    if (type == PRINT_LOCAL) {
	    	printHeader.setTraceSignalValues(
				chromSnapShot.signalValuesShown(),
				chromSnapShot.fifthDyeDataShown(),
				chromSnapShot.getTraceSignalValues());
	    }
	} catch (Exception e) {
	    System.err.println("Error setting printHeader" + e.toString());
	    e.printStackTrace(System.err);
	}
    }

    public TTPrintJob(Component c, ChromSnapShot css, boolean isFirstFile) {
	this(c, css);
	this.isFirstFile = isFirstFile;
    }

    /** Sets the parent component for the message dialogs. */
    public void setParent(Component c) {
	source = c;
    }

    /** Calculates the total number of pages to be printed. */
    protected void calcNumOfPages() {
	if (!chromSnapShot.hasValidSampleAndPhd()) {
	    totalNumOfPages = 0; 
	}
	if (type == PRINT_LOCAL) {
	    totalNumOfPages = 1; 
	} else {
	    // PRINT_SINGLE_GLOBAL or PRINT_ALL_GLOBAL
            double imageableWidth = pageFormat.getImageableWidth();
            double imageableHeight = pageFormat.getImageableHeight();

	    /* calculates the scale on  X axis and Y axis */
	    scaleX = imageableWidth / paRowWidth;
	    scaleY = (imageableHeight - headerHeight - extraSpaceHeight) 
			    / (paRowPerPage * ROW_HEIGHT);

	    /* calculates the total row number and total page number */
	    double traceLen = (double) chromSnapShot.getAnalyzedTraceLength();
            totalRowCounts = (int)Math.ceil(traceLen / (double) paRowWidth);
            totalNumOfPages = (int)Math.ceil((double)totalRowCounts 
				    / (double)paRowPerPage);
	}

    }

    /** Releases the memory resources. */
    public void dispose() {
	if (chromSnapShot != null) {
	    chromSnapShot.dispose();
	    chromSnapShot = null;
	}
	printerJob = null;
	pageFormat = null;
	printHeader = null;
    }

    /** 
     * Gets the total number of pages in the set. Required by 
     * <code>Pageable</code> interface.
     * @return	the total number of pages in the set.
     * @see	Pageable
     */
    public int getNumberOfPages() {
	return totalNumOfPages;
    }


    /** 
     * Gets the PageFormat of the page specified by pageIndex. 
     * Required by <code>Pageable</code> interface.
     * @return	the PageFormat instance.
     * @param	pageIndex	the specified page index.
     * @see	Pageable
     */
    public PageFormat getPageFormat(int pageIndex) {
	return pageFormat;
    }

    /** 
     * Gets the Printable instance responsible for rendering the page
     * specified by pageIndex. Required by <code>Pageable</code> interface.
     * @return	the Printable instance.
     * @param	pageIndex	the specified page index.
     * @see	Pageable
     */
    public Printable getPrintable(int pageIndex) {
	return this;
    }

    /**
     * Creates the printing status frame, does the preparation for the actual
     * printing.  This is the function to be called on a TTPrintJob object
     * to perform a TTPrintJob printing action.
     */
    public void print() {
	printCancelled = false;
	printerJob = PrinterJob.getPrinterJob();
	printerJob.setPageable(this);

	String typeString;
	if (type == PRINT_LOCAL) {
	    typeString = "Print Local View";
	} else if (type == PRINT_SINGLE_GLOBAL) {
	    typeString = "Print Global View of Current File";
	} else {
	    typeString = "Print Global View of All Files";
	}

	JFrame printFrame = new JFrame(typeString + " Status");
	JPanel contentPane = new JPanel();
	printFrame.setContentPane(contentPane);
	contentPane.setPreferredSize(new Dimension(350, 75));
	contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
	contentPane.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));

	printStatusLabel = new JLabel();
        printStatusLabel.setText("");
        printStatusLabel.setAlignmentX(JComponent.CENTER_ALIGNMENT);
        contentPane.add(printStatusLabel);
        contentPane.add(Box.createRigidArea(new Dimension(0,10)));
        JButton cancel = new JButton("Cancel");
        cancel.addActionListener(new ActionHandler());
        cancel.setAlignmentX(JComponent.CENTER_ALIGNMENT);
        contentPane.add(cancel);
        try {
            printFrame.pack();
        } catch (Exception ignore) {
            // Do nothing.
        }

	if (isFirstFile) {
	    if (!printerJob.printDialog()) {
		printCancelled = true;
		return;
	    }
	}
	if (type != PRINT_ALL_GLOBAL) {
	    source.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
	}
	try {
            printFrame.setVisible(true);
            printerJob.print();
        } catch(PrinterAbortException pae) {
            if (!printCancelled) {
                JOptionPane.showMessageDialog(source, pae.getMessage()
					+ "\nPlease check printer settings.",
					"Print Error",
					JOptionPane.INFORMATION_MESSAGE);
                System.err.println("Error printing: " + pae);
		pae.printStackTrace(System.err);
		printCancelled = true;
	    }
	} catch(PrinterIOException pioe) {
            JOptionPane.showMessageDialog(source,
					pioe.getMessage()
					+ "\nPlease check printer settings.",
					"Print Error",
					JOptionPane.INFORMATION_MESSAGE);
            System.err.println("Error printing: " + pioe);
	    pioe.printStackTrace(System.err);
	    printCancelled = true;
        } catch(PrinterException pe) {
            JOptionPane.showMessageDialog(source,
					pe.getMessage()
					+ "\nPlease check printer settings.",
					"Print Error",
					JOptionPane.INFORMATION_MESSAGE);
            System.err.println("Error printing: " + pe);
	    pe.printStackTrace(System.err);
	    printCancelled = true;
        }

        if (type != PRINT_ALL_GLOBAL) {
       	    source.setCursor(
			 Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	}
        printFrame.dispose();
        printFrame = null;

        if ((!printCancelled) && (type != PRINT_ALL_GLOBAL)) {
            JOptionPane.showMessageDialog(source,
                        typeString + " completed", typeString,
                        JOptionPane.INFORMATION_MESSAGE);
	}
    }

    
    /**
     * Prints the page at the specified index into the specified 
     * <code>Graphics</code> context in the specified format.
     * @return	NO_SUCH_PAGE if the pageIndex is too large, or the
     *		requested page does not exist, or print was canceled;
     *		PAGE_EXISTS to signify that the requested page was rendered.
     */
    public int print(Graphics g, PageFormat pageFormat, int pageIndex) {
	if (printCancelled) {
	    return NO_SUCH_PAGE;
	}

	if (pageIndex < 0 || pageIndex >= totalNumOfPages) {
	    return NO_SUCH_PAGE;
	}

        double imageableWidth = pageFormat.getImageableWidth();
        double imageableHeight = pageFormat.getImageableHeight();
        double imageableX = pageFormat.getImageableX();
        double imageableY = pageFormat.getImageableY();
	
	Graphics2D g2d = (Graphics2D)g;
	g2d.setColor(Color.white);
	g2d.fillRect((int) imageableX, (int) imageableY,
		     (int) imageableWidth, (int) imageableHeight);

        int pageNum = pageIndex+1;
	String text;
        if (type == PRINT_ALL_GLOBAL) {
            text = "Now printing file " + chromSnapShot.getCurrentFileIndex()
	    			+ " of " + chromSnapShot.getTotalNumOfFiles();
	} else {
            text = "Now printing page " + pageNum
	    			+ " of " + totalNumOfPages;
	}
	if (printStatusLabel != null) {
            printStatusLabel.setText(text);
	}

	g2d.translate(imageableX, imageableY);

	/* Print the header portion. */
	try {
	    printHeader.setArea(headerHeight, (int) imageableWidth);
	    printHeader.setOrigin(0, 0);
	    printHeader.setPageOrientation(pageFormat.getOrientation());
	    printHeader.setPageNumber(pageNum);
	} catch (Exception e) {
	    System.err.println("Error printing header" + e.toString());
	    e.printStackTrace(System.err);
	}

	printHeader.print(g2d);

	// Print the body
	double shiftX, shiftY;
	if (type == PRINT_LOCAL) {
	    /* Print the trace diagrams. */
            scaleX = imageableWidth 
	    		/ ((double)chromSnapShot.getViewWidth());

            scaleY = (imageableHeight - headerHeight - extraSpaceHeight)
	    		/ ((double)chromSnapShot.getViewHeight());

            shiftX = -(chromSnapShot.getViewPos().getX())*scaleX;
            shiftY = (double)(headerHeight + extraSpaceHeight);

            g2d.translate(shiftX, shiftY);
            g2d.setClip((int)-shiftX, 0, (int) Math.ceil(imageableWidth),
                        (int)Math.ceil(imageableHeight));
            g2d.scale(scaleX, scaleY);

	    chromSnapShot.drawContentPrintView(g2d);

	    g2d.scale(1/scaleX, 1/scaleY);
	    g2d.translate(-shiftX, -shiftY);

	} else {
            // Print the trace diagrams row by row.
            rowIndexOnPage = 0;
            // rowIndex indicates the row index among all the rows
           int rowIndex = pageIndex * paRowPerPage;
           int rowHeight = ((int)imageableHeight - headerHeight
                                - extraSpaceHeight)/paRowPerPage;

            while (rowIndexOnPage < paRowPerPage) {
                if ((rowIndexOnPage >= totalRowCounts)
                            || (rowIndex >= totalRowCounts)) {
                    return PAGE_EXISTS;
		}

                shiftX = -rowIndex * imageableWidth;
                shiftY = (double) (headerHeight + extraSpaceHeight
                                       + rowIndexOnPage * rowHeight);

                g2d.translate(shiftX, shiftY);
                g2d.setClip((int)(imageableWidth * rowIndex),
                            0,
                            (int)Math.ceil(imageableWidth),
                            (int)Math.ceil(imageableHeight));
                g2d.scale(scaleX, scaleY);

                int xPosition = paRowWidth * rowIndex;

	        chromSnapShot.drawContentPrintAll(g2d, xPosition, paRowWidth);

                g2d.scale(1/scaleX, 1/scaleY);
                g2d.translate(-shiftX, -shiftY);

                rowIndexOnPage++;
                rowIndex++;
            }
	}

	return PAGE_EXISTS;
    }

    /**
     * The action handler class to handle the situation when user clicks
     * "Cancel" button in print status window to cancel the current
     * print job.
     */
    public class ActionHandler implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            printerJob.cancel();
            printCancelled = true;
	    String s;
	    if (type == PRINT_LOCAL) {
		s = "Cancelling Print Local View ...";
	    } else if (type == PRINT_SINGLE_GLOBAL) {
	    	s = "Cancelling Print Global View of Current File ...";
	    } else {
	    	s = "Cancelling Print Global View of All Files ...";
	    }
            printStatusLabel.setText(s);
	    JButton b = (JButton) e.getSource();
            b.setEnabled(false);
        }
    }

    public boolean isCanceled() { 
	return printCancelled;
    }
}
