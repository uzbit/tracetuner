/**
	$(CLASS_DIR)/PrintPreview.class \
 * 1.5
 *
 * @(#)PageSetupDialog.java	1.0	2000/11/14
 * 
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.lang.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import javax.swing.*;
import javax.swing.border.*;
import com.paracel.tt.util.*;

/**
 * Creates the page setup dialog, which keeps track of parameters including
 * Paper Size, Orientation, Margins, and two parameters for Print Global View:
 * Approximate Number of Bases per Row, and Number of Rows per Page.
 * The default value of these parameters are: US Letter, Portrait,
 * 1.0", 80, and 6, respectively.
 * <p>The dialog is invoked by PageSetupManager class through its class
 * method, PageSetupManager.pageDialog().  Actually this class is an
 * auxillary class of the PageSetupManager class, which manages page setup
 * parameters application-wise.
 */
public class PageSetupDialog extends JDialog implements ActionListener {

    /** The dialog. */
    private static PageSetupDialog dialog;

    /** The original page format before user's input. */
    private static PageFormat origPageFormat;
   
    /** The original (Print Global View) row num per page and base num per row. */
    private static int origRowNum, origBaseNum;

    /** The (Print Global View) row num per page and base num per row. */
    private static int rowNum, baseNum;

    /** The original page size and the current page size. */
    private static PaperSize origPaperSize, paperSize;

    /** The orientation to print. */
    private static int orientation;

    /** The paper size box. */
    private JComboBox paperSizeBox;

    /** The orientation check boxes. */
    private JRadioButton portrait, landscape;

    /** The text field for the margins. */
    private JTextField left, right, top, bottom, rowNumField, baseNumField;

    /** The "OK", and "Cancel" button. */
    private JButton btnOK, btnCancel;

    /** The margins. */
    private static double leftMargin, rightMargin, topMargin, bottomMargin;

    /** Whether or not the dialog was canceled. */
    private static boolean canceled;

    /** Creates a PageSetupDialog that request inputs from the user.  */
    public PageSetupDialog() {
	super();
	setTitle("Page Setup");
	setModal(true);
	setSize(new Dimension(200, 150));
	setResizable(false);

	JPanel c = (JPanel) getContentPane();
	c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));

	c.add(Box.createVerticalStrut(5));
	c.add(createPaperPanel());
	c.add(Box.createVerticalStrut(5));
	JPanel p = new JPanel();
	p.setLayout(new BoxLayout(p, BoxLayout.X_AXIS));
	p.add(createOrientationPanel());
	p.add(createMarginsPanel());
	c.add(p);
	c.add(Box.createVerticalStrut(5));
	c.add(createPrintGlobalViewPanel());
	c.add(Box.createVerticalStrut(10));
	c.add(createButtonsPanel());

	pack();
    }

    /** Called by PageSetupManager to invoke the page setup dialog. 
        Returns an 3-element Object array containing the PageFormat,
        the Print Global View base num per row, and the Print Global View row
        num per page.  */
    public static Object[] pageDialog(PageFormat pf, int base, int row) {
	if (dialog == null) {
	    dialog = new PageSetupDialog();
	}
	if (pf == null) {
            origPageFormat = defaultPageFormat();
	    origPaperSize = PaperSize.LETTER;
	} else {
	    origPageFormat = (PageFormat) pf.clone();
	    origPaperSize = PaperSize.getPaperSize(origPageFormat);
	}
	origRowNum = (row<=0) ? Constants.PA_DEFAULT_ROW_NUM_PER_PAGE : row;
	origBaseNum = (base<=0) ? Constants.PA_DEFAULT_BASE_NUM_PER_ROW : base;
	dialog.reset();
	dialog.setVisible(true);

	if (canceled) {
	    Object o[] = { origPageFormat, new Integer(origBaseNum), 
		           new Integer(origRowNum) };
	    return o;
	} else {
	    Object o[] = { getPageFormat(), new Integer(baseNum), 
		           new Integer(rowNum) };
	    return o;
	}
    }

    /** Creates the paper size panel. */
    protected JPanel createPaperPanel() {
	JPanel p = new JPanel();
	p.setLayout(new BoxLayout(p, BoxLayout.X_AXIS));

	String osName = Constants.OS_NAME.trim().toLowerCase();
	if (osName.indexOf("windows") >= 0) {
	    PaperSize[] sizes = { PaperSize.LETTER, PaperSize.LEGAL,
				  PaperSize.A3, PaperSize.A4, PaperSize.A5,
				  PaperSize.JIS_B4, PaperSize.JIS_B5 };
	    paperSizeBox = new JComboBox(sizes);
	} else {
	    PaperSize[] sizes = { PaperSize.LETTER, 
				  PaperSize.A4, PaperSize.JIS_B5 };
	    paperSizeBox = new JComboBox(sizes);
	}
	paperSizeBox.setSelectedIndex(0);

	JLabel l = new JLabel("Size:");
	l.setForeground(Color.black);

	p.add(Box.createHorizontalStrut(5));
	p.add(l);
	p.add(Box.createHorizontalStrut(30));
	p.add(paperSizeBox);
	p.add(Box.createHorizontalStrut(5));

        Border b = BorderFactory.createEtchedBorder();
        TitledBorder tb = BorderFactory.createTitledBorder(b, "Paper");
        tb.setTitleJustification(TitledBorder.LEFT);
	p.setBorder(tb);

	return p;
    }

    /** creates the orientation panel. */
    protected JPanel createOrientationPanel() {
	JPanel p = new JPanel(new FlowLayout(FlowLayout.LEFT));
	JPanel p1 = new JPanel(new GridLayout(2, 1));

	ButtonGroup btnGroup = new ButtonGroup();
	portrait = new JRadioButton("Portrait", false);
	landscape = new JRadioButton("Landscape", true);
	btnGroup.add(portrait);
	btnGroup.add(landscape);

	p1.add(portrait);
	p1.add(landscape);

	p.add(Box.createHorizontalStrut(5));
	p.add(p1);

        Border b = BorderFactory.createEtchedBorder();
        TitledBorder tb = BorderFactory.createTitledBorder(b, "Orientation");
        tb.setTitleJustification(TitledBorder.LEFT);
	p.setBorder(tb);

	return p;
    }

    /** Creates the margins panel. */
    protected JPanel createMarginsPanel() {
	JPanel p = new JPanel(new FlowLayout(FlowLayout.LEFT));

	JPanel p1 = new JPanel(new GridLayout(2, 2, 0, 4));
	JPanel p2 = new JPanel(new GridLayout(2, 2, 0, 4));

	left = new JTextField("1.0''", 5);
	left.setMaximumSize(new Dimension(8, 10));
	right = new JTextField("1.0''", 5);
	right.setMaximumSize(new Dimension(8, 10));
	top = new JTextField("1.0''", 5);
	top.setMaximumSize(new Dimension(8, 10));
	bottom = new JTextField("1.0''", 5);
	bottom.setMaximumSize(new Dimension(8, 10));

	JLabel l = new JLabel("Left:");
	l.setForeground(Color.black);
	p1.add(l);
	p1.add(left);
	l = new JLabel("Top:");
	l.setForeground(Color.black);
	p1.add(l);
	p1.add(top);
	l = new JLabel("Right:");
	l.setForeground(Color.black);
	p2.add(l);
	p2.add(right);
	l = new JLabel("Bottom:");
	l.setForeground(Color.black);
	p2.add(l);
	p2.add(bottom);

	p.add(Box.createHorizontalStrut(5));
	p.add(p1);
	p.add(Box.createHorizontalStrut(15));
	p.add(p2);
	p.add(Box.createHorizontalGlue());

        Border b = BorderFactory.createEtchedBorder();
        TitledBorder tb = BorderFactory.createTitledBorder(b,
							   "Margins (inches)");
        tb.setTitleJustification(TitledBorder.LEFT);
	p.setBorder(tb);

	return p;
    }

    /** Creates the Print Global View panel. */
    protected JPanel createPrintGlobalViewPanel() {
	JPanel p = new JPanel(new FlowLayout(FlowLayout.LEFT));
	JPanel p1 = new JPanel(new GridLayout(2, 1, 0, 4));
	JPanel p2 = new JPanel(new GridLayout(2, 1, 0, 4));

	rowNumField = new JTextField(5);
	rowNumField.setMaximumSize(new Dimension(8, 15));

	JLabel l1 = new JLabel("Number of Rows per Page:");
	l1.setForeground(Color.black);

	baseNumField = new JTextField(5);
	baseNumField.setMaximumSize(new Dimension(8, 15));

	JLabel l2 = new JLabel("Approximate Number of Bases per Row:");
	l2.setForeground(Color.black);

	p1.add(l2);
	p1.add(l1);
	p2.add(baseNumField);
	p2.add(rowNumField);

	p.add(Box.createHorizontalStrut(5));
	p.add(p1);
	p.add(Box.createHorizontalStrut(15));
	p.add(p2);
	p.add(Box.createHorizontalGlue());

        Border border = BorderFactory.createEtchedBorder();
        TitledBorder tBorder = BorderFactory.createTitledBorder(
					border, "Print Global View");
        tBorder.setTitleJustification(TitledBorder.LEFT);
	p.setBorder(tBorder);

	return p;
    }

    /** Creates the buttons panel. */
    protected JPanel createButtonsPanel() {
	JPanel bottom = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	btnOK = new JButton("OK");
	btnCancel = new JButton("Cancel");

	btnOK.addActionListener(this);
	btnCancel.addActionListener(this);

	JPanel p = new JPanel(new GridLayout(1, 2, 6, 0));
	p.add(btnOK);
	p.add(btnCancel);

	bottom.add(p);

	return bottom;
    }


    /** Invoked when an action is performed. */
    public void actionPerformed(ActionEvent e) {
	Object source = e.getSource();
	if (source == btnOK) {
	    canceled = false;
	    setVisible(false);
	    OKPressed();
	} else if (source == btnCancel) {
	    canceled = true;
	    setVisible(false);
	}
    }

    /** OK button was pressed. */
    protected void OKPressed() {
	paperSize = (PaperSize) paperSizeBox.getSelectedItem();
	orientation = (portrait.isSelected()) ? PageFormat.PORTRAIT
					      : PageFormat.LANDSCAPE;
	double d;
	d = parseMarginField(left);
	if (d <= 0) {
	    left.setText(leftMargin + "''");
	} else {
	    leftMargin = d * 72;
	    left.setText(d + "''");
	}
	d = parseMarginField(right);
	if (d <= 0) {
	    right.setText(rightMargin + "''");
	} else {
	    rightMargin = d * 72;
	    right.setText(d + "''");
	}
	d = parseMarginField(top);
	if (d <= 0) {
	    top.setText(topMargin + "''");
	} else {
	    topMargin = d * 72;
	    top.setText(d + "''");
	}
	d = parseMarginField(bottom);
	if (d <= 0) {
	    bottom.setText(bottomMargin + "''");
	} else {
	    bottomMargin = d * 72;
	    bottom.setText(d + "''");
	}

	int row = parseIntNumField(rowNumField);
	if (row <= 0) {
	    rowNumField.setText(rowNum + "");
	} else if (row > Constants.PA_MAX_ROW_NUM_PER_PAGE) {
	    rowNumField.setText(Constants.PA_MAX_ROW_NUM_PER_PAGE + "");
	    rowNum = Constants.PA_MAX_ROW_NUM_PER_PAGE;
	} else {
	    rowNum = row;
	}

	int base = parseIntNumField(baseNumField);
	if (base <= 0) {
	    baseNumField.setText(baseNum + "");
	} else {
	    baseNum = base;
	}
    }
    
    /** Parses the specified margin text field. */
    protected double parseMarginField(JTextField tf) {
	double d;
	String s = tf.getText().trim();
	if (s.endsWith("'")) {
	    s = s.substring(0, s.indexOf("'")).trim();
	}
	if (s.length() == 0) {
	    return -1.0;
	}
	try {
	    d = Double.parseDouble(s);
	    if (d >= 4.0 || d < 0) {
		return -1.0;
	    }
	    return d;
	} catch (NumberFormatException e) {
	    return -1.0;
	}
    }

    /** Parses the specified integer number text field. */
    protected int parseIntNumField(JTextField tf) {
	int i;
	String s = tf.getText().trim();
	if (s.length() == 0) {
	    return -1;
	}
	try {
	    i = Integer.parseInt(s);
	    if (i <= 0) {
		return -1;
	    }
	    return i;
	} catch (NumberFormatException e) {
	    return -1;
	}
    }

    /** Returns the default page format.  A wrapper funtion for
	java.awt.PrinterJob.defaultPage().  This function returns
        the java.awt.PrinterJob.getPrinterJob().defaultPage(), if
        the resulted PageFormat is valid; otherwise, it returns a
        PageFormat object with US Letter size, 1 inch margins, and
        LANDSCAPE orientation. */
    public static PageFormat defaultPageFormat() {
	PageFormat pf = PrinterJob.getPrinterJob().defaultPage();
	if (pf.getWidth() <= 0 || pf.getHeight() <= 0
		|| pf.getImageableX() <= 0 || pf.getImageableY() <= 0) {
	    Paper paper = new Paper();
	    paper.setSize(PaperSize.LETTER.width, PaperSize.LETTER.height);
	    paper.setImageableArea(72.0, 72.0,
			           PaperSize.LETTER.width - 72.0 - 72.0,
			           PaperSize.LETTER.height - 72.0 - 72.0);
	    pf.setPaper(paper);
	}
	pf.setOrientation(PageFormat.LANDSCAPE);
	return pf;
    }

    /** Computes the current page format from the parameters entered. */
    protected static PageFormat getPageFormat() {
	double height, width, temp;
	height = paperSize.height;
	width = paperSize.width;
	PageFormat pf = new PageFormat();
	if (orientation == PageFormat.LANDSCAPE) {
	    temp = leftMargin;
	    leftMargin = topMargin;
	    topMargin = temp;
	    temp = rightMargin;
	    rightMargin = bottomMargin;
	    bottomMargin = temp;
	}
	double imgX = leftMargin;
	double imgY = topMargin;
	Paper paper = new Paper();
	paper.setSize(width, height);
	double imgWidth = width - leftMargin - rightMargin;
	double imgHeight = height - topMargin - bottomMargin;
	paper.setImageableArea(imgX, imgY, imgWidth, imgHeight);
	pf.setPaper(paper);
	pf.setOrientation(orientation);
	return (PageFormat) pf.clone();
    }

    protected void reset() {
	PaperSize s;
	paperSizeBox.setSelectedIndex(0);
	for (int i = 0; i < paperSizeBox.getItemCount(); i++) {
	    s = (PaperSize) paperSizeBox.getItemAt(i);
	    if (s.name.equals(origPaperSize.name)) {
		paperSizeBox.setSelectedIndex(i);
		break;
	    }
	}

	if (origPageFormat.getOrientation() == PageFormat.PORTRAIT) {
	    portrait.setSelected(true);
	} else {
	    landscape.setSelected(true);
	}

	double imgX = origPageFormat.getImageableX();
	double imgY = origPageFormat.getImageableY();
	double imgW = origPageFormat.getImageableWidth();
	double imgH = origPageFormat.getImageableHeight();
	double height = origPageFormat.getHeight();
	double width = origPageFormat.getWidth();

	left.setText(pixelToInch(imgX) + "''");
	right.setText(pixelToInch(width - imgX - imgW) + "''");
	top.setText(pixelToInch(imgY) + "''");
	bottom.setText(pixelToInch(height - imgY - imgH) + "''");

	rowNumField.setText(origRowNum + "");
	baseNumField.setText(origBaseNum + "");
    }

    /** Converts the specified value (in pixel) to a 3-letter string 
	representing its value in inch. */
    private static String pixelToInch(double p) {
	double d = p / 72.0 + 0.05;
	String s = d + "000";
	return s.substring(0, 3);
    }

    /** A static nested class of PageSetupDialog class.  It defines paper
	sizes, such as US Letter, US Legal, A4, JIS_B5, etc., along with a
	few methods. */
    static class PaperSize {

	/** US Letter paper size. */
	public static final PaperSize LETTER = new PaperSize("US Letter",
							     612.0, 792.0);

	/** US Legal paper size. */
	public static final PaperSize LEGAL = new PaperSize("US Legal",
							     612.0, 1008.0);

	/** A3 paper size. */
	public static final PaperSize A3 = new PaperSize("A3", 842.0, 1190.0);

	/** A4 paper size. */
	public static final PaperSize A4 = new PaperSize("A4", 595.0, 842.0);

	/** A5 paper size. */
	public static final PaperSize A5 = new PaperSize("A5", 419.0, 595.0);

	/** JIS_B4 paper size. */
	public static final PaperSize JIS_B4 = new PaperSize("B4 [JIS]",
							     728.0, 1031.0);

	/** JIS_B5 paper size. */
	public static final PaperSize JIS_B5 = new PaperSize("B5 [JIS]",
							     499.0, 708.0);

	/** The name of the PaperSize. */
	public String name;

	/** The width and height. */
	public double width, height;
	
	/** Constructor. */
	public PaperSize(String n, double w, double h) {
	    name = n;
	    width = w;
	    height = h;
	}

	/** Returns the name of the PaperSize. */
	public String toString() {
	    String s = (name == null) ? "" : name;
	    return s;
	}

	/** Returns the paper size of the specified page format. */
	public static PaperSize getPaperSize(PageFormat pf) {

	    double height = (pf.getOrientation() == PageFormat.PORTRAIT)
				? pf.getHeight() : pf.getWidth();
	    int h = (int) (height + 0.5);
	    if (h == (int) LETTER.height) {
	    	return LETTER; }
	    else if (h == (int) LEGAL.height) {
	    	return LEGAL;
	    } else if (h == (int) A3.height) {
	    	return A3;
	    } else if (h == (int) A4.height) {
	    	return A4;
	    } else if (h == (int) A5.height) {
	    	return A5;
	    } else if (h == (int) JIS_B4.height) {
	    	return JIS_B4;
	    } else if (h == (int) JIS_B5.height) {
	    	return JIS_B5;
	    } else {
	    	return LETTER;
	    }
    	}
    }
}
