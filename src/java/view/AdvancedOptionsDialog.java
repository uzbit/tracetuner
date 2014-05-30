/*
 * 1.9
 *
 * @(#)AdvancedOptionsDialog.java	2001/02/28
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */
package com.paracel.tt.view;

import java.io.*;
import java.lang.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import com.paracel.tt.util.*;

/**
 * This class creates the Advanced Options dialog window for advanced users
 * to specify advanced ttuner run options, including, trimming, .tip files,
 * .tab files, .tal files, etc.
 */
public class AdvancedOptionsDialog extends JDialog
				   implements ActionListener {
    /** The parent frame (The Options frame). */
    private JFrame frame;

    /** The Advanced Options frame. */
    private static AdvancedOptionsDialog advOptDialog;

    /** This holds the current selections in the Options window. */
    private Options options;

    /** The sub panels included in the window. */
    private JPanel outputPanel, programPanel, buttonPanel;

    /** The check boxes. */
    private JCheckBox tipBox, talBox, tabBox, edBaseBox, consBox, ctBox,
		      trimWindowBox, trimThresholdBox, shiftBox;

    /** The buttons. */
    private JButton okBtn, cancelBtn, resetBtn;

    /** The consensus chooser and the context table chooser. */
    private Filespec consChooser, ctChooser;

    /** The panel has the trim window size text box following its check box. */
    private JPanel trimWindowPanel;

    /** The panel has the trim threshold text box following its check box. */
    private JPanel trimThresholdPanel;

    /** The trim threshold and trim window size text fields. */
    private JTextField trimWindowField, trimThresholdField;

    /** The error message. */
    private String errMesg;

    /** The edited base and shift flag values that was saved
	(temporarily by clicking on the OK button in Advanced Options window
	or permanently into a file by clicking on the save button in Basic
	Options window. */
    private boolean lastSavedEditedBase;

    /** The constructor. */
    protected AdvancedOptionsDialog(JFrame frame) {
	super(frame, "Advanced Options", true);
	this.frame = frame;
	createUI();
    }

    /** Gets the instance of the AdvancedOptionsDialog. */
    public static AdvancedOptionsDialog getInstance(JFrame frame) {
	if (advOptDialog == null) {
	    advOptDialog = new AdvancedOptionsDialog(frame);
	}
	return advOptDialog;
    }

    /** Creates the UI. */
    protected void createUI() {
	createOutputPanel();
	createProgramPanel();
	createButtonPanel();

	JPanel contentPanel = new JPanel();
	if (Constants.runAsDev) {
	    contentPanel.setPreferredSize(new Dimension(555, 420));
	} else {
	    contentPanel.setPreferredSize(new Dimension(555, 390));
	}
	contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

	contentPanel.add(outputPanel);
	contentPanel.add(programPanel);
        contentPanel.add(Box.createVerticalStrut(10));
	contentPanel.add(buttonPanel);

	getContentPane().add(new JScrollPane(contentPanel));
    }

    /** Creates the Output Options panel. */
    protected void createOutputPanel() {
	outputPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
	
	JPanel p1 = new JPanel();
	p1.setLayout(new GridLayout(1, 1, 0, 0));
	JPanel p2 = new JPanel();
	p2.setLayout(new GridLayout(1, 1, 0, 0));
	JPanel p3 = new JPanel();
	p3.setLayout(new GridLayout(1, 1, 0, 0));

	tabBox = new JCheckBox("Tab files");
	talBox = new JCheckBox("Tal files");
	tipBox = new JCheckBox("Tip files");

	p1.add(tabBox);
	p2.add(talBox);
	p3.add(tipBox);

    outputPanel.add(Box.createHorizontalStrut(40));
	outputPanel.add(p1);
    outputPanel.add(Box.createHorizontalStrut(100));
    outputPanel.add(p2);
    outputPanel.add(Box.createHorizontalStrut(160));
    outputPanel.add(p3);

        Border outputBorder = BorderFactory.createEtchedBorder();
        TitledBorder outputTBorder = BorderFactory.createTitledBorder(
					outputBorder, "More output options");
        outputTBorder.setTitleJustification(TitledBorder.LEFT);
        outputTBorder.setTitlePosition(TitledBorder.DEFAULT_POSITION);
	outputPanel.setBorder(outputTBorder);
    }

    /** Creates the Program Options Panel. */
    protected void createProgramPanel() {
	programPanel = new JPanel();
	programPanel.setLayout(new BoxLayout(programPanel, BoxLayout.Y_AXIS));

	consBox = new JCheckBox("Specify consensus/reference sequence file");
	ctBox = new JCheckBox("Specify context table file");
	edBaseBox = new JCheckBox("Edited bases (read and recall " +
                                             "from edited bases)");
	shiftBox = new JCheckBox("Perform mobility shift correction");
	trimWindowBox = new JCheckBox("Specify trim window size "
				    + "( >0, default is "
				    + Constants.DEFAULT_TRIM_WINDOW_SIZE +")");
	trimThresholdBox = new JCheckBox("Specify trim threshold "
				    + "(1-100, default is "
				    + Constants.DEFAULT_TRIM_THRESHOLD + ")");
	consBox.addActionListener(this);
	ctBox.addActionListener(this);
	trimWindowBox.addActionListener(this);
	trimThresholdBox.addActionListener(this);

	JPanel p1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p1.add(Box.createHorizontalStrut(40));
	p1.add(edBaseBox);
	p1.setBorder(new EmptyBorder(-4, 0, -4, 0));

	JPanel p2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p2.add(Box.createHorizontalStrut(40));
	p2.add(shiftBox);
	p2.setBorder(new EmptyBorder(-4, 0, -4, 0));

	// Lower half of program options panel
	trimWindowPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        trimWindowPanel.add(Box.createHorizontalStrut(40));
	trimWindowPanel.add(trimWindowBox);
	trimWindowPanel.setBorder(new EmptyBorder(-4, 0, -4, 0));

	trimThresholdPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        trimThresholdPanel.add(Box.createHorizontalStrut(40));
	trimThresholdPanel.add(trimThresholdBox);
	trimThresholdPanel.setBorder(new EmptyBorder(-4, 0, -4, 0));

	JPanel p4 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p4.add(Box.createHorizontalStrut(40));
	p4.add(ctBox);
	p4.setBorder(new EmptyBorder(-4, 0, -4, 0));

	JPanel p5 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p5.add(Box.createHorizontalStrut(40));
	p5.add(consBox);
	p5.setBorder(new EmptyBorder(-4, 0, -4, 0));

	programPanel.add(p1);
	programPanel.add(p2);
	programPanel.add(trimWindowPanel);
	programPanel.add(trimThresholdPanel);
	if (Constants.runAsDev) {
	    programPanel.add(p4);
	}
	programPanel.add(p5);

        Border progBorder = BorderFactory.createEtchedBorder();
        TitledBorder progTBorder = BorderFactory.createTitledBorder(
					progBorder, "More program options");
        progTBorder.setTitleJustification(TitledBorder.LEFT);
        progTBorder.setTitlePosition(TitledBorder.DEFAULT_POSITION);
	programPanel.setBorder(progTBorder);

	trimWindowField = new JTextField(5);
	trimWindowField.setHorizontalAlignment(JTextField.RIGHT);

	trimThresholdField = new JTextField(5);
	trimThresholdField.setHorizontalAlignment(JTextField.RIGHT);

	consChooser = new Filespec("Consensus/Reference Sequence File",
				   JFileChooser.FILES_ONLY,
				   "Select consensus/reference sequence file");
	consChooser.setBorder(new EmptyBorder(-4, 0, -4, 0));

	ctChooser = new Filespec("Context Table File",
				 JFileChooser.FILES_ONLY,
				 "Select context table file");
	ctChooser.setBorder(new EmptyBorder(-4, 0, -4, 0));
    }

    /** Creates the buttons panel. */
    protected void createButtonPanel() {
	buttonPanel = new JPanel();

	okBtn = new JButton("Ok");
	cancelBtn = new JButton("Cancel");
	resetBtn = new JButton("Reset Defaults");

	okBtn.addActionListener(this);
	cancelBtn.addActionListener(this);
	resetBtn.addActionListener(this);

	buttonPanel.add(okBtn);
	buttonPanel.add(cancelBtn);
	buttonPanel.add(resetBtn);
    }

    /** Loads the advanced options and becomes visible. */
    public void becomeVisible(Options opt) {
	Point p = frame.getLocation();
	setLocation((int) p.getX(), (int) p.getY() + 50);
	options = opt;
	loadAdvOptions();
	pack();
	setVisible(true);
    }

    /** Loads the advanced ptions. */
    protected void loadAdvOptions() {
	reset();

	boolean b;
	b = (options.outputTip()) ? true : false;
	tipBox.setSelected(b);

	b = (options.outputTal()) ? true : false;
	talBox.setSelected(b);

	b = (options.outputTab()) ? true : false;
	tabBox.setSelected(b);

	b = (options.editedBase()) ? true : false;
	edBaseBox.setSelected(b);
	Constants.fromEdited = b;
	lastSavedEditedBase = b;

	b = (options.shift()) ? true : false;
	shiftBox.setSelected(b);

	String s = options.getTrimWindowSize();
	if (s != null) {
	    trimWindowBox.setSelected(true);
	    trimWindowField.setText(s);
	    trimWindowPanel.add(trimWindowField);
	}
	
	s = options.getTrimThreshold();
	if (s != null) {
	    trimThresholdBox.setSelected(true);
	    trimThresholdField.setText(s);
	    trimThresholdPanel.add(trimThresholdField);
	}
	
	if (Constants.runAsDev) {
	    s = options.getContextTable();
	    if (s != null) {
	    	ctBox.setSelected(true);
	    	ctChooser.text1.setText(s);
	    	programPanel.add(ctChooser, 6);
	    }
	}

	s = options.getConsensus();
	if (s != null) {
	    consBox.setSelected(true);
	    consChooser.text1.setText(s);
	    if (Constants.runAsDev) {
	    	int index = (ctBox.isSelected()) ? 8 : 7;
	    	programPanel.add(consChooser, index);
	    } else {
	    	programPanel.add(consChooser);
	    }
	}

	programPanel.revalidate();
	programPanel.repaint();
    }

    /** Validates the advanced options currently chosen on the panels. */
    public boolean validateAdvOptions() {
	String s, e;
	File f;

	if (trimWindowBox.isSelected()) {
	    s = trimWindowField.getText().trim();
	    e = "Please specify a valid trim window size value\n"
		 + "(a positive integer number required)\n"
		 + "in the Trim Window Size text box or turn off\n"
		 + "the Specify Trim Window Size option.\n";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    } else {
		try {
		    int i = Integer.parseInt(s);
		    if (i <= 0) {
			errMesg = e;
			return false;
		    }
		} catch (Exception ee) {
		    errMesg = e;
		    return false;
		}
	    }
	}

	if (trimThresholdBox.isSelected()) {
	    s = trimThresholdField.getText().trim();
	    e = "Please specify a valid trim threshold value\n"
		 + "(an integer number between 1 and 100 required)\n"
		 + "in the Trim Threshold text box or turn off\n"
		 + "the Specify Trim Threshold option.\n";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    } else {
		try {
		    int i = Integer.parseInt(s);
		    if (i < 1 || i > 100) {
			errMesg = e;
			return false;
		    }
		} catch (Exception ee) {
		    errMesg = e;
		    return false;
		}
	    }
	}

	if (consBox.isSelected()) {
	    s = consChooser.text1.getText().trim();
	    e = "Please specify a valid consensus/reference sequence file\n"
		+ "in the Consensus/Reference Sequence File text box or\n"
		+ "turn off the Specify Consensus/Reference Sequence\n"
		+ "File option.";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    }

	    f = new File(s);
	    if (!f.exists()) {
		errMesg = "The specified consensus/reference sequence file"
		    	  + " does not exist.\n" + e;
		return false;
	    }
	    if (!f.canRead()) {
		errMesg = "Read permission is denied to the specified\n"
		          + "consensus/reference sequence file.\n" + e;
		return false;
	    }
	}

	if (Constants.runAsDev) {
	    if (ctBox.isSelected()) {
	    	s = ctChooser.text1.getText().trim();
	    	e = "Please specify a valid context table file\n"
			+ "in the Context Table File text box or\n"
			+ "turn off the Specify Context Table File option.";
	    	if (s.length() == 0) {
		    errMesg = e;
		    return false;
	    	}

	    	f = new File(s);
	    	if (!f.exists()) {
		    errMesg = "The specified context table file"
		    	  	+ " does not exist.\n" + e;
		    return false;
	    	}
	    	if (!f.canRead()) {
		    errMesg = "Read permission is denied to the specified\n"
		          	+ "context table file.\n" + e;
		    return false;
	    	}
	    }
	}

	return true;
    }

    /** Resets the panels to default. */
    protected void reset() {
	//Clear the selections in the window.
	tipBox.setSelected(false);
	talBox.setSelected(false);
	tabBox.setSelected(false);
	edBaseBox.setSelected(false);
	shiftBox.setSelected(false);
	Constants.fromEdited = false;
	if (trimWindowBox.isSelected()) {
	    trimWindowPanel.remove(trimWindowField);
	    trimWindowBox.setSelected(false);
	}
	if (trimThresholdBox.isSelected()) {
	    trimThresholdPanel.remove(trimThresholdField);
	    trimThresholdBox.setSelected(false);
	}
	if (consBox.isSelected()) {
	    programPanel.remove(consChooser);
	    consBox.setSelected(false);
	}
	if (Constants.runAsDev) {
	    if (ctBox.isSelected()) {
	    	programPanel.remove(ctChooser);
	    	ctBox.setSelected(false);
	    }
	}

	if (isVisible()) {
	    programPanel.revalidate();
	    repaint();
	}
    }

    /** Invoked when an action is performed. */
    public void actionPerformed(ActionEvent e) {
	Object source = e.getSource();
	if (source == trimWindowBox ) {
	    if (trimWindowBox.isSelected()) {
		trimWindowPanel.add(trimWindowField);
	    } else {
		trimWindowPanel.remove(trimWindowField);
	    }
	    programPanel.revalidate();
	    programPanel.repaint();
	} else if (source == trimThresholdBox ) {
	    if (trimThresholdBox.isSelected()) {
		trimThresholdPanel.add(trimThresholdField);
	    } else {
		trimThresholdPanel.remove(trimThresholdField);
	    }
	    programPanel.revalidate();
	    programPanel.repaint();
	} else if (source == consBox ) {
	    if (consBox.isSelected()) {
		if (Constants.runAsDev) {
		    int index = (ctBox.isSelected()) ? 8 : 7;
		    programPanel.add(consChooser, index);
		} else {
  		    programPanel.add(consChooser);
		}
	    } else {
		programPanel.remove(consChooser);
	    }
	    programPanel.revalidate();
	    programPanel.repaint();
  	} else if (source == ctBox ) {
	    if (Constants.runAsDev) {
	    	if (ctBox.isSelected()) {
		    programPanel.add(ctChooser, 6);
	    	} else {
		    programPanel.remove(ctChooser);
	    	}
	    	programPanel.revalidate();
	    	programPanel.repaint();
	    }
	} else if (source == okBtn) {
	    if (validateAdvOptions()) {
		/* After the validation, saves the current selections
		   on the Advanced Options window into the Options
		   object passed from the Options window.  Then makes
		   this Advanced Options window invisible. */
		boolean b;
		b = (tipBox.isSelected()) ? true : false;
		options.setOutputTip(b);

		b = (talBox.isSelected()) ? true : false;
		options.setOutputTal(b);

		b = (tabBox.isSelected()) ? true : false;
		options.setOutputTab(b);

		b = (edBaseBox.isSelected()) ? true : false; 
		options.setEdBase(b);
		Constants.fromEdited = b;
		lastSavedEditedBase = b;

		b = (shiftBox.isSelected()) ? true : false;
		options.setShift(b);

		String s = (trimWindowBox.isSelected()) 
	    			? trimWindowField.getText().trim() : null;
		options.setTrimWindowSize(s);

		s = (trimThresholdBox.isSelected()) 
	    			? trimThresholdField.getText().trim() : null;
		options.setTrimThreshold(s);

		s = (consBox.isSelected()) 
		    		? consChooser.text1.getText().trim()
				: null;
		options.setConsensus(s);

		if (Constants.runAsDev) {
		    s = (ctBox.isSelected()) 
		    		? ctChooser.text1.getText().trim()
				: null;
		    options.setContextTable(s);
		}
		setVisible(false);
	    } else {
	        JOptionPane.showMessageDialog(this, errMesg);
	    }
	} else if (source == cancelBtn) {
	    setVisible(false);
	    Constants.fromEdited = lastSavedEditedBase;
	} else if (source == resetBtn) {
	    reset();
	}
    }
}
