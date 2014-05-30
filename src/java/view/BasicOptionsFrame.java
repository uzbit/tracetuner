/*
 * 1.15
 *
 * @(#)BasicOptionsFrame.java	2001/02/28
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
 * This class creates the Options dialog window for the regular user
 * to specify basic ttuner run options.  Also it provides the advanced
 * user the access to more advanced options through the "Advanced
 * Options" button.
 */
public class BasicOptionsFrame extends JFrame
			       implements ActionListener {
    /** Number of default (or built-in) lookup tables. */
    private final static int NUM_DEFAULT_LOOKUP_TABLES = 5;

    /** The parent frame (The Launcher frame). */
    private JFrame frame;

    /** The Basic Options frame. */
    private static BasicOptionsFrame basicOptFrame;

    /** The Advanced Options dialog. */
    private static AdvancedOptionsDialog advOptDialog;

    /** The options. */
    private static Options options;

    /** The first panel on the frame: Output Options panel. */
    private JPanel outputPanel;

    /** The second panel on the frame: Processing Options panel. */
    private JPanel processingPanel;

    /** The third panel on the frame: Program Options panel. */
    private JPanel programPanel;

    /** The fourth panel on the frame: buttons panel. */
    private JPanel buttonPanel;

    /** The panel has the min ratio text box following its check box. */
    private JPanel minRatioPanel;

    /** The output format check boxes. */
    private JCheckBox phdBox, scfBox, qualBox, seqBox, qaBox, saBox;

    /** The run options check boxes. */
    private JCheckBox noCallBox, recallNBox, hetBox, mixBox, minRatioBox, outDirBox;

    /** The lookup table list combo box. */
    private JComboBox tableList;

    /** The min ratio text fields. */
    private JTextField minRatioField;

    /** The qa file name, sa file name text fields. */
    private JTextField qaField, saField;

    /** The panel containing the "File Name" label followed by the
	qa file name textfield; the panel containing the "File Name"
	label followed by the sa file name textfield. */
    private JPanel qaPanel, saPanel;

    /** The left half and right half of the output options panel. */
    private JPanel leftOutputPanel, rightOutputPanel;

    /** The buttons. */
    private JButton saveBtn, cancelBtn, helpBtn, resetBtn, advBtn;

    /** The lookup table chooser and the output directory chooser. */
    private Filespec tableChooser, outDirChooser;

    /** The error message. */
    private String errMesg;

    private Options lastSavedOptions;

    /** The constructor. */
    protected BasicOptionsFrame() {
	super("Basic Options");
	createUI();
	advOptDialog = AdvancedOptionsDialog.getInstance(this);
	options = new Options();
	loadOptions();
    }

    /** Gets the instance of the BasicOptionsFrame. */
    public static BasicOptionsFrame getInstance() {
	if (basicOptFrame == null) {
	    basicOptFrame = new BasicOptionsFrame();
	}
	return basicOptFrame;
    }

    /** Sets the specified frame as the parent frame. */
    public void setParentFrame(JFrame f) {
	frame = f;
    }

    /** Creates the UI. */
    protected void createUI() {
	createOutputPanel();
	createProcessingPanel();
	createProgramPanel();
	createButtonPanel();

	JPanel contentPanel = new JPanel();
	contentPanel.setPreferredSize(new Dimension(540, 480));
	contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

	contentPanel.add(outputPanel);
	contentPanel.add(processingPanel);
	contentPanel.add(programPanel);
        contentPanel.add(Box.createVerticalStrut(10));
	contentPanel.add(buttonPanel);

	getContentPane().add(new JScrollPane(contentPanel));
    }

    /** Creates the Output Options panel. */
    protected void createOutputPanel() {
	outputPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
	
	leftOutputPanel = new JPanel(new GridLayout(4, 1, 0, 2));
	rightOutputPanel = new JPanel(new GridLayout(4, 1, 0, 2));

	phdBox = new JCheckBox("Phd files (required by TTView)");
	scfBox = new JCheckBox("SCF files");
	qualBox = new JCheckBox("Qual files (one file per read)");
	qaBox = new JCheckBox("Qual file (single file for all reads)");
	seqBox = new JCheckBox("Seq files (one file per read)");
	saBox = new JCheckBox("Seq file (single file for all reads)");
	qaBox.addActionListener(this);
	saBox.addActionListener(this);
	phdBox.setSelected(true);

	qaField = new JTextField(10);
	saField = new JTextField(10);

	qaPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        qaPanel.add(Box.createHorizontalStrut(8));
	qaPanel.add(new JLabel("File Name"));
	qaPanel.add(qaField);
	qaPanel.setBorder(new EmptyBorder(-4, 0, -4, 0));

	saPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        saPanel.add(Box.createHorizontalStrut(8));
	saPanel.add(new JLabel("File Name"));
	saPanel.add(saField);
	saPanel.setBorder(new EmptyBorder(-4, 0, -4, 0));

	leftOutputPanel.add(phdBox);
	leftOutputPanel.add(qualBox);
	leftOutputPanel.add(qaBox);
	rightOutputPanel.add(scfBox);
	rightOutputPanel.add(seqBox);
	rightOutputPanel.add(saBox);

        outputPanel.add(Box.createHorizontalStrut(30));
	outputPanel.add(leftOutputPanel);
        outputPanel.add(Box.createHorizontalStrut(10));
	outputPanel.add(rightOutputPanel);
        outputPanel.add(Box.createHorizontalGlue());

        Border outputBorder = BorderFactory.createEtchedBorder();
        TitledBorder outputTBorder = BorderFactory.createTitledBorder(
					outputBorder, "Output options");
        outputTBorder.setTitleJustification(TitledBorder.LEFT);
        outputTBorder.setTitlePosition(TitledBorder.DEFAULT_POSITION);
	outputPanel.setBorder(outputTBorder);
    }

    /** Creates the Processing Options Panel. */
    protected void createProcessingPanel() {
	processingPanel = new JPanel();
	processingPanel.setLayout(new BoxLayout(processingPanel, 
						BoxLayout.Y_AXIS));

	JPanel calibPanel = new JPanel();
	JLabel calibLabel = new JLabel("Calibration tables:");
	calibPanel.add(calibLabel);

	tableList = new JComboBox();
	tableList.addItem(Constants.DEFAULT_CALIBRATION);
    tableList.addItem(Constants.ABI3730_POP7_LOOKUP_TABLE);
	tableList.addItem(Constants.ABI3700_POP5_LOOKUP_TABLE);
	tableList.addItem(Constants.ABI3700_POP6_LOOKUP_TABLE);
	tableList.addItem(Constants.ABI3100_POP6_LOOKUP_TABLE);
    tableList.addItem(Constants.MegaBACE_LOOKUP_TABLE);
        tableList.setSelectedIndex(0);
	File tableDir = new File(Constants.USER_DIR);
	String[] tList = tableDir.list();
	for (int i = 0; i < tList.length; i++) {
	    if (tList[i].endsWith(Constants.LOOKUP_TABLE_SUFFIX)) {
		tableList.addItem(tList[i].substring(0,
		    tList[i].length()-Constants.LOOKUP_TABLE_SUFFIX.length()));
	    }
	}
	tableList.addItem("Others..(specify)");
	tableList.addActionListener(this);
	calibPanel.add(tableList);
	calibPanel.setBorder(new EmptyBorder(-2, 0, -2, 0));
	
	processingPanel.add(calibPanel);

        Border procBorder = BorderFactory.createEtchedBorder();
        TitledBorder procTBorder = BorderFactory.createTitledBorder(
					procBorder, "Processing options");
        procTBorder.setTitleJustification(TitledBorder.LEFT);
        procTBorder.setTitlePosition(TitledBorder.DEFAULT_POSITION);
	processingPanel.setBorder(procTBorder);

	tableChooser  = new Filespec("Specify Lookup Table",
			             JFileChooser.FILES_ONLY,
                                     "Select lookup table file");
	tableChooser.setBorder(new EmptyBorder(-2, 0, -2, 0));
    }

    /** Creates the Program Options Panel. */
    protected void createProgramPanel() {
	programPanel = new JPanel();
	programPanel.setLayout(new BoxLayout(programPanel, BoxLayout.Y_AXIS));

	minRatioBox = new JCheckBox("Specify Minimum Peak Height Ratio "
				    + "(0-1.0, default is "
				    + Constants.DEFAULT_MIN_RATIO + ")");
	hetBox = new JCheckBox("Call heterozygotes");
        mixBox = new JCheckBox("Call mixed bases");
	noCallBox = new JCheckBox("No call (use original base "
                                  + "calls without adjustment)");
	outDirBox = new JCheckBox("Specify output directory "
                                  + "(default goes to input directory)");
	recallNBox = new JCheckBox("Recall N (no change to other "
                                   + "base calls)");
	noCallBox.addActionListener(this);
	recallNBox.addActionListener(this);
	hetBox.addActionListener(this);
        mixBox.addActionListener(this);
	minRatioBox.addActionListener(this);
	outDirBox.addActionListener(this);

	JPanel p0 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p0.add(Box.createHorizontalStrut(30));
	p0.add(noCallBox);
	p0.setBorder(new EmptyBorder(-4, 0, -4, 0));

	JPanel p1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p1.add(Box.createHorizontalStrut(30));
	p1.add(recallNBox);
	p1.setBorder(new EmptyBorder(-4, 0, -4, 0));

	JPanel p2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p2.add(Box.createHorizontalStrut(30));
	p2.add(hetBox);
	p2.setBorder(new EmptyBorder(-4, 0, -4, 0));

        //JPanel p7 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        //p7.add(Box.createHorizontalStrut(30));
        //p7.add(mixBox);
        //p7.setBorder(new EmptyBorder(-4, 0, -4, 0));

	// Lower half of program options panel
	minRatioPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        minRatioPanel.add(Box.createHorizontalStrut(60));
	minRatioPanel.add(minRatioBox);
	minRatioPanel.setBorder(new EmptyBorder(-4, 0, -4, 0));

	JPanel p6 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        p6.add(Box.createHorizontalStrut(30));
	p6.add(outDirBox);
	p6.setBorder(new EmptyBorder(-4, 0, -4, 0));

	programPanel.add(p0);
	programPanel.add(p1);
	programPanel.add(p2);
        //programPanel.add(p7);
	programPanel.add(minRatioPanel);
	programPanel.add(p6);

        Border progBorder = BorderFactory.createEtchedBorder();
        TitledBorder progTBorder = BorderFactory.createTitledBorder(
					progBorder, "Program options");
        progTBorder.setTitleJustification(TitledBorder.LEFT);
        progTBorder.setTitlePosition(TitledBorder.DEFAULT_POSITION);
	programPanel.setBorder(progTBorder);

	minRatioField = new JTextField(5);
	minRatioField.setHorizontalAlignment(JTextField.RIGHT);

	outDirChooser = new Filespec("Output Directory",
				     JFileChooser.DIRECTORIES_ONLY,
				     "Select output directory");
	outDirChooser.setBorder(new EmptyBorder(-4, 0, -4, 0));
    }

    /** Creates the buttons panel. */
    protected void createButtonPanel() {
	buttonPanel = new JPanel();

	saveBtn = new JButton("Save");
	cancelBtn = new JButton("Cancel");
	helpBtn = new JButton(new HelpAction("Help"));
	resetBtn = new JButton("Reset Defaults");
	advBtn = new JButton("Advanced Options");

	saveBtn.addActionListener(this);
	cancelBtn.addActionListener(this);
	resetBtn.addActionListener(this);
	advBtn.addActionListener(this);

	buttonPanel.add(saveBtn);
	buttonPanel.add(cancelBtn);
	buttonPanel.add(helpBtn);
	buttonPanel.add(resetBtn);
	buttonPanel.add(advBtn);
    }

    /** Loads the options and becomes visible. */
    public void becomeVisible() {
	Point p = frame.getLocation();
	Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
	int x = Math.min(p.x + 480, (int) d.getWidth() - getWidth());
	int y = Math.min(p.y, (int) d.getHeight() - getHeight());
	setLocation(x, y);
	refreshLookupTableList();
	loadOptions();
	pack();
	setVisible(true);
    }

    /** Refreshes the combo box which contains the lookup table list. */
    protected void refreshLookupTableList() {
	// List content might be changed because users move .tbl files
	// in or out of the current directory.
	for (int i = tableList.getItemCount() - 1;
	     		i > NUM_DEFAULT_LOOKUP_TABLES; i--) {
	    tableList.removeItemAt(i);
	}
	File tableDir = new File(Constants.USER_DIR);
	String[] tList = tableDir.list();
	for (int i = 0; i < tList.length; i++) {
	    if (tList[i].endsWith(Constants.LOOKUP_TABLE_SUFFIX)) {
		tableList.addItem(tList[i].substring(0,
		    tList[i].length()-Constants.LOOKUP_TABLE_SUFFIX.length()));
	    }
	}
	tableList.addItem("Others..(specify)");
    }

    /** Loads the options. */
    protected void loadOptions() {
	options.load();
	if (options.getStatus() != Options.OK) {
	    options.reset();
	}
	lastSavedOptions = new Options(options);
	Constants.fromEdited = options.editedBase();
	setOptions();
    }

    /** Sets the options on the panels. */
    protected void setOptions() {
	reset();

	boolean b;
	String s;
	b = (options.outputPhd()) ? true : false;
	phdBox.setSelected(b);

	b = (options.outputQual()) ? true : false;
	qualBox.setSelected(b);

	b = (options.outputQa()) ? true : false;
	qaBox.setSelected(b);

	s = options.getQAFileName();
	if (s != null) {
	    qaField.setText(s);
	    leftOutputPanel.add(qaPanel);
	}

	b = (options.outputScf()) ? true : false;
	scfBox.setSelected(b);

	b = (options.outputSeq()) ? true : false;
	seqBox.setSelected(b);

	b = (options.outputSa()) ? true : false;
	saBox.setSelected(b);

	s = options.getSAFileName();
	if (s != null) {
	    saField.setText(s);
	    rightOutputPanel.add(saPanel);
	}

	s = options.getOtherLookupTable();
	if (s != null) {
	    tableList.setSelectedIndex(tableList.getItemCount() - 1);
	    tableChooser.text1.setText(s);
	    processingPanel.add(tableChooser, 1);
	    processingPanel.revalidate();
	    processingPanel.repaint();
	} else {
	    s = options.getLookupTable();
	    if (s.equals(Constants.DEFAULT_CALIBRATION)) {
		tableList.setSelectedIndex(0);
        } else if (s.equals(Constants.ABI3730_POP7_LOOKUP_TABLE)) {
        tableList.setSelectedIndex(1);
	    } else if (s.equals(Constants.ABI3700_POP5_LOOKUP_TABLE)) {
		tableList.setSelectedIndex(2);
	    } else if (s.equals(Constants.ABI3700_POP6_LOOKUP_TABLE)) {
		tableList.setSelectedIndex(3);
	    } else if (s.equals(Constants.ABI3100_POP6_LOOKUP_TABLE)) {
		tableList.setSelectedIndex(4);
	    } else if (s.equals(Constants.MegaBACE_LOOKUP_TABLE)) {
		tableList.setSelectedIndex(5);
	    } else {
	    	s = s.substring(0,
			s.length()-Constants.LOOKUP_TABLE_SUFFIX.length());
	    	for (int i = NUM_DEFAULT_LOOKUP_TABLES + 1;
		     		i < (tableList.getItemCount()-1); i++) {
		    if (s.endsWith((String)tableList.getItemAt(i))) {
		    	tableList.setSelectedIndex(i);
		    	break;
		    }
	    	}
	    }
	}

	b = (options.noCall()) ? true : false;
	noCallBox.setSelected(b);

	b = (options.recallN()) ? true : false;
	recallNBox.setSelected(b);

	b = (options.het()) ? true : false;
	hetBox.setSelected(b);

        b = (options.mix()) ? true : false;
        mixBox.setSelected(b);

	minRatioBox.setEnabled(b);
	minRatioField.setEnabled(b);

	s = options.getMinRatio();
	if (s != null) {
	    minRatioBox.setSelected(true);
	    minRatioField.setText(s);
	    minRatioPanel.add(minRatioField);
	}

	s = options.getOutputDir();
	if (s != null) {
	    outDirBox.setSelected(true);
	    outDirChooser.text1.setText(s);
	    programPanel.add(outDirChooser, 5);
	}
	programPanel.revalidate();
	programPanel.repaint();
    }

    /** Validates the options currently chosen on the Options window
	and the Advanced Options window. */
    public boolean validateOptions() {
	String s, e;
	File f;

	if ( !(phdBox.isSelected() || scfBox.isSelected()
	          || qualBox.isSelected() || seqBox.isSelected()
	          || qaBox.isSelected() || saBox.isSelected()
	          || options.outputTip() || options.outputTal()
	          || options.outputTab()) ) {
	    errMesg = "Please specified at least one output option.";
	    return false;
	}

	if (qaBox.isSelected()) {
	    s = qaField.getText().trim();
	    e = "Please specify a file name in the text box under\n"
		      + "the \"Qual file (single file for all reads)\""
		      + "checkbox.\n";
	    if (s.length() <= 0) {
	    	errMesg = e;
		return false;
	    }
	    if (s.indexOf(Constants.FILE_SEPARATOR) >= 0) {
	    	errMesg = "A file name is required, NOT a pathname.\n" + e;
		return false;
	    }
	}

	if (saBox.isSelected()) {
	    s = saField.getText().trim();
	    e = "Please specify a valid file name in the text box\n"
		 + "under the \"Seq file (single file for all reads)\" "
		 + "check box.\n";
	    if (s.length() <= 0) {
	    	errMesg = e;
		return false;
	    }
	    if (s.indexOf(Constants.FILE_SEPARATOR) >= 0) {
	    	errMesg = "A file name is required, NOT a pathname.\n" + e;
		return false;
	    }
	}

	if (options.outputTab() && (!phdBox.isSelected())) {
	    errMesg = "Phd Files option is required,\n"
		      +"when Tab Files option is selected.";
	    return false;
	}

	if (options.outputTal() && (!phdBox.isSelected())) {
	    errMesg = "Phd Files option is required,\n"
		      +"when Tal Files option is selected.";
	    return false;
	}

	if (options.outputTip() && (!phdBox.isSelected())) {
	    errMesg = "Phd Files option is required,\n"
		      +"when Tip Files option is selected.";
	    return false;
	}

	if (minRatioBox.isSelected() && minRatioBox.isEnabled()) {
	    s = minRatioField.getText().trim();
	    e = "Please specify a valid minimum peak height ratio value\n"
		+ "(a real number between 0.0 and 1.0 required) in the\n"
		+ "Minimum Peak Height Ratio text box or turn off\n"
		+ "the Specify Minimum Peak Height Ratio option.\n";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    } else {
		try {
		    double d = Double.parseDouble(s);
		    if (d < 0 || d > 1.0) {
			errMesg = e;
			return false;
		    }
		} catch (Exception ee) {
		    errMesg = e;
		    return false;
		}
	    }
	}

	if (tableList.getSelectedIndex() == (tableList.getItemCount()-1)) {
	    s = tableChooser.text1.getText().trim();
	    e = "Please specify a valid lookup table in\n"
		          + "the Lookup Table text box or select\n"
		          + "an existing table in the list.";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    }

	    f = null;
	    f = new File(s);
	    if (!f.exists()) {
		errMesg = "The specified lookup table no longer exists.\n" + e;
		return false;
	    }
	    if (!f.canRead()) {
		errMesg = "Read permission is denied to the specified"
		          + " lookup table.\n" + e;
		return false;
	    }
	} else if (tableList.getSelectedIndex() 
		   	>= NUM_DEFAULT_LOOKUP_TABLES + 1) {
	    s = Constants.USER_DIR + Constants.FILE_SEPARATOR
		+ (String)tableList.getSelectedItem() 
		+ Constants.LOOKUP_TABLE_SUFFIX;
	    e = "Please specify a valid lookup table in\n"
		          + "the Lookup Table text box or select\n"
		          + "the default built-in lookup table in the list.";
	    f = null;
	    f = new File(s);
	    if (!f.exists()) {
		errMesg = "The selected lookup table no longer exists.\n" + e;
		return false;
	    }
	    if (!f.canRead()) {
		errMesg = "Read permission is denied to the selected"
		          + " lookup table.\n" + e;
		return false;
	    }
	}

	if (outDirBox.isSelected()) {
	    s = outDirChooser.text1.getText().trim();
	    e = "Please specify a valid output directory in\n"
		          + "the Output Directory text box or turn off\n"
		          + "the Specified Output Directory Option.";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    }

	    f = null;
	    f = new File(s);
	    if (!f.exists()) {
		errMesg = "The specified output directory does not exist.\n"+e;
		return false;
	    }
	    if (!f.canWrite()) {
		errMesg = "Write permission is denied to the specified"
		          + " output directory.\n" + e;
		return false;
	    }
	}

	// Validation of the options in the Advanced Options window.
	s = options.getTrimWindowSize();
	if (s != null) {
	    e = "Please specify a valid trim window size value\n"
		 + "(a positive integer number required)\n"
		 + "in the Advanced Options window by clicking\n"
		 + "on the Advanced Options button.\n";
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

	s = options.getTrimThreshold();
	if (s != null) {
	    e = "Please specify a valid trim threshold value\n"
		 + "(an integer number between 1 and 100 required)\n"
		 + "in the Advanced Options window by clicking\n"
		 + "on the Advanced Options button.\n";
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

	s = options.getContextTable();
	if (s != null) {
	    e = "Please specify a valid context table file\n"
		+ "in the Advanced Optons window by clicking on the Advanced\n"
		+ "Options button.";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    }

	    f = null;
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

	s = options.getConsensus();
	if (s != null) {
	    e = "Please specify a valid consensus/reference sequence file\n"
		+ "in the Advanced Optons window by clicking on the Advanced\n"
		+ "Options button.";
	    if (s.length() == 0) {
		errMesg = e;
		return false;
	    }

	    f = null;
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

	return true;
    }

    /** Creates an Options object from the current options in both 
	Opitions window and the Advanced Options window. */
    protected Options currentOptions() {
	Options temp = new Options();
	boolean b;
	String s;
	b = (phdBox.isSelected()) ? true : false;
	temp.setOutputPhd(b); 

	b = (scfBox.isSelected()) ? true : false; 
	temp.setOutputScf(b);

	b = (qualBox.isSelected()) ? true : false;
	temp.setOutputQual(b);

	b = qaBox.isSelected();
	temp.setOutputQa(b);

	s = b ? qaField.getText().trim() : null;
	temp.setQAFileName(s);

	b = (seqBox.isSelected()) ? true : false;
	temp.setOutputSeq(b);

	b = (saBox.isSelected()) ? true : false;
	temp.setOutputSa(b);

	s = b ? saField.getText().trim() : null;
	temp.setSAFileName(s);

	if (tableList.getSelectedIndex() == (tableList.getItemCount()-1)) {
	    temp.setOtherLookupTable(tableChooser.text1.getText().trim());
	    temp.setLookupTable(null);
	} else {
	    if (tableList.getSelectedIndex() == 0) {
		s = Constants.DEFAULT_CALIBRATION;
        } else if (tableList.getSelectedIndex() == 1) {
        s = Constants.ABI3730_POP7_LOOKUP_TABLE;
	    } else if (tableList.getSelectedIndex() == 2) {
		s = Constants.ABI3700_POP5_LOOKUP_TABLE;
	    } else if (tableList.getSelectedIndex() == 3) {
		s = Constants.ABI3700_POP6_LOOKUP_TABLE;
	    } else if (tableList.getSelectedIndex() == 4) {
		s = Constants.ABI3100_POP6_LOOKUP_TABLE;
	    } else if (tableList.getSelectedIndex() == 5) {
		s = Constants.MegaBACE_LOOKUP_TABLE;
	    } else {
	        s = Constants.USER_DIR + Constants.FILE_SEPARATOR
			+ (String)tableList.getSelectedItem() 
		        + Constants.LOOKUP_TABLE_SUFFIX;
	    }
	    temp.setLookupTable(s);
	    temp.setOtherLookupTable(null);
	}

	b = (noCallBox.isSelected()) ? true : false; 
	temp.setNoCall(b);

	b = (recallNBox.isSelected()) ? true : false; 
	temp.setRecallN(b);

	b = (hetBox.isSelected()) ? true : false; 
	temp.setHet(b);

        b = (mixBox.isSelected()) ? true : false;
        temp.setMix(b);

	s = (minRatioBox.isEnabled() && minRatioBox.isSelected())
	    		? minRatioField.getText().trim() : null;
	temp.setMinRatio(s);

	s = (outDirBox.isSelected()) ? outDirChooser.text1.getText().trim()
	    		             : null;
	temp.setOutputDir(s);

	// The advanced options
	temp.setOutputTip(options.outputTip());
	temp.setOutputTal(options.outputTal());
	temp.setOutputTab(options.outputTab());
	temp.setEdBase(options.editedBase());
	temp.setShift(options.shift());
	s = options.getTrimWindowSize();
	if (s != null) {
	    temp.setTrimWindowSize(new String(s));
	} else {
	    temp.setTrimWindowSize(null);
	}
	s = options.getTrimThreshold();
	if (s != null) {
	    temp.setTrimThreshold(new String(s));
	} else {
	    temp.setTrimThreshold(null);
	}
	s = options.getContextTable();
	if (s != null) {
	    temp.setContextTable(new String(s));
	} else {
	    temp.setContextTable(null);
	}
	s = options.getConsensus();
	if (s != null) {
	    temp.setConsensus(new String(s));
	} else {
	    temp.setConsensus(null);
	}

	return temp;
    }
	
    /**
     * Saves the current options (in both Options window and Advanced
     * Options window) permanently. 
     * @return true if succeeded; false otherwise.
     */
    protected boolean saveOptions() {
	Options temp = currentOptions();
	temp.save();
	if (temp.getStatus() == Options.OK) {
	    options = null;
	    options = temp;
	    Constants.fromEdited = temp.editedBase();
	    return true;
	} else {
	    Constants.fromEdited = options.editedBase();
	    JOptionPane.showMessageDialog(this, temp.getErrorMessage());
	    return false;
	}
    }

    /** Resets the options on the panels to default. */
    protected void reset() {
	phdBox.setSelected(true);
	scfBox.setSelected(false);
	qualBox.setSelected(false);
	seqBox.setSelected(false);

	if (qaBox.isSelected()) {
	    leftOutputPanel.remove(qaPanel);
	    qaBox.setSelected(false);
	}
	if (saBox.isSelected()) {
	    rightOutputPanel.remove(saPanel);
	    saBox.setSelected(false);
	}

	if (tableList.getSelectedIndex() == (tableList.getItemCount()-1)) {
	    processingPanel.remove(tableChooser);
	}
	tableList.setSelectedIndex(0);
	
	noCallBox.setSelected(false);
	recallNBox.setSelected(false);
	hetBox.setSelected(false);
        mixBox.setSelected(false);
	minRatioBox.setEnabled(false);
	minRatioField.setEnabled(false);

	if (minRatioBox.isSelected()) {
	    minRatioPanel.remove(minRatioField);
	    minRatioBox.setSelected(false);
	}

	if (outDirBox.isSelected()) {
	    programPanel.remove(outDirChooser);
	    outDirBox.setSelected(false);
	}

	if (isVisible()) {
	    programPanel.revalidate();
	    repaint();
	}
    }

    /** Gets the output directory.
     * @return "" if no output directory specified; the output directory 
     *         name if the directory is specified, it exists and is
     *	       writable; null, otherwise.
     */
    public String getOutputDir() {
	String s, e, errMsg;
	File f;
	if (isVisible()) {
	    if (outDirBox.isSelected()) {
	    	s = outDirChooser.text1.getText().trim();
	    	e = "Please specify a valid output directory in\n"
		          + "the Output Directory text box or turn off\n"
		          + "the Specified Output Directory Option.";
	    	if (s.length() == 0) {
		    errMsg = e;
	    	    if (getState() == Frame.ICONIFIED) {
			setState(Frame.NORMAL);
		    }
	    	    JOptionPane.showMessageDialog(this, errMsg);
	    	    return null;
	    	}

	    	f = null;
	    	f = new File(s);
	    	if (!f.exists()) {
		    errMsg = "The specified output directory "
				+ "does not exist.\n"+e;
	    	    if (getState() == Frame.ICONIFIED) {
			setState(Frame.NORMAL);
		    }
		    JOptionPane.showMessageDialog(this, errMsg);
		    return null;
	    	}
	    	if (!f.canWrite()) {
		    errMsg = "Write permission is denied to the specified"
		          + " output directory.\n" + e;
	    	    if (getState() == Frame.ICONIFIED) {
			setState(Frame.NORMAL);
		    }
		    JOptionPane.showMessageDialog(this, errMsg);
		    return null;
		}

		return s;
	    } else {
		// outputdir not specified
		return "";
	    }
	} else {
	    String ret = options.getOutputDir();
	    if (ret == null) {
		return "";
	    } else {
	    	return ret;
	    }
	}
    }

    /** 
     * Builds the command array using the specified output directory. 
     * @return the command builder if the options are valid, null otherwise.
     */
    public CommandBuilder buildCommand(String outDir) {
	if (isVisible()) {
	    if (!validateOptions()) {
	    	if (getState() == Frame.ICONIFIED) {
		    setState(Frame.NORMAL);
		}
	    	JOptionPane.showMessageDialog(this, errMesg);
	    	return null;
	    } else {
		Options temp = currentOptions();
	    	return temp.buildCommand(outDir);
	    }
	} else {
	    return options.buildCommand(outDir);
	}
    }

    /** Invoked when an action is performed. */
    public void actionPerformed(ActionEvent e) {
	Object source = e.getSource();
	if (source == tableList) {
	    if (tableList.getSelectedIndex() == (tableList.getItemCount()-1)) {
		processingPanel.add(tableChooser, 1);
	    } else {
		processingPanel.remove(tableChooser);
	    }
	    processingPanel.revalidate();
	    processingPanel.repaint();
	} else if (source == noCallBox) {
	    if (noCallBox.isSelected()) {
	    	recallNBox.setSelected(false);
	    	hetBox.setSelected(false);
                mixBox.setSelected(false);
		minRatioBox.setEnabled(false);
		minRatioField.setEnabled(false);
	    }
	} else if (source == recallNBox) {
	    if (recallNBox.isSelected()) {
	    	noCallBox.setSelected(false);
	    	hetBox.setSelected(false);
                mixBox.setSelected(false); 
		minRatioBox.setEnabled(false);
		minRatioField.setEnabled(false);
	    }
	} else if (source == hetBox) {
	    if (hetBox.isSelected()) {
	    	noCallBox.setSelected(false);
	    	recallNBox.setSelected(false);
                mixBox.setSelected(false);
		minRatioBox.setEnabled(true);
		minRatioField.setEnabled(true);
	    } else {
		minRatioBox.setEnabled(false);
		minRatioField.setEnabled(false);
	    }
        } else if (source == mixBox) {
            if (mixBox.isSelected()) {
                noCallBox.setSelected(false);
                recallNBox.setSelected(false);
                hetBox.setSelected(false);
                minRatioBox.setEnabled(true);
                minRatioField.setEnabled(true);
            } else {
                minRatioBox.setEnabled(false);
                minRatioField.setEnabled(false);
            }
	} else if (source == qaBox) {
	    if (qaBox.isSelected()) {
		leftOutputPanel.add(qaPanel);
	    } else {
		leftOutputPanel.remove(qaPanel);
	    }
	    outputPanel.revalidate();
	    outputPanel.repaint();
	} else if (source == saBox) {
	    if (saBox.isSelected()) {
		rightOutputPanel.add(saPanel);
	    } else {
		rightOutputPanel.remove(saPanel);
	    }
	    outputPanel.revalidate();
	    outputPanel.repaint();
	} else if (source == minRatioBox) {
	    if (minRatioBox.isSelected()) {
		minRatioPanel.add(minRatioField);
	    } else {
		minRatioPanel.remove(minRatioField);
	    }
	    programPanel.revalidate();
	    programPanel.repaint();
	} else if (source == outDirBox ) {
	    if (outDirBox.isSelected()) {
		programPanel.add(outDirChooser, 5);
	    } else {
		programPanel.remove(outDirChooser);
	    }
	    programPanel.revalidate();
	    programPanel.repaint();
	} else if (source == saveBtn) {
	    if (validateOptions()) {
		if (saveOptions()) {
	    	    setVisible(false);
		} 
	    } else {
	        JOptionPane.showMessageDialog(this, errMesg);
	    }
	} else if (source == cancelBtn) {
	    // remove the unsaved changes and restore the original options
	    options = null;
	    options = lastSavedOptions;
	    Constants.fromEdited = options.editedBase();
	    setVisible(false);
	} else if (source == resetBtn) {
	    reset();
	} else if (source == advBtn) {
	    // Makes the Advanced Options window visible and passes
	    // it the Options object.
	    advOptDialog.becomeVisible(options);
	}
    }
}
