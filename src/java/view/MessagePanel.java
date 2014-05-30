/*
 * 1.11
 *
 * @(#)MessagePanel.java       1.0     2000/10/25
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.table.*;
import com.paracel.tt.chrom.*;
import com.paracel.tt.io.*;
import com.paracel.tt.util.*;
import com.paracel.tt.event.*;

/**
 * This class creates the lower panel of TTView window, which contains
 * information such as file index, file names, file status, signal values.
 */
public class MessagePanel extends JSplitPane 
			  implements DisplaySignalValueEventListener,
				     FileNavigateEventListener,
				     SearchConstants {

    static final Font titleFont = new Font("Dialog", Font.BOLD, 11);
    static final Font labelFont = new Font("Arial", Font.PLAIN, 11);

    /* The TTViewFrame this panel resides in. */
    private JFrame frame;

    /* The sub panel containing file related information. */
    private JPanel fileInfoPanel;

    /* The sub panel containing SNPs related information. */
    private JScrollPane mixScrollPane;

    /* The text area containing SNPs related information. */
    private JTextArea mixSummaryArea;

    /* The table that lists the indices, codes, positions and QVs
	of all the detected hets/mixed bases in the current shown phd file. */
    private JTable mixTable;

    /* The table model used by the detected het/mixed bases table. */
    private SNPsTableModel mixTableModel;

    /* The sub panel containing trace signal related information. */
    private JScrollPane signalValuePanel;

    /* The splitpanel containing the file info panel and the het/mixed bases
	table panel. */
    private JSplitPane leftPanel;

    /* Labels used in the fileInfoPanel. */
    private JLabel fileIndex, chromName, phdName, phd2Name, tabName, talName, 
		   tipName, chromStatus, phdStatus, phd2Status, tabStatus,
		   talStatus, tipStatus;

    /* Labels used in the signalValuePanel. */
    private JLabel location, ASignal, CSignal, GSignal, TSignal,
		   fifthDyeSignal, fifthDyeLabel;

    /* The array containing trace signal values. */
    private int[] signalValues;

    /* Whether signal value sub panel is shown. */
    private boolean signalValuePanelShown;

    /* Whether signal values exist. */
    private boolean signalValuesExisting;

    /* Whether fifth dye data signal value exists. */
    private boolean fifthDyeDataExisting;

    /* The other phd file name if 2 phd files are shown in one viewer. */
    private String phd2FileName;

    /* The other phd file status if 2 phd files are shown in one viewer. */
    private int phd2StatusValue;

    /* Whether to display the other phd file. */
    private boolean displayPhd2;

    /* SearchEvent listeners. */
    private Vector searchListeners = new Vector();

    /*
     * Creates a MessagePanel displaying file and signal values information,
     * which is related to the files currently shown in the viewer.
     * @param	f	the TTViewFrame this panel resides in.
     * @param	phd2Name	the other phd file name.
     * @param	phd2Status	the other phd file status.
     */
    public MessagePanel(JFrame f, String phd2Name, int phd2Status) {
	super(JSplitPane.HORIZONTAL_SPLIT);
	setBackground(Color.white);
	frame = f;
	signalValuesExisting = false;
	fifthDyeDataExisting = false;
	phd2FileName = phd2Name;
	if (phd2Name == null) {
	    displayPhd2 = false; 
	} else {
	    displayPhd2 = true; 
	}
	phd2StatusValue = phd2Status;
	createUI();
    }

    /** Creates the UI of both sub panels and add them into a split panel. */
    protected void createUI() {
	createLeftPanel();
	createSignalValuePanel();

	setLeftComponent(leftPanel);
	setDividerSize(4);
	setResizeWeight(0.8);

	setPreferredSize(new Dimension(frame.getSize().width, 80));

	signalValuePanelShown = false;
    }

    /* Creates the left panel, which is a split pane with the file info
     * panel at left and the detected het/mixed bases table panel at right. 
     */
    protected void createLeftPanel() {
	createFileInfoPanel();
	createSNPsTablePanel();
	JScrollPane s = new JScrollPane(fileInfoPanel);
	s.getHorizontalScrollBar().setUnitIncrement(10);
	s.getVerticalScrollBar().setUnitIncrement(5);
	leftPanel = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
	leftPanel.setLeftComponent(s);
	leftPanel.setDividerSize(4);
	leftPanel.setDividerLocation(0.7);
	leftPanel.setResizeWeight(0.8);
    }

    /** Creates the file info sub panel. */
    protected void createFileInfoPanel() {
	int row = Constants.runAsDev ? 6 : 5;

	if (displayPhd2) {
	    row++;
	}
	
	fileInfoPanel = new JPanel();
	fileInfoPanel.setBackground(Color.white);
	fileInfoPanel.setBorder(new EmptyBorder(-6, -2, -6, -2));
	fileInfoPanel.setLayout(new FlowLayout(FlowLayout.LEFT)); 

	JPanel left = new JPanel();
	left.setBackground(Color.white);
	left.setLayout(new GridLayout(row, 1, 0, -1));
	JPanel middle = new JPanel();
	middle.setBackground(Color.white);
	middle.setLayout(new GridLayout(row, 1, 0, -1));
	JPanel right = new JPanel();
	right.setBackground(Color.white);
	right.setLayout(new GridLayout(row, 1, 0, -1));

	JLabel label0 = makeLabel("File Index:", null);
	JLabel label1 = makeLabel("Sample File:", null);
	JLabel label2 = makeLabel("Phd File:", null);
	JLabel label6 = makeLabel("Tab File:", null);

	left.add(label0);
	left.add(label1);
	left.add(label2);
	left.add(label6);
	
	fileIndex = makeLabel(null, Color.black);
	chromName = makeLabel(null, Color.black);
	phdName = makeLabel(null, Color.black);
	tabName = makeLabel(null, Color.black);

	middle.add(fileIndex);
	middle.add(chromName);
	middle.add(phdName);
	middle.add(tabName);

	JLabel statusLabel = makeLabel(" ", null);
	chromStatus = makeLabel(null, null, titleFont);
	phdStatus = makeLabel(null, null, titleFont);
	tabStatus = makeLabel(null, null, titleFont);

	right.add(statusLabel);
	right.add(chromStatus);
	right.add(phdStatus);
	right.add(tabStatus);

	if (displayPhd2) {
	    JLabel label3 = makeLabel("2nd Phd File:", null);
	    phd2Name = makeLabel(null, Color.black);
	    phd2Status = makeLabel(null, null, titleFont);

	    left.add(label3);
	    middle.add(phd2Name);
	    right.add(phd2Status);
	}

	JLabel label4 = makeLabel("Tal File:", null);
	talName = makeLabel(null, Color.black);
	talStatus = makeLabel(null, null, titleFont);

	left.add(label4);
	middle.add(talName);
	right.add(talStatus);

	JLabel label5 = makeLabel("Tip File:", null);
	tipName = makeLabel(null, Color.black);
	tipStatus = makeLabel(null, null, titleFont);

        left.add(label5);
        middle.add(tipName);
        right.add(tipStatus);

	fileInfoPanel.add(left);	
	fileInfoPanel.add(Box.createHorizontalStrut(2));	
	fileInfoPanel.add(middle);	
	fileInfoPanel.add(Box.createHorizontalStrut(2));	
	fileInfoPanel.add(right);	
    }

    /* Creates the scrollpanel that wraps the detected het/mixed table. */
    protected void createSNPsTablePanel() {
	// Constructs the talbe and set attributes.
	mixTableModel = new SNPsTableModel();
	mixTable = new JTable();
	mixTable.setAutoCreateColumnsFromModel(false);
	mixTable.setModel(mixTableModel);
	mixTable.setRowSelectionAllowed(true);
	mixTable.setShowHorizontalLines(false);
	mixTable.setShowVerticalLines(false);
	mixTable.setShowGrid(false);
	mixTable.getSelectionModel().setSelectionMode(
					ListSelectionModel.SINGLE_SELECTION);

	// Customizes the table columns.
	for (int i = 0; i < SNPsTableModel.columns.length; i++) {
	    DefaultTableCellRenderer renderer = new DefaultTableCellRenderer();
	    renderer.setHorizontalAlignment(
					SNPsTableModel.columns[i].alignment);
	    TableColumn column = new TableColumn(i,
					SNPsTableModel.columns[i].width, 
					renderer, null);
	    mixTable.addColumn(column);
	}

	// Sets the header attributes and adds listener for sorting purposes.
	JTableHeader header = mixTable.getTableHeader();
	header.setUpdateTableInRealTime(true);
	header.addMouseListener(mixTableModel.new ColumnListener(mixTable));

	// Adds mouse listener to the SNP table.
	mixTable.addMouseListener(new MouseAdapter() {
	    public void mouseClicked(MouseEvent me) {
		/* If the user double-clicks in the SNP table, the view
		   will be scrolled to show the clicked SNP. */
		if (me.getClickCount() == 2) {
		    Point p = me.getPoint();
		    int row = mixTable.rowAtPoint(p);
		    if (row < 0) {
		    	return;
		    }
		    /* The reason of " - 1" is that the indexes in the table
		       use 1-offset whereas the java arrays use zero-offset. */
		    int index = ((Integer)mixTable.getValueAt(row, 0))
							.intValue() - 1;
		    String base = ((Character)mixTable.getValueAt(row, 1))
								.toString();
		    /* Borrows the SearchEvent facility to draw a black solid
		       arrow above the clicked SNP in the chrom view. */
		    SearchEvent se = new SearchEvent(MessagePanel.this,
						     TT, base, FORWARD, index);
		    for (int i = 0; i < searchListeners.size(); i++) {
		    	((SearchEventListener) searchListeners.elementAt(i))
			    		.search(se);
		    }
		}
	    }
	});

	mixScrollPane = new JScrollPane(mixTable);
	mixScrollPane.getViewport().setBackground(Color.white);
	mixScrollPane.getHorizontalScrollBar().setUnitIncrement(10);
	mixScrollPane.getVerticalScrollBar().setUnitIncrement(5);
    }

    /** Make a label with the specified text, color, and font. */
    protected JLabel makeLabel(String text, Color c, Font f) {
	JLabel l = makeLabel(text, c);
	l.setFont(f);
	return l;
    }

    /** Make a label with the specified text, color, and "labelFont". */
    protected JLabel makeLabel(String text, Color c) {
	JLabel l;
	if (text != null) {
	    l = new JLabel(text);
	} else {
	    l = new JLabel();
	}
	if (c != null) {
	    l.setForeground(c);
	}
	l.setFont(labelFont);
	return l;
    }

    /** Creates the signal value sub panel. */
    protected void createSignalValuePanel() {
	JPanel panel = new JPanel();
	panel.setBackground(Color.white);
	panel.setBorder(new EmptyBorder(-6, -2, -6, -2));
	panel.setLayout(new FlowLayout(FlowLayout.LEFT));

	JLabel label3 = makeLabel("X Loc:", Color.black,titleFont);
	JLabel label4 = makeLabel("A:", Trace.findTraceColor("A"), titleFont);
	JLabel label5 = makeLabel("C:", Trace.findTraceColor("C"), titleFont);
	JLabel label6 = makeLabel("G:", Trace.findTraceColor("G"), titleFont);
	JLabel label7 = makeLabel("T:", Trace.findTraceColor("T"), titleFont);
	fifthDyeLabel = makeLabel("5th-Dye:", Color.white, titleFont);

	location = makeLabel(null, Color.black);
	ASignal = makeLabel(null, Color.black);
	CSignal = makeLabel(null, Color.black);
	GSignal = makeLabel(null, Color.black);
	TSignal = makeLabel(null, Color.black);
	fifthDyeSignal = makeLabel(null, Color.black);

	JPanel p1 = new JPanel();
	JPanel p2 = new JPanel();

	p1.setBackground(Color.white);
	p1.setLayout(new GridLayout(6, 1, 0, -1));
	p2.setBackground(Color.white);
	p2.setLayout(new GridLayout(6, 1, 0, 0));

	p1.add(label3);
	p1.add(label4);
	p1.add(label5);
	p1.add(label6);
	p1.add(label7);
	p1.add(fifthDyeLabel);

	p2.add(location);
	p2.add(ASignal);
	p2.add(CSignal);
	p2.add(GSignal);
	p2.add(TSignal);
	p2.add(fifthDyeSignal);

	panel.add(p1);
	panel.add(Box.createHorizontalStrut(2));	
	panel.add(p2);
	signalValuePanel = new JScrollPane(panel);
	signalValuePanel.getHorizontalScrollBar().setUnitIncrement(10);
	signalValuePanel.getVerticalScrollBar().setUnitIncrement(5);
    }
    
    /** Attempts to release the resouces. */
    public void releaseResources() {
	fileInfoPanel = null;
	mixScrollPane = null;
	mixSummaryArea = null;
	mixTable = null;
	mixTableModel = null;
	signalValuePanel = null;
	leftPanel = null;
	searchListeners = null;
	System.gc();
    }
    
    /** Invoked when the file index changed, i.e., navigating to 
	another file.  Sets the file index, file names, and status
        on the file info sub panel. */
    public void navigateTo(FileNavigateEvent e) {
	int ps;
	String pName;
	String temp = e.getFilesStatus().getIndex()
		     	+ " of " + e.getNumOfFiles();
	fileIndex.setText(temp);
	
	temp = e.getFilesStatus().getSampleFileName();
	chromName.setText(temp);
	
	pName = e.getFilesStatus().getPhdFileName();
	phdName.setText(pName);

	temp = e.getFilesStatus().getTabFileName();
	tabName.setText(temp);

	chromStatus.setText(FilesStatus.getFileStatusString(e.getFilesStatus().
					        getSampleFileStatus()));
	chromStatus.setForeground(FilesStatus.getFileStatusStringColor(
				    e.getFilesStatus().getSampleFileStatus()));

	ps = e.getFilesStatus().getPhdFileStatus();
	phdStatus.setText(FilesStatus.getFileStatusString(ps));
	phdStatus.setForeground(FilesStatus.getFileStatusStringColor(ps));

	tabStatus.setText(FilesStatus.getFileStatusString(e.getFilesStatus().
					      getTabFileStatus()));
	tabStatus.setForeground(FilesStatus.getFileStatusStringColor(
				    e.getFilesStatus().getTabFileStatus()));

	temp = e.getFilesStatus().getTalFileName();
	talName.setText(temp);
	talStatus.setText(FilesStatus.getFileStatusString(
				  e.getFilesStatus().getTalFileStatus()));
	talStatus.setForeground(FilesStatus.getFileStatusStringColor(
				       e.getFilesStatus().getTalFileStatus()));
	
	temp = e.getFilesStatus().getTipFileName();
	tipName.setText(temp);
	tipStatus.setText(FilesStatus.getFileStatusString(
				  e.getFilesStatus().getTipFileStatus()));
	tipStatus.setForeground(FilesStatus.getFileStatusStringColor(
				       e.getFilesStatus().getTipFileStatus()));
	
	if (displayPhd2) {
	    temp = phd2FileName;
	    phd2Name.setText(temp);
	    phd2Status.setText(FilesStatus.getFileStatusString(
						   phd2StatusValue));
	    phd2Status.setForeground(FilesStatus.getFileStatusStringColor(
						   phd2StatusValue));
	}

	if (ps == FilesStatus.VALID) {
	    setSNPsTableData(pName);
	} else {
	    clearSNPsTableData();
	}
	clearSignalValuePanel();
    }

    /** Registers the specified SearchEventListener. */
    public void addSearchEventListener(SearchEventListener sel) {
	searchListeners.add(sel);
    }

    /** Adds the signal value sub panel to the right side. */
    public void addSignalValuePanel() {
	if (getRightComponent() == null) {
	    setRightComponent(signalValuePanel);
	    setDividerLocation(0.85);
	    signalValuePanelShown = true;
	}
    }

    /** Removes the signal value sub panel. */
    public void removeSignalValuePanel() {
	remove(signalValuePanel);
	signalValuePanelShown = false;
    }

    /* Adds the detected hets/mixed bases table panel to the right of the
     * file info panel. 
     */
    public void addSNPsTablePanel() {
	if (leftPanel.getRightComponent() == null) {
	    leftPanel.setRightComponent(mixScrollPane);
	}
    }

    /** Removes the detected hets/mixed bases table panel. */
    public void removeSNPsTablePanel() {
	if (leftPanel.getRightComponent() != null) {
	    leftPanel.remove(mixScrollPane);
	}
    }

    /* Sets the detected hets/mixed bases table with the data in the
     * specified phd file. 
     */
    protected void setSNPsTableData(String pName) {
	try {
	    PhdFile p = new PhdFile(pName);
		
	    int num = p.getNumOfSNPs();
	    mixTableModel.setData(num, p.getSnpIndices(),
				   p.getSnpBases(), p.getSnpPositions(),
				   p.getSnpQVs());
	} catch (Exception e) {
	    clearSNPsTableData();
	}
    }

    /* Clears the detected hets/mixed bases table data. */
    protected void clearSNPsTableData() {
    	mixTableModel.setData(0, null, null, null, null);
    }
    
    /** Invoked when the signal values are set.  Sets the trace signal
        values on the signal value sub panel. */
    public void signalValueSet(DisplaySignalValueEvent e) {
	signalValuesExisting = true;
	signalValues = e.getSignalValues();

	location.setText(String.valueOf(signalValues[0]));
	ASignal.setText(String.valueOf(signalValues[1]));
	CSignal.setText(String.valueOf(signalValues[2]));
	GSignal.setText(String.valueOf(signalValues[3]));
	TSignal.setText(String.valueOf(signalValues[4]));
	if (e.fifthDyeDataExists()) {
	    fifthDyeDataExisting = true;
	    fifthDyeLabel.setForeground(Trace.findTraceColor("Fifth Dye"));
	    fifthDyeSignal.setText(String.valueOf(signalValues[5]));
	} else {
	    fifthDyeLabel.setForeground(Color.white);
	    fifthDyeSignal.setText(" ");
	}
    }

    /** Clears the signal value panel. */
    protected void clearSignalValuePanel() {
	fifthDyeLabel.setForeground(Color.white);
	location.setText(" ");
	ASignal.setText(" ");
	CSignal.setText(" ");
	GSignal.setText(" ");
	TSignal.setText(" ");
	fifthDyeSignal.setText(" ");
	signalValuesExisting = false;
	fifthDyeDataExisting = false;
    }

    /** Returns true if signal values are shown; false otherwise. */
    public boolean signalValuesShown() { 
	return signalValuePanelShown && signalValuesExisting; 
    }

    /** Returns true if fifth dye data is shown; false otherwise. */
    public boolean fifthDyeDataShown() { 
	return signalValuePanelShown && fifthDyeDataExisting; 
    }

    /** Returns the trace signal values. */
    public int[] getSignalValues() { return signalValues; }
}
