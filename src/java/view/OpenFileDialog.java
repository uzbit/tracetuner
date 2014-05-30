/**
 * 1.6
 *
 * @(#)OpenFileDialog.java	1.0	2000/11/14
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
import javax.swing.*;
import javax.swing.border.*;
import com.paracel.tt.util.*;

/**
 * Creates the "Open file" dialog, which is invoked by selecting the
 * "Open" option in the "file" menu of TTViewFrame.
 */
public class OpenFileDialog extends JDialog implements ActionListener {

    /** The text field to get the sample file name from the user. */
    private JTextField cNameField;

    /** The text field to get the phd file name from the user. */
    private JTextField pNameField;

    /** The "browse" buttom for browsing the sample file. */
    private JButton    cNameBrowse;

    /** The "browse" buttom for browsing the phd file. */
    private JButton    pNameBrowse;

    /** The "OK", "Clear", and "Cancel" buttons. */
    private JButton    btnOK, btnClear, btnCancel;

    /** The frame this dialog is originated. */
    private JFrame frame;
    
    /** The file choosers for sample file and phd file. */
    private JFileChooser sampleChooser, phdChooser;

    /**
     * Creates a OpenFileDialog that request input from the user, indicating
     * which sample file and phd file the user wants to open.
     * @param	frame	The source TTViewFrame.
     */
    public OpenFileDialog(JFrame frame) {
	super(frame, "Open ...", true);

	this.frame = frame;

	JPanel c = (JPanel) getContentPane();
	c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));

	JPanel top = new JPanel();

	JPanel p1 = new JPanel();
	JPanel p2 = new JPanel();
	JPanel p3 = new JPanel();
	p1.setLayout(new GridLayout(2, 1, 2, 8));
	p2.setLayout(new GridLayout(2, 1, 2, 8));
	p3.setLayout(new GridLayout(2, 1, 2, 8));

	JLabel label1 = new JLabel("Sample file:");
	JLabel label2 = new JLabel("Phd file (Optional):");
	cNameField = new JTextField(20);
	pNameField = new JTextField(20);
	cNameBrowse = new JButton("Browse");
	pNameBrowse = new JButton("Browse");
	cNameBrowse.setMargin(new Insets(0, 2, 0, 2));
	pNameBrowse.setMargin(new Insets(0, 2, 0, 2));

	cNameBrowse.addActionListener(this);
	pNameBrowse.addActionListener(this);

	p1.add(label1);
	p2.add(cNameField);
	p3.add(cNameBrowse);
	p1.add(label2);
	p2.add(pNameField);
	p3.add(pNameBrowse);

	top.add(p1);
	top.add(Box.createHorizontalStrut(4));
	top.add(p2);
	top.add(Box.createHorizontalStrut(4));
	top.add(p3);
	JPanel bottom = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	btnOK = new JButton("OK");
	btnClear = new JButton("Clear");
	btnCancel = new JButton("Cancel");

	btnOK.addActionListener(this);
	btnClear.addActionListener(this);
	btnCancel.addActionListener(this);

	JPanel p = new JPanel(new GridLayout(1, 3, 6, 0));
	p.add(btnOK);
	p.add(btnClear);
	p.add(btnCancel);
	bottom.add(p);

	c.add(Box.createVerticalStrut(5));
	c.add(top);
	c.add(Box.createVerticalStrut(10));
	c.add(bottom);

	pack();

	sampleChooser = new JFileChooser(Constants.USER_DIR);
	sampleChooser.setDialogTitle("Select Sample File");
	Filter filterScf = new Filter(".scf");
        Filter filterAbi = new Filter(".abi");
        Filter filterAbd = new Filter(".abd");
        Filter filterAb1 = new Filter(".ab1");
        sampleChooser.addChoosableFileFilter(filterScf);
        sampleChooser.addChoosableFileFilter(filterAbi);
        sampleChooser.addChoosableFileFilter(filterAbi);
        sampleChooser.addChoosableFileFilter(filterAb1);
	sampleChooser.setFileFilter(filterAb1);
	sampleChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	
	phdChooser = new JFileChooser(Constants.USER_DIR);
	phdChooser.setDialogTitle("Select PHD file");
	Filter filterPhd = new Filter(Constants.PHD_FILE_SUFFIX);
	phdChooser.addChoosableFileFilter(filterPhd);
	phdChooser.setFileFilter(filterPhd);
	phdChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	
	int x = (int) frame.getLocation().getX();
	int y = (int) frame.getLocation().getY();
	setLocation(x + 50, y + 50);
	setVisible(true);
    }

    /** Invoked when an action was performed. */
    public void actionPerformed(ActionEvent e) {
	Object source = e.getSource();
	if (source == cNameBrowse) {
	    int returnVal = sampleChooser.showDialog(this, "Open");
	    if (returnVal == JFileChooser.APPROVE_OPTION) {
		cNameField.setText(sampleChooser
				   .getSelectedFile().getAbsolutePath());
	    }
	} else if (source == pNameBrowse) {
	    int returnVal = phdChooser.showDialog(this, "Open");
	    if (returnVal == JFileChooser.APPROVE_OPTION) {
		pNameField.setText(phdChooser
				   .getSelectedFile().getAbsolutePath());
	    }
	} else if (source == btnOK) {
	    String cName = cNameField.getText();
	    if (cName == null || (cName != null && cName.length() == 0)) {
		JOptionPane.showMessageDialog(this, 
					      "Sample File not selected.");
		return;
	    }

	    String pName = pNameField.getText();
	    if (pName == null || (pName != null && pName.length() == 0)) {
		pName = cName + Constants.PHD_FILE_SUFFIX;
	    }

	    setVisible(false);
	    String[] c = { cName };
	    String[] p = { pName };
	    new TTViewFrame(c, p, null, false);
	    dispose();
	} else if (source == btnClear) {
	    cNameField.setText("");
	    pNameField.setText("");
	} else if (source == btnCancel) {
	    setVisible(false);
	    dispose();
	}
    }
}

