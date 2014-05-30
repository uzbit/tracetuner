/**
 * 1.2
 *
 * @(#)GoToDialog.java	1.0	2000/11/14
 * 
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 *
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
 * Creates the "Go to ..." dialog, which is invoked by selecting the
 * "Go to ..." option in the "file" menu of TTViewFrame.
 */
public class GoToDialog extends JDialog implements ActionListener {

    /** The reference to the TTViewFrame this dialog is originated from. */
    private JFrame frame;

    /** The text field to accept the user input. */
    private JTextField textField;

    /** The "OK", and "Cancel" button. */
    private JButton btnOK, btnCancel;

    /** The total number of files. */
    private int numOfFiles;

    /** The file index to go to. */
    private int gotoIndex;

    /** Whether or not the dialog was canceled. */
    private boolean canceled;

    /**
     * Creates a GoToDialog that request an input from the user, indicating
     * which file the user wants to go to.
     * @param	frame	The source TTViewFrame.
     * @param	numOfFiles	The total number of files.
     */
    public GoToDialog(JFrame frame, int numOfFiles) {
	super(frame, "Go To ...", true);
	this.frame = frame;
	this.numOfFiles = numOfFiles;

	JPanel c = (JPanel) getContentPane();
	c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));

	JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT));
	JLabel label1 = new JLabel("Go to file:");
	JLabel label2 = new JLabel(" of " + numOfFiles);
	textField = new JTextField(5);
	textField.setBorder(BorderFactory.createLoweredBevelBorder());
	textField.setHorizontalAlignment(JTextField.RIGHT);
	textField.addActionListener(this);

	Dimension d = new Dimension(30, label1.getSize().height);
	top.add(Box.createRigidArea(d));
	top.add(label1);
	top.add(textField);
	top.add(label2);
	top.add(Box.createRigidArea(d));

	JPanel bottom = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	btnOK = new JButton("OK");
	btnCancel = new JButton("Cancel");
	btnOK.addActionListener(this);
	btnCancel.addActionListener(this);

	JPanel p = new JPanel(new GridLayout(1, 2, 6, 0));
	p.add(btnOK);
	p.add(btnCancel);
	bottom.add(p);

	c.add(Box.createVerticalStrut(5));
	c.add(top);
	c.add(Box.createVerticalStrut(15));
	c.add(bottom);

	addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {
		canceled = true;
	    }
	});
	
	pack();
    }

    public void becomeVisible() {
	int x = (int) frame.getLocation().getX();
	int y = (int) frame.getLocation().getY();
	setLocation(x + 50, y + 50);
	setVisible(true);
	textField.selectAll();
    }

    /** Invoked when an action is performed. */
    public void actionPerformed(ActionEvent e) {
	Object source = e.getSource();
	if (source == btnOK || source == textField) {
	    canceled = false;
	    OKPressed();
	} else if (source == btnCancel) {
	    canceled = true;
	    setVisible(false);
	}
    }

    /** The "OK" button was pressed. */
    protected void OKPressed() {
	String temp = textField.getText();
	if (temp == null) {
	    invalidInput();
	} else {
	    temp = temp.trim();
	    try {
		gotoIndex = Integer.parseInt(temp) - 1;
		if (gotoIndex < 0 || gotoIndex >= numOfFiles) {
		    invalidInput();
		} else {
		    setVisible(false);
		}
	    } catch (NumberFormatException ne) {
		invalidInput();
	    }
	}
    }
    
    /** 
     * Gets the file index to go to.
     * @return the file index entered by the user.
     */
    public int getGotoIndex() { return gotoIndex; }

    /** Shows an error message if the input was invalid. */
    protected void invalidInput() {
	textField.selectAll();
	JOptionPane.showMessageDialog(this, "Invalid input.");
    }

    /**
     * Indicates whether this dialog was canceled.
     * @return true if the dialog was canceled; false otherwise.
     */
    public boolean isCanceled() { return canceled; }
}

