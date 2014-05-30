/**
 * 1.4
 *
 * @(#)CustomZoomDialog.java	1.0	2000/11/14
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
 * Creates the "Custom Zoom" dialog, which is invoked by selecting the
 * "Custom Zoom" options in the "view" menu of TTViewFrame.
 */
public class CustomZoomDialog extends JDialog implements ActionListener {

    /** The reference to the TTViewFrame where this dialog is originated. */
    private JFrame frame;

    /** The text field where the user can type in zoom values. */
    private JTextField textField;

    /** The "OK", and "Cancel" button. */
    private JButton btnOK, btnCancel;

    /** The minimum acceptable value. */
    private int min;

    /** The maximum acceptable value. */
    private int max;

    /** The zoom factor value calculated from the number the user entered. */
    private float zoomValue;

    /** The axis. */
    private char zoomAxis;

    /** Whether or not the dialog was canceled. */
    private boolean canceled;

    /**
     * Creates a CustomZoomDialog that request a zoom value input from 
     * the user, 
     * @param	frame	the source TTViewFrame.
     * @param	axis	the axis.
     * @param	min	the minimum acceptable value.
     * @param	max	the maximum acceptable value.
     */
    public CustomZoomDialog(JFrame frame, char axis, int min, int max) {
	super(frame, "Custom Zoom", true);
	this.frame = frame;
	this.zoomAxis = axis;
	this.min = min;
	this.max = max;

	JPanel c = (JPanel) getContentPane();
	c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));

	JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT));
	JLabel label = new JLabel("Enter a zoom " + axis + " value in %:");
	Dimension d = new Dimension(5, label.getSize().height);
	textField = new JTextField(6);
	textField.setBorder(BorderFactory.createLoweredBevelBorder());
	textField.addActionListener(this);

	top.add(Box.createRigidArea(d));
	top.add(label);
	top.add(Box.createRigidArea(d));
	top.add(textField);

	JPanel bottom = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	JPanel p = new JPanel(new GridLayout(1, 2, 6, 0));
	btnOK = new JButton("OK");
	btnCancel = new JButton("Cancel");

	btnOK.addActionListener(this);
	btnCancel.addActionListener(this);

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
	textField.selectAll();
	setVisible(true);
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
	    float value;
	    try {
		value = Float.parseFloat(temp);
		if (value < min || value > max) {
		    invalidInput();
		} else {
		    zoomValue = value/100.0f;
		    setVisible(false);
		}
	    } catch (NumberFormatException ne) {
		invalidInput();
	    }
	}
    }
    
    /**
     * Gets the calculated zoom factor.
     * @return	the zoom factor calculated from the input number.
     */
    public float getZoomValue() { return zoomValue; }

    /**
     * Indicates whether this dialog was canceled.
     * @return true if the dialog was canceled; false otherwise.
     */
    public boolean isCanceled() { return canceled; }

    /** Shows an error message if the input was invalid. */
    protected void invalidInput() {
	textField.selectAll();
	JOptionPane.showMessageDialog(this, "Invalid input.\n Expecting a "
				      + "value between " + min + " and "
				      + max);
    }
}

