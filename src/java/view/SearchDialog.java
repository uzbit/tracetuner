/**
 * 1.5
 *
 * @(#)SearchDialog.java	1.0	2001/03/21
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
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import com.paracel.tt.event.*;
import com.paracel.tt.util.*;

/**
 * Creates the "Open file" dialog, which is invoked by selecting the
 * "Open" option in the "file" menu of TTViewFrame.
 */
public class SearchDialog extends JDialog implements ActionListener,
						     SearchConstants {

    /** The parent frame. */
    private JFrame frame;

    /** The text field to get the search pattern. */
    private JTextField patternField;

    /** The combo box for the search subject selection. */
    private JComboBox subjectBox;

    /** The "Next Heterozygote", "Specify a sequence", "Right-to-left",
	and "Left-to-right" radio buttons. */
    private JRadioButton nextSNPButton, specifyButton, leftButton, rightButton;

    /** The "Find Next", "Close", and "Help" buttons. */
    private JButton btnFind, btnClose, btnHelp;

    /** The search pattern specified by the user. */
    private String pattern;

    /** The search subject type and the search direction. */
    private int subject, direction;

    /** The search event listeners vector. */
    private Vector listeners;

    /** The search event. */
    private SearchEvent searchEvent;

    /**
     * Creates a SearchDialog that requests input from the user to construct
     * a SearchEvent, and allows the user to execute such a SearchEvent.
     */
    public SearchDialog(JFrame frame) {
	super(frame);
	String s = frame.getTitle();
	setTitle("Find ... " + s.substring(s.lastIndexOf('['), s.length()));
	this.frame = frame;

	JPanel c = (JPanel) getContentPane();
	c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));

	// "top" panel contains "Find what" info
	JPanel top = new JPanel(new FlowLayout(FlowLayout.LEFT));

	JPanel p1 = new JPanel(new GridLayout(3, 1));
	p1.setPreferredSize(new Dimension(78, 75));
	JPanel p2 = new JPanel(new GridLayout(3, 1));
	p2.setPreferredSize(new Dimension(150, 75));
	JPanel p3 = new JPanel(new GridLayout(3, 1));
	p3.setPreferredSize(new Dimension(200, 75));

	nextSNPButton = new JRadioButton("Next het/mixed base", true);
	specifyButton = new JRadioButton("Specify sequence(s)", false);
	nextSNPButton.addActionListener(this);
	specifyButton.addActionListener(this);
	ButtonGroup bg1 = new ButtonGroup();
	bg1.add(nextSNPButton);
	bg1.add(specifyButton);
	patternField = new JTextField();
	patternField.setEnabled(false);
	patternField.addActionListener(this);

	p1.add(new JLabel("Find what:"));
	p1.add(new JLabel(""));
	p1.add(new JLabel(""));
	p2.add(nextSNPButton);
	p2.add(new JLabel("or"));
	p2.add(specifyButton);
	p3.add(new JLabel(""));
	p3.add(new JLabel(""));
	p3.add(patternField);

	top.add(Box.createHorizontalStrut(10));
	top.add(p1);
	top.add(p2);
	top.add(p3);

	// "middle" panel contains search attributes, such as where to
	// search, which direction.
	JPanel middle = new JPanel();
	middle.setLayout(new BoxLayout(middle, BoxLayout.Y_AXIS));
	JPanel p4 = new JPanel(new FlowLayout(FlowLayout.LEFT));
	p4.add(Box.createHorizontalStrut(5));
	p4.add(new JLabel("Search:"));
	p4.add(Box.createHorizontalStrut(24));
	String[] subjects = {"TraceTuner base calls",
			    "Original (ABI) base calls",
			    "Consensus/reference sequence"};
	subjectBox = new JComboBox(subjects);
	p4.add(subjectBox);

	p4.add(Box.createHorizontalStrut(15));
	p4.add(new JLabel("Direction:"));
	JPanel directionPanel = new JPanel(new GridLayout(2, 1));
	rightButton = new JRadioButton("Left-to-right", true);
	leftButton = new JRadioButton("Right-to-left", false);
	ButtonGroup bg = new ButtonGroup();
	bg.add(rightButton);
	bg.add(leftButton);
	directionPanel.add(rightButton);
	directionPanel.add(leftButton);
	p4.add(Box.createHorizontalStrut(5));
	p4.add(directionPanel);

	middle.add(p4);
        Border b = BorderFactory.createEtchedBorder();
        TitledBorder tb = BorderFactory.createTitledBorder(
					b, "Search options");
        tb.setTitleJustification(TitledBorder.LEFT);
	middle.setBorder(tb);

	//"bottom" panel contains the buttons.
	JPanel bottom = new JPanel(new FlowLayout(FlowLayout.RIGHT));
	btnFind = new JButton("Find Next");
	btnClose = new JButton("Close");
	btnHelp = new JButton(new HelpAction("Help"));

	btnFind.addActionListener(this);
	btnClose.addActionListener(this);

	JPanel p = new JPanel(new GridLayout(1, 3, 6, 0));
	p.add(btnFind);
	p.add(btnClose);
	p.add(btnHelp);
	bottom.add(p);

	c.add(Box.createVerticalStrut(5));
	c.add(top);
	c.add(Box.createVerticalStrut(5));
	c.add(middle);
	c.add(Box.createVerticalStrut(10));
	c.add(new JSeparator());
	c.add(Box.createVerticalStrut(10));
	c.add(bottom);

	pack();
	listeners = new Vector();
	searchEvent = null;
    }

    /** Invoked when an action was performed. */
    public void actionPerformed(ActionEvent e) {
	Object objSource = e.getSource();
	if (objSource == btnClose) {
	    setVisible(false);
	} else if (objSource == nextSNPButton || objSource == specifyButton) {
	    patternField.setEnabled(specifyButton.isSelected());
	} else if (objSource == patternField || objSource == btnFind) {
	    // search-for-what
	    if (nextSNPButton.isSelected()) {
		// the space between the bases means "or"
		pattern = "R Y K M S W";
	    } else {
		pattern = patternField.getText().trim();
		if (pattern.length() <= 0) {
		    JOptionPane.showMessageDialog(frame,
						  "Nothing to search for.");
		    patternField.requestFocus();
		    return;
		}
		// validate pattern
		char c;
		int numOfBases = 0;
		for (int i = 0; i < pattern.length(); i++) {
		    c = pattern.charAt(i);
		    // ; , whitespaces can be used to delimit multiple
		    // sequences
		    if (c == ';' || c == ',' || Character.isWhitespace(c)) {
			continue;
		    } 
		    numOfBases++;
		    if (!IUBCode.isValid(c)) {
			int result = JOptionPane.showConfirmDialog(frame,
				"Invalid base code \"" + c
				+ "\" found in the specified "
				+ "sequence(s).  Start searching anyway?",
				"Question",
				JOptionPane.YES_NO_OPTION);
			if (result == JOptionPane.YES_OPTION) {
			    break;
			} else {
			    patternField.requestFocus();
			    return;
			}
		    }
		}
		if (numOfBases == 0) {
		    JOptionPane.showMessageDialog(frame,
						  "Nothing to search for.");
		    patternField.requestFocus();
		    return;
		}
	    }

	    // where-to-search
	    int ind = subjectBox.getSelectedIndex();
	    switch (ind) {
		case 0: subject = TT; break;
		case 1: subject = ABI; break;
		case 2: subject = REF; break;
		default: subject = TT; break;
	    }

	    direction = (rightButton.isSelected()) ? FORWARD : BACKWARD;

	    searchEvent = new SearchEvent(frame, subject, pattern,
					     direction);
	    for (int i = 0; i < listeners.size(); i++) {
		((SearchEventListener)listeners.elementAt(i))
		    				.search(searchEvent);
	    }
	}
    }

    /** Registers the specified SearchEvent Listener. */
    public void addSearchListener(SearchEventListener l) {
	listeners.add(l);
    }

    /** Centers the search dialog (relative to the Viewer window), and
	makes it visible. */
    public void becomeVisible() {
	int x = frame.getLocation().x;
	int y = frame.getLocation().y;
	int w1 = frame.getWidth();
	int h1 = frame.getHeight();
	int w2 = getWidth();
	int h2 = getHeight();
	setLocation(x + (w1-w2)/2, y + (h1-h2)/2);
	setVisible(true);
    }

    /** Returns the current search event. */
    public SearchEvent getSearchEvent() { return searchEvent; }
}

