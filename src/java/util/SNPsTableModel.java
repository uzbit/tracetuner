/*
 * 1.2
 *
 * @(#)SNPsTableModel.java	1.0	2001/01/26
 *
 * Copyright (c) 2001 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.util;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.table.*;
import javax.swing.event.*;
import java.util.*;

/**
 * This class creates the table model for the sortable detected het/mix bases   
 * table, which is included in the viewer window.
 */
public class SNPsTableModel extends AbstractTableModel {

    /** This inner class encapsulates information about each TableColumn,
	including the title, width, and alignment.*/
    public static class ColumnData {
	public String title;
	public int width;
	public int alignment;

	public ColumnData(String t, int w, int a) {
	    title = t;
	    width = w;
	    alignment = a;
	}
    }

    /** The static list of column names, widths, and alignments. */
    public static final ColumnData columns[] = {
	new ColumnData("Index", 35, JLabel.CENTER),
	new ColumnData("Het/Mixed Base", 75, JLabel.CENTER),
	new ColumnData("Position", 50, JLabel.CENTER),
	new ColumnData("QV", 25, JLabel.CENTER),
    };

    /** The vector which holds the data. */
    private Vector vector;

    /** The index of the column which performs the sorting. */
    private int sortColumn = 0;

    /** The sorting style: true of ascending and false for descending. */
    private boolean sortAsc = true;

    /** Constructor. */
    public SNPsTableModel() {
	vector =new Vector();
    }

    /** 
     * Constructor.
     * @param	num	total number of rows in the table.
     * @param	i	the indices of all detected het/mix bases 
     * @param	b	the base codes of all detected /mix bases    
     * @param	p	the positions of all detected het/mix bases 
     * @param	q	the QVs of all detected het/mix bases 
     */
    public SNPsTableModel(int num, int[] i, char[] b, int[] p, int[] q) {
	vector =new Vector();
	setData(num, i, b, p, q);
    }

    /** Returns the total number of rows in the model. */
    public int getRowCount() {
	return (vector==null) ? 0 : vector.size();
    }

    /** Returns the total number of columns. */
    public int getColumnCount() {
	return columns.length;
    }

    /** Returns the name of the column at the specified index. */
    public String getColumnName(int column) {
	return columns[column].title;
    }

    /** Returns the value for the cell at nrow and ncolumn. */
    public Object getValueAt(int nrow, int ncolumn) {
	if (vector == null) {
	    return "";
	}
	if (nrow < 0 || nrow >= vector.size()) {
	    return "";
	}
	SNPData row = (SNPData) (vector.elementAt(nrow));
	switch (ncolumn) {
	case 0: return new Integer(row.index);
	case 1: return new Character(row.base);
	case 2: return new Integer(row.position);
	case 3: return new Integer(row.qv);
	default: return "";
	}
    }

    /** Returns true of the cell at nrow and ncolumn is editable. */
    public boolean isCellEditable(int nrow, int ncolumn) {
	return false;
    }
    
    /** 
     * Sets the table with the specified het/mix bases      
     * @param	num	total number of rows in the table.
     * @param	ind	the indices of all detected het/mix bases   
     * @param	b	the base codes of all detected het/mix bases   
     * @param	p	the positions of all detected het/mix bases   
     * @param	q	the QVs of all detected het/mix bases   
     */
    public void setData(int num, int[] ind, char[] b, int[] p, int[] q) {
	vector.removeAllElements();
	for (int i = 0; i < num; i++) {
	    vector.addElement(new SNPData(ind[i], b[i], p[i], q[i]));
	}
	if (vector.size() > 1) {
	    Collections.sort(vector, new SNPComparator(sortColumn, sortAsc));
	}
	fireTableDataChanged();
    }

    /** The inner class that listens to mouse clicks.  If mouse clicked on
	a column header, sort on that column or reverse order of sort on
	that column. */
    public class ColumnListener extends MouseAdapter {
	protected JTable table;

	public ColumnListener(JTable t) {
	    table = t;
	}
	
	public void mouseClicked(MouseEvent e) {
	    TableColumnModel colModel = table.getColumnModel();
	    int colModelIndex = colModel.getColumnIndexAtX(e.getX());
	    int modelIndex = colModel.getColumn(colModelIndex).getModelIndex();

	    if (modelIndex < 0) {
		return;
	    }
	    if (sortColumn == modelIndex) {
		sortAsc = !sortAsc;
	    } else {
		sortColumn = modelIndex;
	    }

	    Collections.sort(vector, new SNPComparator(modelIndex, sortAsc));
	    table.tableChanged(new TableModelEvent(SNPsTableModel.this));
	    table.repaint();
	}
    }

    /** An inner class that encapsulates a het/mix bases.  It contains
        information about the index, code, position, and quality value
        of the het/mix base */
    class SNPData {
	public int index, position, qv;
	public char base;

	public SNPData(int i, char b, int p, int q) {
	    index = i;
	    base = b;
	    position = p;
	    qv = q;
	}
    }

    /** The comparator for SNPData objects, uses current sort column and
	direction to determine what to compare. */
    class SNPComparator implements Comparator {
	private int sortCol;
	private boolean sortAsc;

	public SNPComparator(int sortCol, boolean sortAsc) {
	    this.sortCol = sortCol;
	    this.sortAsc = sortAsc;
	}

	public int compare(Object o1, Object o2) {
	    if (!(o1 instanceof SNPData) || !(o2 instanceof SNPData)) {
		return 0;
	    }
	    SNPData s1 = (SNPData) o1;
	    SNPData s2 = (SNPData) o2;
	    int result = 0;

	    switch (sortCol) {
	    case 0:
		result = compare(s1.index, s2.index);
		break;
	    case 1:
		result = compare((int) s1.base, (int) s2.base);
		break;
	    case 2:
		result = compare(s1.position, s2.position);
		break;
	    case 3:
		result = compare(s1.qv, s2.qv);
		break;
	    }
	    if (!sortAsc) {
		result = -result;
	    }
	    return result;
	}

	public boolean equals(Object obj) {
	    if (obj instanceof SNPComparator) {
		SNPComparator compObj = (SNPComparator) obj;
		return (compObj.sortCol == sortCol)
		       && (compObj.sortAsc == sortAsc);
	    }
	    return false;
	}

	protected int compare (int i1, int i2) {
	    if (i1 < i2) {
		return -1;
	    } else if (i1 == i2) {
		return 0;
	    } else {
		return 1;
	    }
	}
    }
}
