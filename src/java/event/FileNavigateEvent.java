/*
 * 1.1
 *
 * @(#)FileNavigateEvent.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.event;

import java.awt.*;
import java.awt.event.*;
import com.paracel.tt.util.*;

/**
 * This class creates a FileNavigateEvent object. Such an event is
 * fired when the user navigate to another file in the file list.
 */
public class FileNavigateEvent extends AWTEvent {

    /** The <code>FilesStatus</code> of the targeted sample file. */
    private FilesStatus filesStatus;

    /** The total number of files in the file list. */
    private int numOfFiles;

    /**
     * Creates a FileNavigateEvent object.
     * @param	o	the object where the event originated
     * @param	fs	the files status. 
     * @param	n	the total number of files in the sample file list. 
     */
    public FileNavigateEvent(Object o, FilesStatus fs, int n) {
	super(o, AWTEvent.RESERVED_ID_MAX + 3);
	filesStatus = fs;
	numOfFiles = n;
    }

    /**
     * Gets the FilesStatus of the targeted sample file.
     * @return the FilesStatus of the targeted sample file.
     */
    public FilesStatus getFilesStatus() { return filesStatus; }

    /**
     * Gets the total number of files in the sample files list.
     * @return the total number of files in the sample files list.
     */
    public int getNumOfFiles() { return numOfFiles; }
}
