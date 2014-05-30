/*  1.1
 *  Filter.java - Simple JFileChooser FileFilter 
 *              - MFM Jan-00
 */

package com.paracel.tt.util;

import java.io.*;
import javax.swing.filechooser.FileFilter;

public class Filter extends FileFilter 
{
    String extension, description, filterspec;

    public Filter(String extension) {
	if(extension.charAt(0) == '.')
	    this.extension = extension;
	else
	    this.extension = "." + extension;
	//  Provide simple default description
	this.description = "(*" + this.extension + ")";
        this.filterspec = "*" + this.extension;
    }

    public Filter(String extension, String description) {
	this(extension);
	this.description = description;  // Over-ride default description
    }

    public boolean accept(File f) {
	File sf;

	if(f.isDirectory()) 
	    return true;

	if(f.getPath().endsWith(extension) &&
	   f.canRead())
	    return true;
	return false;
    }

     public String getDescription() {
	 return description;
     }

     public String getFilterspec() {
	 return filterspec;
     }

}
