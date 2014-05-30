/*
 * 1.1
 *
 * @(#)LangTools.java	1.0	2000/10/24
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.util;

import java.lang.*;
import java.io.*;
import java.util.*;
import java.util.jar.*;
import java.util.zip.*;

/** This class contains convenient language tools.  */
public class LangTools {

    /**
     * Creates an upper-case version of the specified char array.
     * @return 	the upper-case version of the specified array.
     * @param	chars	the specified char array.
     */
    public static char[] toUpperCase(char[] chars) {
	if (chars == null) {
	    return null;
	}
	char[] result = new char[chars.length];
	for  (int i = 0; i < result.length; i++) {
	    result[i] = Character.toUpperCase(chars[i]);
	}
	return result;
    }

    /**
     * Gets the InputStream to read the specified jar entry from the
     * speified jar file.
     */
    public static InputStream getZipEntryInputStream(String jarFileName,
						     String jarEntryName)
				throws IOException {
        String classPath = System.getProperty("java.class.path");
        StringTokenizer st = new StringTokenizer(classPath, ":");

        while (st.hasMoreTokens()) {
              String pathElement = st.nextToken().trim();
              if (pathElement.endsWith(jarFileName)) {
                  JarFile jarFile = new JarFile(pathElement);
                  ZipEntry zipEntry = jarFile.getEntry(jarEntryName);
                  return jarFile.getInputStream(zipEntry);
              }
        }
        return null;
    }
}
