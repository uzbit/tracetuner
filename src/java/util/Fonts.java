/*
 * 1.1
 *
 * @(#)Fonts.java	1.0	2000/12/21
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.util;

import java.awt.*;

/** This class creates the fonts used by TraceTuner launcher and viewer. */
public class Fonts {

    public final static Font BASE_FONT = new Font("Courier", Font.PLAIN, 12);
    public final static Font INDEX_FONT = new Font("Arial", Font.PLAIN, 10);
    public final static Font QV_FONT = new Font("Arial", Font.PLAIN, 8);
    public final static Font LABEL_FONT1 = new Font("Dialog", Font.BOLD, 11);
    public final static Font LABEL_FONT2 = new Font("Arial", Font.PLAIN, 11);
    public final static Font TITLE_FONT = new Font("Arial", Font.BOLD, 8);
    public final static Font CONTENT_FONT = new Font("Arial", Font.PLAIN, 8);
    public final static Font MENU_FONT = new Font("Dialog", Font.BOLD, 11);
    public final static Font TEXT_FONT = new Font("Courier", Font.PLAIN, 12);

    /**
     * Gets the <code>FontMetrics</code> of the specified font.
     * @return 	the fontMetrics.
     * @param	f	the specified font.
     */
    public static FontMetrics getFontMetrics(Font f) {
	Container c = new Container();
	FontMetrics fm = c.getFontMetrics(f);
	c = null;
	return fm;
    }
}
