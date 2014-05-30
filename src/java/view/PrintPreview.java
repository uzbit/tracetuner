/*
 * 1.4
 *
 * @(#)PrintPreview.java       1.0     2000/11/22
 */

package com.paracel.tt.view;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.util.*;
import java.awt.print.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

/** This class creates a Print Preview window. */
public class PrintPreview extends JFrame {
        protected int m_wPage;
        protected int m_hPage;
        protected Printable m_target;
        protected JComboBox m_cbScale;
        protected PreviewContainer m_preview;
	protected PageFormat m_pageFormat;

	// Maximum scale value. Set to 100% to avoid OutOfMemory Error.
	private static int MAX_SCALE = 100;

        public PrintPreview(Printable target, String title) {
            super(title);
            setSize(350, 300);
            m_target = target;
	    m_pageFormat = ((TTPrintJob)m_target).getPageFormat(0);
	    ((TTPrintJob)m_target).setParent(this);
            JToolBar tb = new JToolBar();
            JButton bt = new JButton("Print");

            ActionListener lst = new ActionListener() { 
                public void actionPerformed(ActionEvent e) { 
		    Thread runner = new Thread() {
			public void run() {
			    ((TTPrintJob)m_target).print();
		    	    PrintPreview.this.setVisible(false);
		    	    PrintPreview.this.releaseResources();
		    	    PrintPreview.this.dispose();
			    System.gc();
			}
	    	    };
		    runner.start();
                }
            };
        
            bt.addActionListener(lst);
            bt.setAlignmentY(0.5f);
            bt.setMargin(new Insets(2,6,2,6));
            tb.add(bt);
            bt = new JButton("Close");
        
            lst = new ActionListener() { 
                public void actionPerformed(ActionEvent e) { 
		    releaseResources();
		    dispose();
		    System.gc();
                }
            };
        
            bt.addActionListener(lst);
            bt.setAlignmentY(0.5f);
            bt.setMargin(new Insets(2,6,2,6));
            tb.add(bt);
        
            String[] scales = { "10 %", "25 %", "50 %", "100 %" };
        
            m_cbScale = new JComboBox(scales);
            lst = new ActionListener() { 
                public void actionPerformed(ActionEvent e) { 
                    Thread runner = new Thread() {
                        public void run() {
                            String str = m_cbScale.getSelectedItem()
			    			.toString();
                            if (str.endsWith("%"))
                                str = str.substring(0, str.length()-1);
                            str = str.trim();
                            int scale = 0;
                            try { 
                                scale = Integer.parseInt(str); 
                            } catch (NumberFormatException ex) { return; }

			    if (scale > MAX_SCALE)
				scale = MAX_SCALE;
    
                            int w = (int)(m_wPage*scale/100);
                            int h = (int)(m_hPage*scale/100);
                            Component[] comps = m_preview.getComponents();
                            for (int k=0; k<comps.length; k++) {
                                if (!(comps[k] instanceof PagePreview))
                                    continue;
                                PagePreview pp = (PagePreview)comps[k];
                                pp.setScaledSize(w, h);
                            }
            
                            m_preview.doLayout();
                            m_preview.getParent().getParent().validate();
                        }
                    };
                    runner.start();
                }
            };

            m_cbScale.addActionListener(lst);
            m_cbScale.setMaximumSize(m_cbScale.getPreferredSize());
            m_cbScale.setEditable(true);

            tb.addSeparator();
            tb.add(m_cbScale);

            getContentPane().add(tb, BorderLayout.NORTH);
            m_preview = new PreviewContainer();
            if (m_pageFormat.getHeight()==0 || m_pageFormat.getWidth()==0) {
                System.err.println("Unable to determine default page size");
                return;
            }

            m_wPage = (int)(m_pageFormat.getWidth());
            m_hPage = (int)(m_pageFormat.getHeight());
            int scale = 10;
            int w = (int)(m_wPage*scale/100);
            int h = (int)(m_hPage*scale/100);
            int pageIndex = 0;
	    ((TTPrintJob)m_target).getNumberOfPages();
            try {
                while (true) {
                    BufferedImage img = new BufferedImage(m_wPage, m_hPage, 
						BufferedImage.TYPE_INT_RGB);
                    Graphics2D g = (Graphics2D)(img.getGraphics());
                    g.setColor(Color.white);
                    g.fillRect(0, 0, m_wPage, m_hPage);
                    if (m_target.print(g, m_pageFormat, pageIndex) 
			    != Printable.PAGE_EXISTS)
                        break;
                    PagePreview pp = new PagePreview(w, h, img);
                    m_preview.add(pp);
                    pageIndex++;
                }
            } catch (PrinterException e) {
                e.printStackTrace();
                System.err.println("Printing error: "+e.toString());
            }
    
            JScrollPane ps = new JScrollPane(m_preview);
            getContentPane().add(ps, BorderLayout.CENTER);
            setDefaultCloseOperation(DISPOSE_ON_CLOSE);
            setVisible(true);
        }

        /** Attempts to release resources. */
        public void releaseResources() {
	    if (m_preview == null) {
		return;
	    }
            Component[] comps = m_preview.getComponents();
            for (int k=0; k<comps.length; k++) {
                if (!(comps[k] instanceof PagePreview)) {
	    	    comps[k] = null;
                    continue;
		}
                PagePreview pp = (PagePreview)comps[k];
		pp.releaseResources();
		pp = null;
            }
	    m_preview = null;
	    ((TTPrintJob)m_target).dispose();
	    m_target = null;
	}

        class PreviewContainer extends JPanel {
                protected int H_GAP = 16;
                protected int V_GAP = 10;

                public Dimension getPreferredSize() {

                    int n = getComponentCount();
                    if (n == 0)
                        return new Dimension(H_GAP, V_GAP);
    
                    Component comp = getComponent(0);
                    Dimension dc = comp.getPreferredSize();
                    int w = dc.width;
                    int h = dc.height;

                    Dimension dp = getParent().getSize();
                    int nCol = Math.max((dp.width-H_GAP)/(w+H_GAP), 1);
                    int nRow = n/nCol;
                    if (nRow*nCol < n)
                        nRow++;

                    int ww = nCol*(w+H_GAP) + H_GAP;
                    int hh = nRow*(h+V_GAP) + V_GAP;
                    Insets ins = getInsets();

                    return new Dimension(ww+ins.left+ins.right, hh+ins.top+ins.bottom);
                }
                
                public Dimension getMaximumSize() {
                    return getPreferredSize();
                }
    
                public Dimension getMinimumSize() {
                    return getPreferredSize();
                }
            
                public void doLayout() {
                    Insets ins = getInsets();
                    int x = ins.left + H_GAP;
                    int y = ins.top + V_GAP;
                    int n = getComponentCount();
                    if (n == 0)
                        return;

                    Component comp = getComponent(0);
                    Dimension dc = comp.getPreferredSize();
                    int w = dc.width;
                    int h = dc.height;

                    Dimension dp = getParent().getSize();
                    int nCol = Math.max((dp.width-H_GAP)/(w+H_GAP), 1);
                    int nRow = n/nCol;
                    if (nRow*nCol < n)
                        nRow++;

                    int index = 0;
                    for (int k = 0; k<nRow; k++) {
                        for (int m = 0; m<nCol; m++) {
                            if (index >= n)
                                return;
                            comp = getComponent(index++);
                            comp.setBounds(x, y, w, h);
                            x += w+H_GAP;
                        }
                        y += h+V_GAP;
                        x = ins.left + H_GAP;
                    }
                }
        }

        class PagePreview extends JPanel {
                protected int m_w;
                protected int m_h;
                protected Image m_source;
                protected Image m_img;

                public PagePreview(int w, int h, Image source) {
                    m_w = w;
                    m_h = h;
                    m_source= source;
                    m_img = m_source.getScaledInstance(m_w, m_h, 
						       Image.SCALE_SMOOTH);
                    m_img.flush();
                    setBackground(Color.white);
                    setBorder(new MatteBorder(1, 1, 2, 2, Color.black));
                }

                public void setScaledSize(int w, int h) {
                    m_w = w;
                    m_h = h;
                    m_img = m_source.getScaledInstance(m_w, m_h, 
						       Image.SCALE_SMOOTH);
                    repaint();
                }

                public Dimension getPreferredSize() {
                    Insets ins = getInsets();
                    return new Dimension(m_w+ins.left+ins.right, 
					 m_h+ins.top+ins.bottom);
                }

                public Dimension getMaximumSize() {
                    return getPreferredSize();
                }
        
                public Dimension getMinimumSize() {
                    return getPreferredSize();
                }
                
                public void paint(Graphics g) {
                    g.setColor(getBackground());
                    g.fillRect(0, 0, getWidth(), getHeight());
                    g.drawImage(m_img, 0, 0, this);
                    paintBorder(g);
                }
	    
	        public void releaseResources() {
		    m_source = null;
		    m_img = null;
		}
        }
}
