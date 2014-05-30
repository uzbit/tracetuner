/*  1.3
 *  Filespec.java - Utility object to accept input of a filename
 *                - MFM Oct-99
 *                - Utility object to accept input of a filename, including
 *                  some handy multiple selection functions to work around
 *                  bugs in the Java 1.2 JFileChooser.  We also hacked in a
 *                  "Select All Files" and "Clear" buttons as a 
 *                  JFileChooser accessory.  Suspect we've violated swing's
 *                  encapsulation so that this code will most likely not run
 *                  properly on future releases of java.
 *                - MFM Jan-00
 */

package com.paracel.tt.util;

import java.beans.*;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

public class Filespec extends JPanel implements ActionListener, 
    PropertyChangeListener {
    
    private JFileChooser chooser;
    public JTextField text1;
    private JTextField jt = new JTextField(); // "File Name:" JTextField
					      // object inside the JFileChooser
    private JList  jlist;  // JList object inside the JFileChooser
    private JScrollPane js; // JScrollPane object inside the JFileChooser
    private String lastentry;
    private boolean selectall=false;
    private boolean cancel = false;
    private boolean foundAMatch = false;
    private Object[] getentries;
    private File[] selectedfiles;
    private int[] index;
    private JFrame jframe;

    public Filespec(String s, int filemode, String dialogTitle) {
	File f = new 
	    File("d:/Perkin-Elmer/Abi/3700/DataExtractor/Completed Runs");
          
	if (f.exists())
	    chooser = new JFileChooser(f);
	else {
	    f = new File("c:/");
	    if (f.exists())
		chooser = new JFileChooser(f);
	    else
		chooser = new JFileChooser(System.getProperty("user.dir"));
	}

	JLabel     label1 = new JLabel(s,JLabel.CENTER);
	JButton    butt1  = new JButton("Browse");

	text1  = new JTextField(20);
	add(label1);
	add(text1); 
	add(butt1); 
	butt1.addActionListener(this);
	chooser.setFileSelectionMode(filemode);
 	setSingleSelection(chooser);
        chooser.setDialogTitle(dialogTitle);

	//  To support multiple file select, grab JList for the file chooser.
	//  To fix the File Chooser bugs from JFileChooser, grab JScrollPane
	//  for the file chooser.
        Container c1 = (Container) chooser.getComponent(3);
        while (c1 != null) {
            Container c = (Container)c1.getComponent(0);
            if (c instanceof JScrollPane)
                js = (JScrollPane)c;
            if (c instanceof JList) {
                jlist = (JList)c;
                break;
            }
            c1 = c;
        }

	chooser.addPropertyChangeListener(this);

	Container c2 = (Container)chooser.getComponent(5);
	Container c3 = (Container)c2.getComponent(3);
	Component c = c3.getComponent(1);
	if ((c != null) && (c instanceof JTextField))
	    jt = (JTextField)c;

	if (filemode == JFileChooser.DIRECTORIES_ONLY) {
	    // Change the label text from "File name:" to "Directory name:"
	    Container c4 = (Container)c2.getComponent(1);
	    c = c4.getComponent(0);
	    if ( c!= null && c instanceof JLabel) {
		((JLabel) c).setText("Directory name:");
	    }
	}

	/* Work-around for the following JFileChooser bugs: 
	   1. Consider having a list of files in the File Chooser with multiple
	   selection turned on.  First you click on one file to select it.
	   Then you scroll down the file list until the first file you selected
	   is not shown in the window, hold down Ctrl or Shift, and click on 
	   another file.  Then the file list is instantly scrolled up to 
	   show the file you first selected. 

	   2. Consider renaming a file, the first single-click on a list element
	   selects it, the second single-click on it brings up a text field in
	   order to change the name of the list element (either a directory or
	   a file).  The text field length does not change with the length of 
	   the name.

	   3. While multi-selection is enabled, hold down shift or control key
	   to select a list of files.  Then 2 consecutive single-clicks on a
	   file brings up the filename-renaming text field containing the name 
	   of the first file in the selected files list.

	   4. With multiple files selection enabled, selection failed when
	   choose item above the already selected ones using either shift or
	   control key.

	   For example: from a 10 items list, select item 5 then hold down the
	   shift key and go backward and click on item 2.  Instead of items 2
	   through 5 all get selected, only item 2 is selected.

	   For example: from a 10 items list, select a item 5 then hold down
	   the control key and go backward and select item 2.  Instead of
	   both items 2 and 5 get selected, only item 2 is selected.  

	   5. With multiple files selection enabled, deselection failed if the
	   lowest numbered selected index on the list gets deselected when
	   there are still selected items remained on the list.

	   For example: from a 10 items list, items 2 through 6 are selected.
	   Hold down control key and deselect item 3, item 2 and items 4
	   through 6 are still selected.  Again hold down control key and now
	   deselect item 2.  Instead of items 4 through 6 remained selected,
 	   only item 4 is selected. 

	   6. While renaming an item in the file list, the current directory
	   gets rescanned.  But the content of the "File Name:" text field
	   is not updated.  Additionally, if the renamed item is repositioned
           in the list (because the item list is possibly reordered when the
           current directory is rescanned), the item at the old position gets
           highlighted, instead of the renamed item. 

	   7. While renaming an item in the file list, the file chooser allows
	   renaming when the renaming could overwrite an existing item, on
	   Solaris paltform.  (Windows platform JFileChooser does not seem
	   to have this bug.)

	 */
	jlist.addMouseListener(new MouseHandler());
        jlist.getModel().addListDataListener(new ListDataHandler());

	/* end of work-around code */
    }

    /* Work-around code */
    String newName = "";  // the new name the user types in while renaming
			  //  a file or a directory in File Chooser.
    String oldName = "";  // the old name before the user types in a new name


    public class ListDataHandler implements ListDataListener {

        public void contentsChanged(ListDataEvent e) {
	    // The following code takes effect on file or directory
	    // renaming action.
            boolean matchFound = false;
	    int modelSize = jlist.getModel().getSize();
            if (modelSize > 0
                    && jlist.getSelectedIndex() >= 0
                    && newName.length() > 0 ) {
                int currIndex = jlist.getSelectedIndex();

		// To find the index of the new item in the list
                for (int i=0; i< modelSize; i++) {
                    File f = (File)(jlist.getModel().getElementAt(i));
                    if ((f.getName().trim()).equals(newName)) {
                        if (i == currIndex) {
			    // If the position of the renamed item did not
			    // change in the list, still, we need to clear
			    // the list first, then set the selection index
			    // in order to force the "File Name:" text field
			    // to update.
                            jlist.clearSelection();
                        }
                        jlist.setSelectedIndex(i);
			jlist.ensureIndexIsVisible(i);
			jt.setText(newName);
                        matchFound = true;
                        break;
                    }
                }

                if (!matchFound && jt != null) {
		    // If the renamed item does not appear in the list anymore,
		    // clear the "File Name:" text field.
                    jt.setText("");
		}
            }
        }

        public void intervalChanged(ListDataEvent e) {
        }

        public void intervalRemoved(ListDataEvent e) {
        }

        public void intervalAdded(ListDataEvent e) {
	    // new method since JDK 1.3
        }

    }//public class ListDataHandler implements ListDataListener


    public class MouseHandler extends MouseAdapter{
	// A variable used to remember the viewport position when the
	// mouse was pressed down.
	Point p;
	int lastClickedIndex = -1;
	int minimumIndex = -1;
	int anchorIndex = -1;
	int[] lastSelectedIndices;

	public void mousePressed(MouseEvent e) {
	    if (!SwingUtilities.isLeftMouseButton(e))
		return;
	    p = js.getViewport().getViewPosition();
	}

	public void mouseReleased(MouseEvent e) {
	    if (!SwingUtilities.isLeftMouseButton(e))
		return;

	    // after the mouse is released, sets the viewport back to 
	    // the value when the mouse was pressed down, because somehow
	    // the viewport was set to a wrong value when the mouse was
	    // released.
	    Point pp = js.getViewport().getViewPosition();
	    if ((!pp.equals(p)) && (e.isControlDown() || e.isShiftDown()))
	    	js.getViewport().setViewPosition(p);
	}

	public void mouseClicked(MouseEvent e) {
	    // currentIndex: the index of the file just clicked.
	    // lastClickedIndex: the index of the file clicked just before
	    //  		 the "currentIndex" file was clicked.
	
	    int newSelectedLen;
	    int currentIndex = jlist.locationToIndex(e.getPoint());	

	    Rectangle rc = jlist.getCellBounds(currentIndex, currentIndex);
	    
	    if (!e.isControlDown() && !e.isShiftDown()) {
		jlist.addSelectionInterval(currentIndex, currentIndex);
		minimumIndex = jlist.getMinSelectionIndex();
		anchorIndex = jlist.getAnchorSelectionIndex();
		lastSelectedIndices = new int[1];
		lastSelectedIndices[0] = minimumIndex;
		jlist.setSelectedIndices(lastSelectedIndices);
	    }

	    if (e.isControlDown()) {
		// Keep track of the anchorIndex if shift key is used 
		// for the next file selection. 
		if (currentIndex == minimumIndex) {
		    // Deselect an item from the list
		    newSelectedLen = lastSelectedIndices.length-1;
		    int[] newSelectedIndices = new int[newSelectedLen]; 
		    System.arraycopy(lastSelectedIndices, 1,
				     newSelectedIndices, 0, newSelectedLen);
		    jlist.setSelectedIndices(newSelectedIndices);
		    anchorIndex = currentIndex;

		} else if (currentIndex < minimumIndex) {
		    // Choose an item above the already selected one,
		    // add this item to the list. 
		    newSelectedLen = lastSelectedIndices.length+1;
		    int[] newSelectedIndices = new int[newSelectedLen];
		    newSelectedIndices[0] = currentIndex;
		    System.arraycopy(lastSelectedIndices, 0,
				     newSelectedIndices, 1, newSelectedLen-1);
		    anchorIndex = currentIndex;
		    jlist.addSelectionInterval(anchorIndex, currentIndex);
		    jlist.setSelectedIndices(newSelectedIndices);

		} else 
		    // JFileChooser already taken care of this case
		    anchorIndex = jlist.getAnchorSelectionIndex();
	    }

	    if (e.isShiftDown()) {
		if (currentIndex < anchorIndex) {
		    // Choose file(s) above the already selected one 
		    jlist.clearSelection();
		    jlist.addSelectionInterval(anchorIndex, currentIndex);
		} else {
		    // Control key was used for the previous file selection
	            jlist.addSelectionInterval(currentIndex, anchorIndex);
	            jlist.addSelectionInterval(anchorIndex, currentIndex);
		}
	    }

	    if (currentIndex == lastClickedIndex) {
		// if the last 2 clicks all happened on one file,
		// fetch the component at the position where the mouse click
		// happened.  "25" is used to skip the image at the front of
		// the list element.
        	Container c = (Container)jlist.getComponentAt(25, 
						(int)e.getPoint().getY());
        	if (c != null) {
            	    if (c instanceof JTextField) {
			// JFileChooser has brought up a JTextField for
			// list element renaming purpose.
			if (e.isControlDown() || e.isShiftDown() ) {
			    // if the control or shift key is down, remove
			    // the renaming JTextField.
                    	    jlist.remove(c);
			} else {
			    // otherwise, readjust the JTextField width to
			    // fit the entire list element name string.
		    	    String nameText = ((JTextField)c).getText();
			    if (nameText != null) {
				// To remember the old name before renaming
                                oldName = nameText.trim();
			    }
		    	    Font nameFont = c.getFont();
		    	    if (nameFont != null && nameText != null) {
		    		int width = c.getFontMetrics(
					    nameFont).stringWidth(nameText)+20;
		    		c.setSize(new Dimension(width, 
						(int)rc.getHeight()));
                                ((JTextField)c).addActionListener(new
							EditActionHandler());
		    	    }
			}
			jlist.repaint();
            	    }
        	}
	    }
	    lastClickedIndex = currentIndex;
	    lastSelectedIndices = jlist.getSelectedIndices();
	    if (!jlist.isSelectionEmpty()) 
	        minimumIndex = jlist.getMinSelectionIndex();
	}

    }//public class MouseHandler extends MouseAdapter


    public class EditActionHandler implements ActionListener {

	public void actionPerformed(ActionEvent e) {
            JTextField jtf = (JTextField)e.getSource();
	    String absoluteName = "";
	    if (jtf != null && jtf.getText() != null) {
            	newName = jtf.getText().trim();
		absoluteName = chooser.getCurrentDirectory() 
					+ File.separator + newName;
	    }

	    // For windows platform, JFileChooser is already doing the
	    // the right job, in terms of disallowing overwriting existing
	    // file (or directory) while renaming.
	    String osName = System.getProperty("os.name").trim().toLowerCase();
	    if (osName.indexOf("windows") >= 0)
		return;

	    // The user entered nothing for the new name.
	    if (newName.length() == 0) return;

	    if (!newName.equals(oldName)) {
		// The user entered a different name.  Check if an
		// item with the same name already exists. 
	    	File f = new File(absoluteName);
	    	if (f.exists()) {
		    // To avoid overwriting, show the error message and
		    // revoke the renaming action, by setting the new
		    // name back to the old name.
		    JOptionPane.showMessageDialog(jframe,
				absoluteName + " already exists.");
		   jtf.setText(oldName);
	        }
	    }
	}
    }//public class EditActionHandler implements ActionListener

    /* end of work-around code */
    

    //  Add an filename extension filter and use as default
    //  for this Filespec
    public void addFilter(Filter filt) {
	chooser.addChoosableFileFilter(filt);
        chooser.setFileFilter(filt);
    }

    
    public void propertyChange(PropertyChangeEvent pe) {
	if (pe.getPropertyName()
	    .equals(JFileChooser.FILE_FILTER_CHANGED_PROPERTY)) {

	    // Clear text in File name text box if Files of type is changed
	    jlist.clearSelection();
	    if (jt != null) 
		chooser.setSelectedFile(new File(""));

	}
    }

    
    //  This function turns multiple selection on in the JFileChooser's
    //  list and adds a select all and select none button as an accessory
    //  to help manipulate the list.  We also add a private action listern
    //  to catch the button presses in the context of this Filespec.
    public void setMultipleSelection() {
	jlist.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

	JButton b1 = new JButton("Clear");
	JButton b2 = new JButton("Select All Files");

	JPanel extra_buttons = new JPanel();
	extra_buttons.setLayout( new BoxLayout(extra_buttons,  
					       BoxLayout.Y_AXIS));
	extra_buttons.add(b1);
	extra_buttons.add(Box.createVerticalStrut(8));
	extra_buttons.add(b2);
	chooser.setAccessory(extra_buttons);

	ActionListener button_listener = new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Clear")) {
		    jlist.clearSelection();
		    chooser.setSelectedFile(new File(""));
		}
		if (e.getActionCommand().equals("Select All Files")) {
		    // Determine the index of the first JList item that
		    // represents a file since we're trying to select all
		    // files and not include directories
		    int first = 0;
		    ListModel list = jlist.getModel();
		    while(first < list.getSize() &&
			  ((File)list.getElementAt(first)).isDirectory())
			first++;

		    //  Catch case of a directory of just directories
		    if(first >= list.getSize()) {
			chooser.setSelectedFile(new File(""));
			return;
		    }

		    // Set all actual files as selected.  Need to repeat
		    // this twice to make it take the first time.
		    jlist.setSelectionInterval(first, list.getSize()-1);
		    jlist.setSelectionInterval(first, list.getSize()-1);
                    selectall=true;
		}
	    }
	};
	b1.addActionListener(button_listener);
	b2.addActionListener(button_listener);
    }

    public File[] getSelectedFiles() {

	// If cancel button is pushed, use the previously selected files.
	// The value of lastentry is retained in the text box.
	if (cancel) {
	   return selectedfiles; 
	}

        //  If the text box is empty, nothing has been selected, return
        //  a File list of length 0.

        if(text1.getText().equals(""))
            return null;

        //  If the text box matches the last jchooser selection, return the
        //  whole list of what was last selected in the jchooser.
        if (text1.getText().equals(lastentry)) {
            Object[] entries = jlist.getSelectedValues();
            File[] files = new File[entries.length];
            for (int k=0; k<entries.length; k++) {
                if (entries[k] instanceof File)
                    files[k] = (File)entries[k];
            }
            return files;
        }

        //  If the text box has different contents than we last wrote into it,
        //  the user has been mucking with it, just go with what he entered.
        File[] files = new File[1];
        files[0] = new File(text1.getText());
        return files;
    }


    //  The following three procedures are work arounds to bug in the
    //  Java 1.2 JFileChooser.  They reachs into the file chooser and
    //	find the JList and set its mode or retrieve its contents Seems
    //	to work.  This is a variation on a bug fix in the book "Swing" by
    //	Robinson and Vorobiev found at http://manning.spindoczine.com/sbe.
    //
    static public void setSingleSelection(JFileChooser chooser) {
	Container c1 = (Container) chooser.getComponent(3);
        JList list = null;
        while (c1 != null) {
            Container c = (Container)c1.getComponent(0);
            if (c instanceof JList) {
                list = (JList)c;
                break;
            }
            c1 = c;
        }
	if(list != null) {
	    list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
	    /* work-around code to fix the following JFileChooser bug:

	       Renaming a file (or a directory): the first single-click on it
               selects it, the second single-click on it brings up a text 
	       field in order to change the name of the file (or the directory).
               The text field length does not change with the length of the 
	       name.
	     */
	    list.addMouseListener(new SSMouseHandler(list)); 
	    /* end of work-around code */
	}
    }

    static public void setMultipleSelection(JFileChooser chooser) {
	chooser.setMultiSelectionEnabled(true);
	Container c1 = (Container) chooser.getComponent(3);
        JList list = null;
        while (c1 != null) {
            Container c = (Container)c1.getComponent(0);
            if (c instanceof JList) {
                list = (JList)c;
                break;
            }
            c1 = c;
        }
	if(list != null)
	    list.setSelectionMode(
		 ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    }

    static public File[] getSelectedFiles(JFileChooser chooser) {
	// Although JFileChooser won't give us this information,
	// we need it...
	Container c1 = (Container)chooser.getComponent(3);
	JList list = null;
	while (c1 != null) {
	    Container c = (Container)c1.getComponent(0);
	    if (c instanceof JList) {
		list = (JList)c;
		break;
	    }
	    c1 = c;
	}
	Object[] entries = list.getSelectedValues();
	File[] files = new File[entries.length];
	for (int k=0; k<entries.length; k++) {
	    if (entries[k] instanceof File)
		files[k] = (File)entries[k];
	}
	return files;
    }
    
    static public class SSMouseHandler extends MouseAdapter {
	JList jlist;

	SSMouseHandler(JList jl) {
	    this.jlist = jl;
	}

	public void mouseClicked(MouseEvent e) {
	    int index = jlist.locationToIndex(e.getPoint());
	    Rectangle rc = jlist.getCellBounds(index, index);
	    Container c = (Container)jlist.getComponentAt(25,
				(int)e.getPoint().getY());
	    if (c != null && c instanceof JTextField) {
		// if the file-rename text box exists
		if (e.isControlDown() || e.isShiftDown()) {
		    // if Control or Shift key is held down, remove the
		    // file-rename text box.
		    jlist.remove(c);
		} else {
		    // otherwise, adjust the text box size according to
		    // the length of the file or directory name string.
		    String nameText = ((JTextField)c).getText();
		    Font nameFont = c.getFont();
		    if (nameFont != null && nameText != null) {
			int width = c.getFontMetrics(
					nameFont).stringWidth(nameText) + 20;
			c.setSize(new Dimension(width,
					(int)rc.getHeight()));
		    }
		}
		jlist.repaint();
	    }
	}
    }

    private int getNumOfFiles(JList jl) {
	int num = 0;
	ListModel list = jl.getModel();
	for (int i = 0; i<list.getSize(); i++) {
	    if (((File)list.getElementAt(i)).isFile())
	    	num ++;
	}
	return num;
    }

    public void actionPerformed(java.awt.event.ActionEvent e) {
	if (e.getActionCommand().equals("Browse")) {
	    newName = "";    // clear the "newName" variable
            chooser.rescanCurrentDirectory();

	    // Clear text in File name text box if selection was cancelled
	    if (cancel && selectedfiles == null) {
		jlist.clearSelection();
		chooser.setSelectedFile(new File(""));
	    }
		
	    if (cancel && selectedfiles != null && selectedfiles.length >= 1) {
		String parent = selectedfiles[0].getParent();
		if (parent != null) {
		    File parentF = 
                         chooser.getFileSystemView().createFileObject(parent);
		    // Clear text in File name box if directory was switched 
		    if (!parentF.equals(chooser.getCurrentDirectory())) { 
		        chooser.setCurrentDirectory(parentF);
			jlist.clearSelection();
			chooser.setSelectedFile(new File(""));
		    }
		}
	    }
            try {
		Component c = (this.getParent() != null) ? this.getParent()
							 : this;
	        int returnVal = chooser.showDialog(c, "OK");
	        if(returnVal == JFileChooser.APPROVE_OPTION) {
		    Object[] entries = jlist.getSelectedValues();
		    int totalFileNum = getNumOfFiles(jlist);
		    if (totalFileNum == entries.length
				&& totalFileNum > 1 )
			selectall = true;

		    // Keep a copy of the selected files 
		    getentries = entries;
		    selectedfiles = new File[getentries.length];

                    javax.swing.filechooser.FileFilter f 
                   	     = chooser.getFileFilter();

                    if (entries.length > 1) { 
                        if (f instanceof Filter) {
                            // Files of type is 'ab1' 
                            if (selectall) {
                                // select all files of 'ab1'
		                lastentry=
                                    (chooser.getSelectedFile().getParent() +
			            File.separatorChar + 
                                    ((Filter)chooser.getFileFilter())
                                    .getFilterspec());
                            } else {
                            // select some files of 'ab1' from a directory
                                String fileext = 
                                ((Filter)chooser.getFileFilter()).
                                getFilterspec();
                                fileext = fileext.substring(1,5);
                                lastentry=
                                    (chooser.getSelectedFile().getParent() +
                                    File.separatorChar + "<LIST>" + fileext);
                            }
                        }
                        else {
                            // Files of type is 'All Files' 
                            lastentry=(chooser.getSelectedFile().getParent() +
                            File.separatorChar + "<LIST>");
                        }
		    } else {
		        lastentry=
                            (chooser.getSelectedFile().getAbsolutePath());
                    }
                    text1.setText(lastentry);
                    selectall=false;
		    cancel = false;
	        }
		if (returnVal == JFileChooser.CANCEL_OPTION) {
		    cancel = true;
		    jlist.clearSelection();

		    if (getentries != null) {
		        // Retrieved a copy of the previously saved files
                        for (int k = 0; k < getentries.length; k++) {
                            if (getentries[k] instanceof File) 
                                selectedfiles[k] = (File)getentries[k];
                        }
			
			/* JDK1.3 still has not implemented multiSelection
			 * properly, more work around is required in the
			 * case of setSelectedFiles(). 
			 * Use the indices from jlist and try to match with
			 * what was previously selected.
			 */
		        index = new int[selectedfiles.length];
		        for (int i = 0, j=0; i < jlist.getModel().getSize(); i++) {
			    for (int k = 0; k < selectedfiles.length; k++) {
			        File f = (File)jlist.getModel().getElementAt(i);
			        if (selectedfiles[k].equals(f)) {
			            index[j++] = i;
				    // Need this if switch directory
				    foundAMatch = true;
				}
			    }
		        }
			if (foundAMatch) {
		            if (selectedfiles.length > 1) {
		                chooser.setSelectedFiles(selectedfiles);
			        jlist.setSelectedIndices(index);
		            } else {
			        chooser.setSelectedFile(selectedfiles[0]);
		                jlist.setSelectedIndices(index); 
		            }
			}
		        foundAMatch = false;
		    }
		}
            }
            catch(ArrayIndexOutOfBoundsException excp) {
		JOptionPane.showMessageDialog(jframe, 
                "Error occurred in the number of files selected\n"+
                "Please redo the selection");
                jlist.clearSelection();
                return;
            }
	}
    }


    public static void main(String s[]) {   // Used for debugging only
	JFrame   frame = new JFrame("Filespec debug");
	final Filespec panel = new Filespec("Filespec Test", 
				      JFileChooser.FILES_ONLY,
				      "Select Sample files(s)");
	panel.setMultipleSelection();
	panel.addFilter(new Filter(".ab1"));

	frame.addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e) {
		    int i;
                    System.out.println("lastentry "+panel.lastentry);
		    File[] files = panel.getSelectedFiles();
		    if (files != null && files.length >= 1) {
                        System.out.println("number of files selected: " +
                                            files.length);
		        for(i=0; i < files.length; i++) {
			    System.out.println("Selected: " + files[i]);
		        }
		    }
                    System.exit(0);}
	    });
	 frame.setSize(500,80);
	 frame.getContentPane().add(panel);
	 frame.setVisible(true);
    }
}
