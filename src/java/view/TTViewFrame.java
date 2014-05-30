/*
 * 1.34
 *
 * @(#)TTViewFrame.java       1.0     2000/10/25
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.view;

import java.io.*;
import java.lang.*;
import java.util.*;
import java.net.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.print.*;
import javax.swing.*;
import javax.swing.border.*;
import com.paracel.tt.chrom.*;
import com.paracel.tt.event.*;
import com.paracel.tt.util.*;
import com.paracel.tt.io.*;

/**
 * This class creates the main frame of TTView window, which contains the
 * menus, buttons, etc.
 */
public class TTViewFrame extends JFrame 
    			 implements ActionListener,
				    SearchEventListener {

    final static Font menuFont = new Font("Dialog", Font.BOLD, 11);
    final static Font labelFont = new Font("Arial", Font.PLAIN, 11);

    /** The stale time.  If .tab, .tal, or .tip file is 1 minute older
	than the .phd.1 file, it is considered as stale. */
    final static int STALE_TIME = 60000;

    /** The chromatogram file names. */
    private String[] chromatFileNames;

    /** The phd file names. */
    private String[] phdFileNames;

    /** The total number of chromat files. */
    private int numOfFiles;

    /** The file index. */
    private int fileIndex;

    /** The TTView panel resides in this frame. */
    private ChromViewPanel chromViewPanel;

    /** The message panel. */
    private MessagePanel msgPanel;

    /** The menu bar of this frame. */
    private JMenuBar menuBar;

    /** The tool bar of this frame. */
    private JToolBar toolBar;

    /** The flag indicating whether or not the application will exit when
        closing this viewer. */
    private boolean exitOnClose;

    /** The FilesStatus of the list of files and the related auxillary files
	(such as .tal file and .tip file) to be viewed. */
    private FilesStatus filesStatus;

    /** The ZoomEvent listeners. */
    private Vector zoomEventListeners; 

    /** The DisplayFeatureEvent listeners. */
    private Vector displayFeatureEventListeners;

    /** The DisplaySignalValueEvent listeners. */
    private Vector displaySigValEventListeners;

    /** The FileNavigateEvent listeners. */
    private Vector fileNavEventListeners;

    /** The SearchEvent listeners. */
    private Vector searchEventListeners;

    /** The menus. */
    private JMenu fileMenu, viewMenu, featureMenu, helpMenu,
		  zoomXMenu, zoomYMenu, printPreviewMenu, advMenu;

    /** The menu items. */
    private JMenuItem openItem, gotoItem, quitItem, infoItem,
    		      pageSetupItem, printViewItem, printAllItem,
    		      previewPAItem, previewPVItem, refreshItem, 
		      findItem, findAgainItem, printAllFilesItem;

    /** The menu items. */
    private JCheckBoxMenuItem origItem, rawDataItem, signalItem,
			      hetItem, trimItem, tipItem, 
                              totalISItem, alnItem, abcItem,
			      selectedZoomXItem, selectedZoomYItem;

    /** The buttons. */
    private ButtonGroup btnGroup1, btnGroup2;
    private JButton openBtn, pgSetupBtn, prnViewBtn, refreshBtn,
		    firstBtn, backBtn, forwardBtn, lastBtn, infoBtn,
		    findBtn, findAgainBtn;

    /** The trim threshold label. */
    private JLabel trimLabel;
    
    /** The trim threshold text field. */
    private JTextField trimField;

    /** The trim threshold. */
    private int trimThreshold;
    
    /** The zoom factors. */
    private float zoomX, oldZoomX, oldZoomY;

    /** Whether "Find Again" option/button should be enabled. */
    private boolean enableFindAgain;

    /** Whether "Align with consensus" option should be enabled. */
    private boolean enableDisplayAln;

    /** Whether "Display Alternative Base Calls" option should be enabled. */
    private boolean enableDisplayAbc;

    /** Whether "Display Intrinisic Peaks" option should be enabled. */
    private boolean enableDisplayTip;

    /** Whether "Display Trimmed Regions" option should be enabled. */
    private boolean enableDisplayTrim;

    /** A vector to hold the references to all the opened preview windows. */
    private Vector previewWindows;

    /** The file navigation dialog. */
    private GoToDialog gotoDialog;

    /** The custom zoom dialogs. */
    private CustomZoomDialog zoomXDialog, zoomYDialog;

    /** The search dialog. */
    private SearchDialog searchDialog;

    /** A static variable to count the total number of TTView frames. */
    private static int frameCount;

    /** Whether the application is running on Windows OSes. */
    private static boolean isWindows;

    private int numPrintAllFilesJobs;

    /**
     * Creates a TTViewFrame displaying chromatograms from the specified
     * chromat files and information from the specified phd file and other
     * auxillary files.
     * @param	cFileNames   A String array of the names of all chromat files.
     * @param   pFileNames   A String array of the names of the phd files.
     * @param   secondPhdFileName   The second phd file name.
     * @param	exitOnCls	Whether the app exits on closing the viewer.
     */
    public TTViewFrame(String[] cFileNames,
		       String[] pFileNames,
		       String   secondPhdFileName,
		       boolean exitOnCls) {
	super("Paracel TraceTuner " + Constants.VERSION 
	      	+ " Viewer  [ " + (++frameCount) + " ]");
	setSize(new Dimension(800, 600));

	numOfFiles = cFileNames.length;
	if (numOfFiles <= 0) {
	    dispose();
	    return;
	}

	this.exitOnClose = exitOnCls;

	String osName = Constants.OS_NAME.toLowerCase().trim();
	isWindows = (osName.indexOf("windows") >= 0) ? true : false;
	previewWindows = new Vector();
	zoomEventListeners = new Vector(1);
	displayFeatureEventListeners = new Vector(1);
	displaySigValEventListeners = new Vector(1);
	fileNavEventListeners = new Vector(1);
	searchEventListeners = new Vector(1);
	
	zoomX = 1.0f;
	oldZoomX = 1.0f;
	oldZoomY = 1.0f;

        chromatFileNames = new String[numOfFiles];
        phdFileNames = new String[numOfFiles];
	System.arraycopy(cFileNames, 0, chromatFileNames, 0, numOfFiles);
	System.arraycopy(pFileNames, 0, phdFileNames, 0, numOfFiles);

	// sps: second phd file status
	int sps = FilesStatus.UNKNOWN;
	if (secondPhdFileName != null && secondPhdFileName.length() == 0) {
	    secondPhdFileName = null;
	}

	if (secondPhdFileName != null) {
	    File f = new File(secondPhdFileName);
	    if (!f.exists()) {
		sps = FilesStatus.DOES_NOT_EXIST;
	    } else if (!f.canRead()) {
		sps = FilesStatus.CAN_NOT_READ;
	    } else {
		sps = FilesStatus.VALID;
	    } 
	}
	
	createMenuBar();
	setJMenuBar(menuBar);

	createToolBar();
	getContentPane().add(toolBar, BorderLayout.NORTH);

	chromViewPanel = new ChromViewPanel(this, secondPhdFileName, sps);
	chromViewPanel.setBorder(
		BorderFactory.createBevelBorder(BevelBorder.LOWERED));
	zoomEventListeners.add(chromViewPanel);
	displayFeatureEventListeners.add(chromViewPanel);
	fileNavEventListeners.add(chromViewPanel);
	searchEventListeners.add(chromViewPanel);

	msgPanel = new MessagePanel(this, secondPhdFileName, sps);
	msgPanel.addSearchEventListener(chromViewPanel);
	fileNavEventListeners.add(msgPanel);
	chromViewPanel.addDisplaySignalValueEventListener(msgPanel);

	JSplitPane splitPanel = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
					  chromViewPanel.getWrapperComponent(),
					  msgPanel);
	splitPanel.setDividerSize(2);
	splitPanel.setResizeWeight(0.9);
	getContentPane().add(splitPanel);

	ToolTipManager.sharedInstance().setInitialDelay(100);
	pack();
	centerWindow();

	setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
	addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {
		if (numPrintAllFilesJobs > 0) {
		    String message = "Closing this window will terminate"
			+ " the ongoing Print-Global-View-of-All-Files job.\n"
			+ "Do you want to close it?";
		    int result = JOptionPane.showConfirmDialog(
						TTViewFrame.this, message,
						"Close Confirmation",
						JOptionPane.YES_NO_OPTION,
						JOptionPane.QUESTION_MESSAGE);
		    if (result == JOptionPane.NO_OPTION) {
			return;
		    }
		}
		if (!exitOnClose) {
		    TTViewFrame.this.dispose();
		    TTViewFrame.this.releaseResources();
		} else {
		    System.exit(1);
		}
	    }
	});

	setVisible(true);

	fileIndexChanged();
    }

    /** Creates the menu bar. */
    protected void createMenuBar() {

	menuBar = new JMenuBar();

	menuBar.setFont(menuFont);
	menuBar.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));

	createFileMenu();
	createViewMenu();
	createFeatureMenu();
	createHelpMenu();

	menuBar.add(fileMenu);
	menuBar.add(viewMenu);
	menuBar.add(featureMenu);
	menuBar.add(helpMenu);
    }
		       
    /** Creates the file menu. */
    protected void createFileMenu() {
	openItem = makeMenuItem("Open...");
	gotoItem = makeMenuItem("Go to...");
	if (numOfFiles <= 1) {
	    gotoItem.setEnabled(false);
	}

	pageSetupItem = makeMenuItem("Page Setup...");
	printPreviewMenu = makeMenu("Print Preview");
	previewPVItem = makeMenuItem("Print Local View Preview");
	previewPAItem = makeMenuItem("Print Global View of Current File Preview");

	printPreviewMenu.add(previewPVItem);
	printPreviewMenu.add(previewPAItem);

	printViewItem = makeMenuItem("Print Local View...");
	printAllItem = makeMenuItem("Print Global View of Current File...");
	printAllFilesItem = makeMenuItem("Print Global View of All Files...");
	if (numOfFiles <= 1) {
	    printAllFilesItem.setEnabled(false);
	}

	quitItem = makeMenuItem("Quit");

	fileMenu = makeMenu("File");

	fileMenu.add(openItem);
	fileMenu.addSeparator();
	fileMenu.add(gotoItem);
	fileMenu.addSeparator();
	fileMenu.add(pageSetupItem);
	fileMenu.add(printPreviewMenu);
	fileMenu.add(printViewItem);
	fileMenu.add(printAllItem);
	fileMenu.add(printAllFilesItem);
	fileMenu.addSeparator();
	fileMenu.add(quitItem);
    }
		       
    /** Creates the View menu. */
    protected void createViewMenu() {
	btnGroup1 = new ButtonGroup();
	btnGroup2 = new ButtonGroup();

	zoomXMenu = makeMenu("Zoom X");

	createZoomItems("default", 'X');
	createZoomItems("6.25", 'X');
	createZoomItems("300", 'X');
	createZoomItems("600", 'X');
	createZoomItems("1000", 'X');
	createZoomItems("Other...", 'X');

	zoomYMenu = makeMenu("Zoom Y");

	createZoomItems("default", 'Y');
	createZoomItems("200", 'Y');
	createZoomItems("400", 'Y');
	createZoomItems("600", 'Y');
	createZoomItems("800", 'Y');
	createZoomItems("Other...", 'Y');

	signalItem = makeCheckBoxMenuItem("Show signal value");
	hetItem = makeCheckBoxMenuItem("List detected hets/mixed bases");
	refreshItem = makeMenuItem("Refresh");
	findItem = makeMenuItem("Find...");
	findAgainItem = makeMenuItem("Find again");
	findAgainItem.setEnabled(false);

	viewMenu = makeMenu("View");
	viewMenu.add(zoomXMenu);
	viewMenu.add(zoomYMenu);
	viewMenu.addSeparator();
	viewMenu.add(findItem);
	viewMenu.add(findAgainItem);
	viewMenu.add(refreshItem);
	viewMenu.addSeparator();
	viewMenu.add(signalItem);
	viewMenu.add(hetItem);
    }

    /** Creates the Feature menu. */
    protected void createFeatureMenu() {
	origItem = makeCheckBoxMenuItem("Display original/ABI base calls");

	rawDataItem = makeCheckBoxMenuItem("Display raw data");

	alnItem = makeCheckBoxMenuItem(
				"Align with consensus/reference sequence");

	abcItem = makeCheckBoxMenuItem("Display alternative base calls");

	tipItem = makeCheckBoxMenuItem("Display intrinsic peaks");

	totalISItem = makeCheckBoxMenuItem("Display total intrinsic signal");

	trimItem = makeCheckBoxMenuItem("Display trimmed regions");

	advMenu = makeMenu("Advanced");

	featureMenu = makeMenu("Feature");
	
	advMenu.add(alnItem);
	advMenu.add(abcItem);
        advMenu.add(tipItem);
	if (Constants.runAsDev) {
	    advMenu.add(totalISItem);
	}
	advMenu.add(trimItem);

	featureMenu.add(origItem);
	featureMenu.add(rawDataItem);
	featureMenu.addSeparator();
	featureMenu.add(advMenu);
    }

    /** Creates the Help menu. */
    protected void createHelpMenu() {
	infoItem = makeMenuItem(new HelpAction("Info"));
	helpMenu = makeMenu("Help");
	helpMenu.add(infoItem);
    }

    protected JMenu makeMenu(String text) {
	JMenu m = new JMenu(text);
	if (!isWindows) {
	    m.setFont(menuFont);
	}
	return m;
    }

    protected JMenuItem makeMenuItem(String text) {
	JMenuItem m = new JMenuItem(text);
	if (!isWindows) {
	    m.setFont(menuFont);
	}
	m.addActionListener(this);
	return m;
    }

    public JMenuItem makeMenuItem(Action action) {
        JMenuItem m = new JMenuItem();
	if (!isWindows) {
	    m.setFont(menuFont);
	}
        m.setAction(action);
        return m;
    }

    protected JCheckBoxMenuItem makeCheckBoxMenuItem(String text) {
	JCheckBoxMenuItem m = new JCheckBoxMenuItem(text);
	if (!isWindows) {
	    m.setFont(menuFont);
	}
	m.addActionListener(this);
	return m;
    }

    /** 
     * Creates the zoom menu items.
     * @param	value	the zoom value
     * @param	type	zoom type (or axis).
     */
    protected void createZoomItems(String value, char type) {
	float dbValue = 1.0f;
	JCheckBoxMenuItem mi;
	boolean isCustom;
	if (value.equals("Other...")) {
	    mi = new JCheckBoxMenuItem("Other...");
	    isCustom = true;
	    dbValue = 0.0f;
	} else if (value.equals("default")){
	    mi = new JCheckBoxMenuItem("default");
	    isCustom = false;
	    mi.setState(true);
	    if (type == 'X') {
		selectedZoomXItem = mi;
	    } else {
		selectedZoomYItem = mi;
	    }
	    dbValue = 1.0f;
	} else {
	    mi = new JCheckBoxMenuItem(value + "%");
	    isCustom = false;
	    if (value.equals("6.25")) {
	    	dbValue = 0.0625f;
	    } else {
	    	int intValue = Integer.parseInt(value);
	    	dbValue = ((float) intValue)/100.0f;
	    }
	}
	if (!isWindows) {
	    mi.setFont(menuFont);
	}
	if (type == 'X') {
	    mi.addActionListener(new ActionHandler(this, 'X', 
						   dbValue, isCustom));
	    zoomXMenu.add(mi);
	    btnGroup1.add(mi);
	} else {
	    mi.addActionListener(new ActionHandler(this, 'Y',
						   dbValue, isCustom));
	    zoomYMenu.add(mi);
	    btnGroup2.add(mi);
	}
    }

    /** inner Action Handler class for the zoom menu items. */
    public class ActionHandler implements ActionListener {
	JFrame f;
	char type;
	float oldValue, newValue;
	boolean isCustom;

	public ActionHandler(JFrame f, char type, 
			     float newValue, boolean isCustom) {
	    this.f = f;
	    this.type = type;
	    this.newValue = newValue;
	    this.isCustom = isCustom;
	}

        public void actionPerformed(ActionEvent e) {
	    float oldValue = (type == 'X') ? oldZoomX : oldZoomY;
	    
	    if (!isCustom) {
		if (type == 'X') {
		    selectedZoomXItem = (JCheckBoxMenuItem) e.getSource();
		} else {
		    selectedZoomYItem = (JCheckBoxMenuItem) e.getSource();
		}
	    } else {
		boolean isCanceled;
		if (type == 'X') {
		    if (zoomXDialog == null) {
	            	zoomXDialog = new CustomZoomDialog(f, type, 5, 5000);
		    }
		    zoomXDialog.becomeVisible();
	            isCanceled = zoomXDialog.isCanceled();
		} else {
		    if (zoomYDialog == null) {
	            	zoomYDialog = new CustomZoomDialog(f, type, 5, 5000);
		    }
		    zoomYDialog.becomeVisible();
	            isCanceled = zoomYDialog.isCanceled();
		}
	        if (isCanceled) {
		    newValue = oldValue;
		    if (type == 'X') {
			selectedZoomXItem.setSelected(true);
		    } else {
			selectedZoomYItem.setSelected(true);
		    }
	        } else {
		    if (type == 'X') {
	                newValue = zoomXDialog.getZoomValue();
		        selectedZoomXItem = (JCheckBoxMenuItem) e.getSource();
		    } else {
	                newValue = zoomYDialog.getZoomValue();
		        selectedZoomYItem = (JCheckBoxMenuItem) e.getSource();
		    }
		}
	    }

	    if (oldValue == newValue) {
		return;
	    }

	    if (type == 'X') {
	        oldZoomX = newValue;
		zoomX = newValue;
		int tipStatus = filesStatus.getTipFileStatus();
		if (oldValue <= 0.5 && newValue > 0.5
		    	&& (tipStatus == FilesStatus.VALID
			    || tipStatus == FilesStatus.VALID_BUT_STALE)) {
		    enableDisplayTip = true;
		    if ((!rawDataItem.isEnabled()) 
			    || (!rawDataItem.getState())) {
	            	tipItem.setEnabled(true);
	            	totalISItem.setEnabled(true);
	    	    	fireDisplayFeatureEvent(false);
		    }
		}
		if (oldValue > 0.5 && newValue <= 0.5 
			&& enableDisplayTip) {
		    enableDisplayTip = false;
	            tipItem.setEnabled(false);
	            totalISItem.setEnabled(false);
	    	    fireDisplayFeatureEvent(false);
	        }
	    } else {
	        oldZoomY = newValue;
	    }

	    ZoomEvent ze = new ZoomEvent(f, type, oldValue, newValue);
            for (int i = 0; i < zoomEventListeners.size(); i++) {
                ((ZoomEventListener)(zoomEventListeners.elementAt(i)))
                				.zoomValueChanged(ze);
            }
        }
    }

    /** Creates the tool bar. */
    protected void createToolBar() {
	toolBar = new JToolBar();

	toolBar.setMargin(new Insets(0,0,0,0));
	toolBar.setSize(new Dimension(this.getSize().width, 18));
	toolBar.setFloatable(false);

	Dimension d = new Dimension(25, toolBar.getSize().height);
	Dimension d1 = new Dimension(4, toolBar.getSize().height);
	Dimension d2 = new Dimension(4, toolBar.getSize().height);

	toolBar.addSeparator(d1);
	openBtn = makeToolButton("Open16.gif", "Open File", "Open");
	toolBar.add(openBtn);
	toolBar.addSeparator(d);

	pgSetupBtn = makeToolButton("PageSetup16.gif", "Page Setup", "PageSetup");
	toolBar.add(pgSetupBtn);
	toolBar.addSeparator(d1);
	prnViewBtn = makeToolButton("Print16.gif", "Print Local View", "PrintLV");
	toolBar.add(prnViewBtn);
	toolBar.addSeparator(d);

	findBtn = makeToolButton("Find16.gif", "Find", "Find");
	toolBar.add(findBtn);
	toolBar.addSeparator(d1);
	findAgainBtn = makeToolButton("FindAgain16.gif", "Find Again", "Find Again");
	findAgainBtn.setEnabled(false);
	toolBar.add(findAgainBtn);
	toolBar.addSeparator(d);

	refreshBtn = makeToolButton("Refresh16.gif", "Refresh Screen", "Refresh");
	toolBar.add(refreshBtn);
	toolBar.addSeparator(d);

	firstBtn = makeToolButton("StepBack16.gif", "First Sample File", "First");
	toolBar.add(firstBtn);
	toolBar.addSeparator(d1);
	backBtn = makeToolButton("back.gif", "Previous Sample File", "Previous");
	toolBar.add(backBtn);
	toolBar.addSeparator(d1);
	forwardBtn = makeToolButton("forward.gif", "Next Sample File", "Next");
	toolBar.add(forwardBtn);
	toolBar.addSeparator(d1);
	lastBtn = makeToolButton("StepForward16.gif", "Last Sample File","Last");
	toolBar.add(lastBtn);
	toolBar.addSeparator(d);

	infoBtn = makeToolButton("Help16.gif", "Help Info", "Help");
	toolBar.add(infoBtn);
	toolBar.add(Box.createHorizontalGlue());

	trimLabel = new JLabel("Trim Threshold");
	trimLabel.setFont(labelFont);
	trimLabel.setForeground(Color.black);
	toolBar.add(trimLabel);
	toolBar.addSeparator(d2);

	trimField = new JTextField(2);
	trimField.setMaximumSize(new Dimension(8, 18));
	trimField.setEditable(false);
	trimField.setHorizontalAlignment(SwingConstants.CENTER);
	trimField.setBackground(Color.white);
	toolBar.add(trimField);
	toolBar.addSeparator(d1);
    }

    /**
     * Creates the tool button.
     * @param	imgSrc	image file name.
     * @param	toolTipText	the tool tip text.
     */
    protected JButton makeToolButton(String imgSrc, String toolTipText,
				     String text) {
	JButton btn;
	byte[] imageData = null;
	try {
	    InputStream in = LangTools.getZipEntryInputStream(
					"ttuner_tools.jar", "images/"+imgSrc);
	    imageData = new byte[in.available()];
	    in.read(imageData);
	} catch (Exception e) {}
	if (imageData == null) {
	    btn = new JButton(text);
	    btn.setFont(menuFont);
	    btn.setForeground(Color.blue);
	} else {
	    ImageIcon icon = new ImageIcon(imageData);
	    btn = new JButton(icon);
	}

	if (toolTipText != null) {
	    btn.setToolTipText(toolTipText);
	}

	btn.addActionListener(this);
	if (imageData == null) {
	    btn.setBorder(BorderFactory.createRaisedBevelBorder());
	} else {
	    btn.setBorder(new EmptyBorder(2,2,2,2));
	    btn.addMouseListener(new MouseAdapter() {
	    	public void mouseEntered(MouseEvent e) {
		    JButton b = (JButton) e.getSource();
		    b.setBorder(BorderFactory.createRaisedBevelBorder());
		    toolBar.repaint();
	    	}
	    	public void mouseExited(MouseEvent e) {
		    JButton b = (JButton) e.getSource();
	    	    b.setBorder(new EmptyBorder(2,2,2,2));
		    toolBar.repaint();
	    	}
	    });
	}

	return btn;
    }

    /** Attempts to release resources. */
    public void releaseResources() {
	chromatFileNames = null;
	phdFileNames = null;
	if (gotoDialog != null) {
	    gotoDialog.dispose();
	    gotoDialog = null;
	}
	if (zoomXDialog != null) {
	    zoomXDialog.dispose();
	    zoomXDialog = null;
	}
	if (zoomYDialog != null) {
	    zoomYDialog.dispose();
	    zoomYDialog = null;
	}
	if (searchDialog != null) {
	    searchDialog.dispose();
	    searchDialog = null;
	}
        disposePreviewWindows();
	if (chromViewPanel != null) {
	    chromViewPanel.releaseResources();
	    chromViewPanel = null;
	}
	if (msgPanel != null) {
	    msgPanel.releaseResources();
	    msgPanel = null;
	}
	menuBar = null;
	toolBar = null;
	if (filesStatus != null) {
	    filesStatus.releaseResources();
	    filesStatus = null;
	}
	zoomEventListeners = null;
	displayFeatureEventListeners = null;
	displaySigValEventListeners = null;
	fileNavEventListeners = null;
	searchEventListeners = null;
        dispose();
        System.gc();
    }

    /** Invoked when an action was performed. */
    public void actionPerformed(ActionEvent e) {
	final Object source = e.getSource();

	if (source == quitItem) {
	    if (numPrintAllFilesJobs > 0) {
		String message = "Quiting this window will terminate"
			+ " the ongoing Print-Global-View-of-All-Files job.\n"
			+ "Do you want to quit?";
		int result = JOptionPane.showConfirmDialog(
						TTViewFrame.this, message,
						"Quit Confirmation",
						JOptionPane.YES_NO_OPTION,
						JOptionPane.QUESTION_MESSAGE);
		if (result == JOptionPane.NO_OPTION) {
		    return;
		}
	    }
            if (!exitOnClose) {
		releaseResources();
            } else {
                System.exit(1);
            }
	} else if (source == infoBtn) {
	    // a tmp solution
	    infoItem.doClick();
	} else if (source == openItem || source == openBtn) {
	    new OpenFileDialog(this);
	} else if (source == findItem || source == findBtn) {
	    if (searchDialog == null) {
		searchDialog = new SearchDialog(TTViewFrame.this);
		searchDialog.addSearchListener(chromViewPanel);
		searchDialog.addSearchListener(TTViewFrame.this);
	    }
	    searchDialog.becomeVisible();
	} else if (source == findAgainItem || source == findAgainBtn) {
	    if (searchDialog == null) {
		searchDialog = new SearchDialog(TTViewFrame.this);
		searchDialog.addSearchListener(chromViewPanel);
		searchDialog.addSearchListener(TTViewFrame.this);
	    }
	    SearchEvent se = searchDialog.getSearchEvent();
	    if (se != null) {
		for (int i = 0; i < searchEventListeners.size(); i++) {
		    ((SearchEventListener)
				searchEventListeners.elementAt(i)).search(se);
		}
	    }
	} else if (source == signalItem) {
	    if (signalItem.getState()) {
	    	msgPanel.addSignalValuePanel();
	    } else {
		msgPanel.removeSignalValuePanel();
	    }
	} else if (source == hetItem) {
	    if (hetItem.getState()) {
	    	msgPanel.addSNPsTablePanel();
	    } else {
		msgPanel.removeSNPsTablePanel();
	    }
	} else if (source == origItem) {
	    fireDisplayFeatureEvent(false);
	} else if (source == abcItem) {
	    fireDisplayFeatureEvent(false);
	} else if (source == alnItem) {
	    fireDisplayFeatureEvent(false);
	} else if (source == tipItem) {
	    fireDisplayFeatureEvent(false);
	} else if (source == totalISItem) {
	    fireDisplayFeatureEvent(false);
	} else if (source == trimItem) {
	    setTrimThreshold();
	    fireDisplayFeatureEvent(false);
	} else if (source == rawDataItem) {
	    if (rawDataItem.getState()) {
		origItem.setEnabled(false);
		abcItem.setEnabled(false);
		alnItem.setEnabled(false);
		tipItem.setEnabled(false);
		totalISItem.setEnabled(false);
		trimItem.setEnabled(false);
	    } else {
		origItem.setEnabled(true);
		abcItem.setEnabled(enableDisplayAbc);
		alnItem.setEnabled(enableDisplayAln);
		tipItem.setEnabled(enableDisplayTip);
		totalISItem.setEnabled(enableDisplayTip);
		trimItem.setEnabled(enableDisplayTrim);
	    }
	    setTrimThreshold();
	    msgPanel.clearSignalValuePanel();
	    fireDisplayFeatureEvent(false);
	} else if (source == gotoItem) {
	    if (gotoDialog == null) {
	        gotoDialog = new GoToDialog(this, numOfFiles);
	    }
	    gotoDialog.becomeVisible();
	    if (gotoDialog.isCanceled()) {
		return;
	    }
	    int temp = gotoDialog.getGotoIndex();
	    if (temp == fileIndex) {
		return;
	    }
	    fileIndex = temp;
	    fileIndexChanged();
	} else if (source == refreshItem || source == refreshBtn) {
	    /* The file index is actually not changed. Only using
	       this fileIndexChanged() method to refresh the screen. */
	    fileIndexChanged();
	} else if (source == firstBtn) {
	    fileIndex = 0;
	    fileIndexChanged();
	} else if (source == backBtn) {
	    fileIndex -= 1;
	    fileIndexChanged();
	} else if (source == forwardBtn) {
	    fileIndex += 1;
	    fileIndexChanged();
	} else if (source == lastBtn) {
	    fileIndex = numOfFiles - 1;
	    fileIndexChanged();
	} else if (source == pageSetupItem || source == pgSetupBtn) {
	    PageSetupManager.pageDialog();
	} else if (source == printViewItem 
		   || source == prnViewBtn
		   || source == previewPVItem ) {
	    Thread runner = new Thread() {
		public void run() {
		    ChromSnapShot css = new ChromSnapShot(
						      TTPrintJob.PRINT_LOCAL,
						      chromViewPanel,
						      msgPanel);
	    	    TTPrintJob tpJob = new TTPrintJob(chromViewPanel, css);
		    if (source == previewPVItem) {
			PrintPreview pp = new PrintPreview(tpJob,
						 "Print Local View Preview");
			previewWindows.add(pp);
		    } else {
		    	tpJob.print(); 
			tpJob.dispose();
	    	        tpJob = null;
		    }
		}
	    };
	    runner.start();
	} else if (source == printAllItem || source == previewPAItem) {
	    Thread runner = new Thread() {
		public void run() {
		    ChromSnapShot css = new ChromSnapShot(
						TTPrintJob.PRINT_SINGLE_GLOBAL,
						chromViewPanel, msgPanel);
	    	    TTPrintJob tpJob = new TTPrintJob(chromViewPanel, css);
		    if (source == previewPAItem) {
			PrintPreview pp = new PrintPreview(tpJob,
						  "Print Global View Preview");
			previewWindows.add(pp);
		    } else {
		    	tpJob.print(); 
			tpJob.dispose();
	    	    	tpJob = null;
		    }
		}
	    };
	    runner.start();
	} else if (source == printAllFilesItem) {
	    Thread runner = new Thread() {
		public void run() {
		    numPrintAllFilesJobs++;
		    ChromSnapShot css;
		    TTPrintJob tpJob;
		    int numOfInvalidFiles = 0;
		    boolean jobCanceled = false;
		    StringBuffer invalidFileIndexes = new StringBuffer();
		    invalidFileIndexes.append("File:\n    ");
		    /* Prints the global view of each file in the batch
		       selected in this Viewer. */
		    for (int i = 0; i < numOfFiles; i++) {
		    	css  = new ChromSnapShot(TTPrintJob.PRINT_ALL_GLOBAL,
						getSampleAndPhdFilesStatus(i));
			/* Both sample file and phd file need to be valid. */
		    	if (!css.hasValidSampleAndPhd()) {
			    css.dispose();
			    css = null;
			    numOfInvalidFiles++;
			    invalidFileIndexes.append(i+1);
			    invalidFileIndexes.append(" ");
			    if (numOfInvalidFiles%10 == 0) {
			    	invalidFileIndexes.append("\n    ");
			    }
			    continue;
		    	}
			/* The sample file and the phd file are valid,
			   create the print job for this file pair and print
			   them. */
	    	    	tpJob = new TTPrintJob(chromViewPanel, css, i==0);
		    	tpJob.print(); 
			if (tpJob.isCanceled()) {
			    jobCanceled = true;
			    tpJob.dispose();
			    tpJob = null;
			    break;
			}
			tpJob.dispose();
	    	    	tpJob = null;
		    }
		    if (!jobCanceled) {
			/* Show the report of how many files were printed. */
		    	StringBuffer sb = new StringBuffer();
			sb.append("Print Global View of All Files completed.\n");
		    	if (numOfInvalidFiles > 0) {
			    int numValidFiles = numOfFiles - numOfInvalidFiles;
			    sb.append(numValidFiles);
			    sb.append(" out of ");
			    sb.append(numOfFiles);
			    sb.append(" files ");
			    if (numValidFiles == 1) {
				sb.append("was ");
			    } else {
				sb.append("were ");
			    }
			    sb.append("printed.\n");
			    sb.append(invalidFileIndexes);
			    sb.append("\n");
			    if (numOfInvalidFiles == 1) {
				sb.append("was ");
			    } else {
				sb.append("were ");
			    }
			    sb.append("not printed due to any of ");
			    sb.append("the following possible reasons:\n");

			    sb.append(" * The .phd.1 file does not exist,\n");
			    sb.append(" * The input file contains invalid data, or\n");
			    sb.append(" * The .phd.1 file contains invalid data.\n");
		    	} else {
			    sb.append("All files were printed.");
		    	}
			JOptionPane.showMessageDialog(TTViewFrame.this,
					sb.toString(), "Print Summary",
					JOptionPane.INFORMATION_MESSAGE);
		    }
		    numPrintAllFilesJobs--;
		}
	    };
	    runner.start();
	}
    }

    /** Returns the FilesStatus object which contains the information of
	the sample file and its phd file at the specified index. */
    public FilesStatus getSampleAndPhdFilesStatus(int index) {
	int cs, ps;
	File f;
	ABIFile chromatFile = null;
	PhdFile phdFile = null;

	String chromatFileName = chromatFileNames[index];
	String phdFileName = phdFileNames[index];

	String s = phdFileName.substring(0,
			phdFileName.lastIndexOf(Constants.PHD_FILE_SUFFIX));

	f = new File(chromatFileName);
	if (!f.exists()) {
	    cs = FilesStatus.DOES_NOT_EXIST;
	} else if (!f.canRead()) {
	    cs = FilesStatus.CAN_NOT_READ;
	} else {
	    cs = FilesStatus.VALID;
	    try {
		chromatFile = new ABIFile(chromatFileName);
	    } catch (Exception e) {
		cs = FilesStatus.UNKNOWN;
	    }
	}

	f = null;
	f = new File(phdFileName);
	if (!f.exists()) {
	    ps = FilesStatus.DOES_NOT_EXIST;
	} else if (!f.canRead()) {
	    ps = FilesStatus.CAN_NOT_READ;
	} else {
	    ps = FilesStatus.VALID;
	    try {
	    	phdFile = new PhdFile(phdFileName);
	    } catch (Exception e) {
		System.err.println("TTViewFrame--fileIndexChanged:"
				   + e.getMessage());
		e.printStackTrace(System.err);
		ps = FilesStatus.UNKNOWN;
	    }
	}
	return (new FilesStatus(index + 1, numOfFiles,
				chromatFile, phdFile, null, null,
				chromatFileName, phdFileName, null, null, null,
				cs, ps, 0, 0, 0));
    }

    /** Invoked when file index was changed. */
    protected void fileIndexChanged() {
        setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

	int cs, ps, as, bs, ts;
	long phdModTime;
	File f;
	ABIFile chromatFile = null;
	PhdFile phdFile = null;
	TabFile tabFile = null;
	TalFile talFile = null;

	String chromatFileName = chromatFileNames[fileIndex];
	String phdFileName = phdFileNames[fileIndex];

	// .tip and .tal files are stored together with phd files.
	String s = phdFileName.substring(0,
			phdFileName.lastIndexOf(Constants.PHD_FILE_SUFFIX));

	String tipFileName = s + Constants.TIP_FILE_SUFFIX;
	String talFileName = s + Constants.TAL_FILE_SUFFIX;
	String tabFileName = s + Constants.TAB_FILE_SUFFIX;

	f = new File(chromatFileName);
	if (!f.exists()) {
	    cs = FilesStatus.DOES_NOT_EXIST;
	} else if (!f.canRead()) {
	    cs = FilesStatus.CAN_NOT_READ;
	} else {
	    cs = FilesStatus.VALID;
	    try {
		chromatFile = new ABIFile(chromatFileName);
	    } catch (Exception e) {
		cs = FilesStatus.UNKNOWN;
	    }
	}

	f = null;
	f = new File(phdFileName);
	phdModTime = f.lastModified();
	if (!f.exists()) {
	    ps = FilesStatus.DOES_NOT_EXIST;
	    enableDisplayTrim = false;
	} else if (!f.canRead()) {
	    ps = FilesStatus.CAN_NOT_READ;
	    enableDisplayTrim = false;
	} else {
	    ps = FilesStatus.VALID;
	    try {
		phdFile = new PhdFile(phdFileName);
		int trim[] = phdFile.getTrimData();
		if (trim == null) {
		    enableDisplayTrim = false;
		} else {
		    enableDisplayTrim = true;
		    trimThreshold = trim[2];
		}
	    } catch (Exception e) {
		System.err.println("TTViewFrame--fileIndexChanged:"
				   + e.getMessage());
		e.printStackTrace(System.err);
		ps = FilesStatus.UNKNOWN;
		enableDisplayTrim = false;
	    }
	} 

	f = null;
	f = new File(tabFileName);
	if (!f.exists()) {
	    bs = FilesStatus.DOES_NOT_EXIST;
	} else if (!f.canRead()) {
	    bs = FilesStatus.CAN_NOT_READ;
	} else {
	    try {
		tabFile = new TabFile(tabFileName);
		if (tabFile.hasInvalidFileData()) {
		    bs = FilesStatus.INVALID_DATA;
		} else {
		    if (tabFile.getNumOfAbcs() == 0) {
			bs = FilesStatus.NO_ABC;
		    } else if (tabFile.hasInvalidAbcData()) {
			bs = FilesStatus.INVALID_ABC_DATA;
		    } else {
			/* If the last modification time of the .tab file is
			   >1 minute older than the .phd.1 file, the .tab
			   file is considered to be stale. */
			long temp = Math.abs(f.lastModified() - phdModTime);
			if (temp > STALE_TIME) {
			    bs = FilesStatus.VALID_BUT_STALE;
			} else {
			    bs = FilesStatus.VALID;
			}
		    }
		}
	    } catch (Exception e) {
		System.err.println("TTViewFrame--fileIndexChanged:"
				   + e.getMessage());
		e.printStackTrace(System.err);
		bs = FilesStatus.UNKNOWN;
	    }
	}

        f = null;
        f = new File(tipFileName);
	if (!f.exists()) {
	    ts = FilesStatus.DOES_NOT_EXIST;
	} else if (!f.canRead()) {
	    ts = FilesStatus.CAN_NOT_READ;
	} else {
	    /* If the last modification time of the .tal file is >1 minute
	     * older than the .phd.1 file, the .tal file is considered to
	     * be stale. 
             */
	    long temp = Math.abs(f.lastModified() - phdModTime);
	    if (temp > STALE_TIME) {
	        ts = FilesStatus.VALID_BUT_STALE;
	    } else {
	        ts = FilesStatus.VALID;
	    }
	} 

	f = null;
	f = new File(talFileName);
	if (!f.exists()) {
	    as = FilesStatus.DOES_NOT_EXIST;
	} else if (!f.canRead()) {
	    as = FilesStatus.CAN_NOT_READ;
	} else {
	    try {
	    	talFile = new TalFile(talFileName);
		if (talFile.hasInvalidFileData()) {
		    as = FilesStatus.INVALID_DATA;
		} else {
	    	    int alignStatus = talFile.getAlignStatus();
		    if (alignStatus != TalFile.OK) {
		    	switch (alignStatus) {
                    	case TalFile.NO_GOOD_ALIGN:
			    as = FilesStatus.NO_GOOD_ALIGN;
                            break;
                    	case TalFile.POSSIBLE_REPEATS:
			    as = FilesStatus.POSSIBLE_REPEATS;
                            break;
                    	case TalFile.BAD_PROCESSING:
			    as = FilesStatus.BAD_PROCESSING;
                            break;
                    	case TalFile.UNKNOWN:
                    	default:
			    as = FilesStatus.UNKNOWN;
                            break;
                    	}
		    } else if (talFile.hasInvalidAlignData()) {
			as = FilesStatus.INVALID_ALIGN_DATA;
		    } else {
			/* If the last modification time of the .tip file is
			   >1 minute older than the .phd.1 file, the .tip
			   file is considered to be stale. */
	    		long temp = Math.abs(f.lastModified() - phdModTime);
	    		if (temp > STALE_TIME) {
			    as = FilesStatus.VALID_BUT_STALE;
			} else {
			    as = FilesStatus.VALID;
			}
		    }
		}
	    } catch (Exception e) {
		System.err.println("TTViewFrame--fileIndexChanged:"
				   + e.getMessage());
		e.printStackTrace(System.err);
		as = FilesStatus.UNKNOWN;
	    }
	} 

	filesStatus = new FilesStatus(fileIndex + 1, numOfFiles,
				     chromatFile, phdFile, tabFile, talFile,
				     chromatFileName, phdFileName,
				     tabFileName, talFileName, tipFileName,
				     cs, ps, bs, as, ts);

	setNavigationButtons();
	setMenuStates();
	fireDisplayFeatureEvent(true);

	FileNavigateEvent e = 
	 	new FileNavigateEvent(this, filesStatus, numOfFiles);
        for (int i = 0; i < fileNavEventListeners.size(); i++) {
            ((FileNavigateEventListener)
                      (fileNavEventListeners.elementAt(i))).navigateTo(e);
	}
    }

    /** Sets the states of the navigation tool buttons. */
    protected void setNavigationButtons() {
	if (numOfFiles <= 1) {
	    firstBtn.setEnabled(false);
	    backBtn.setEnabled(false);
	    forwardBtn.setEnabled(false);
	    lastBtn.setEnabled(false);
	    return;
	}

	// more than 1 files
	if (fileIndex == 0) {
	    firstBtn.setEnabled(false);
	    backBtn.setEnabled(false);
	    forwardBtn.setEnabled(true);
	    lastBtn.setEnabled(true);
	    return;
	}
	if ( fileIndex == (numOfFiles - 1) ) {
	    firstBtn.setEnabled(true);
	    backBtn.setEnabled(true);
	    forwardBtn.setEnabled(false);
	    lastBtn.setEnabled(false);
	    return;
	}

	firstBtn.setEnabled(true);
	backBtn.setEnabled(true);
	forwardBtn.setEnabled(true);
	lastBtn.setEnabled(true);
    }

    /** Sets the states of the menus and menu items. */
    protected void setMenuStates() {
	int sampleFileStatus = filesStatus.getSampleFileStatus();
	int phdFileStatus = filesStatus.getPhdFileStatus();
	int tabFileStatus = filesStatus.getTabFileStatus();
	int tipFileStatus = filesStatus.getTipFileStatus();
	int talFileStatus = filesStatus.getTalFileStatus();

	if ((sampleFileStatus != FilesStatus.VALID)
	        || (phdFileStatus != FilesStatus.VALID)) {
	    featureMenu.setEnabled(false);
	    zoomXMenu.setEnabled(false);
	    zoomYMenu.setEnabled(false);
	    findItem.setEnabled(false);
	    findAgainItem.setEnabled(false);
	    findBtn.setEnabled(false);
	    findAgainBtn.setEnabled(false);
	    signalItem.setEnabled(false);
	    hetItem.setEnabled(false);
	    printPreviewMenu.setEnabled(false);
	    printViewItem.setEnabled(false);
	    printAllItem.setEnabled(false);
	    pgSetupBtn.setEnabled(false);
	    prnViewBtn.setEnabled(false);
	    return;
	} else {
	    featureMenu.setEnabled(true);
	    zoomXMenu.setEnabled(true);
	    zoomYMenu.setEnabled(true);
	    findItem.setEnabled(true);
	    findBtn.setEnabled(true);
	    findAgainItem.setEnabled(enableFindAgain);
	    findAgainBtn.setEnabled(enableFindAgain);
	    signalItem.setEnabled(true);
	    hetItem.setEnabled(true);
	    printPreviewMenu.setEnabled(true);
	    printViewItem.setEnabled(true);
	    printAllItem.setEnabled(true);
	    if (numOfFiles > 1) {
	    	printAllFilesItem.setEnabled(true);
	    }
	    pgSetupBtn.setEnabled(true);
	    prnViewBtn.setEnabled(true);

	    String sampleFileName = filesStatus.getSampleFileName();
	    int chromType = 0;
	    try {
	    	chromType = AbstractChromatogram.getChromType(sampleFileName);
            } catch (Exception e) {
	    	System.err.println("TTViewFrame--setMenuStates: " 
					   + e.getMessage());
	    	e.printStackTrace(System.err);
	    }
	
	    if (chromType == AbstractChromatogram.ABI) {
	    	rawDataItem.setEnabled(true);
	    } else {
	    	rawDataItem.setEnabled(false);
	    }
	}

	enableDisplayAbc = (tabFileStatus == FilesStatus.VALID
			       || tabFileStatus == FilesStatus.VALID_BUT_STALE)
	    		   ? true : false;
	enableDisplayAln = (talFileStatus == FilesStatus.VALID
			       || talFileStatus == FilesStatus.VALID_BUT_STALE)
	    		   ? true : false;
	enableDisplayTip = ((tipFileStatus == FilesStatus.VALID
			       || tipFileStatus == FilesStatus.VALID_BUT_STALE)
			        && zoomX > 0.5)
	    		   ? true : false;

	if (rawDataItem.isEnabled() && rawDataItem.getState() == true) {
	    alnItem.setEnabled(false);
	    origItem.setEnabled(false);
	    abcItem.setEnabled(false);
	    tipItem.setEnabled(false);
	    totalISItem.setEnabled(false);
	    trimItem.setEnabled(false);
	} else {
	    origItem.setEnabled(true);
	    abcItem.setEnabled(enableDisplayAbc);
	    tipItem.setEnabled(enableDisplayTip);
	    totalISItem.setEnabled(enableDisplayTip);
	    alnItem.setEnabled(enableDisplayAln);
	    trimItem.setEnabled(enableDisplayTrim);
	}

	setTrimThreshold();
    }

    /** Sets the states of trim threshold label, field, and the value. */
    protected void setTrimThreshold() {
	if (trimItem.isEnabled() && trimItem.getState()) {
	    trimLabel.setEnabled(true);
	    trimField.setEnabled(true);
	    trimField.setText(trimThreshold+"");
	} else {
	    trimLabel.setEnabled(false);
	    trimField.setEnabled(false);
	    trimField.setText("");
	}
    }

    /** Fires a DisplayFeatureEvent. */
    protected void fireDisplayFeatureEvent(boolean fileIndexChanged) {
	boolean displayRaw = rawDataItem.isEnabled() && rawDataItem.getState();
	boolean displayOrig = origItem.isEnabled() && origItem.getState();
	boolean displayAbc = abcItem.isEnabled() && abcItem.getState();
	boolean displayAln = alnItem.isEnabled() && alnItem.getState();
	boolean displayTip = tipItem.isEnabled() && tipItem.getState();
	boolean displayTotalIS = Constants.runAsDev && totalISItem.isEnabled()
	    				&& totalISItem.getState();
	boolean displayTrim = trimItem.isEnabled() && trimItem.getState();

	DisplayFeatureEvent de = new DisplayFeatureEvent(this, displayRaw, 
						    displayOrig, displayAbc,
						    displayAln, displayTip,
						    displayTotalIS,
						    displayTrim,
						    fileIndexChanged);
        for (int i = 0; i < displayFeatureEventListeners.size(); i++) {
            ((DisplayFeatureEventListener)
                      (displayFeatureEventListeners.elementAt(i)))
	    				.displayFeatureChanged(de);
	}
    }

    /** Disposes all existing preview windows, if any. */
    public void disposePreviewWindows() {
	if (previewWindows == null) {
	    return;
	}
	Object o;
	for (int i = previewWindows.size() - 1; i >= 0; i--) {
	    o = previewWindows.elementAt(i);
	    if (o == null) {
	        previewWindows.removeElementAt(i);
		continue; 
	    }
	    ((PrintPreview)o).releaseResources();
	    ((Window)o).dispose();
	    o = null;

	    previewWindows.removeElementAt(i);
	}
	previewWindows = null;
    }

    /** Invoked when a SearchEvent is fired.  Does nothing but enable the
        "Find Again" option and button when the first search event is fired.
	If none search event is ever fired, "Find Again" should be disabled. */
    public void search(SearchEvent se) {
	if (!enableFindAgain) {
	    enableFindAgain = true;
	    findAgainItem.setEnabled(true);
	    findAgainBtn.setEnabled(true);
	}
    }

    /** Locates the View window at the center of the screen. */
    protected void centerWindow() {
	Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
	Dimension size = getPreferredSize();
	int x = (screen.width - size.width) / 2;
	int y = (screen.height - size.height) / 2;
	if (x < 0) {
	    x = 72;
	}
	if (y < 0) {
	    y = 72;
	}
	setLocation(x, y);
    }
}
