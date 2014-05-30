/*
 * 1.14
 *
 * @(#)TTrun.java            2000/11/22
 *
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
 *
 * This software is the proprietary information of Paracel, Inc.
 * Use is subject to license terms.
 */

package com.paracel.tt.run;

import java.io.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import java.util.*;
import java.text.*;
import com.paracel.tt.util.*;
import com.paracel.tt.view.*;

/** This class creates the TraceTuner Launcher.  */
public class TTrun extends JPanel implements ActionListener, Runnable {

    /** The frame this panel resides in. */
    static JFrame frame;

    /** The thread for executing ttuner command. */
    private Thread runThread = null;

    /** The command builder. */
    private CommandBuilder cmdBuilder;

    /** The sample file chooser or browser. */
    private Filespec inputFileChooser;

    /** The buttons. */
    private JButton chgOptsBtn, runBtn, viewBtn, helpBtn, exitBtn;

    /** The text area for displaying the log messages. */
    JTextArea   logText;

    /** The scroll panel that the log message text area resides in. */
    JScrollPane logPanel;

    /** The writer to output the log messages to a file. */
    private BufferedWriter logWriter = null;

    /** The File contains quality value report. */
    private File qrFile = null;

    /** File contains name of sample files */
    private File ifFile = null;         

    /** The Options frame. */
    private BasicOptionsFrame optsFrame;

    /** The error message. */
    private String errMesg;

    private Process process;

    /** The constructor which sets up the Launcher frame. */
    public TTrun() {
    	inputFileChooser  = new Filespec("Input File(s)",
                                           JFileChooser.FILES_ONLY,
                                           "Select one or more sample files"); 
	inputFileChooser.setMultipleSelection();
        inputFileChooser.addFilter(new Filter(".scf"));
        inputFileChooser.addFilter(new Filter(".abi"));
        inputFileChooser.addFilter(new Filter(".abd"));
        inputFileChooser.addFilter(new Filter(".ab1"));

	setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

	chgOptsBtn  = new JButton("Change Options...");
	runBtn  = new JButton("Run");
	viewBtn  = new JButton("View");
	helpBtn  = new JButton(new HelpAction("Help"));
	exitBtn  = new JButton("Exit");


        JLabel labelmain = new JLabel("Select input sample files for " + 
                                      "analysis", JLabel.CENTER);
	JLabel labelopt = new JLabel("Advanced users may wish to change " +
                                     "output and processing " +
                                     "options", JLabel.CENTER);
        JLabel labelrun  = new JLabel("Run TraceTuner or view results for " + 
                                      "selected files", JLabel.CENTER);
	JPanel checkmain  = new JPanel();
	JPanel checkopt   = new JPanel();
	JPanel checkrun   = new JPanel();
        JPanel checkarea1 = new JPanel();
	JPanel checkarea2 = new JPanel();
	JPanel checkarea3 = new JPanel();

        checkmain.setLayout( new BoxLayout(checkmain,  BoxLayout.X_AXIS));

	JPanel buttonarea1 = new JPanel();
	JPanel buttonarea2 = new JPanel(new GridLayout(2, 2, 6, 16));

	Box logtitlearea = Box.createHorizontalBox();
	logtitlearea.add(new JLabel(" Log Messages:"));
	logtitlearea.add(Box.createHorizontalGlue());

	//logPanel.getViewport().add(logText);
	logText = new JTextArea("");
	logPanel = new JScrollPane(logText);
	logPanel.setVerticalScrollBarPolicy(
			JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
	logPanel.getViewport().setScrollMode(JViewport.SIMPLE_SCROLL_MODE);

	logText.setFont(Fonts.TEXT_FONT);
	logText.setLineWrap(true);
        logText.setEditable(false);
	logPanel.setBorder(new LineBorder(Color.black,2));
	logPanel.setPreferredSize(new Dimension(420, 250));

	chgOptsBtn.addActionListener(this);
	runBtn.addActionListener(this);
	viewBtn.addActionListener(this);
	//helpBtn.addActionListener(this);
	exitBtn.addActionListener(this);

        checkmain.add(labelmain);
	checkarea1.add(checkmain);
	add(checkarea1);

	add(inputFileChooser);
	checkopt.add(labelopt);
	checkarea2.add(checkopt);
	add(checkarea2);
        buttonarea1.add(chgOptsBtn);
        add(buttonarea1); 
        checkrun.add(labelrun);
        checkarea3.add(checkrun);
        add(checkarea3);
	buttonarea2.add(runBtn);
	buttonarea2.add(viewBtn);
	buttonarea2.add(helpBtn);
	buttonarea2.add(exitBtn);
	add(buttonarea2);

	//  Put on the log area title and log area
	add(Box.createRigidArea(new Dimension(0,10)));
	add(logtitlearea);
	add(Box.createRigidArea(new Dimension(0,4)));
	add(logPanel);

	//  Want the scrolling region to absorb any extra space 
	//  when the window resizes.  Best way to do this seems 
	//  to be to set a maximum size on the other components.
	//  Boxlayout in particular will obey these maximum sizes.
	checkarea1.setMaximumSize(new Dimension(420, 60));
	checkarea2.setMaximumSize(new Dimension(720, 60));
	checkarea3.setMaximumSize(new Dimension(420, 60));

	inputFileChooser.setMaximumSize(new Dimension(500, 80));
	buttonarea1.setMaximumSize(new Dimension(420, 60));
	buttonarea2.setMaximumSize(new Dimension(140, 260));

	optsFrame = BasicOptionsFrame.getInstance();
	optsFrame.setParentFrame(frame);
    }

    /** Invoked when an action is performed. */
    public void actionPerformed(ActionEvent e) {
	if (e.getActionCommand().equals("Run")){
	    //  Only allow one ttuner thread at a time.  
	    if(runThread == null) {
		/* Changes the button to Stop button, which allows the user
		   to terminate the current run. */
		runBtn.setText("Stop");
		runBtn.setForeground(Color.red);
		runTT();
	    }
	}

	if (e.getActionCommand().equals("Stop")){
	    if(runThread != null) {
	        DateFormat fmt = DateFormat.getDateTimeInstance();
	    	printLog(Constants.NEW_LINE
		      + "Stopped running TTuner at " + fmt.format(new Date()) 
		      + Constants.NEW_LINE + Constants.NEW_LINE
		      + Constants.NEW_LINE);
		cleanUpRun();
	    }
	}

	if (e.getActionCommand().equals("View")) {
	    File[] chromatFiles = inputFileChooser.getSelectedFiles();
	    if (chromatFiles == null) {
		JOptionPane.showMessageDialog(frame,
			"Please specify input sample file(s) for viewing.\n");
		return;
	    }

	    if (chromatFiles.length <= 0) {
		JOptionPane.showMessageDialog(frame,
		   "Please specify valid input sample file(s) for viewing.\n");
		return;
	    }

	    String outputDir = optsFrame.getOutputDir();

	    if (outputDir == null) {
		// invalid output dir specified in the Options frame
		return;
	    }

	    String[] chromatFileNames = new String[chromatFiles.length];
	    String[] phdFileNames = new String[chromatFiles.length];

	    for (int i = 0; i < chromatFiles.length; i++) {
	     
		chromatFileNames[i] = chromatFiles[i].getAbsolutePath();
	   	if (outputDir.length() == 0) {
		    // no output dir specified in the Options frame
		    phdFileNames[i] = chromatFiles[i].getAbsolutePath()
					+ Constants.PHD_FILE_SUFFIX;
	    	} else {
		    phdFileNames[i] = outputDir + Constants.FILE_SEPARATOR
					+ chromatFiles[i].getName()
					+ Constants.PHD_FILE_SUFFIX;
		}
	    }

	    new TTViewFrame(chromatFileNames, phdFileNames, null, false);
	}

	if (e.getActionCommand().equals("Exit")) {
	    System.exit(0);
	}

        // Use default setting if no previously saved configuration exists
        if (e.getSource() == chgOptsBtn) {
	    if (!optsFrame.isVisible()) {
	    	optsFrame.becomeVisible();
	    } else if (optsFrame.getState() == Frame.ICONIFIED) {
		optsFrame.setState(Frame.NORMAL);
	    }
	    optsFrame.toFront();
	}		
    }

    /** Validates the inputs and starts an asyncronous thread to run 
        the TraceTuner command compiled from the Options. */
    private void runTT() {
	runThread = new Thread(this, "run ttuner");
        File[] inputFiles = inputFileChooser.getSelectedFiles();
	if (inputFiles == null) {
	    JOptionPane.showMessageDialog(frame,
			"Please specify input sample file(s)\n"
			+ "to run TraceTuner.\n");
	    cleanUpRun();
	    return;
	}

	if (inputFiles.length <= 0) {
	    JOptionPane.showMessageDialog(frame,
			"Please specify valid input sample file(s)\n"
			+ "to run TraceTuner.\n");
	    cleanUpRun();
	    return;
	}

	String outputDir = optsFrame.getOutputDir();
	if (outputDir == null) {
	    cleanUpRun();
	    return;
	}

	if (outputDir.length() == 0) {
	    outputDir = inputFiles[0].getAbsolutePath();
	    outputDir = outputDir.substring(0, outputDir.lastIndexOf(
						Constants.FILE_SEPARATOR));
	}

	cmdBuilder = optsFrame.buildCommand(outputDir);
	if (cmdBuilder == null || cmdBuilder.isEmpty()) {
	    cleanUpRun();
	    return;
	}

	//  Prepare file-of-files for ttuner command line.  This will be the
	//  first place where we'll run into permission or other problems
	//  with the output directory.
        ifFile = new File(outputDir, Constants.FILES_NAME_FILE);
        try {
            ifFile.createNewFile();
            BufferedWriter ifWriter = new BufferedWriter(
                                    new FileWriter(ifFile.getPath()));
            for (int j=0; j < inputFiles.length; j++) {
                 ifWriter.write(inputFiles[j].getPath());
                 ifWriter.write("\n");
            }
            ifWriter.close();
        } catch (IOException e) {
	    JOptionPane.showMessageDialog(frame, 
	        "Error writing sample file names to file:\n" + 
		ifFile.getPath() + "\n" + e.getMessage());
	    cleanUpRun();
	    return;
        } 

        // Generate quality report if running with more than one sample file
        if (inputFiles.length > 1) {
            qrFile = new File(outputDir, Constants.QV_REPORT_FILE);
            try {
                qrFile.createNewFile();
                BufferedWriter qrWriter = new BufferedWriter(
                                        new FileWriter(qrFile.getPath()));
                qrWriter.close();
            } catch (IOException e) {
	        JOptionPane.showMessageDialog(frame, 
	            "Error writing quality report to file: \n" 
		    + qrFile.getPath() + e.getMessage());
	    	cleanUpRun();
                return;
            }         
        }

	if (cmdBuilder.contains("-pd")) {
	    if (phdFilesExist(inputFiles, outputDir)) {
		Object[] options = {"Yes", "No"};
		int n = JOptionPane.showOptionDialog(frame, errMesg,
						 "Question",
						 JOptionPane.YES_NO_OPTION,
						 JOptionPane.QUESTION_MESSAGE,
						 null, options, options[0]);
		if (n == 1) {
	    	    cleanUpRun();
		    return;
		}
	    }
	}

	try {
	    logWriter = null;
	    File logFile = new File(outputDir, Constants.LOG_FILE);
	    logWriter = new BufferedWriter(
			     new FileWriter(logFile.getPath(),true));
	} catch (Exception excp) {
	    System.out.println("Error occurred with logfile: " 
				   + excp.getMessage());
	    System.err.println("Error occurred with logfile: " 
				   + excp.getMessage());
	    logWriter = null;
	}

	// Build the command string we'll actually use to run ttuner
	cmdBuilder.append("-if");
	cmdBuilder.append(ifFile.getPath());

	int logLength = logText.getText().length();
	if (logLength >= 1500000) {
	        JOptionPane.showMessageDialog(frame, 
	            "Log Messages window will now be cleared to recycle memory it occupies.\n" 
		    + "The record of the previous runs can be found in ttlog.txt file located in the output directory.");
		logText.setText("");
	}
        printLog(cmdBuilder.getCommandString() 
		 + Constants.NEW_LINE + Constants.NEW_LINE); 
            
	//  Okay, everything is setup.  Run ttuner is a separate thread.
	runThread.start();
	return;
    }

    /** Checks whether phd file(s) of the specified sample file(s) 
	exist(s) in the specified output directory. */
    protected boolean phdFilesExist(File[] inputFiles, String outDir) {
	if (inputFiles.length == 1) {
	    String phdName = outDir + Constants.FILE_SEPARATOR
		    		 + inputFiles[0].getName() 
		    		 + Constants.PHD_FILE_SUFFIX;
	    File f = new File(phdName);
	    if (f.exists()) {
		errMesg = "File " + phdName + "\n exists. Overwrite?";
		return true;
	    } else {
		return false;
	    }
	} else if (inputFiles.length > 1) {
	    File f = new File(outDir);
	    String[] fileList = f.list();
	    for (int i = 0; i < fileList.length; i++) {
		if (fileList[i].length() > 6
		        && fileList[i].endsWith(Constants.PHD_FILE_SUFFIX)) {
		    errMesg = "The directory " + outDir
			      + "\ncontains phd files.  Overwrite?";
		    return true;
		}
	    }
	    return false;
	} else {
	    return false;
	}
    }

    /**
     *  Adds to the bottom of the attached log message text
     *  area.  Note that if it's called from within an event handler
     *  java/swing will not actually draw any of the new text until the
     *  event handling thread is completed.
     *
     *  We also use this function to copy all of the log messages out to
     *  an external file in the output directory.  If the output directory
     *  changes we need to re-open the file in a new location.
     */
    public void printLog(String s) {
	logText.append(s);    // Append s to the log messages textarea

	try {
	    // Scroll to the bottom of the log message window
	    int actualHeight = logText.getSize().height;
	    int windowHeight = logPanel.getSize().height;
	    int windowWidth = logPanel.getSize().width;
	    if (actualHeight > windowHeight) {
		logPanel.getViewport().scrollRectToVisible(
				new Rectangle(0, actualHeight - windowHeight, 
					      windowWidth, windowHeight));
	    }

	    if(logWriter != null) {
		// Append s to the text logfile
		logWriter.write(s);
	    }
	}
	catch (Exception e) {
	    System.err.println("TTrun: " + e.toString());
	    e.printStackTrace(System.err);
	}
    }

    /** Runs in a separate thread to execute the ttuner command. */
    public void run() {
	Runtime runTime = Runtime.getRuntime();
	try {
	    process = runTime.exec(cmdBuilder.getCommandArray());  
	    InputStream inputStream = process.getErrorStream();
	    byte[]      buffer = new byte[1024];
	    int		length = inputStream.read(buffer);

	    while(length != -1) {
		printLog(new String(buffer, 0, length));
		length = inputStream.read(buffer);
	    }  

	    // Trash varaibles that will not be used any more.
	    buffer = null;

            // Output content of the quality report to
            // the logText area and to the ttlog file
            if (qrFile != null) {
                BufferedReader qrReader = new BufferedReader(
						new FileReader(qrFile));
 
		String outString = "";
		String line = qrReader.readLine();
                while (line != null) {
		    outString += line + Constants.NEW_LINE;
		    if (outString.length() >= 400) {
                    	printLog(outString);
			outString = "";
		    }
                    line = qrReader.readLine();
                }
		if (outString.length() > 0) {
		    printLog(outString);
		}
                qrReader.close();
		qrReader = null;
                printLog(Constants.NEW_LINE);
                qrFile.delete();
		qrFile = null;
            }

	    // runTime.exec should now have finished
	    DateFormat fmt = DateFormat.getDateTimeInstance();
	    if (process != null) {
	    	printLog("Done running TTuner at " + fmt.format(new Date()) 
		      	 + Constants.NEW_LINE + Constants.NEW_LINE
		      	 + Constants.NEW_LINE);
	    }
		      

	} catch (IOException e) {
	    JOptionPane.showMessageDialog(frame,
                   "Create run ttuner process failed\n"+
                   "Configuration error occurred\n"+
                   "Check for ttuner executable in current directory\n"+
                   e.getMessage());
	} finally {
	    cleanUpRun();
	}
    }

    /** Performs the clean-up work once a run is done or stopped. */
    protected void cleanUpRun() {
        if (ifFile != null) {
            ifFile.delete();
	    ifFile = null;
	}
        if (qrFile != null) {
            qrFile.delete();
	    qrFile = null;
	}
	// Close the output logfile
	if(logWriter != null) {
	    try {
	    	logWriter.close();
	    } catch(Exception e) {}
	    logWriter = null;
	}
	if (process != null) {
	    process.destroy();
	    process = null;
	}
	runThread = null;
	runBtn.setText("Run");
	runBtn.setForeground(Color.black);
    }

    /** The main function. */
    public static void main(String args[]) {
	ErrorLog.turnOn(Constants.ERROR_LOG_FILE);
	ErrorLog.printDate();

	// Plug in the Windows Look&Feel if running on windows
	try {
	    if (Constants.OS_NAME.toLowerCase().trim().indexOf("windows")>=0) {
		UIManager.setLookAndFeel(
			"com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
	    }
	} catch (Exception e) {
	    System.err.println("Couldn't load Window look and feel: " + e);
	}

	if (args.length > 0 && args[0].equalsIgnoreCase("-dev")) {
	    Constants.runAsDev = true;
	}
	frame = new JFrame("Paracel TraceTuner " + Constants.VERSION 
			   + " Launcher");

	frame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {
	    	System.exit(0);
	    }
	});
	frame.getContentPane().add(new TTrun());
	frame.setLocation(36, 36);
	frame.pack();
	frame.setVisible(true);
    }
}
