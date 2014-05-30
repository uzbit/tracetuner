/**
 * 1.22
 *
 * @(#)HelpFrames.java	1.0	2000/12/18
 * 
 * Copyright (c) 2000 Paracel, Inc.  All rights reserved.
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
import com.paracel.tt.util.*;

/**
 * Creates the help info frames, which are invoked by clicking the
 * "Help" buttons.  There are totally three different help buttons
 * across the application: one for the Launcher window, one for the
 * Basic Options window and Advance Options window, one for the Viewer window.
 */
public class HelpFrames {

    /** The help text for the Viewer window. */
    public final static String TTVIEW_HELP =
"----------------------------------                            ------------\n"+
"Paracel TraceTuner " + Constants.VERSION + " Viewer Help                               June-2001\n"+
"----------------------------------                            ------------\n"+
"\n"+
"The TraceTuner Viewer displays chromatograms and quality values.\n"+
"\n"+
"The four colored plots are the chromatograms from the ABI sample files.\n"+
"Color maps are GREEN for A, BLUE for C, BLACK for G, RED for T.  ABI's\n"+
"base calls contains some Ns for undertermined peaks and they are written\n"+
"in LIGHT GRAY.  TraceTuner's base calls may contain some het/mixed bases\n"+
"(R, Y, K, M, S, or W), and they are also written in LIGHT GRAY.\n"+
"\n"+
"In the upper part of the window, the TraceTuner called bases and their \n"+
"indices are displayed.  The quality value bars and the numerical quality\n"+
"values are displayed under the called bases.  If users choose to display\n"+
"trimmed regions, the TraceTuner called bases within the trimmed regions\n"+
"will be shown with a yellow background.  Also users can choose to display\n"+ 
"the ABI called bases above the TraceTuner called bases and to display\n"+
"the best alignment between the consensus or reference sequence and the\n"+
"TraceTuner called bases, with the indices of the consensus or reference\n"+
"sequence displayed above the indices of the TraceTuner called bases.\n"+
"The consensus or reference sequence indices are written in gray, and the\n"+
"TraceTuner base indices in black.  The sources of the base calls are\n"+
"displayed at the left side of the window utilizing five labels: TT/PH,\n"+
"ALT, ABI, and REF, where TT means TraceTuner base calls, PH means Phred\n"+
"base calls, ABI means TraceTuner alternative base calls, ABI means\n"+
"original ABI base calls, and REF means consensus/reference sequence bases.\n"+
"In addition, three other labels: tt, ph, and ref, are used to mark the\n"+
"indexes of TraceTuner called bases, the indexes of Phred called bases,\n"+
"and the indexes of consensus/reference sequence, respectviely.\n"+
"\n"+
"Users can find out the exact peak locations of the TT base calls by double\n"+
"clicking at the neighborhood of a peak or the TT bases. A vertical line\n"+
"will be shown that passes through the called base, the quality value bar\n"+
"and the peak of the trace. The peak location is indicated at the bottom of\n"+
"the line. If ABI called bases are displayed, one can look at the ABI peak\n"+
"locations by clicking at the ABI bases.  Similarly one can also look at\n"+
"the alternative base call locations by clicking at the alternative base\n"+
"calls, if they are displayed.\n"+
"\n"+
"Users can find out the signal strength for each of the traces by selecting\n"+
"the Show Signal Value menu option in the View menu, then clicking at\n"+
"any point on a trace.  The location where the mouse click occurred and\n"+
"the intensity of each of the traces are shown in a color coded scheme at\n"+
"the bottom of the Viewer window.  The same information will be printed as\n"+
"part of the header information if user decides to do a screen shot.\n"+
"\n"+
"File names of the PHD and ABI sample files and other auxiliary files are\n"+
"displayed at the bottom of the Viewer window, as well as their statuses .\n"+
"Users can specify the PHD file, .tab, and .tal file directory by\n"+
"clicking on the Change Options... button in the Launcher window.  The\n"+
"Basic Options window will appear, in which the users can specify the\n"+
"Output Directory.\n"+
"\n"+
"The text on the title bar indicates the TraceTuner version and the viewer\n"+
"window number.\n"+
"\n"+
"Users can view a list of input files and their PHD files along with their\n"+
"available auxiliary files.  The First Sample File, Previous Sample\n"+
"File, Next Sample File, and Last Sample File tool buttons on the tool bar\n"+
"and the Go to... option in the File menu can be used to navigate through\n"+
"the list of the input files selected for the Viewer.\n\n"+
"When the content of the currently shown files changes, the Refresh\n"+
"option in the View menu or the Refresh button on the tool bar can be\n"+
"used to update the display content of the viewer window.\n\n"+
"Viewer has the following menu options:\n\n"+
"    o  File - -\n"+
"           *   Open  - A browse window will be given to choose a sample\n"+
"                       file and a PHD file (file with .phd.1 as suffix).\n"+
"                       If the PHD file is not chosen, it will be chosen\n"+
"                       automatically from the sample file directory.  This\n"+
"                       new set of PHD and sample file and their available\n"+
"                       auxiliary files will be displayed in a new Viewer.\n"+
"\n"+
"           *   Go to  - A window will be given to choose the index of an\n"+
"                       input file in the batch for display in the Viewer.\n"+
"\n"+
"           *   Page Setup  - Displays the page setup that will be used by\n"+
"                       the printer.  The user may change any of the\n"+
"                       available options.  It is recommended that the\n"+
"                       Landscape parameter be selected in the Page Setup\n"+
"                       option.\n"+
"\n"+
"           *   Print Preview - A window will be given to display what the\n"+
"                       printout looks like on screen.  The submenu items\n"+
"                       are as follow:\n"+
"                       *   Print Local View\n"+
"                       *   Print Global View of Current File\n"+
"                       Each of the submenu will have the same features.  A\n"+
"                       pull down selection list with a set of four zoom\n"+
"                       values to choose from.  The numerical value can\n"+
"                       be overwritten by inputting a more suitable zoom\n"+
"                       factor of user's choice.  Use the Print button if\n"+
"                       the printout looks desirable on the print preview\n"+
"                       screen.  Use the Close button to quit the print\n"+
"                       preview window.\n"+
"\n"+
"           *   Print Local View - Prints a screen shot providing a view of\n"+
"                       what is currently displayed in the Viewer window.\n"+
"                       This will be influenced by user's selection of\n"+
"                       features at the time when the print button is\n"+
"                       pushed.\n"+
"\n"+
"           *   Print Global View of Current File - Prints an overview of\n"+
"                       the four colored analyazed traces, TT base calls,\n"+
"                       the quality value bars and the actual numerical\n"+
"                       quality values.  The print out may not reflect what\n"+
"                       is currently displayed in the Viewer window.\n"+
"\n"+
"           *   Print Global View of All Files - Prints overviews of\n"+
"                       the four colored analyazed traces, TT base calls,\n"+
"                       the quality value bars and the actual numerical\n"+
"                       quality values of all files selected in the Viewer.\n"+
"                       This option is disabled if only one sample file is\n"+
"                       selected for this Viewer window.\n\n"+
"                       NOTE: In the current release, this print option\n"+
"                       does not function correctly if the user chooses to\n"+
"                       print to a file or a non-system-default printer.\n"+
"                       Please make sure the print destination is set to\n"+
"                       the system-default printer while using this print\n"+
"                       option.\n"+
"\n"+
"           *   Quit  - Closes the Viewer window, but not the Launcher.\n"+
"\n"+
"    o  View - -\n"+
"           *   Zoom x - Displays data at different scales along the x axis - \n"+
"                        default: the default scale is one pixel per datum.\n"+
"                        x 6.25%: 1/16 of the size of default\n"+
"                        x 300 %: 300 % of the default size\n"+
"                        x 600 %: 600 % of the default size\n"+
"                        x 1000%: 1000 % of the default size\n"+
"                        other  : display a dialog window for input by user\n"+
"\n"+
"           *   Zoom y - Displays data at different scales along the y axis - \n"+
"                        default: the default scale is one pixel per datum.\n"+
"                        x 200 %: 200 % of the default size\n"+
"                        x 400 %: 400 % of the default size\n"+
"                        x 600 %: 600 % of the default size\n"+
"                        x 800 %: 800 % of the default size\n"+
"                        other  : display a dialog window for input by user\n"+
"\n"+
"           *   Find... - Displays a Find... dialog for the user to perform\n"+
"                        search actions in the upper frame of the viewer\n"+
"                        window that search for any sequences that the user\n"+
"                        specifies.\n\n"+
"           *   Find Again - Repeats the previous search action.\n\n"+
"           *   Refresh - Re-reads the files related to the current screen\n"+
"                        display and refresh the viewer window.\n\n"+
"           *   Show signal values - Displays the trace signal values at\n"+
"                        the bottom of the Viewer window.  The location of\n"+
"                        the mouse click and the strength of the trace\n"+
"                        signals at that X location are displayed.\n\n"+
"           *   List detected het/mixed bases - Displays the het/mixed bases bases\n"+
"                        detected by TraceTuner at the bottom of the viewer\n"+
"                        window.  The data are displayed in table format\n"+
"                        with four columns: index, het/mixed base, position,\n"+
"                        and QV, which means the index of the het/mixed\n"+
"                        base call, its code (R means G and A, Y means T\n"+
"                        and C, K means G and T, M means A and C, S means\n"+
"                        G and C, W means A and T.), its position, and its\n"+
"                        quality value, respectively.  Users can sort the\n"+
"                        data within the table by clicking on any of the\n"+
"                        four column headers.  By double-clicking on a data\n"+
"                        row, the user will be able to see the clicked \n"+
"                        het/mixed base marked with a black solid arrow in\n"+
"                        the display area of the Viewer window.\n"+
"\n" +
"    o  Feature - -\n"+
"           *   Display original (ABI) base calls - \n"+
"               Displays the original ABI base calls on top of the\n"+
"               TraceTuner base calls.  This feature is disabled if\n"+
"               'Display raw data' is selected since the ABI base calls\n"+
"               are pertained to analyzed data only.\n\n"+
"           *   Display raw data - \n"+
"               Displays the raw data from the sample file.  Selecting this\n"+
"               feature will automatically disable the selections of 'Align\n"+
"               with consensus/reference sequence', 'Display original(ABI)\n"+
"               base calls', 'Display alternative base calls', and 'Display\n"+
"               intrinsic peaks', since they are irrelevant to the raw\n"+
"               data.  This feature will be automatically disabled if the\n"+
"               sample file currently being displayed is in SCF format.\n\n"+
"           *   Advanced - \n"+
"               Displays a sub-menu that contains the following items:\n\n"+
"               *   Align with consensus/reference sequence - \n"+
"               Displays the alignment between the TraceTuner base calls\n"+
"               and the consensus.  The consensus sequence is shown at\n"+
"               the top with colored dots for exact literal matches between\n"+
"               the consensus base and the TT called base; base codes\n"+
"               for mismatches.  The index of the consensus bases is shown\n"+
"               above the index of the called bases.  The alignment data\n"+
"               is stored in the .tal file of the currently shown sample\n"+
"               file.  If this .tal file does not exist, or it does not\n"+
"               contain good alignment data, or 'Display raw data' is\n"+
"               selected, this feature is disabled.\n\n"+
"               *   Display alternative base calls - \n"+
"               Displays the alternative base calls along with their\n"+
"               quality values (displayed beneath the bases) above the\n"+
"               corresponding TraceTuner base calls which are wrapped in\n"+
"               orange rectangle boxes.  The alternative base calls are\n"+
"               stored in the .tab file of the currently shown sample file.\n"+
"               If this .tab file does not exist, or it contains no\n"+
"               TraceTuner's alternative base calls, or the 'Display raw\n"+
"               data' option is selected, this feature is disabled.\n\n"+
//  "               *   Display intrinsic peaks - \n"+
//  "               Displays the intrinsic peaks in dashed lines utilizing the\n"+
//  "               same color scheme as the chromatogram.  The intrinsic peaks\n"+
//  "               data is stored in the .tip file of the currently shown\n"+
//  "               sample file.  If this .tip file does not exist, or the\n"+
//  "               'Display raw data' is selected, this feature is disabled.\n\n"+
//  "               *   Display total intrinsic signal - \n"+
//  "               Displays the total intrinsic signal in dashed lines in the\n"+
//  "               same color scheme as the chromatogram.  The total intrinsic\n"+
//  "               signal data is calculated from the intrinsic peaks data\n"+
//  "               stored in the .tip file of the currently shown sample file.\n"+
//  "               If this .tip file does not exist, or the 'Display raw data'\n"+
//  "               is selected, this feature is disabled.\n\n"+
"               *   Display trimmed regions - \n"+
"               Displays the trimmed/bad quality regions in yellow back-\n"+
"               ground.  The trim threshold is shown at the right end of\n"+
"               the toolbar of the viewer window.  This feature is disabled\n"+
"               if no valid trim data found in the PHD file, or the \n"+
"               'Display raw data' option is selected.\n"+
"\n"+
"    o  Help - -\n"+
"           *   Info - Displays information about TraceTuner Viewer.\n\n" +
"\n"+
"Viewer has the following tool buttons (from the leftmost to the rightmost):\n"+
"\n"+
"    o  Open - -\n"+
"           *   The same effect as the Open option in the File menu.\n"+
"\n"+
"    o  Page Setup - -\n"+
"           *   The same effect as the Page Setup option in the File menu.\n"+
"\n"+
"    o  Print View - -\n"+
"           *   The same effect as the Print View option in the File menu.\n"+
"\n"+
"    o  Find... - -\n"+
"           *   The same effect as the Find... option in the View menu.\n"+
"\n"+
"    o  Find Again - -\n"+
"           *   The same effect as the Find Again option in the View menu.\n"+
"\n"+
"    o  Refresh - -\n"+
"           *   The same effect as the Refresh option in the View menu.\n"+
"\n"+
"    o  First Sample File - -\n"+
"           *   Displays the first sample file / PHD file pair along with\n"+
"               the available auxiliary files among the batch of the files\n"+
"               opened in the Viewer.\n"+
"\n"+
"    o  Previous Sample File - -\n"+
"           *   Displays the previous sample file / PHD file pair along\n"+
"               with the available auxiliary files among the batch of files\n"+
"               opened in the Viewer.\n"+
"\n"+
"    o  Next Sample File - -\n"+
"           *   Displays the next sample file / PHD file pair along with\n"+
"               the available auxiliary files among the batch of the files\n"+
"               opened in the Viewer.\n"+
"\n"+
"    o  Last Sample File - -\n"+
"           *   Displays the last sample file / PHD file pair along with\n"+
"               the available auxiliary files among the batch of the files\n"+
"               opened in the Viewer.\n"+
"\n"+
"    o  Info - -\n"+
"           *   The same effect as the Info option in the Help menu.\n\n\n";


    /** The help text for the Launcher window. */
    public final static String TTRUN_HELP =
"-----------------------------------------                     ------------\n"+
"Paracel TraceTuner " + Constants.VERSION + " Launcher Help                             June-2001\n"+
"-----------------------------------------                     ------------\n"+
"\n"+
"The TraceTuner Launcher is a Graphical User Interface (GUI) that allows a\n"+
"user to run TraceTuner without needing to use the MS-DOS command line.\n"+
"TraceTuner itself is a program that assigns accurate, calibrated quality\n"+
"values to a sample file from an ABI 3700 DNA sequencer.\n"+
"\n"+
"TraceTuner by default also recalls Ns and may adjust certain other base\n"+ 
"calls from the original ABI calls.\n"+
"\n"+
"The program input may be a single sample file or a list of sample files in\n"+
"a directory.  The selected files will be displayed in the Input Files box\n"+
"showing the directory path with the file extension.  TraceTuner will\n"+
"process all the selected sample files.  By default, TraceTuner outputs a\n"+
"PHD-format (.phd.1) file containing a list of updated base calls, the\n"+
"predicted quality values and the peak location associated with the base\n"+
"call.\n"+
"\n"+
"TraceTuner does not alter the original information contained in the ABI\n"+
"sample file.\n"+
"\n"+
"The .phd.1 files must be created or already available in order to view\n"+
"them using the Viewer.\n"+
"\n"+
"Some processing options are available which alter the behavior of\n"+
"TraceTuner.  Although it is generally recommended to use the defaults,\n"+
"advanced users can change output and processing options by clicking on the\n"+
"Change Options button.\n"+
"\n"+
"TraceTuner has 4 main operations, each invoked by a corresponding button.\n"+
"\n"+
"     o  Run -  Once an input file or a list of files has been specified,\n"+
"               click on the Run button to launch TraceTuner.  A list of\n"+
"               the names of the selected files will be saved in the output\n"+
"               directory in a file named 'ttfiles.txt'.  This list will be\n"+
"               used to run TraceTuner.  This file will automatically be\n"+
"               removed from the directory after the run is completed.  A\n"+
"               log message window will display messages from the run.  The\n"+
"               same log information is also saved in a file called\n"+ 
"               'ttlog.txt' in the output directory.  The log message from\n"+
"               each run are appended to 'ttlog.txt' file if the same\n"+
"               output directory is selected.\n\n"+
"               Once the run gets started, the Run button will be changed\n"+
"               to \"Stop\".  The user can click on the \"Stop\" button to\n"+
"               terminate the current run.  Once the \"Stop\" button is\n"+
"               clicked, the button will be switched back to \"Run\".\n"+
"\n"+
"     o  View - Launch TraceTuner's Viewer to look at the combined PHD file\n"+
"               output and traces from the original sample file and other.\n"+
"               available auxiliary files.  You can only view TraceTuner's\n"+
"               .phd.1 file output.\n"+
"\n"+
"     o  Help - Display this screen\n"+
"\n"+
"     o  Exit - Exit the launcher, any currently running TraceTuner jobs,\n"+
"               and any viewer displays.\n"+
"\n"+
"\n";

    /** The help text for the Basic Options window and the Advanced Options
	window. */
    public final static String OPT_HELP =
"----------------------                                        ------------\n"+
"Options Help                                                     June-2001\n"+
"----------------------                                        ------------\n"+
"\n"+
"Changing these options can alter the processing behavior and output formats\n"+
"used by TraceTuner.  In general, best results are assured by accepting the\n"+ 
"default settings.  However, some applications and environments may benefit\n"+
"from different settings.  The Basic Options window provides access to the\n"+
"basic options.\n\n"+
"TraceTuner can output six well known formats:\n"+
"\n"+
"    o   Phd files - PHD-format files with the names obtained by appending\n"+
"                    \".phd.1\" to the names of the input files.\n"+
"\n"+
"    o   SCF files - SCF-format files with the names obtained by replacing\n"+
"                    the \".ab1\" suffix of the names of the input files\n"+
"                    with the \".scf\" suffix.  If the input files have\n"+
"                    the \".scf\" suffix, the output file names are\n"+
"                    obtained by appending \".scf\" to the names of the\n"+
"                    input files.\n"+
"\n"+
"    o   Qual files (one file per read) - .qual-formatted quality value\n"+
"                    files with the names obtained by appending \".qual\"\n"+
"                    to the names of the input files.  (One .qual file is\n"+
"                    output for each input file.)\n"+
"\n"+
"    o   Qual file (single file for all reads) - A quality value file with\n"+
"                    a user-specified name.  The file contains the quality\n"+
"                    values of all the reads processed in this TraceTuner\n"+
"                    run.  Selecting this option will make visible an input\n"+
"                    field where the user can specify the file name.\n"+
"\n"+
"    o   Seq files (one file per read) - FastA-formatted sequence output\n"+
"                    files with the names obtained by appending \".seq\"\n"+
"                    to the names of the input files.  (One .seq file is\n"+
"                    output for each input file.)\n"+
"\n"+
"    o   Seq file (single file for all reads) - A sequence output file with\n"+
"                    a user-specified name.  The file contains the base\n"+
"                    calls of all the reads processed in this TraceTuner\n"+
"                    run.  Selecting this option will make visible an input\n"+
"                    field where the user can specify the file name.\n"+
"\n"+
"Also TraceTuner can output two additional formats, which regular users\n"+
"may not need to know at all.  Any combination of these files may be output\n"+
"although typically the PHD files are all that's needed.  The PHD files are\n"+
"required by the Viewer.  At least one output type must be selected.\n"+
"\n"+
"A customized calibration tool is available:\n"+
"\n"+
"    o  Calibration - The lookup table is used for TraceTuner in assigning\n"+
"                quality values.  TraceTuner has four generic built-in\n"+
"                lookup tables: ABI3730_POP7_BigDye, ABI3700_POP5_BigDye, \n"+
"                ABI3700_POP6_BigDye and MegaBACE. The ABI3730_POP7_BigDye was\n"+
"                calibrated for the ABI 3730 with Pop7 separation polymer\n"+
"                and BigDye terminator chemistry; The ABI3700_POP5_BigDye was\n"+
"                calibrated for the ABI 3700 with Pop5 and BigDye terminator\n"+
"                chemistry; the ABI3700_POP6_BigDye was calibrated for the\n"+
"                ABI 3700 with Pop6 and BigDye terminator chemistry; the\n"+
"                the MegaBACE was calibrated for Amersham's MegaBace data.\n"+
"                The first item listed in the pull down menu ('Automatic \n"+
"                Selection') allows TraceTuner to automatically select among\n"+
"                these four built-in lookup tables based on the type of\n"+
"                chemistry used to generate the sample file.  If the type\n"+
"                of chemistry can not be determined from the sample file,\n"+
"                the ABI3700_POP5_BigDye is used by default.\n\n"+
"                The menu then lists these four built-in lookup tables.\n"+
"                By selecting any one of them, the automatic lookup table\n"+
"                selection will be turned off and the selected built-in\n"+
"                lookup table will be used.\n\n"+
"                The menu also lists other customized or alternate ones, if\n"+
"                they are stored in the current directory with '.tbl' as the\n"+
"                file extension.  If the lookup table file is not in the\n"+
"                current directory or does not have '.tbl' as an extension,\n"+
"                select the last item 'Other...(specify)' from the menu.  This\n"+
"                will make visible another file selection line where you may\n"+
"                specify your selection.  Customized or alternate lookup tables\n"+
"                are available from Paracel at an additional charge.\n"+
"\n"+
"TraceTuner provides various program options:\n"+
"\n"+
"    o  No call - Check this box to turn off TraceTuner's own base calling\n"+
"                adjustments to the original ABI base calls.  When unchecked,\n"+
"                by default, TraceTuner will recall Ns using its best\n"+
"                estimate of what the base at that locations is, and it may\n"+
"                adjust base calls at other locations if its algorithms\n"+
"                suggest that this would improve the overall calling of the\n"+
"                sequence.  This option can not be used in conjunction\n"+
"                with Recall N or Call het/mixed bases.\n"+
"\n"+
"    o  Recall N - Check this box will force TraceTuner to recall Ns and do\n"+
"                not change any other base calls which it reads from sample\n"+
"                file.  This option can not be used in conjunction with\n"+
"                Nocall or Call het/mixed bases.\n"+
"\n"+
"    o  Call het/mixed bases - Check this box will force TraceTuner to\n"+
"                call het/mixed bases.  This optons can not be used in\n"+
"                conjunction with No Call or Recall N.\n"+
"\n"+
"    o  Specify Minimum Peak Height Ratio - Check this box will make visible\n"+
"                a text box where you may specify a value between 0 and 1.0\n"+
"                to force TraceTuner to use the specified value as the\n"+
"                threshold ratio of the current peak height to the called\n"+
"                peak height.  If the ratio is under this specified value,\n"+
"                the current peak will be considered as noise.  The default\n"+
"                value is "+Constants.DEFAULT_MIN_RATIO+".  This option is\n"+
"                used in combination with Call het/mixed bases.  Therefore\n"+
"                will only be enabled if Call het/mixed bases option\n"+
"                is checked.\n"+
"\n"+
"    o  Specify Output Directory - Checking this box will make visible\n"+
"                another file selection line where you may specify a\n"+
"                particular directory for TraceTuner's output to be sent.\n"+
"                The default is to output to the same directory containing\n"+
"                the sample files.\n"+
"\n"+
"Selections in the option window can be saved by clicking on the Save button.\n"+
"The saved setting will be displayed next time when Change Option is selected.\n"+
"The Cancel button will revert back to what is previously saved.  Click on the\n"+
"Reset Defaults button will put the setting back to default configuration.\n"+
"\n"+
"Click on the Advanced Options button will open an Advanced Options window\n"+
"which allows advanced user to select more output formats:\n"+
"\n"+
"    o   Tab (TraceTuner alternative base calls) files with a .tab suffix.\n"+
"\n"+
"    o   Tal (TraceTuner alignment) files with a .tal suffix.\n"+
"\n"+
//  "    o   Tip (TraceTuner intrinsic peaks) files with a .tip suffix.\n"+
//  "\n"+
"Any combination of these files must be output along with .phd.1 format\n"+
"and any combination of the other six output formats included in the Basic\n"+
"Options window.\n\n"+
"Also six additional program options are available in the Advanced Options\n"+
"window:\n"+
"\n"+
"    o  Edited bases - Check this box will force TraceTuner to read edited\n"+
"                base calls and locations from sample file and start from\n"+
"                them when recalling bases.  The default is to read and\n"+
"                start from called bases and locations.\n"+
"\n"+
"    o  Perform mobility shift correction - Check this box will force\n"+
"                TraceTuner to perform mobility shift corrections when\n"+
"                recalling bases.  The default is to perform no mobility\n"+
"                shift correction.\n"+
"\n"+
"    o  Specify trim window size - Check this box will make visible a text \n"+
"                box where you may specify a positive integer value to\n"+
"                force TraceTuner to use the specified value as the size\n"+
"                of the sliding window of bases used to trim the beginning\n"+
"                and the end of a sequence.  The default value is "+
Constants.DEFAULT_TRIM_WINDOW_SIZE +".\n"+
"                The trimming stops when the average quality value of\n"+
"                bases in the window is equal to or greater than the trim\n"+
"                threshold value.\n"+
"\n"+
"    o  Specify trim threshold - Check this box will make visible a text \n"+
"                box where you may specify an integer value between 1 and\n"+
"                100 to force TraceTuner to use the specified value as the\n"+
"                threshold for average quality value of bases in the\n"+
"                sliding window.  The default value is "+
Constants.DEFAULT_TRIM_THRESHOLD+".\n"+
"\n"+
//  "    o  Specify Context Table File - Checking this box will make visible\n"+
//  "                a file selection line where you may specify a context\n"+
//  "                table file for TraceTuner to call bases taking the context\n"+
//  "                effect into account.\n"+
//  "\n"+
"    o  Specify Consensus/Reference Sequence File - Checking this box will\n"+
"                make visible a file selection line where you may specify\n"+
"                a consensus/reference sequence file for TraceTuner to\n"+
"                align with the TraceTuner base calls and generate .tal\n"+
"                files.\n"+
"\n"+
"The Advanced Options window also contains three buttons.  By clicking on\n"+
"the Ok button, selections in the Advanced Options window will be validated\n"+
"and remembered.  (Note: These selections will not be saved permanently,\n"+
"until the user clicks the Save button in the Basic Options window.)  The\n"+
"Cancel button will revert back to what was previously saved.  Clicking the\n"+
"Reset Defaults button will put the setting back to default configuration.\n"+
"\n\n";

    /** The help text for the Find... window. */
    public final static String FIND_HELP =
"----------------------                                        ------------\n"+
"Find... Help                                                     June-2001\n"+
"----------------------                                        ------------\n"+
"\n"+
"The user can specify the following attributes to perform a search task in\n"+
"the upper frame of the current viewer window:\n\n"+
"    o Find What - The user can either select the \"Next het/mixed\"\n"+
"                  option, or select the \"Specify a sequence\" and type\n"+
"                  one or more base sequences that he/she want to find\n"+
"                  in the text field.  Find performs literal match.  It is\n"+
"                  NOT case sensitive.  The sequences can be one or more\n"+
"                  bases long.  More than one sequences can be specified\n"+
"                  by separating them with white spaces, ',', or ';'.\n\n"+
"    o Search -    Drop-down box in which to select the Find scope/subject.\n"+
"                  Three options are provided:\n"+
"                  * TraceTuner base calls\n"+
"                  * Original (ABI) base calls\n"+
"                  * Consensus/reference sequence.\n\n"+
"    o Direction - Buttons to select the search direction.\n\n"+
"If any of the specified sequences is found, each one of the bases will be\n"+
"marked with a black solid arrow above it in the upper frame of the current\n"+
"viewer window.\n\n";

    /** The Options help window. */
    private static JFrame optHelp;

    /** The Launcher help window. */
    private static JFrame ttrunHelp;

    /** The Viewer help window. */
    private static JFrame ttviewHelp;

    /** The Viewer Find... help window. */
    private static JFrame findHelp;

    /** Returns the Options help window. */
    public static JFrame getOptionsHelp() {
	if (optHelp == null) {
	    optHelp = createHelpFrame("Options Help", OPT_HELP);
	}
//  	System.out.println(OPT_HELP);
	return optHelp;
    }

    /** Returns the Launcher help window. */
    public static JFrame getTTRunHelp() {
	if (ttrunHelp == null) {
	    ttrunHelp = createHelpFrame("TraceTuner Launcher Help", 
					TTRUN_HELP);
	}
//  	System.out.println(TTRUN_HELP);
	return ttrunHelp;
    }

    /** Returns the Viewer help window. */
    public static JFrame getTTViewHelp() {
	if (ttviewHelp == null) {
	    ttviewHelp = createHelpFrame("TraceTuner Viewer Help", 
					 TTVIEW_HELP);
	}
//  	System.out.println(TTVIEW_HELP);
	return ttviewHelp;
    }

    /** Returns the Find... help window. */
    public static JFrame getFindHelp() {
	if (findHelp == null) {
	    findHelp = createHelpFrame("TraceTuner Find... Help", 
				       FIND_HELP);
	}
//  	System.out.println(FIND_HELP);
	return findHelp;
    }

    /** Creates a frame with the specified title and the specified 
	content text. */
    protected static JFrame createHelpFrame(String title, String content) {
	final JFrame f = new JFrame(title);
	JTextArea text = new JTextArea(content);
	text.setEditable(false);
	text.setFont(Fonts.TEXT_FONT);

	JScrollPane sp = new JScrollPane(text);
	sp.setPreferredSize(new Dimension(600, 300));
	
	f.getContentPane().add(sp);
	f.pack();

	f.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {
		f.setVisible(false);
	    }
	});

	return f;
    }
}
