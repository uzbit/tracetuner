
/* This program creates data file for plotting           
// the color data using xgraph 
*/
#include <jni.h>
#include "Visualize.h"

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ABI_Toolkit.h"
#include "FileHandler.h"
#include "Btk_qv.h"
#include "util.h"

#define MAX(a,b)  (a>b)?a:b
#define MIN(a,b)  (a<b)?a:b
#define DEBUG 1
#define INF 1000000
#define NUMG 0
#define NUMA 1
#define NUMT 2
#define NUMC 3
#define MAXPATHLEN 100
#define BUFLEN 1000


/* Function prototypes:  
   --------------------
*/
int process(char* filename, BtkMessage* message, int InputPhD, int OutputCaml, int OutputXgraph);

int read_nums(int* num_called_bases, int* num_datapoints);

int read_arrays(int num_called_bases, char* called_bases, int* called_peak_locs,
    int num_datapoints, int** chromatogram, char* color2base, BtkMessage* message);

int print_to_caml  (char* file_name, int data_begin, int data_end,
                          int num_datapoints, int** chromatogram, 
                          int num_called_bases, char* called_bases,
                          int* qualities, int* called_peak_locs,
                          char* color2base, BtkMessage* message);

int print_to_xgraph(char* file_name, int data_begin, int data_end,
                          int num_datapoints, int** chromatogram, 
			  int num_called_bases, char* called_bases, int* called_peak_locs,
                          char* color2base, BtkMessage* message);

int release(char* called_bases, int* called_peak_locs, int* chromatogram[NUM_COLORS],
                 int* quality_values);

int Btk_input_phd_file( char *file_name, char *path,
    char *called_bases, int *qualities, int *called_locs, int num_bases, 
    BtkMessage *message);


int process(char* filename , BtkMessage* message, int InputPHD, int OutputCaml, int OutputXgraph)
{
   int            r;
   void*          p;
   long           fileSize;
   int            num_called_bases; 
   char*          called_bases; 
   int*           called_peak_locs;
   int            num_datapoints; 
   int*           chromatogram[NUM_COLORS]; 
   int*           quality_values;
   char           color2base[NUM_COLORS];
   int            data_begin, data_end;

// Store command line arguments as integer variables
   data_begin=0;
   data_end  =INF;

   if (filename > 0)
   {
     printf("Mark 1\n");
      r = F_Open(filename, &p, &fileSize);
          if (r == ERROR) return r;
     printf("Mark 2\n");
      r = read_nums(&num_called_bases, &num_datapoints);
          if (r == ERROR) return r;

	fprintf(stderr, "File - %s\n", filename);
	fprintf(stderr, "Number of called bases=%d\n Number of datapoints=%d\n",
	        num_called_bases, num_datapoints); 

// Allocate memory for data arrays
   called_bases = (char*) malloc(sizeof(char)*num_called_bases);
   MEM_ERROR(called_bases);
 
   called_peak_locs = (int*) malloc(sizeof(int)*num_called_bases);
   MEM_ERROR(called_peak_locs);
 
   quality_values = (int*) malloc(sizeof(int)*num_called_bases);
   MEM_ERROR(quality_values);
    (void)memset(quality_values, 0, sizeof(int)*num_called_bases);

   r = read_arrays(num_called_bases, called_bases, called_peak_locs,
       num_datapoints, chromatogram, color2base, message);
       if (r == ERROR) return r;       

    if( InputPHD ) {
	r = Btk_input_phd_file( filename, NULL, called_bases,
	    quality_values, called_peak_locs, num_called_bases, message);
	if (r == ERROR) return r;       
    }

    if( OutputCaml ) {
	r = print_to_caml(filename, data_begin, data_end, 
	    num_datapoints, chromatogram, num_called_bases, called_bases,
	    quality_values, called_peak_locs, color2base, message);
	if (r == ERROR) return r;       
    }

    if( OutputXgraph ) {
	r = print_to_xgraph(filename, data_begin, data_end, 
	    num_datapoints, chromatogram, num_called_bases, called_bases,
	    called_peak_locs, color2base, message);
	if (r == ERROR) return r;       
    }

   r = release(called_bases, called_peak_locs, chromatogram, quality_values);
       if (r == ERROR) return r;          

   F_Close(p);
   }

error:

   return SUCCESS;
}
// ---------------------------------------------------------
int read_nums(int* num_called_bases, int* num_datapoints)
{
   ABIError   r;
   int        j;
   long       num_bases;
   long       num_peaks;
   long       num_points;
   short      dye_number;
   short      lane=0;

   *num_datapoints = 100000;
 
// Read the number of called bases
   r = ABI_NumCalledBases(&num_bases);
       if (r != kNoError) return (int)r;
   *num_called_bases = (int)num_bases;
 
// Read the number of called peak locations
   r = ABI_NumCalledPeakLocations(&num_peaks);
       if (r != kNoError) return (int)r;
   if (num_peaks != num_bases)
   {
      printf("Number of bases != number of peaks\n");
//    return ERROR;
   }
// printf("Number of called bases = %d\n", (int)num_bases);

// Choose any color number
   for (j=0; j < NUM_COLORS; j++)
   {
      dye_number = j+1;
 
// Read the number of data points for this color
      r = ABI_NumAnalyzedData(lane, dye_number, &num_points);
      if (r != kNoError) return (int)r;
      if (*num_datapoints > (int)num_points)
          *num_datapoints = (int)num_points;
   }
// printf("Number of datapoints = %d\n", (int)num_points);
   
   return SUCCESS;
}   
// --------------------------------------------------
int read_arrays(int num_called_bases, char* called_bases, int* called_peak_locs,
                int num_datapoints, int** chromatogram,
                    char* color2base, BtkMessage * message)
{
   short*     peak_locs;
   short*     analyzed_data;
   char       base;
   int i,j;
   short      dye_number;
   short      lane=0;
   ABIError   r;

// Read the called bases
// r = ABI_CalledBases(called_bases);
   r = ABI_EditedBases(called_bases);
       if (r != kNoError) return (int)r;
   for (i=0; i<num_called_bases; i++)
   {
/*    if (i<2 || i>num_called_bases-3)
      {
         base = called_bases[i];
         printf("In main: called_bases[%d]=",i);
         putchar(base);
         printf("\n");
      }
*/
   }
   
// Create an auxiliary array of called peak locations
   peak_locs = (short*) malloc(sizeof(short)*num_called_bases);
   MEM_ERROR(peak_locs);
 
// Read the called peak locations and free the auxiliary array
// r = ABI_CalledPeakLocations(peak_locs);
   r = ABI_EditedPeakLocations(peak_locs);
       if (r != kNoError) return (int)r;
   for (i=0; i<num_called_bases; i++)
   {
      called_peak_locs[i] = (int)peak_locs[(long)i];
/*    if (i<2 || i>num_called_bases-3)
         printf("In main: called_peak_locs[%d]=%d\n",i,called_peak_locs[i]);
*/
   }
   free(peak_locs);
   peak_locs = NULL;

// Choose any color number
   for (j=0; j < NUM_COLORS; j++)
   {
      dye_number = j+1;

// Which base corresponds to the selected dye number?
      r = ABI_DyeIndexToBase(dye_number, &base);
      if (r != kNoError) return (int)r;
      color2base[j] = base;
 
// Create the major and auxiliary array of data points
      analyzed_data = (short*) malloc(sizeof(short)*num_datapoints);
      MEM_ERROR(analyzed_data);
      chromatogram[j] = (int*)malloc(sizeof(int)*(int)num_datapoints);
      MEM_ERROR(chromatogram);

// Read the chromatograms and free the auxiliary array
      r = ABI_AnalyzedData(lane, dye_number, analyzed_data);
          if (r != kNoError) return (int)r;
      for (i=0; i<num_datapoints; i++)
      {
         chromatogram[j][i] = (int)analyzed_data[(long)i];
//       if (i<20 || i>num_datapoints-3)
//       printf("In main: chromatogram[%d][%d]=%d\n",j,i,chromatogram[j][i]);
  
      }
      free (analyzed_data);
      analyzed_data = NULL;
   }

error:

   return SUCCESS;

}

// ------------------------------------------------------------------
int print_to_xgraph(char* file_name, int data_begin, int data_end, 
        int num_datapoints, int** chromatogram, 
	int num_called_bases, char* called_bases, int* called_peak_locs,
        char* color2base, BtkMessage* message)
{
   int i, j;
   int x, y;
   char suffix[10] = ".data";
   char *out_name;                        /* name of output file */
   FILE* xgr_out;
   int name_length;
   char  base;
 
// Create the name of xgraph file by appending
// suffix ".data" to the sample file name
   name_length = (int)strlen(file_name)+(int)strlen(suffix);
   out_name = (char*)malloc(name_length+1);
   out_name = strcpy(out_name, file_name);
   out_name = strcat(out_name,suffix);
   xgr_out = fopen(out_name,"w");
 
// Output the data    
   fprintf(xgr_out,"TitleText: colordata for sample %s\n",
           file_name);      
   data_end = MIN(data_end,num_datapoints-1);
   for (j=0; j<4; j++)
   {
      base = color2base[j];
      fprintf(xgr_out,"\n\"Color %d (%s)\n",j,&base);
      for (i = data_begin; i <= data_end; i++)
      {
         fprintf(xgr_out,"%d   %d   \n",i,chromatogram[j][i]);
      }
   }   

   fprintf(xgr_out,"\n\"Called Peaks\n");

   for (i = 0; i < num_called_bases; i++)
   {
	base = called_bases[i];
	switch (base) {
		case 'A': j=NUMA; break;
		case 'C': j=NUMC; break;
		case 'T': j=NUMT; break;
		case 'G': j=NUMG; break;
		default: j=-1;
	}
	x= called_peak_locs[i];
	if ( x > data_end ) x = data_end;
	if ( x < data_begin ) x = data_begin;
	if ( j > -1 )
		y= chromatogram[j][x];
	else
		y=0;
	fprintf(xgr_out, "move %d   %d   \n", x-1, y-5);
	fprintf(xgr_out, "%d   %d   \n", x-1, y+5);
	fprintf(xgr_out, "%d   %d   \n", x+1, y+5);
	fprintf(xgr_out, "%d   %d   \n", x+1, y-5);
	fprintf(xgr_out, "%d   %d   \n", x-1, y-5);
   }

   /* Put color information at bottom of file. */
    for(j=0; j<4 && color2base[j]!='G'; j++)
    { }
   fprintf (xgr_out, "%d.Color: black\n", j);
    for(j=0; j<4 && color2base[j]!='A'; j++)
    { }
   fprintf (xgr_out, "%d.Color: green\n", j);
    for(j=0; j<4 && color2base[j]!='T'; j++)
    { }
   fprintf (xgr_out, "%d.Color: red\n",   j);
    for(j=0; j<4 && color2base[j]!='C'; j++)
    { }
   fprintf (xgr_out, "%d.Color: blue\n",  j);

    free(out_name);
    fclose(xgr_out);

   return SUCCESS;
}

// ------------------------------------------------------------------
int print_to_caml(char* file_name, int data_begin, int data_end, 
        int num_datapoints, int** chromatogram, 
	int num_called_bases, char* called_bases, int* qualities,
        int* called_peak_locs, char* color2base, BtkMessage* message)
{
   int i, j;
   char suffix[10] = ".caml";
   char *out_name;                        /* name of output file */
   FILE* xgr_out;
   int name_length;
   char *no_path;
 
// Create the name of xgraph file by appending
// suffix ".data" to the sample file name
   name_length = (int)strlen(file_name)+(int)strlen(suffix);
   out_name = (char*)malloc(name_length+1);
   out_name = strcpy(out_name, file_name);
   out_name = strcat(out_name,suffix);
   xgr_out = fopen(out_name,"w");
   if( (no_path = strrchr(file_name, '/')) == NULL )
   {
       no_path=file_name;
   } else {
       no_path++;
   }
 
// Output the data    
/* Header information. */
    fprintf(xgr_out, "<?xml version=\"1.0\" encoding=\"UTF-8\" ");
    fprintf(xgr_out, "standalone=\"no\"?>\n");
    fprintf(xgr_out, "<CAML>\n");
    fprintf(xgr_out, "<ALIGNMENT_SET SOURCE=\"qv\" ALIGNMENTS=\"1\">\n");
    fprintf(xgr_out, "<SET_SUMMARY>\n");
    fprintf(xgr_out, "<TITLES>%s</TITLES>\n", no_path);
    fprintf(xgr_out, "<NUM_ROWS>1</NUM_ROWS>\n");
    fprintf(xgr_out, "</SET_SUMMARY>\n");
    fprintf(xgr_out, "\n");

    fprintf(xgr_out, "<MULTI_ALIGNMENT TITLE=\"%s\" ROWS=\"1\">\n", no_path);
    fprintf(xgr_out, "\n");

    fprintf(xgr_out, "<ROW TITLE=\"%s\" ORIENT=\"FORWARD\" OFFSET=\"0\">\n",
            no_path);
    fprintf(xgr_out, "\n");

/* Sequence. */
    fprintf(xgr_out, "<MA_SEQUENCE LENGTH=\"%d\">", num_called_bases);
    for (i = 0; i < num_called_bases; i++)
    {
        fprintf(xgr_out, "%c", called_bases[i]);
    }
    fprintf(xgr_out, "</MA_SEQUENCE>\n");
    fprintf(xgr_out, "\n");

/* Qualities. */
    fprintf(xgr_out, "<MA_QUALITY LENGTH=\"%d\">", num_called_bases);
    for (i = 0; i < num_called_bases-1; i++)
    {
        fprintf(xgr_out, "%03d ", qualities[i]);
    }
    fprintf(xgr_out, "%03d", qualities[num_called_bases-1]);
    /* No trailing space.*/

    fprintf(xgr_out, "</MA_SEQUENCE>\n");
    fprintf(xgr_out, "\n");

/* X-Coordinates. */
    fprintf(xgr_out, "<CHROMATOGRAM NUM_TRACES=\"4\" ");
    fprintf(xgr_out, "COLOR_MAP=\"ACGT\" MAX_Y_VALUE=\"1600\">\n");
    data_end = MIN(data_end,num_datapoints-1);
    fprintf(xgr_out, "<COORDINATES LENGTH=\"%d\">", num_called_bases);
    for (i = 0; i < num_called_bases-1; i++)
    {
        fprintf(xgr_out,"%04d ", called_peak_locs[i]);
    }
    fprintf(xgr_out,"%04d", called_peak_locs[num_called_bases-1]);
    /* No trailing space.*/

    fprintf(xgr_out, "</COORDINATES>\n");
    fprintf(xgr_out, "\n");

/* Y-Coordinates. */
    for(j=0; j<4 && color2base[j]!='A'; j++)
    { }
    fprintf(xgr_out, "<TRACE COLOR_NUMBER=\"9\" LENGTH=\"%d\">", 
            data_end-data_begin+1);
    for (i = data_begin; i < data_end; i++)
    {
        fprintf(xgr_out,"%04d ", chromatogram[j][i]);
    }
    fprintf(xgr_out,"%04d", chromatogram[j][data_end]); /* No trailing space.*/
    fprintf(xgr_out, "</TRACE>\n");
    fprintf(xgr_out, "\n");

    for(j=0; j<4 && color2base[j]!='C'; j++)
    { }
    fprintf(xgr_out, "<TRACE COLOR_NUMBER=\"10\" LENGTH=\"%d\">", 
            data_end-data_begin+1);
    for (i = data_begin; i < data_end; i++)
    {
        fprintf(xgr_out,"%04d ", chromatogram[j][i]);
    }
    fprintf(xgr_out,"%04d", chromatogram[j][data_end]); /* No trailing space.*/
    fprintf(xgr_out, "</TRACE>\n");
    fprintf(xgr_out, "\n");

    for(j=0; j<4 && color2base[j]!='G'; j++)
    { }
    fprintf(xgr_out, "<TRACE COLOR_NUMBER=\"11\" LENGTH=\"%d\">", 
            data_end-data_begin+1);
    for (i = data_begin; i < data_end; i++)
    {
        fprintf(xgr_out,"%04d ", chromatogram[j][i]);
    }
    fprintf(xgr_out,"%04d", chromatogram[j][data_end]); /* No trailing space.*/
    fprintf(xgr_out, "</TRACE>\n");
    fprintf(xgr_out, "\n");

    for(j=0; j<4 && color2base[j]!='T'; j++)
    { }
    fprintf(xgr_out, "<TRACE COLOR_NUMBER=\"12\" LENGTH=\"%d\">", 
            data_end-data_begin+1);
    for (i = data_begin; i < data_end; i++)
    {
        fprintf(xgr_out,"%04d ", chromatogram[j][i]);
    }
    fprintf(xgr_out,"%04d", chromatogram[j][data_end]); /* No trailing space.*/
    fprintf(xgr_out, "</TRACE>\n");
    fprintf(xgr_out, "\n");

/* Ending information. */
    fprintf(xgr_out, "</CHROMATOGRAM>\n");
    fprintf(xgr_out, "</ROW>\n");
    fprintf(xgr_out, "</MULTI_ALIGNMENT>\n");
    fprintf(xgr_out, "</ALIGNMENT_SET>\n");
    fprintf(xgr_out, "</CAML>\n");
    fprintf(xgr_out, "\n");

    /* Clean up. */
    free (out_name);
    fclose(xgr_out);

    return SUCCESS;
}


// ------------------------------------------------------------------
int release(char* called_bases, int* called_peak_locs, int** chromatogram,
                 int* quality_values)
{
   int i;
 
   free(called_bases);
   called_bases = NULL;
   free(called_peak_locs);
   called_peak_locs = NULL;
   free(quality_values);
   quality_values = NULL;
 
   for (i=0; i<NUM_COLORS; i++)
   {
      free(chromatogram[i]);
      chromatogram[i] = NULL;
   }
 
return SUCCESS;
}


int
Btk_input_phd_file(
    char *file_name,
    char *path,
    char *called_bases, 
    int *qualities,
    int *called_locs,
    int num_bases,
    BtkMessage *message
    )
{
    char *seq_name, phd_file_name[MAXPATHLEN];
    FILE *phd_out;
    int i; 
    char buffer[BUFLEN];
    char c;
    int qv, pos;


    /* Use the name of the sample file, sans path, as the sequence name */
    if ((seq_name = strrchr(file_name, '/')) != NULL) {
	seq_name++;
    }
    else {
	seq_name = file_name;
    }

    /*
     * Create the name of phd file by appending suffix ".phd.1" to the sample
     * file name.
     */
    if (path != NULL) {
	(void)sprintf(phd_file_name, "%s/%s.phd.1", path, seq_name);
    }
    else {
	/* put it in the current dir */
	(void)sprintf(phd_file_name, "%s.phd.1", seq_name);
    }
    if ((phd_out = fopen(phd_file_name, "r")) == NULL) {
	(void)sprintf(message->text, "No such file: %s", phd_file_name);
	return ERROR;
    }

    while ( ( fgets(buffer, BUFLEN, phd_out) != NULL) 
	    && (strcmp(buffer, "BEGIN_DNA\n")!= 0 ) ) { 
    }

    i=0;
    while ( ( fgets(buffer, BUFLEN, phd_out) != NULL) 
	    && (strcmp(buffer, "END_DNA\n")!= 0 ) && (i < num_bases) ) { 
        sscanf( buffer, "%c %d %d\n", &c, &qv, &pos );
	if(islower((int)c)) c=toupper((int)c); /* always use upper */
	called_bases[i]=c;
	qualities[i]=qv;
        called_locs[i]=pos;
	i++;
    }

    (void)fclose(phd_out);
    return SUCCESS;
}


JNIEXPORT void JNICALL Java_Visualize_visualize(JNIEnv *env, jobject obj, jstring filename)
{
    int r, InputPHD, OutputCaml, OutputXgraph;
    BtkMessage message;
    const char* mystring;
    jboolean isCopy;
    mystring=(*env)->GetStringUTFChars(env,filename, &isCopy);

    InputPHD=1;
    OutputCaml=1;
    OutputXgraph=0;

    printf("In C, Filename=%s\n",mystring);

    r = process((char*) mystring, &message, InputPHD, OutputCaml, OutputXgraph);
        if(r==ERROR){
	    fprintf(stderr,"error:%s\n",message.text);
	    exit(-1);
        }
	if (isCopy==JNI_TRUE){
	  (*env)->ReleaseStringUTFChars(env,filename, mystring);
	}
}
