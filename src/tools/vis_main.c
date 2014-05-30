/* main.c */

/* This program creates data file for plotting           
// the color data using xgraph 
*/

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
int process(char* filename, BtkMessage* message, int InputPHD, int OutputCaml, int  OutputXgraph);

static void usage(char** argv)
{
    fprintf(stderr, "usage: %s -[pcx] sample_file_name \n",argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "p = input phd.1 file\n");
    fprintf(stderr, "c = output caml file\n");
    fprintf(stderr, "x = output xgraph file\n");
}




int main(int argc, char** argv)
{

static int InputPHD;
static int OutputCaml;
static int OutputXgraph;

    int r, i;
    BtkMessage message;

    if(argc<2){
	usage(argv);
	exit(2);
    }

    /* Set defaults. */
    InputPHD     = 0;
    OutputCaml   = 0;
    OutputXgraph = 0;

    while ( (i=getopt(argc, argv, "pcx")) != EOF) {
	switch(i) {
	case 'p':
	    InputPHD++;
	    break;
	case 'c':
	    OutputCaml++;
	    break;
	case 'x':
	    OutputXgraph++;
	    break;
	default:
	    usage(argv);
	    exit(2);
	}
    }

    if ( !InputPHD && !OutputCaml && !OutputXgraph ) {
	(void)fprintf(stderr, "Error: No input or output specified.\n");
	usage(argv);
	exit(2);
    }

    /* Turn line buffering off, so that someone watching can see progress. */
    (void) setbuf(stderr, NULL);

    /* Loop for each sample file. */
    for(i = optind; i< argc; i++)
    {
        r = process(argv[i], &message, InputPHD, OutputCaml, OutputXgraph);
        if(r==ERROR){
	    fprintf(stderr,"error:%s\n",message.text);
	    exit(-1);
        }
    }

    exit(0);
}

