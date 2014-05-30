
/**************************************************************************
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
 *
 * Copyright (c) 1999 Paracel Inc.  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

/*
 * CheckLicense.c    2.1
 */

#include "LicenseKey.h"
#include <string.h>
#include "Btk_qv.h"

const char *USAGE = "
Usage: chklicense licensefile

Opens specified license file and prints out the information contained
therein. Returns 0 if the file was read correctly, 1 otherwise.

";


int main(int argc, char *argv[])
{
    const char *fname = argv[1];
    FILE *licF = fopen(fname, "r");
    int ret;
    int noexpire;
    time_t timestamp;
    char product[PRODNAMELEN];

    if ((argc != 2) || !strncmp(argv[1], "-h", 2))
    {
        fputs(USAGE, stderr);
        exit(1);
    }
    else if ((argc == 2) && (!strcmp(argv[1], "-v")))
    {
        fprintf(stderr, "chklicense %s\n", TT_VERSION);
        exit(1);
    }

    /* open up the specified input license file */
    if (!licF)
    {
        fprintf(stderr, "ERROR: Unable to open license file %s. \n", 
            fname);
        exit(1);
    }

    /* read the license file */
    if (ReadLicenseFile(product, &noexpire, &timestamp, licF))
    {
        printf("License file:      %s\n", fname);
        if (noexpire)
        {
            printf("Product:           %s\n", product);
            printf("License expires:   never\n");
        }
        else
        {
            const int timebufsize = 200;
            char timebuf[200];
      
            strftime(timebuf, timebufsize, "%b %d, %Y %H:%M:%S", 
	       localtime(&timestamp));

            printf("Product:           %s\n", product);
            printf("License expires:   %s\n", timebuf);
        }
        ret = 0;
    }
    else
    {
        printf("ERROR: Invalid or corrupt license file: %s\n", fname);
        ret = 1;
    }
    fclose(licF);
    exit(ret);
}
