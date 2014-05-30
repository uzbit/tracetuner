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

static const char *CVSREV = "CreateLicense.c : 2.1";
static void dummyfunc(const char *a) { dummyfunc(CVSREV); }

#include <string.h>
#include <unistd.h>   // for getpid()

#include "Btk_qv.h"
#include "LicenseKey.h"

const char *USAGE = "
Usage: crlicense productname never|unlimited|MM/DD/YYYY|+days

Creates a license key with the specified parameters and writes the 
key to STDOUT. The license key will be created for the specified
product name, which must match the product name in the source code that
checks the license file. 

Licenses can be created to expire:
    o  never (unlimited)
    o  on a specific date (MM/DD/YYYY)
    o  after a specified number of days (+days)

";


/*
  Reads in an expiration string in the form of MM/DD/YYYY or +days
  and sets the expiretime argument to be the time_t representation
  of that expiration date. Returns 1 on success, 0 on error (e.g.
  invalid format). 
*/
int ProcessExpireString(const char *str, time_t *expiretime)
{
  int ret = 1;
  int a, b, c;
  time_t tmp;

  // read in MM/DD/YYYY. Set expiretime to be the last second of that 
  // date. This has to do error checking to make sure that the input 
  // values are valid.
  if (sscanf(str, "%d/%d/%d", &a, &b, &c) == 3)
  {
    struct tm exp_tm;
    memset(&exp_tm, 0, sizeof(struct tm));
    exp_tm.tm_sec = 59;
    exp_tm.tm_min = 59;
    exp_tm.tm_hour = 23;
    exp_tm.tm_mday = b;
    exp_tm.tm_mon = a - 1;
    exp_tm.tm_year = c - 1900;

    tmp = mktime(&exp_tm);
    if (tmp > 0)
    {
      *expiretime = tmp;
    }
    else
    {
      ret = 0;
    }
  }
  // read in +days. Set expiretime to be current time plus that number 
  // of days 
  else if (sscanf(str, "+%d", &a) == 1)
  {
    *expiretime = (time(NULL) + (86400 * a));
  }
  // invalid format
  else
  {
    ret = 0;
  }
  
  return ret;
}


int main(int argc, char *argv[])
{
  int pid = getpid();
  long seed;
  char product[PRODNAMELEN];
  const char *expirestr = argv[2];
  int ret = 0;

  pid += (pid << 15);
  seed = (time(NULL) * pid);
  srand48(seed);

  if (argc != 3)
  {
    if ((argc == 2) && (!strcmp(argv[1], "-v")))
    {
      fprintf(stderr, "crlicense %s\n", TT_VERSION);
      exit(1);
    }
    else
    {
      fputs(USAGE, stderr);
      exit(1);
    }
  }

  strncpy(product, argv[1], PRODNAMELEN);

  // if this is an expiring license...
  if (strcmp(expirestr, "unlimited") && strcmp(expirestr, "never"))
  {
    time_t expiretime;

    // parse the expiration time
    if (!ProcessExpireString(expirestr, &expiretime))
    {
      fprintf(stderr, "ERROR: Could not process expiration string: "
	      "%s\n", expirestr);
      ret = 1;
    }
    else
    {
      // write out the expiring license file
      if (WriteLicenseFile(product, 0, expiretime, stdout))
      {
	const int timebufsize = 200;
	char timebuf[200];
      
	strftime(timebuf, timebufsize, "%b %d, %Y %H:%M:%S", 
		 localtime(&expiretime));
	
	fprintf(stderr, "Created license for %s that expires on %s.\n",
	       product, timebuf);
	ret = 0;
      }
      else
      {
	fprintf(stderr, "ERROR: Failed to create license for %s.\n", 
		product);
 	ret = 1;
      }
    }
  }
  // else this is a non-expiring license
  else
  {
    // write out the non-expiring license
    if (WriteLicenseFile(product, 1, time(NULL), stdout))
    {
      fprintf(stderr, "Created unlimited license for %s.\n", product);
      ret = 0;
    }
    else
    {
      fprintf(stderr, "ERROR: Failed to create unlimited license for %s.\n", 
	      product);
      ret = 1;
    }
  }

  exit(ret);
}
