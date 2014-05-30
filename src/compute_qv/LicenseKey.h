/**************************************************************************
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
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
 * Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
 *
 * 2.4 2003/11/06 18:57:08
 */

#ifndef _LICSENSEKEY_H_
#define _LICSENSEKEY_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>

#if 1
#define NEW_LICENSE_FORMAT
#endif

#define LICENSE_VERSION "2.1"

/* maximum length of product name field */
#define PRODNAMELEN (200)

/* maximum length of product version field */
#define VERSIONLEN (20)

/* maximum length of hardware MAC address field */
#define MACLEN (20)

/* maximum length of OS field */
#define OSLEN (10)

/* Return codes for license functions that use them. */
#define LICENSE_SUCCESS          (0)
#define LICENSE_EXCEPTION_THROWN (-1)


typedef int LICENSE_RESULT;     /* Result of function calls. */

/* Definitions for license key fields. */

#define LICENSE_EXCEPTION_MESSAGELENGTH   250

/* Exception field for return information */

typedef enum {
    /* These are equal to the corresponding BTK values, where possible. */
    LICENSE_EXCEPTION_NONE             = 0x0000,
    LICENSE_EXCEPTION_INVALIDFILE      = 0x0101,
    LICENSE_EXCEPTION_MISMATCH         = 0x0102,
    LICENSE_EXCEPTION_EXPIRED          = 0x0103,
    /*
     * Add new exceptions above here!!!!!
     * (Instructions at top of this enum.)
     */
    LICENSE_EXCEPTION_MAX
} LICENSE_EXCEPTION_CODE;

/* Exception structure used to convey return information. */

typedef struct sLicenseException {
    LICENSE_EXCEPTION_CODE code;
    char                text[LICENSE_EXCEPTION_MESSAGELENGTH];
} LICENSE_EXCEPTION;

/*
  Explanation of license key file format


  The license key file consists of a single line containing the
  following fields:
  PRODUCT NOEXPIRE TIME1 TIME2 CHKSUM

  PRODUCT     This field contains a product name, which must not have
              any whitespace. E.g., PTA, PGA. This field is not
              encrypted (so it is easy to tell which product the
              license is for). The product name *is* case-sensitive.

  NOEXPIRE    This field is set to the value of the NOEXPIRE
              pre-processor macro if the license is non-expiring.
              Otherwise, a random value is chosen for the field.

  TIME1       These fields contain encrypted time stamps of when the
  TIME2       license expires. The (non-encrypted) time stamps are in
              the format returned by the time() function. The
              encryption has two phases:
              (1) Flipping of bits in the time stamp
              (2) XOR'ing of the time stamp with pre-generated random
              numbers.
              Both TIME1 and TIME2 are encrypted this way, but using
              different bit flips and different random numbers. When
              TIME1 and TIME2 are decrypted, they must be the exact
              same time stamp.

  CHKSUM      This is a checksum of the preceeding fields. The checksum
              is calculated using a series of bit shifts and
              multiplications on the preceeding fields. This makes it
              so changing one of the previous fields will invalidate
              the checksum and the file will not be able to be read.
              The checksum calculation is not bullet-proof, but should
              be appropriate for the level of security needed for
              the license file.
*/



#ifndef __WIN32		// we don't need this functionality on Windows
			// and the drand48() in the implementation is not
			// available in MinGW at this point.
/*
  Writes out license file for specified product with specified time stamp
  as the expiration time. Returns 1 on success, 0 otherwise.

  Arguments:
  product   Product name to write to the license file (single word).
  noexpire  Set to 1 if this should be a non-expiring license, 0
            otherwise.
  timestamp Time at which license should expire. It doesn't matter
            what you set this to if noexpire==1.
  licF      FILE to which the license is written.
*/

#ifdef NEW_LICENSE_FORMAT  // for the new license format

int WriteLicenseFile(const char product[PRODNAMELEN],
                     const char version[VERSIONLEN],
                     int noexpire,
                     const char *expiretime,
                     time_t timestamp,
                     const char mac[MACLEN],
                     int port,
                     int numCPUs,
                     char OS[OSLEN],
                     int some_integer,
                     char some_string[PRODNAMELEN],
                     FILE *licF);

#else // for old license format

int WriteLicenseFile(const char product[PRODNAMELEN],
		     int noexpire,
		     time_t timestamp,
		     FILE *licF);
#endif
#endif

/*
  Reads in license file and sets arguments for read-in values.
  Returns LICENSE_SUCCESS on success (valid file and read OK),
  LICENSE_EXCEPTION_THROWN otherwise.

  Arguments:
  product   Product name as read from the license file. The caller
            should check this value to make sure it has an expected
	    value.
  noexpire  Set to 1 if the read file is a non-expiring license, 0
            otherwise. If set to 1, then the caller should disregard
	    the value of timestamp.
  timestamp Time at which read license expires (not applicable if
            noexpire was set to 1).
  licF      FILE from which the license is read.
  exception Return information, error code and message
*/


#ifdef NEW_LICENSE_FORMAT  // for the new license format

LICENSE_RESULT
ReadLicenseFile(const char *product_in,
                int *noexpire,
                char *expiretime,
                time_t *timestamp,
                const char *version_in,
                int *port,
                int *numCPUs,
                int *some_integer,
                char *some_string,
                FILE *licF,
                LICENSE_EXCEPTION *lexception);

#else

LICENSE_RESULT
ReadLicenseFile(const char *product_in,
		int *noexpire,
		time_t *timestamp,
		FILE *licF,
		LICENSE_EXCEPTION *lexception);
#endif

/*
  Prints a license message for specified program name, which expires
  at specified expire_time. If expire_time is before the current time,
  then this function will print out an "already expired" message and
  will return 0. If expire_time is after the current time, then this
  function will print out a "will expire at..." message and will return 1.
  Note that this function does not take into account the noexpire flag,
  which must be checked before this function is called.
*/
int PrintLicenseMessage(const char *product, int noexp, time_t expire_time);


/*
  Checks the license file for specified product. What this does is take
  prodname and looks for the corresponding license file at
  $prodnameHOME/license.prodname. If the file doesn't exist, this function
  also checks ~/.license.prodname. If neither file exists, the license
  is corrupted, or if the license is expired, then the function prints
  a message and calls exit(2). If the license hasn't exipired yet, a
  message is printed that indicates the expiration date. If it is a
  non-expiring license, then nothing is printed.
*/
void CheckLicense(const char *prodname, const char *version);


#endif
