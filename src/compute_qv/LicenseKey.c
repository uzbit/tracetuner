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
 * 2.6 2003/11/14 23:31:02
 *
 * Source code needed to create and check license
 */

#include "LicenseKey.h"


/* magic number for non-expiration field. */
#define NOEXPIRE (228769551)

const int NUM_RAND_INT = 20;
const unsigned int RANDOMINT[20] =
{
  91066118u,
  2805734754u,
  1882820942u,
  4014223571u,
  1144696427u,
  2526554599u,
  1212808621u,
  2508456494u,
  288355551u,
  765440394u,
  3118137097u,
  353309801u,
  62061601u,
  802286842u,
  1220643822u,
  438044638u,
  1173472177u,
  180855140u,
  1542851914u,
  544274012u
};

const int NUM_RAND_INT2 = 20;
const unsigned int RANDOMINT2[20] =
{
  3940420493u,
  627266379u,
  1127871766u,
  2356436219u,
  3729353586u,
  2685654607u,
  2394412533u,
  351914199u,
  4189092454u,
  2914040945u,
  3469284305u,
  2105122193u,
  1361662504u,
  3152002657u,
  2426873141u,
  3170171781u,
  3592057401u,
  3460833611u,
  645329683u,
  1814330499u,
};

const int NUM_RAND_BIT = 40;
const unsigned int RANDOMBIT[40] =
{
  0, 27,
  1, 25,
  2, 31,
  3, 23,
  4, 28,
  5, 22,
  6, 30,
  7, 26,
  8, 24,
  9, 20,
  10, 13,
  11, 29,
  12, 24,
  13, 17,
  14, 15,
  15, 0,
  16, 5,
  17, 3,
  18, 19,
  20, 24
};

const int NUM_RAND_BIT2 = 360;
const unsigned int RANDOMBIT2[360] =
{
   8, 16, 31, 20,  5, 27, 21,  6, 20, 10,  9,  1, 31, 10,  3, 25,  7,  9,
   7,  1, 24, 12,  3,  2, 27, 19, 18,  6, 23,  3, 13,  0, 12,  1, 27, 18,
  15, 16, 23,  2, 14, 16, 17, 12,  3, 28,  9, 22, 12,  7, 10, 18, 22, 23,
   9, 30,  5, 12, 18, 20, 29, 16, 12, 13, 30,  8,  7,  4, 22, 20, 18, 18,
   0,  3, 16,  2, 21,  4, 19, 22, 12,  7,  0,  5,  4, 18, 30, 22, 24, 18,
  17,  4,  3,  1, 16,  7, 30, 11,  3, 10,  0,  7, 29, 21,  9, 15,  0, 11,
  10, 23, 30, 20, 14, 15, 22, 21, 14,  4,  6, 25, 22, 29,  3, 10, 24, 24,
   3, 21,  6,  2, 17,  0,  5,  4, 29, 30,  1,  8, 10, 18, 15, 14,  2,  6,
  18, 13, 20, 17, 12, 24, 13, 25,  8, 30,  4, 10, 30, 29, 28, 21,  6, 21,
   9, 28, 24, 10, 27, 19, 13, 16, 21,  9,  0,  0,  2, 11,  1, 10, 24, 21,
  12,  5, 28, 10, 23, 29, 14, 29,  1,  1, 24, 25, 24,  1,  6, 24,  7,  6,
   1, 27, 21, 25, 25, 20, 26,  2,  2, 20, 29,  0,  7,  7,  6, 28,  4, 11,
   0, 11,  0, 16, 20, 26,  0,  3,  4,  4, 23, 23, 26,  7, 10, 29,  4,  0,
  23, 28, 15, 25, 26,  6,  4, 30, 26,  3, 14, 11, 20,  6, 11,  4, 22, 31,
  19, 17, 18,  8, 19,  2, 18, 14,  9, 22, 24, 31, 12,  3, 22, 26, 17,  9,
   2, 18, 17, 29,  6, 14, 14,  7, 18,  5, 17, 30, 11, 16, 14, 26, 22,  8,
  17,  2, 13, 11, 12, 24, 26,  9, 10,  5, 12,  0, 10, 30,  9,  1, 14, 26,
  26, 14, 27, 14, 24, 10, 16, 11,  7,  6, 25, 12, 27, 12,  9, 25, 17,  2,
  14, 13, 26,  1, 18, 24,  1,  7,  0, 16,  4, 26,  8, 31, 24, 18, 21, 10,
  14, 31, 14,  7,  3, 25,  3, 25, 22,  4,  2, 10, 19, 21, 12,  0, 23, 31
};


void PrintBits(unsigned int i)
{
  int j;
  for (j = 31; j >= 0; j--)
  {
    printf("%d", (i >> j) & 0x01);
  }
  printf("\n");
}

/*
  Switchs the two specified bits in the input, returning the modified
  number. Assumes bit1 and bit2 are in the range 0 <= x <= 31.
*/
unsigned int SwitchBits(unsigned int input,
			unsigned int bit1,
			unsigned int bit2)
{
  unsigned int mask1, mask2;
  unsigned int setmask1, setmask2;

  /* masks to clear the bits from input */
  mask1 = ~(1 << bit1);
  mask2 = ~(1 << bit2);

  /*
     Masks to set the switched bits. We get the bit from input and
     shift it to its new location.
  */
  setmask1 = (((input >> bit1) & 0x01) << bit2);
  setmask2 = (((input >> bit2) & 0x01) << bit1);

  return ((input & mask1 & mask2) | setmask1 | setmask2);
}


/*
   Assume timestamp is at least 4 bytes. We apply some random bit
   shuffling and xor-ing to make it look more random.
*/
unsigned int EncryptTimeStamp1(unsigned int timestamp)
{
  unsigned int ret;
  int i;
  ret = (timestamp & 0xFFFFFFFF);

  for (i = 0; i < NUM_RAND_BIT; i += 2)
  {
    ret = SwitchBits(ret, RANDOMBIT[i], RANDOMBIT[i+1]);
  }
  for (i = 0; i < NUM_RAND_INT; i++)
  {
    ret ^= RANDOMINT[i];
  }
  return ret;
}


unsigned int EncryptTimeStamp2(unsigned int timestamp)
{
  unsigned int ret;
  int i;
  ret = (timestamp & 0xFFFFFFFF);

  for (i = 0; i < NUM_RAND_BIT2; i += 2)
  {
    ret = SwitchBits(ret, RANDOMBIT2[i], RANDOMBIT2[i+1]);
  }
  for (i = 0; i < NUM_RAND_INT2; i++)
  {
    ret ^= RANDOMINT2[i];
  }
  return ret;
}


unsigned int DecryptTimeStamp1(unsigned int timestamp)
{
  unsigned int ret;
  int i;
  ret = (timestamp & 0xFFFFFFFF);

  for (i = NUM_RAND_INT - 1; i >= 0; i--)
  {
    ret ^= RANDOMINT[i];
  }
  for (i = NUM_RAND_BIT - 1; i >= 0; i -= 2)
  {
    ret = SwitchBits(ret, RANDOMBIT[i-1], RANDOMBIT[i]);
  }
  return ret;
}


unsigned int DecryptTimeStamp2(unsigned int timestamp)
{
  unsigned int ret;
  int i;
  ret = (timestamp & 0xFFFFFFFF);

  for (i = NUM_RAND_INT2 - 1; i >= 0; i--)
  {
    ret ^= RANDOMINT2[i];
  }
  for (i = NUM_RAND_BIT2 - 1; i >= 0; i -= 2)
  {
    ret = SwitchBits(ret, RANDOMBIT2[i-1], RANDOMBIT2[i]);
  }
  return ret;
}


/* calculates checksum from specified arguments. this is nothing fancy,
   just some randomly chosen multiplications and hashing to get a
   somewhat unique checksum.
*/

#ifdef NEW_LICENSE_FORMAT  // for the new license format

unsigned int CheckSum(const char *s,
                      const char *s2,
		      unsigned int a,
		      unsigned int b,
		      unsigned int c,
                      const char *s3,
                      unsigned int d,
                      unsigned int e,
                      const char *s4,
                      unsigned int f,
                      const char *s5)
{
  unsigned int cksum = 0;
  const char *p = s;
  while (*p)
  {
    cksum <<= 1;
    cksum ^= *p++;
  }

  p = s2;
  while (*p)
  {
    cksum <<= 1;
    cksum ^= *p++;
  }

  cksum = (cksum << 1) * (a >> 5);
  cksum = (cksum >> 4) * (b >> 13);
  cksum = (cksum << 3) * (c >> 11);

  p = s3;
  while (*p)
  {
    cksum <<= 1;
    cksum ^= *p++;
  }

  cksum = cksum + 2*d + 13;
  cksum = cksum + 4*e + 117;

  p = s4;
  while (*p)
  {
    cksum ^= *p++;
  }

  cksum = cksum + 3*f + 24;

  p = s5;
  while (*p)
  {
    cksum ^= *p++;
  }

  return cksum;
}

/* returns 1 on success, 0 on failure */

#ifndef __WIN32
int
WriteLicenseFile(
  const char product[PRODNAMELEN],
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
  FILE *licF)
{
  int expire_flag;
  unsigned int enc_timestamp1, enc_timestamp2, checksum;

  if (noexpire)
  {
    expire_flag = NOEXPIRE;
  }
  else
  {
    /* choose random value for noexpire flag but make sure it isn't
       randomly set to NOEXPIRE. */
    do
    {
      expire_flag = ((int) (drand48() * (double) INT_MAX));
    } while (expire_flag == NOEXPIRE);
  }
  enc_timestamp1 = EncryptTimeStamp1((unsigned int) timestamp);
  enc_timestamp2 = EncryptTimeStamp2((unsigned int) timestamp);

  /* calculate checksum for these encrypted values */
  checksum = CheckSum(product, version, expire_flag, enc_timestamp1, 
                      enc_timestamp2, mac, port, numCPUs, OS, some_integer,
                      some_string);

  if (!fprintf(licF, "%s %s \"%s\" %u %u %u %s %u %u %s  %u %s %u\n",
	       product, version, expiretime, expire_flag, 
               enc_timestamp1, enc_timestamp2, mac, port, numCPUs, OS, 
               some_integer, some_string, checksum))
  {
    return 0;
  }
  else
  {
    return 1;
  }
}
#endif

/* returns 1 on success, 0 on failure */

int 
isHexDigit(int c) {

    if ((c >= '0' && c <= '9') 
    ||  (c >= 'a' && c <= 'f')
    ||  (c >= 'A' && c <= 'F'))
    {
        return 1;
    }
    else {
        return 0;
    }

}

/* returns 1 on success, 0 on failure */

int
parseForMAC(FILE *fp, char ***mac, int *numMAC) {

    int c;
    int numGroups;
    int macNumber;
    char possibleMAC[MACLEN];
    int validMAC;
    int i; /* counter for MAC characters from start of string; 
              between 11 and 17 */
    int j; /* counter for characters per group; either one or two */
    int k;

    i = 0;
    j = 0;
    numGroups = 0;
    macNumber = 0;
    validMAC = 0;
    memset(possibleMAC, 0, MACLEN);
    
    while ((c = fgetc(fp)) != EOF) {
	if (numGroups == 6 && validMAC == 1) {
	    macNumber++;
	    if (macNumber > 1) {
	        *mac = (char **) realloc(*mac, macNumber * sizeof (char *));
	    } else {
	        *mac = (char **) malloc(sizeof (char *));
	    }
	    if (*mac == NULL) {
	      fprintf(stderr, "out of memory\n");
	      exit(2);
	    }
	    (*mac)[macNumber-1] = (char *) malloc(i + 1);
	    if ((*mac)[macNumber - 1] == NULL) {
	      fprintf(stderr, "out of memory\n");
	      exit(2);
	    }
	    for (k = 0; k < i; k++) {
		(*mac)[macNumber-1][k] = possibleMAC[k];
	    }
	    (*mac)[macNumber-1][i] = '\0';
	    *numMAC = macNumber;
	    memset(possibleMAC, 0, MACLEN);
	    numGroups = 0;
	    i = 0;
	    j = 0;
	    validMAC = 0;
	}
        if (!isspace(c)) {
            if (isHexDigit(c) && j < 2 && numGroups < 6) {
	        possibleMAC[i] = c;
	        i++;
                j++;
		if (numGroups == 5 && j == 2) {
		    numGroups++;
		    validMAC = 1;
		}
	    }
	    else if (c == ':' && (j == 1 || j == 2)) {
	        possibleMAC[i] = c;
	        i++;
	        j = 0;
		numGroups++;
	    }
#ifdef __WIN32
	    else if (c == '-' && (j == 1 || j == 2)) {
	        possibleMAC[i] = c;
	        i++;
	        j = 0;
		numGroups++;
	    }
#endif
            else {
		memset(possibleMAC, 0, MACLEN);
	        i = 0;
	        j = 0;
		numGroups = 0;
		validMAC = 0;
	    }
        }
	else if (numGroups == 5 && j > 0) {
	    numGroups++;
            validMAC = 1;
	}
	else {
	    memset(possibleMAC, 0, MACLEN);
	    i = 0;
	    j = 0;
	    numGroups = 0;
	    validMAC = 0;
	}
    }

    return 1; 
}


int
getMACList(char ***mac, int *mac_count) {
    char *OSNAME = OS_NAME;  /* defined in Makefile */
    FILE *fp;

    /* Cache the results so we don't have to call the external program
       more than once. */
    static char **mac_address = NULL;
    static int mac_address_count = -1;

    if (mac_address_count < 0) {
	mac_address = NULL;
	mac_address_count = 0;

	if (strcmp(OSNAME, "SunOS") == 0) {
	    if ((fp = popen("/usr/sbin/arp `hostname`", "r")) == NULL) {
		fprintf(stderr, "RetrieveMAC: unable to get MAC address.\n");
		return 0;
	    }
	}
	else if (strcmp(OSNAME, "Linux") == 0) {
	    if ((fp = popen("/sbin/ifconfig", "r")) == NULL) {
		fprintf(stderr, "RetrieveMAC: unable to get MAC address.\n");
		return 0;
	    }
	}
	else if (strcmp(OSNAME, "OSF1") == 0) {
	    if ((fp = popen("/usr/sbin/netstat -i", "r")) == NULL) {
		fprintf(stderr, "RetrieveMAC: unable to get MAC address.\n");
		return 0;
	    }
	}
	else if (strcmp(OSNAME, "IRIX") == 0) {
	    if ((fp = popen("/usr/etc/netstat -ia", "r")) == NULL) {
		fprintf(stderr, "RetrieveMAC: unable to get MAC address.\n");
		return 0;
	    }
	}
	else if (strcmp(OSNAME, "CYGWIN_NT-4.0") == 0) {
	    if ((fp = popen("ipconfig /all", "r")) == NULL) {
		fprintf(stderr, "RetrieveMAC: unable to get MAC address.\n");
		return 0;
	    }
	}
	else if (strcmp(OSNAME, "CYGWIN_NT-5.0") == 0) {
	    if ((fp = popen("ipconfig /all", "r")) == NULL) {
		fprintf(stderr, "RetrieveMAC: unable to get MAC address.\n");
		return 0;
	    }
	}
	else {
	    fprintf(stderr, "RetrieveMAC: unsupported system.\n");
	    return 0;
	}
    
	if (parseForMAC(fp, &mac_address, &mac_address_count) == 0) {
	    pclose(fp);
	    return 0;
	}

	pclose(fp);
    }

    *mac = mac_address;
    *mac_count = mac_address_count;
    return 1;
}

/* returns 1 on success, 0 on failure */

int
isMACValid(const char *licensedMAC) {
    char **mac = NULL;
    int mac_count = 0;
    int k;

    getMACList(&mac, &mac_count);

    for (k = 0; k < mac_count; k++) {
	if (strcmp(licensedMAC, mac[k]) == 0) {
	    // MAC in license file matches one of the machine's MAC addresses
	    return 1;
	}
    }
    
    // MAC in license file didn't match any of the machine's MAC addresses    
    return 0; 
}



/* returns 1 on success, 0 on failure */

int
isOSValid(const char *licensedOS) {

    char *OSNAME = OS_NAME;  /* defined in Makefile */

    if (strcmp(OSNAME, licensedOS)) {
	// OS in license file does not match that of machine
	return 0;
    }
    
    // OS in license file matches that of machine   
    return 1; 
}

LICENSE_RESULT
ReadLicenseFile(
  const char *product_in,
  int *noexpire,
  char *expiretime,
  time_t *timestamp,
  const char *version_in,
  int *port,
  int *numCPUs,
  int *some_integer,
  char *some_string,
  FILE *licF,
  LICENSE_EXCEPTION *lexception)
{
  char product[PRODNAMELEN];
  char licensedVersion[VERSIONLEN];
  char mac[MACLEN];
  char OS[OSLEN];
  time_t decrypt1, decrypt2;
  int expire_flag;
  unsigned int enc_timestamp1, enc_timestamp2, checksum_read, checksum_calc;
  int k;

  /* Read the contents of the license file. */
  if (fscanf(licF, "%s %s %s %u %u %u %s %u %u %s %u %s %u\n",
	     product, licensedVersion, expiretime, &expire_flag, 
             &enc_timestamp1, &enc_timestamp2, mac, port, numCPUs, OS,
             some_integer, some_string, &checksum_read)
      != 13)
  {
    sprintf(lexception->text, "License key corrupt");
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  /* calculate the checksum for the values we read in */
  checksum_calc = CheckSum(product, licensedVersion, expire_flag, 
                           enc_timestamp1, enc_timestamp2, mac, *port, 
                           *numCPUs, OS, *some_integer, some_string);

  /*
    if the read checksum does not agree with the calculated one, then
    the file is corrupt.
  */
  if (checksum_read != checksum_calc)
  {
    sprintf(lexception->text, "License key checksum mismatch");
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  decrypt1 = (time_t) DecryptTimeStamp1(enc_timestamp1);
  decrypt2 = (time_t) DecryptTimeStamp2(enc_timestamp2);

  /* make sure the two decrypted timestamps agree */
  if (decrypt1 == decrypt2)
  {
    *timestamp = decrypt1;
    *noexpire = (expire_flag == NOEXPIRE) ? 1 : 0;
  }
  else
  {
    /* If not, return an error. */
    sprintf(lexception->text, "Decrypted timestamp mismatch");
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }
  
  /* Make sure the product names match. */
  if ( strcmp(product, product_in))
  {
    sprintf(lexception->text, "Product name mismatch");
    lexception->code = LICENSE_EXCEPTION_MISMATCH;
    return LICENSE_EXCEPTION_THROWN;
  }

  /* Make sure there is a product version in the license */
  if (licensedVersion == NULL) {
    sprintf(lexception->text, "Invalid product version");
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  /* Make sure the versions match */
  k = 0;
  while (licensedVersion[k] != '\0') {
    if (licensedVersion[k] != version_in[k]) {
      sprintf(lexception->text, "Product version mismatch");
      lexception->code = LICENSE_EXCEPTION_MISMATCH;
      return LICENSE_EXCEPTION_THROWN;
    }
    k++;
  }

  /* Make sure that the OS matches that in the license. */
  if (!isOSValid(OS)) {
    sprintf(lexception->text, "License was issued for a different OS");
    lexception->code = LICENSE_EXCEPTION_MISMATCH;
    return LICENSE_EXCEPTION_THROWN;
  }
  
  /* 
   * Make sure that the MAC address in the license file matches one of the 
   * hardware MAC addresses.
   */
  if (!isMACValid(mac)) {
    sprintf(lexception->text, "License was issued for a different machine");
    lexception->code = LICENSE_EXCEPTION_MISMATCH;
    return LICENSE_EXCEPTION_THROWN;
  }

  return LICENSE_SUCCESS;
}

#else  // for the old license format

unsigned int CheckSum(const char *s,
		      unsigned int a,
		      unsigned int b,
		      unsigned int c)
{
  unsigned int cksum = 0;
  const char *p = s;
  while (*p)
  {
    cksum <<= 1;
    cksum ^= *p++;
  }
  cksum = (cksum << 1) * (a >> 5);
  cksum = (cksum >> 4) * (b >> 13);
  cksum = (cksum << 3) * (c >> 11);

  return cksum;
}

#ifndef __WIN32
int
WriteLicenseFile(
  const char product[PRODNAMELEN],
  int noexpire,
  time_t timestamp,
  FILE *licF)
{
  int expire_flag;
  unsigned int enc_timestamp1, enc_timestamp2, checksum;

  if (noexpire)
  {
    expire_flag = NOEXPIRE;
  }
  else
  {
    /* choose random value for noexpire flag but make sure it isn't
       randomly set to NOEXPIRE. */
    do
    {
      expire_flag = ((int) (drand48() * (double) INT_MAX));
    } while (expire_flag == NOEXPIRE);
  }
  enc_timestamp1 = EncryptTimeStamp1((unsigned int) timestamp);
  enc_timestamp2 = EncryptTimeStamp2((unsigned int) timestamp);

  /* calculate checksum for these encrypted values */
  checksum = CheckSum(product, expire_flag, enc_timestamp1, enc_timestamp2);

  if (!fprintf(licF, "%s %u %u %u %u\n",
	       product, expire_flag, enc_timestamp1, enc_timestamp2,
	       checksum))
  {
    return 0;
  }
  else
  {
    return 1;
  }
}
#endif


LICENSE_RESULT
ReadLicenseFile(
  const char *product_in,
  int *noexpire,
  time_t *timestamp,
  FILE *licF,
  LICENSE_EXCEPTION *lexception)
{
  char product[PRODNAMELEN];
  time_t decrypt1, decrypt2;
  int expire_flag;
  unsigned int enc_timestamp1, enc_timestamp2, checksum_read, checksum_calc;

  /* Read the contents of the license file. */
  if (fscanf(licF, "%s %u %u %u %u\n",
	     product, &expire_flag, &enc_timestamp1, &enc_timestamp2,
	     &checksum_read)
      != 5)
  {
    sprintf(lexception->text, "License key corrupt");
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  /* calculate the checksum for the values we read in */
  checksum_calc = CheckSum(product, expire_flag, enc_timestamp1,
			   enc_timestamp2);

  /*
    if the read checksum does not agree with the calculated one, then
    the file is corrupt.
  */
  if (checksum_read != checksum_calc)
  {
    sprintf(lexception->text, "License key checksum mismatch");
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  decrypt1 = (time_t) DecryptTimeStamp1(enc_timestamp1);
  decrypt2 = (time_t) DecryptTimeStamp2(enc_timestamp2);

  /* make sure the two decrypted timestamps agree */
  if (decrypt1 == decrypt2)
  {
    *timestamp = decrypt1;
    *noexpire = (expire_flag == NOEXPIRE) ? 1 : 0;
  }
  else
  {
    /* If not, return an error. */
    sprintf(lexception->text, "Decrypted timestamp mismatch");
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  /* Make sure the product names match. */
  if ( strcmp(product, product_in))
  {
    sprintf(lexception->text, "Product name mismatch");
    lexception->code = LICENSE_EXCEPTION_MISMATCH;
    return LICENSE_EXCEPTION_THROWN;
  }


  return LICENSE_SUCCESS;
}

#endif

/*
  This function checks the expiration of a license at the current
  time.

  product is the name of the product
  noexp is the license non-expiring?
  expire_time is the expiration time of the license
  exception contains return information.

  product is the name of the product, passed in.
  noexp is true if the license has no expiration.
  expire_time is the time when the license expires.

  Returns LICENSE_SUCCESS on success,
  LICENSE_EXCEPTION_THROWN otherwise.

 */
LICENSE_RESULT
CheckLicenseExpiration(
  const char *product,
  int noexp,
  time_t expire_time,
  LICENSE_EXCEPTION *lexception)
{
  time_t currtime;
  const int timebufsize = 200;
  char timebuf[200];

    if(noexp) {
	/* If the license didn't expire, it's fine. */
	return LICENSE_SUCCESS;
    } else {
	/* Otherwise check the expiration. */
	currtime = time(NULL);

	if (currtime >= expire_time) {
	    /* If it expired, tell them so. */
	    strftime(timebuf, timebufsize, "%b %d, %Y %H:%M:%S",
		     localtime(&expire_time));

	    sprintf(lexception->text, "License for %s expired on %s",
		    product, timebuf);
	    lexception->code = LICENSE_EXCEPTION_EXPIRED;
	    return  LICENSE_EXCEPTION_THROWN;
	}


    }
    return LICENSE_SUCCESS;

}

int PrintLicenseMessage(const char *product, int noexp, time_t expire_time)
{
  time_t currtime;
  const int timebufsize = 200;
  char timebuf[200];
  int ret;

  if (noexp) {
    printf("The license for %s is unlimited.\n", product);
    ret = 0;
    return ret;
  }


  currtime = time(NULL);

  strftime(timebuf, timebufsize, "%b %d, %Y %H:%M:%S",
	   localtime(&expire_time));

  if (currtime >= expire_time)
  {
    printf("The license for %s expired on %s.\n", product, timebuf);
    ret = 0;
    return ret;
  }
  else
  {
    int numdays;

    ret = 1;
    /* number of integral days left on license */
    numdays = (expire_time - currtime) / 86400;

    printf("The license for %s will expire on %s.\n", product, timebuf);

    /*
      If expiration time is less than a week away, then we print extra
      messages telling the user how to renew.
    */
    if (numdays < 7)
    {
      printf("You have less than %d day%s remaining to renew the license!\n",
	     numdays+1, numdays ? "s" : "");
    }
  }

  printf("Please contact Paracel customer support (support@paracel.com)\n");
  printf("to purchase an extended license.\n");
  return ret;
}

#ifdef NEW_LICENSE_FORMAT  // for the new license format

/*
  This function checks a license given a product name and license file that
  was written in the new license format (productName productVersion no_expire 
  timeStamp1 timeStamp2 MACAdress portNumber numLicensedCPUs checksum).

  product is the name of the product
  filename is the file with the license
  exception contains return information

  Returns LICENSE_SUCCESS on success,
  LICENSE_EXCEPTION_THROWN otherwise.

 */
LICENSE_RESULT
LicenseCheck(
  const char *prodname,
  const char *filename,
  const char *version,
  int *port,
  int *numCPUs,
  LICENSE_EXCEPTION *lexception)
{
  FILE *F;
  int noexp;
  time_t exptime;
  LICENSE_RESULT r;
  int some_integer;
  char some_string[PRODNAMELEN];
  char expiretime[PRODNAMELEN];


  /* Open the file. */
  F = fopen(filename, "r");

  if (!F) {
    sprintf(lexception->text, "Unable to open file %s", filename);
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  /* Read the license file, aborting on corrupted file. */
  r = ReadLicenseFile(prodname, &noexp, expiretime, &exptime, version, 
                      port, numCPUs, &some_integer, some_string, F, 
                      lexception);
  fclose(F);
  
  if (r != LICENSE_SUCCESS)
  {
    return r;
  }

  /* Make sure the license hasn't expired. */

  r = CheckLicenseExpiration(prodname, noexp, exptime, lexception);
  if (r != LICENSE_SUCCESS) {
    return r;
  }

  /* Everything checks out OK. */
  return LICENSE_SUCCESS;
}

#else  // for the old license format

/*
Checks a license file, with some environment variable stuff in there.
Also does some stuff with environment variables.
Prints messages to stderr.
Exits on error.

prodname is the product name passed in.
filename_in is the filename to check.
noexp is set to true if it's a non-expiring license.
exptime is set to the expiration time of the license.
 */

void
CheckLicenseFile(
  const char *prodname,
  const char *filename_in,
  int *noexp,
  time_t *exptime)
{
  FILE *F;
  char filename[1024];
  LICENSE_RESULT r;
  LICENSE_EXCEPTION lexception;

  memset(&lexception, 0, sizeof(LICENSE_EXCEPTION));

  strcpy(filename, filename_in);

  F = fopen(filename, "r");
  if (!F)
  {
    /* If license file is not under envhome, then we look for .license.*
       under the user's home directory. */
    char filename2[1024];
    const char *homedir = getenv("HOME");

    /* If HOME is set... */
    if (homedir)
    {
#ifdef __WIN32
      sprintf(filename2, "%s\\.license.%s", homedir, prodname);
#else
      sprintf(filename2, "%s/.license.%s", homedir, prodname);
#endif
      F = fopen(filename2, "r");
    }

    /* If we opened up the license file successfully, then we remember
       the file name so we can reference it later. */
    if (F)
    {
      strcpy(filename, filename2);
    }
    /* Else we print a message about not being able to find the license
       file under prodhome (we don't tell the user that we also looked
       under their home directory). */
    else
    {
      fprintf(stderr, "ERROR: Unable to open license file %s.  Abort.\n",
	      filename);
      exit(2);
    }
  }

  /* Read the license file, aborting on corrupted file. */
  r = ReadLicenseFile(prodname, noexp, exptime, F, &lexception);
  fclose(F);
  if (r != LICENSE_SUCCESS)
  {
    fprintf(stderr, "ERROR reading license file %s:%s.  Abort\n",
	    filename, lexception.text);
    exit(2);
  }

  (void) PrintLicenseMessage(prodname, *noexp, *exptime);
}


/*
  This function checks a license given a product name and license file that
  written in the old license format (productName no_expire timestamp1 
  timestamp2 checksum).

  product is the name of the product
  filename is the file with the license
  lexception contains return information

  Returns LICENSE_SUCCESS on success,
  LICENSE_EXCEPTION_THROWN otherwise.

 */
LICENSE_RESULT
LicenseCheck(
  const char *prodname,
  const char *filename,
  LICENSE_EXCEPTION *lexception)
{
  FILE *F;
  int noexp;
  time_t exptime;
  LICENSE_RESULT r;


  /* Open the file. */
  F = fopen(filename, "r");

  if (!F) {
    sprintf(lexception->text, "Unable to open file %s", filename);
    lexception->code = LICENSE_EXCEPTION_INVALIDFILE;
    return LICENSE_EXCEPTION_THROWN;
  }

  /* Read the license file, aborting on corrupted file. */
  r = ReadLicenseFile(prodname, &noexp, &exptime, F, lexception);
  fclose(F);

  if (r != LICENSE_SUCCESS)
  {
    return r;
  }

  /* Make sure the license hasn't expired. */

  r = CheckLicenseExpiration(prodname, noexp, exptime, lexception);
  if (r != LICENSE_SUCCESS) {
    return r;
  }

  /* Everything checks out OK. */
  return LICENSE_SUCCESS;
}

#endif

void CheckLicense(const char *prodname, const char *version)
{
    char filename[1024];
    char filename_public[1024];
    char envname[256];
    char *envhome;
    struct stat statbuf;
    int r;
    int i;
#ifdef NEW_LICENSE_FORMAT
    char **mac;
    int mac_count;
#endif
    const char *home = getenv("HOME");

#ifndef CHECK_DEVELOPMENTAL_LICENSE
    /* Check for "DEVELOPMENTAL" without being too obvious.  If this is
       a developmental version, return success.  */
    if (version[0] == 'D' && version[2] == 'V' && version[4] == 'L' &&
	version[1] == 'E' && version[3] == 'E') {
	PrintLicenseMessage(prodname, 1, 0);
	return;
    }
#endif

    /* Check for a license file under the following names:
       $PRDHOME/license.PRD-VER-MAC
       $HOME/.license.PRD-VER-MAC
       $PRDHOME/license.PRD
       $HOME/.license.PRD */

    /* Get the $PRDHOME environment variable. */
    sprintf(envname, "%sHOME", prodname);
    envhome = getenv(envname);
    if (envhome == NULL) {
	fprintf(stderr, "ERROR: Environment variable %s not set.  Abort.\n",
		envname);
	exit(2);
    }

#ifdef NEW_LICENSE_FORMAT
    getMACList(&mac, &mac_count);

    /* Try $PRDHOME/license.PRD.MAC. */
    for (i = 0; i < mac_count; i++) {
	sprintf(filename, "%s/license.%s-%s-%s",
		envhome, prodname, version, mac[i]);
	r = stat(filename, &statbuf);
	if (r == 0) goto license_found;
    }

    /* Try $HOME/.license.PRD.MAC. */
    if (home != NULL) {
	for (i = 0; i < mac_count; i++) {
	    sprintf(filename, "%s/.license.%s-%s-%s",
		    home, prodname, version, mac[i]);
	    r = stat(filename, &statbuf);
	    if (r == 0) goto license_found;
	}
    }
#endif

    /* Try $PRDHOME/license.PRD. */
#ifdef __WIN32
    sprintf(filename, "%s\\license.%s", envhome, prodname);
#else
    sprintf(filename, "%s/license.%s", envhome, prodname);
#endif
    r = stat(filename, &statbuf);
    if (r == 0) goto license_found;

    /* Save this last filename to report as the one we couldn't find. */
    strcpy(filename_public, filename);

    /* Try $HOME/.license.PRD. */
    if (home != NULL) {
	sprintf(filename, "%s/.license.%s", home, prodname);
	r = stat(filename, &statbuf);
	if (r == 0) goto license_found;
    }

    /* Didn't find a license file under any name. */
    fprintf(stderr, "ERROR: Unable to open license file %s.  Abort.\n",
	    filename_public);
    exit(2);

license_found:

#ifdef NEW_LICENSE_FORMAT  // for the new license format
    {
	int port;
	int numCPUs;
	LICENSE_EXCEPTION lexception;
	LICENSE_RESULT lr;

	lr = LicenseCheck(prodname, filename, version, &port, &numCPUs,
			  &lexception);
	if (lr != LICENSE_SUCCESS) {
	    fprintf(stderr, "ERROR: License error: %s.  Abort.\n",
		    lexception.text);
	    exit(2);
	}
    }
#else
    {
	int noexp;
	time_t timestamp;

	CheckLicenseFile(prodname, filename, &noexp, &timestamp);
    }
#endif
}


// Local Variables:
// c-basic-offset: 4
// End:

