/**************************************************************************
 * This file is part of TraceTuner, a DNA sequencing quality value,
 * base calling and trace processing program.
 * Copyright (C) 2007-2008, J. Craig Venter Institute
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

#define TT_SFF_MAGIC 0x2E736666
#define MAX_BASES_LEN 1000
#define SFF_VERSION "\0\0\0\1"
#define MANIFEST_INDEX_MAGIC_NUMBER ( (((unsigned int) '.') << 24) + (((unsigned int) 'm') << 16) + (((unsigned int) 'f') << 8) + ((unsigned int) 't') )
#define SORT_INDEX_MAGIC_NUMBER ( (((unsigned int) '.') << 24) + (((unsigned int) 's') << 16) + (((unsigned int) 'r') << 8) + ((unsigned int) 't') )
#define INDEX_VERSION "1.00"

extern void getBytes(FILE *fp, unsigned int size, char *buf, char *filename);
extern int32_t getuint32(char *b);
extern unsigned int lookupAccno(char *accno, char *indexBuf, int indexSize);
extern int32_t getuint8( char *b);
extern int32_t getuint16(char *b);
extern int32_t getuint32(char *b);
extern int64_t getuint64(char *b);
extern void parseManifest(char *buf);
extern int getUAccnoInfo(char *accno, char *timeOut, int *regionOut, char *xyOut);
extern int getManifestInfo(char *accno, char *runNameOut,
                    char *analysisNameOut, char *runPathOut);
typedef struct {
  unsigned int  score  : 30;
  unsigned int  action : 2;
} dpCell;

typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32_t   magic_number;
  char       version[4];
  uint64_t   index_offset;
  uint32_t   index_length;
  uint32_t   number_of_reads;
  uint16_t   header_length;
  uint16_t   key_length;
  uint16_t   number_of_flows_per_read;
  uint8_t    flowgram_format_code;

  char    *flow_chars;        //  h->number_of_flows_per_read
  char    *key_sequence;      //  h->key_length

  void    *data_block;
  uint32_t   data_block_len;

  uint32_t   swap_endianess;

  //  Mate-finding DP related storage.
  //
  //  These really don't belong in here, but we know
  //  this is allocated.
  //
  char     alignA[MAX_BASES_LEN + MAX_BASES_LEN + 2];
  char     alignB[MAX_BASES_LEN + MAX_BASES_LEN + 2];
  dpCell   matrix[MAX_BASES_LEN][MAX_BASES_LEN];
} sffHeader;


typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32_t    magic_number;
  char      version[4];
  uint32_t    manifest_length;
  uint32_t    nothing;

  char     *manifest;
} sffManifest;

typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint16_t   read_header_length;
  uint16_t   name_length;
  uint32_t   number_of_bases;
  uint16_t   clip_quality_left;
  uint16_t   clip_quality_right;
  uint16_t   clip_adapter_left;
  uint16_t   clip_adapter_right;

  char      *name;                 //  r->name_length
  uint16_t  *flowgram_values;      //  h->number_of_flows_per_read
  uint8_t   *flow_index_per_base;  //  r->number_of_bases
  char      *bases;                //  r->number_of_bases
  uint8_t   *quality_values;       //  r->number_of_bases
  char      *quality;              //  quality_scores converted to CA-format qv

  int        final_length;         //  trimmed, processed read, ready for
  char      *final_bases;          //  loading.  NOT zero terminated.
  char      *final_quality;        //  DO NOT ZERO TERMINATE.

  void      *data_block;
  uint32_t   data_block_len;
} sffRead;

extern void readsff_manifest(FILE *sff, sffHeader *h, sffManifest *m);
extern int  readsff_header(FILE *sff, sffHeader *h, sffManifest *m);
extern void readsff_read(FILE *sff, sffHeader *h, sffRead *r);
