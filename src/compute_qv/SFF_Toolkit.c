
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

//"$Id: SFF_Toolkit.c,v 1.4 2009/01/08 20:45:03 gdenisov Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include "Btk_qv.h"
#include "SFF_Toolkit.h"
#include "util.h"

inline
uint64_t
uint64Swap(uint64_t x) {
  x = ((x >>  8) & 0x00ff00ff00ff00ffLLU) | ((x <<  8) & 0xff00ff00ff00ff00LLU);
  x = ((x >> 16) & 0x0000ffff0000ffffLLU) | ((x << 16) & 0xffff0000ffff0000LLU);
  x = ((x >> 32) & 0x00000000ffffffffLLU) | ((x << 32) & 0xffffffff00000000LLU);
  return(x);
}

inline
uint32_t
uint32Swap(uint32_t x) {
  x = ((x >>  8) & 0x00ff00ff) | ((x <<  8) & 0xff00ff00);
  x = ((x >> 16) & 0x0000ffff) | ((x << 16) & 0xffff0000);
  return(x);
}

inline
uint16_t
uint16Swap(uint16_t x) {
  x = ((x >>  8) & 0x000000ff) | ((x <<  8) & 0x0000ff00);
  return(x);
}

off_t
TT_ftell(FILE *stream) {
  off_t  pos = 0;
  errno = 0;
  pos = ftello(stream);
  if (errno == ESPIPE)
    //  Not a seekable stream.  Return some goofy big number.
    return(((off_t)1) < 42);
  if (errno) {
    fprintf(stderr, "AS_UTL_ftell()--  Failed with %s.\n", strerror(errno));
    assert(errno == 0);
  }
  return(pos);
}

// Read nobj objects of size size from buffer
size_t
TT_safeRead(FILE *file, void *buffer, char *desc, size_t size, size_t nobj) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024 / size;
  size_t  toread   = 0;
  size_t  written  = 0;  //  readen?

  while (position < nobj) {
    toread = length;
    if (position + toread > nobj)
      toread = nobj - position;

    errno = 0;
    written = fread(((char *)buffer) + position * size, size, toread, file);
    position += written;

    if (feof(file) || (written == 0))
      goto finish;

    if ((errno) && (errno != EINTR)) {
      fprintf(stderr, "sffRead()-- Read failure on %s: %s.\n", desc, strerror(errno));
//    fprintf(stderr, "sffRead()-- Wanted to read "F_SIZE_T" objects (size="F_SIZE_T"), read "F_SIZE_T".\n",
//            toread, size, written);
      assert(errno == 0);
    }
  }

 finish:
  //  Just annoys developers.  Stop it.
  //if (position != nobj)
  //  fprintf(stderr, "AS_UTL_sffRead()--  Short read; wanted "F_SIZE_T" objects, read "F_SIZE_T" instead.\n",
  //          nobj, position);
  return(position);
}

void
readsff_manifest(FILE *sff, sffHeader *h, sffManifest *m) {

  if (h->index_length == 0)
    //  No manifest.
    return;

  if (TT_ftell(sff) != h->index_offset)
    //  Not at the manifest.
    return;

  if (m->manifest)
    //  Already got it?!
    return;

  TT_safeRead(sff, m, "readsff_manifest", sizeof(char), 16);

  if (h->swap_endianess) {
    m->magic_number    = uint32Swap(m->magic_number);
    m->manifest_length = uint32Swap(m->manifest_length);
  }

  m->manifest = CALLOC(char, m->manifest_length + 1);
  TT_safeRead(sff, m->manifest, "readsff_manifest_text", sizeof(char), m->manifest_length);

  m->manifest[m->manifest_length] = 0;

  //  We only read the manifest.  There is still an index in there.

  uint64_t  padding_length = h->index_length - 16 - m->manifest_length;
  if (padding_length > 0) {
    //fprintf(stderr, "manifest pad "F_U64"\n", padding_length);
    char *junk = CALLOC(char, padding_length);
    TT_safeRead(sff, junk, "readsff_manifest_pad", sizeof(char), padding_length);
    FREE(junk);
  }
}


int 
readsff_header(FILE *sff, sffHeader *h, sffManifest *m) {

  TT_safeRead(sff, h, "readsff_header_1", 31, 1);

  if (h->magic_number != 0x2e736666) {
    h->swap_endianess           = 1;
    h->magic_number             = uint32Swap(h->magic_number);
    h->index_offset             = uint64Swap(h->index_offset);
    h->index_length             = uint32Swap(h->index_length);
    h->number_of_reads          = uint32Swap(h->number_of_reads);
    h->header_length            = uint16Swap(h->header_length);
    h->key_length               = uint16Swap(h->key_length);
    h->number_of_flows_per_read = uint16Swap(h->number_of_flows_per_read);
  }

  if (h->magic_number != 0x2e736666) return -8;

  uint32_t newlen = h->number_of_flows_per_read + h->key_length + 2;
  if (h->data_block_len < newlen) {
    h->data_block_len = newlen;
    h->data_block     = REALLOC(h->data_block, void, h->data_block_len);
  }

  memset(h->data_block, 0, h->data_block_len);

  h->flow_chars   = h->data_block;
  h->key_sequence = h->data_block + (h->number_of_flows_per_read + 1) * sizeof(char);

  TT_safeRead(sff,  h->flow_chars,   "readsff_header_2", sizeof(char), h->number_of_flows_per_read);
  TT_safeRead(sff,  h->key_sequence, "readsff_header_3", sizeof(char), h->key_length);

  uint64_t  padding_length = h->header_length - 31 - h->number_of_flows_per_read - h->key_length;
  if (padding_length > 0) {
    //fprintf(stderr, "header pad "F_U64"\n", padding_length);
    char *junk = CALLOC(char, padding_length);
    TT_safeRead(sff, junk, "readsff_header_4", sizeof(char), padding_length);
    FREE(junk);
  }

  //  The spec says the index might be here, however, all files I've
  //  seen have the index at the end of the file.
  //
  readsff_manifest(sff, h, m);
  return SUCCESS;
}


void
readsff_read(FILE *sff, sffHeader *h, sffRead *r) 
{
  int i;
  unsigned char b[2];

  TT_safeRead(sff, r, "readsff_read_1", 16, 1);

  if (h->swap_endianess) {
    r->read_header_length = uint16Swap(r->read_header_length);
    r->name_length        = uint16Swap(r->name_length);
    r->number_of_bases    = uint32Swap(r->number_of_bases);
    r->clip_quality_left  = uint16Swap(r->clip_quality_left);
    r->clip_quality_right = uint16Swap(r->clip_quality_right);
    r->clip_adapter_left  = uint16Swap(r->clip_adapter_left);
    r->clip_adapter_right = uint16Swap(r->clip_adapter_right);
  }

  //  Can you say UGLY?  Hey, it's a lot better than what I originally came up with.

  uint32_t ss[6];
  ss[0] = (r->name_length + 1)          * sizeof(char);
  ss[1] = (h->number_of_flows_per_read) * sizeof(uint16_t) + ss[0];
  ss[2] = (r->number_of_bases)          * sizeof(uint8_t)  + ss[1];
  ss[3] = (r->number_of_bases + 1)      * sizeof(char)     + ss[2];
  ss[4] = (r->number_of_bases + 1)      * sizeof(uint8_t)  + ss[3];
  ss[5] = (r->number_of_bases + 1)      * sizeof(char)     + ss[4];

  if (r->data_block_len < ss[5]) {
    r->data_block_len = ss[5];
    r->data_block     = REALLOC(r->data_block, void, r->data_block_len);
  }

  memset(r->data_block, 0, r->data_block_len);

  r->name                 = r->data_block;
  r->flowgram_values      = r->data_block + ss[0];
  r->flow_index_per_base  = r->data_block + ss[1];
  r->bases                = r->data_block + ss[2];
  r->quality_values       = r->data_block + ss[3];
  r->quality              = r->data_block + ss[4];

  TT_safeRead(sff, r->name, "readsff_read_2", sizeof(char), r->name_length);
  r->name[r->name_length] = 0;

  uint64_t  padding_length = r->read_header_length - 16 - r->name_length;
  if (padding_length > 0) {
    //fprintf(stderr, "read pad 1 "F_U64"\n", padding_length);
    uint64_t  junk;
    TT_safeRead(sff, &junk, "readsff_read_3", sizeof(char), padding_length);
  }

  for (i=0; i < h->number_of_flows_per_read; i++) {
      fread(b, 1, 2, sff);
      r->flowgram_values[i] = (uint32_t) b[0] * 256 + (uint32_t) b[1];
  }

  TT_safeRead(sff, r->flow_index_per_base, "readsff_read_5", sizeof(uint8_t),  r->number_of_bases);

  TT_safeRead(sff, r->bases,               "readsff_read_6", sizeof(char),     r->number_of_bases);
  TT_safeRead(sff, r->quality_values,      "readsff_read_7", sizeof(uint8_t),  r->number_of_bases);

  for (i=0; i<r->number_of_bases; i++)
    r->quality[i] = r->quality_values[i] + '0';

  r->bases[r->number_of_bases] = 0;
  r->quality[r->number_of_bases] = 0;

  //  The padding_length is the number of bytes to make the above four
  //  chunks of data be of size that is divisible by 8.  The
  //  padding_length we compute directly below is the number of bytes
  //  we read past the last multiple of 8, and if that is non-zero, we
  //  need to read 8-padding_length bytes.
  //
  padding_length = (h->number_of_flows_per_read * sizeof(uint16_t) +
                    r->number_of_bases * sizeof(uint8_t) +
                    r->number_of_bases * sizeof(char) +
                    r->number_of_bases * sizeof(uint8_t)) % 8;
  if (padding_length > 0) {
    //fprintf(stderr, "read pad 2 "F_U64"\n", 8-padding_length);
    char *junk = CALLOC(char, padding_length > 8 ? padding_length : 8);
    TT_safeRead(sff, junk, "readsff_read_8", sizeof(char), 8 - padding_length);
    FREE(junk);
  }
}

void
TT_fseek(FILE *stream, off_t offset, int whence) 
{
  off_t   beginpos = TT_ftell(stream);

  //  If the stream is already at the correct position, just return.
  //
  //  Unless we're on FreeBSD.  For unknown reasons, FreeBSD fails
  //  updating the gkpStore with mate links.  It seems to misplace the
  //  file pointer, and ends up writing the record to the wrong
  //  location.  ftell() is returning the correct current location,
  //  and so AS_PER_genericStore doesn't seek() and just writes to the
  //  current position.  At the end of the write, we're off by 4096
  //  bytes.
  //
  //  LINK 498318175,1538 <-> 498318174,1537
  //  TT_fseek()--  seek to 159904 (whence=0); already there
  //  sffWrite()-- write nobj=1x104 = 104 bytes at position 159904
  //  sffWrite()-- wrote nobj=1x104 = 104 bytes position now 164000
  //  sffWrite()-- EXPECTED 160008, ended up at 164000
  //
#if !defined __FreeBSD__ && !defined __osf__
  if ((whence == SEEK_SET) && (beginpos == offset)) {
#ifdef DEBUG_SEEK
    //  This isn't terribly informative, and adds a lot of clutter.
    //fprintf(stderr, "TT_fseek()--  seek to "F_OFF_T" (whence=%d); already there\n", offset, whence);
#endif
    return;
  }
#endif  //  __FreeBSD__

  errno = 0;
  fseeko(stream, offset, whence);
  if (errno) {
    fprintf(stderr, "TT_fseek()--  Failed with %s.\n", strerror(errno));
    assert(errno == 0);
  }

#ifdef DEBUG_SEEK
  fprintf(stderr, "TT_fseek()--  seek to "F_OFF_T" (requested "F_OFF_T", whence=%d) from "F_OFF_T"\n",
          TT_ftell(stream), offset, whence, beginpos);
#endif

  if (whence == SEEK_SET)
    assert(TT_ftell(stream) == offset);
}

