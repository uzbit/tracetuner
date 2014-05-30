#!/usr/local/gnu/bin/perl   

#
#  lut_to_default_lut.perl $Revision: 1.2 $
#

# This Perl program converts a lookup table from ASCII
# format into a fragment of C code 

printf("\/**    Copyright (c) 1999 Paracel Inc.  All rights reserved.\n");
printf(" ** \n");
printf(" **    Btk_default_table.c \n");
printf(" ** \n");
printf(" **    \$ Revision:  \$ \n");
printf(" **\/  \n\n");

print("\#include \"Btk_lookup_table.h\"\n");
print("static BtkLookupEntry DefaultTableEntries\[\] = \{\n\n");
$line = 0;

while ($_ = <>)
{
   $line++;
   chop;
   @F = split;
   $words = @F;

# Print all
  if ($words == 5)
  {
      print("\{$F[1],  $F[2],  $F[3],  $F[4],  $F[0]\}, \n");
  }

}
print("\};\n\n");

print("static BtkLookupTable DefaultTable = \{ \
    sizeof\(DefaultTableEntries\) \/ sizeof\(DefaultTableEntries\[0\]\), \
    DefaultTableEntries \
\};\n\n \
BtkLookupTable *  \
Btk_get_default_table\(void\) \
\{ \
return \&DefaultTable; \
\}\n");

