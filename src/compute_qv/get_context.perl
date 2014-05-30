#!/usr/local/gnu/bin/perl   

# This Perl program reads the output file generated from
# the "context" program (context.c).
# It generates a an ascii file with the appropriate table

# get_context.perl  $Revision: 1.3 $

use strict;

sub usage
{
    print( " usage: $0 -dim <dimension> infile > outfile\n" );
}

use Getopt::Long;
my($dim) = 0;
my($ret) = GetOptions ( 'dim=i', \$dim );
if( $dim==0 ) { usage; die( "no dimension argument given\n" ); }

# Search for the beginning of the desired information

# (The desired context length, as specified on commmand line to this script)
my($ok) = 0;
while ($_ = <>)
{
    chomp;
    my(@F) = split;
    my($words);
    $words = @F;
    if ($words > 1 && $F[0] eq 'Context:' && $F[1]==$dim ) 
    {
	# Skip one line (header info)
    	$_ = <>;
    	chomp;
    	my(@H) = split;
    	my($numHeaders);
    	$numHeaders = @H;
    	#debug print("substr($H[0],0,1)=",substr($H[0],0,1),"=\n");
	# Sanity check to see if it is the header
	# (1st non-header line will begin with 'A')
    	if ($numHeaders > 1 && substr($H[0],0,1) ne 'A') 
    	{
		$ok = 1;
		last;
	}
	else
	{
		die("Error:  Column headers are missing from input file.\n");
	}
    }
}
if( !$ok ) { die( "cannot find entries for dimension = ", $dim, "\n" ) };
	     
my($total) = 4**$dim;

################################################################################


my($line) = 0;
while ($_ = <>)
{
    chop;
    my(@F) = split;
    my($words);
    $words = @F;    
    if( $words-$dim < 8 ) { die( "not enough entries\n" ) }
    
    my($entry) = $F[$dim+2];
    
    if( $line == $total-1 ) {
	print( $entry );
	last;
    } else {
	if( $line%16 == 0 ) {
	    print( "\n" );
	}
	print( $entry, " " );
    }
    $line++;
}

