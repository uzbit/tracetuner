#!/usr/local/bin/perl

# Purpose: generate a consensus sequence from the reference sequence and alignment
use strict;
use warnings;

(@ARGV == 2) or 
die "Usage: GetVARReport.pl <reference_fasta> <trainphd_out>\n"; 

# Parse the input string:
my $ref_file  = $ARGV[0];
my $trainphd_out  = $ARGV[1];

my $header;
my $right_donor = 0;
my $direction = 0;
my $read = "";

my $prev_header = "#";
my $alignment_len = 0;
my @lines;
my @rclines;
my $base;
my $QV;
my $prev_QV;
my $sind = 0; # string index
my $MAX_QV = 100;

# Read the var file and create new consenssus
my $reference_seq = read_fasta_file($ref_file);
#

my @ref_array = split "", $reference_seq;
my $ref_len = length($reference_seq);
my @variants;  # array of references to arrays
               # $variants[$i]->[k] is the reference to a hash 
               # which stores a sum of QVs of different bases 
               # at k-th position of the variant consensus sequence
               # corresponding to i-th base of the reference sequence
my @consensus; # array of references to hashes

# Initialize consensus
for (my $i=0; $i<$ref_len; $i++)
{
    initialize_hash_ref($variants[$i]->[0]);
    $variants[$i]->[0]->{$ref_array[$i]} = 0;
    $consensus[$i]->{"len"} = 1;
    $consensus[$i]->{"string"} = $ref_array[$i];
}

# Process trainphd file
open(TPHD, "$trainphd_out");
while ($_ = <TPHD>)
{
    chop;
    my @T = split;
    my $dim = @T;
 
    if ($dim > 0) { $header = $T[0]; }
    if ($dim > 1 && $header eq "#")
    {
        if ($T[1] eq "File:") { $read = $T[2]; }
        if ($dim > 2 && ($T[1] eq "Num" || $T[1] eq "No"))
        {
            if ($alignment_len > 0)
            {
                # Process previously stored alignment
                my @first= split " ", $lines[0];
                my @last = split " ", $lines[$alignment_len-1];
                if ($first[0] < $last[0]) { $direction =  1; }
                else                      { $direction = -1; }
#                 print "read= $read direction = $direction align_len= $alignment_len\n";
                if ($direction > 0)
                {
                    store_alignment($alignment_len, \@lines, \@variants);
                }
                else
                {
                    reverse_complement_alignment($alignment_len, \@lines, \@rclines);
                    store_alignment($alignment_len, \@rclines, \@variants);
                }
            }

            $right_donor = 1;
            $alignment_len = 0;
        }
    }
    else
    {
        if ($dim > 1)
        {
            $lines[$alignment_len] = $_;
            $alignment_len++;
        }
    }
}
close(TPHD);

# Determine and report the consensus sequence and quality
# print "Determining the consensus sequence\n";
print "Ref.pos.\tRef.base\tCons.base\tCons.QV\t\tAlt.cons.base\n";

for (my $i=0; $i<$ref_len; $i++)
{
    my $ip1 = $i+1;
    my @array = ("A", "C", "G", "T", "M", "S", "K", "R", "W", "Y");
    my $clen = $consensus[$i]->{"len"};
    my $QVstring = "";
    my $minQV = 100;
    
    my ($cbase, $cQV, $cbase2) = get_consensus_base_and_QV($variants[$i]->[0]);

    $consensus[$i]->{"string"} = $cbase;
    $QVstring .= "$cQV ";
    if ($minQV > $cQV) { $minQV = $cQV; }
    if ($clen > 1)
    { 
        for ( my $j=1; $j < $clen; $j++)
        {
            ($cbase, $cQV, $cbase2) = get_consensus_base_and_QV($variants[$i]->[$j]);
            if ($cbase ne "") { $consensus[$i]->{"string"} .= $cbase; }
            if ($cQV   >  0 ) { $QVstring .= "$cQV ";  }
        }
    }
    my $cstring = $consensus[$i]->{"string"};
#   print "    i= $ip1 clen= $clen  ref_base= $ref_array[$i] consensus_string= $cstring\n";
    if ($consensus[$i]->{"string"} ne $ref_array[$i] || ($minQV < 20 && $minQV > 0))
    {
       $cbase = $consensus[$i]->{"string"};
       if ($minQV >= 20) {
           print "$ip1\t\t$ref_array[$i]       \t$cbase   \t\t$QVstring\n";
       }
       elsif ($minQV > 0 && $minQV < 20)
       {
           if (length($cbase2) > 0) {
               print "$ip1\t\t$ref_array[$i]       \t$cbase   \t\t$QVstring\t\t$cbase2\n";
           }
           else 
           {
               print "$ip1\t\t$ref_array[$i]       \t$cbase   \t\t$QVstring\n";
           }
       }
    }
}

sub mixed_base
{
    my ($b1, $b2) = @_;
    if ($b1 eq $b2) { return $b1; }
    if ($b1 eq "A" && $b2 eq "C") { return "M"; }
    if ($b1 eq "A" && $b2 eq "G") { return "R"; }
    if ($b1 eq "A" && $b2 eq "T") { return "W"; }
    if ($b1 eq "C" && $b2 eq "G") { return "S"; }
    if ($b1 eq "C" && $b2 eq "T") { return "Y"; }
    if ($b1 eq "G" && $b2 eq "T") { return "K"; }
    return "N";
}

sub complement_base
{
    my $b = shift;
    
    if ($b eq "A") { return "T"; }
    if ($b eq "C") { return "G"; }
    if ($b eq "G") { return "C"; }
    if ($b eq "T") { return "A"; }
    if ($b eq "K") { return "M"; }
    if ($b eq "M") { return "K"; }
    if ($b eq "R") { return "Y"; }
    if ($b eq "S") { return "W"; }
    if ($b eq "W") { return "S"; }
    if ($b eq "Y") { return "R"; }
    if ($b eq "-") { return "-"; }
#   if ($b eq "N") { return "N"; }
    return "N";
}

sub get_position
{
    my $string = shift;
    my @array = split "", $string;
    my $len = @array;
    my $pos = 0;
    my $i = 0;
    while ($array[$i] eq "0") { $i++; }
    $pos = substr($string, $i, $len - $i);
    return $pos;
}

sub read_fasta_file
{
    my $file_name = shift;
    my $sequence = "";
    open(RFILE, $file_name);
    while ($_=<RFILE>)
    {
        chomp;
        if (substr($_, 0, 1) eq ">") {}
        else {
           $_ =~ tr/a-z/A-Z/;
           $sequence .= $_;
        }
    }
    close(RFILE);
    return $sequence;
}

sub initialize_hash_ref
{
    my $ref = shift;
    $ref->{"A"} = 0;
    $ref->{"C"} = 0;
    $ref->{"G"} = 0;
    $ref->{"T"} = 0;
    $ref->{"K"} = 0;
    $ref->{"M"} = 0;
    $ref->{"R"} = 0;
    $ref->{"S"} = 0;
    $ref->{"W"} = 0;
    $ref->{"Y"} = 0;
    $ref->{"-"} = 0;
    $ref->{"N"} = 0;
}

sub get_consensus_base_and_QV 
{
    my $ref = shift;
    my %hash = %$ref;
    my $cbase = "";
    my $cbase2 = "";
    my $cQV = 0;
    my $bestQVsum  = 0;
    my $best2QVsum = 0;

    foreach my $b (keys %hash)
    {
        if ($bestQVsum < $hash{$b})
        {
            $best2QVsum = $bestQVsum;
            $bestQVsum  = $hash{$b};
            $cbase2  = $cbase;
            $cbase   = $b;
            $cQV     = $bestQVsum - $best2QVsum;
        }
        elsif ($bestQVsum >= $hash{$b} && $best2QVsum < $hash{$b}) 
        {
            $best2QVsum = $hash{$b};
            $cbase2  = $b;
            $cQV     = $bestQVsum - $best2QVsum;
        }
    }
    if ($cQV > $MAX_QV) { $cQV = $MAX_QV; }
    return ($cbase, $cQV, $cbase2);
}

sub min
{
    my ($a, $b) = @_;
    if ($a < $b) { return $a; }
    return $b;
}

sub store_alignment
{
    my ($alen, $ref_lines, $ref_variants) = @_;
    for (my $i=0; $i<$alen; $i++)
    {
        my @T1 = split "\t", $ref_lines->[$i];
        my $dim1 = @T1;
        my $rind = $T1[0] - 1; # base index in the reference seq
        my $sumQVs;
        if ($T1[1] ne "-" && $T1[4] ne "-") # match or substitution
        {
            $sind = 0;
            $base = $T1[4];
            $QV   = $T1[5];
            $ref_variants->[$rind]->[$sind]->{$base} += $QV;
            if ($T1[1] ne $T1[4])
            {
                $sumQVs = $ref_variants->[$rind]->[$sind]->{$base};
            }
            if ($dim1 > 6) {
                for (my $j=6; $j<$dim1; $j += 2)
                {
                    my $alt_base = $T1[$j];
                    my $alt_QV   = $T1[$j+1];
                    if (!defined $ref_variants->[$rind]->[$sind]->{$alt_base})
                    {
                        $ref_variants->[$rind]->[$sind]->{$alt_base} = 0;
                    }
                    $ref_variants->[$rind]->[$sind]->{$alt_base} += $alt_QV;
                }
            }
            $prev_QV = $QV;
        }
        else {
            if ($T1[4] eq "-") # deletion in a read
            {
                $sind = 0;
                $QV = $prev_QV;
                my $next_QV = -1;
                my $j = 0;
                while ($next_QV < 0 && $i + $j < $alen)
                {
                    my @T2 = split "\t", $ref_lines->[$i + $j];
                    if ($T2[4] ne "-") { $next_QV = $T2[5]; }
                    $j++;
                }
                if ($next_QV < 0) { $next_QV = $prev_QV; }
                $ref_variants->[$rind]->[$sind]->{"-"} += min($prev_QV, $next_QV);
            }
            else {             # insertion in a read; ref. base is '-'
                $sind++;
                if ($consensus[$rind]->{"len"} < $sind + 1)
                {
                    $consensus[$rind]->{"len"} = $sind + 1;
#                   print "   pos= $T1[0]  sind = $sind\n";   
                }
                $base = $T1[4];
                $QV   = $T1[5];
                if (!defined $ref_variants->[$rind]->[$sind]) {
                    initialize_hash_ref($ref_variants->[$rind]->[$sind]);
                }
                $ref_variants->[$rind]->[$sind]->{$base} += $QV;
                if ($dim1 > 6) {
                    for (my $j=6; $j<$dim1; $j += 2)
                    {
                        my $alt_base = $T1[$j];
                        my $alt_QV   = $T1[$j+1];
                        $ref_variants->[$rind]->[$sind]->{$alt_base} += $alt_QV;
                    }
                }
                $prev_QV = $QV;
            }
        }
    }
}

sub reverse_complement_alignment
{
    my ($alen, $ref_lines, $ref_rclines) = @_;
    my $i;
    my $j;
    for ($i=0; $i<$alen; $i++)
    {
        $j= $alen-1-$i;
        my @T = split "\t", $ref_lines->[$i];
        my $dim = @T;
        my $rpos  = ($T[1] eq "-") ? $T[0]-1 : $T[0];
        my $rbase = complement_base($T[1]);
        my $qbase = complement_base($T[4]);
        if (defined $T[5]) {
            $ref_rclines->[$j] = "$rpos\t$rbase\t$T[2]\t$T[3]\t$qbase\t$T[5]"; 
        }
        else {
            $ref_rclines->[$j] = "$rpos\t$rbase\t$T[2]\t$T[3]\t$qbase";
        }
        if ($dim > 6)
        {
            for (my $k=6; $k<$dim; $k += 2)
            {
                $qbase = complement_base($T[$k]);
                $ref_rclines->[$j] .= "\t$qbase\t$T[$k + 1]";
            }
        }
    }
}
