#!/usr/local/gnu/bin/perl

#
# Purpose: extract the names of files containing base calls
#          with extremely high trace parameters 
#          or extremely low iheight from training file(s)
# Usage: 
#   extract_bad_files.perl < training_file > out_file
# or
#   cat <list_of_training_files> | extract_bad_files.perl > out_file        
#

$line = 0;
$phr_max = 1600.;
$psr_max = 100;
$pres_max= 10.;
$iheight_min = 5;

while ($_ = <>) {
    $line++;
    chop;
    @F = split;
    $words = @F;

    if (($F[0] eq '#') && ($F[1] eq 'Version')) {
        $version_line = $line;
    }

    if (($F[0] eq '#') && ($F[1] eq 'File:')) {
        $file_name = $F[2];
    }

    if (($F[0] eq '#') && ($F[1] eq 'Consensus') && ($line == $version_line +7)) {
        print("Consensus == $F[3]\n"); 
    }
    
    if (($words == 10) && ($F[5] > $phr_max)) {
        print("$file_name:\nbase=$F[1]/$F[0] match=$F[4] phr3 = $F[5]\n");
    }
    
    if (($words == 10) && ($F[6] > $phr_max)) {
        print("$file_name:\nbase=$F[1]/$F[0] match=$F[4] phr7 = $F[6]\n");
    }

    if (($words == 10) && ($F[7] > $psr_max)) {
        print("$file_name:\nbase=$F[1]/$F[0] match=$F[4] psr7 = $F[7]\n");
    }

    if (($words == 10) && ($F[8] > $pres_max)) {
        print("$file_name:\nbase=$F[1]/$F[0] match=$F[4] pres = $F[8]\n");
    }

    if (($words == 10) && ($F[9] < $iheight_min)) {
        print("$file_name:\nbase=$F[1]/$F[0] match=$F[4] iheight = $F[9]\n");
    }

}
