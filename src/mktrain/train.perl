#!/usr/local/gnu/bin/perl -w

# 2.6

# For documentation of "Getopt::Std" 
#   http://www.perl.com/pub/doc/manual/html/lib/Getopt/Std.html
# Example:
#
# use Getopt::Std;
# getopt('cd');
# print( "c=", $opt_c, " d=", $opt_d, "\n" );
 
# For documentation of "Getopt::Long"
#   http://www.perl.com/pub/doc/manual/html/lib/Getopt/Long.html
#
# Example:
# use Getopt::Long;
# $ret = GetOptions ('foo=s', \$foo, 'bar=i', 'ar=s', \@ar);
#
# With command line options ``-foo blech -bar 24 -ar xx -ar yy'' 
#   this will result in: 
#   $foo = 'blech'
#   $opt_bar = 24
#   @ar = ('xx','yy')
# print( "foo=", $foo, " opt_bar=", $opt_bar, " ar=", @ar, "\n" );

# Parse the command line and set default command line values -------------------

use Getopt::Long;
#$ret = 
GetOptions ( 'opts=s', \$opts, 'outtag=s', \$outtag, 'doit', \$doit );

@keywords = @ARGV;
$num_keywords = @keywords;

$default_outtag =  "TT_1_2_2_beta";
if( !$outtag ) { $outtag = $default_outtag; }



# sanity check on arguments ----------------------------------------------------

@legal_keywords = ( 
 "ALL",         "ALL_pop6",    "ALL_pop5",    "ALL_377_dt",
 "baylor",      "celera.ATKS", "celera.CPKH", "drosophila",
 "celera.pop6", "staff.pop6",  "tigr",        "tigr.bac2",  "whitehead",
 "celera.pop5", "staff.pop5",  "tigr.pop5",   "whitehead.pop5",
 "Project_377_dt", "tim-data", "438",         "AL",         "B27"
		    );
# Eliminated: "staff.pop5.2"

if( !&keywords_are_legal( \@keywords, \@legal_keywords ) ) 
{ print("\n"); &usage; }

if( $num_keywords==0 )  {  
    print( "\n ERROR:  No dataset keyword entered.\n\n");  
    &usage;
}


sub usage
{
    print(" usage: $0 [-opts \"<options>\"] " );
    print("[-outtag <output file extension>] keyword#1 ... keyword#n\n" );
    print( "\t -opts           e.g., \"-n\" \n" );
    print( "\t -outtag         default is $default_outtag\n" );
    print( "\t Dataset keywords are:");
    my($i) = 0;
    my($key);
    foreach $key (@legal_keywords) 
    {
	if( $i%5 == 0 ) { print("\n"); }
	print( "\t", $key );
        $i++;
    }
    print( 
   "\n Example: $0 -opts \"-n\" -outtag test.TT_1_2_3 baylor celera.ATKS\n\n");
    exit;
}



#print( "opts=", $opts, " outtag=", $outtag, "\n" );
#print( "ARGV ", @keywords, "\n" );

# Some default values ----------------------------------------------------------
$outdir="/home/spiderman3/train_files";         
#$outdir = "out";

#===============================================================================
# Loop over "scripts" twice:, 1st time just to check stuff out, 
#                             2nd to really do it
#for( $iter=0; $iter<2; $iter++ )
{

    $iter = $doit;
    
    if( !$iter )  {  
	print( "\n\nRe-enter command with -doit to execute the following:\n");
	print("\t(DONE below means command will NOT be executed" );
	print(" because output file already exists):\n");
    }

#  if( $num_keywords==0 && $opts ne "" )
#    {
#      mysystem( $iter, "train $opts", "$outdir/default.$outtag" );
#    }

    # ALL_pop6 -----------------------------------------------------------------

    if( intersection_exists( \@keywords,  "ALL", "ALL_pop6", "baylor" ) )    # 1
    {
	my($projectfile) =
      "/home/spiderman1/customer-data/project_files/baylor_3700_dt.projectfile";

	mysystem( $iter, 
		  "train $opts -p $projectfile", 
		  "$outdir/baylor.3700_dt.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "celera.ATKS" ) )# 2
    {
	my($consensus) =
	    "/home/spiderman1/implementation/celera.ATKS/ATKS.consensus.txt";
	    
	mysystem( $iter, 
		  "train $opts -C $consensus -d " .
		      "/home/spiderman1/implementation/celera.ATKS/ATKS", 
		  "$outdir/celera.ATKS.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "celera.CPKH" ) )# 3
    {
	my($vector)  = "/home/blackstone16/gena/factura2.2b3/vectors/PUC19";
	#my($site)    = "/home/blackstone16/gena/factura2.2b3/sites/BstXI";
	my($datadir) = "/home/spiderman1/implementation/celera.CPKH";
		
	mysystem( $iter,
		"train $opts -C /home/blackstone16/gena/abi/CO*/CPKHfinal.txt" .
		  " -d $datadir/CPKH_1",
		"train $opts -C /home/blackstone16/gena/abi/CO*/CPKHfinal.txt" .
		  "-d $datadir/CPKH_2",
		"$outdir/celera.CPKH.$outtag" );
    }
    
    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "drosophila" ) ) # 4
    {
	my($projectfile) = 
          "/home/spiderman1/customer-data/project_files/drosophila.projectfile";

	mysystem( $iter,
		  "train $opts -p $projectfile",
		  "$outdir/celera.drosophila.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "celera.pop6" ) )# 6
    {
	my($consendir) = "/home/blackstone16/gena/abi/CONSENSUSes";
	# NOTICE: pop5 vs. pop6 ????????
	my($datadir)   = "/home/spiderman2/implementation/celera.pop5";

	mysystem( $iter,
		  "train $opts -C $consendir/POP5_A259H10final.txt" .
		      " $datadir/Pop6*/*ab1",
		  "$outdir/celera.pop6.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "staff.pop6" ) ) # 9
    {
	my($projectfile_F) = 
        "/home/spiderman1/customer-data/project_files/staff.pop6_F.projectfile";
	my($projectfile_R) = 
        "/home/spiderman1/customer-data/project_files/staff.pop6_R.projectfile";
	my($vector)        = 
                      "/home/blackstone16/gena/factura2.2b3/vectors/pUC18staff";
	my($primer_F)      = 
                      "/home/blackstone16/gena/factura2.2b3/primers/PACf";
	my($primer_R)      = 
                      "/home/blackstone16/gena/factura2.2b3/primers/PACr";
	my($site)          = 
                      "/home/blackstone16/gena/factura2.2b3/sites/SmaI";
	
	mysystem( $iter,
		  "train -V $vector -P $primer_F -S $site $opts" .
		       " -p $projectfile_F",
		  "train -V $vector -P $primer_R -S $site $opts" .
		       " -p $projectfile_R",
		  "$outdir/staff.pop6.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "tigr" ) )	    # 10
    {
	my($projectfile) = 
	    "/home/spiderman1/customer-data/project_files/tigr.projectfile";

	mysystem( $iter,
		  "train $opts -p $projectfile",
		  "$outdir/tigr.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "tigr.bac2" ) ) # 11
    {
	my($projectfile) = 
	   "/home/spiderman1/customer-data/project_files/tigr.bac2.projectfile";
	
	mysystem( $iter, 
		  "train $opts -p $projectfile",
		  "$outdir/tigr.bac2.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop6", "whitehead" ) ) # 13
    {
	my($consendir) = "/home/spiderman1/customer-data/CONSENSUSes";
	my($datadir)   = "/home/spiderman1/implementation/whitehead";
	my($vps_dir)   = "/home/blackstone16/gena/factura2.2b3";
	my($vectors)   = "$vps_dir/vectors/PUC18";
	my($primers)   = "$vps_dir/primers/M13-21";
	my($sites)     = "$vps_dir/sites/SmaI";

	# I (S. Tyler) changed argument order to be consistent with usage msg.
	# I removed the "-r" option (I believe this has changed to -recalln)
	# I talked this over with David Ho and we mutually agreed to remove it
	mysystem( $iter,	
		  "train -V $vectors -P $primers -S $sites $opts" .
		      " -C $consendir/L690final.txt -d $datadir/L690_1",
		  "train -V $vectors -P $primers -S $sites $opts" .
		      " -C $consendir/L690final.txt -d $datadir/L690_2",
		  "train -V $vectors -P $primers -S $sites $opts" .
		      " -C $consendir/L723final.txt -d $datadir/L723_1",
		  "train -V $vectors -P $primers -S $sites $opts" .
		      " -C $consendir/L723final.txt -d $datadir/L723_2",
		  "$outdir/whitehead.$outtag" ); 
	# NOTICE: was "$TF3/whitehead.$outtag"
    }

    # ALL_pop5 -----------------------------------------------------------------

    if( intersection_exists( \@keywords, "ALL", "ALL_pop5", "celera.pop5" ) )# 5
    {
	my($consendir)    = "/home/blackstone16/gena/abi/CONSENSUSes";
	my($datadir)      = "/home/spiderman2/implementation/celera.pop5";
	my($short_vector) =
	 "/home/spiderman2/implementation/celera.pop5/celera.pop5.short_vector";
	
	mysystem( $iter,
		  "train " . #-V $short_vector $opts" .
		      " -C $consendir/POP5_A259H10final.txt" .
		      " -d $datadir/pop5_707_BC-POP5LR",
		  "train " . #-V $short_vector $opts" .
		      " -C $consendir/POP5_A259H10final.txt" .
		      " -d $datadir/pop5_719_BC-POP5LR",
		  "$outdir/celera.pop5.POP5LR.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL", "ALL_pop5", "staff.pop5" ) ) # 7
    {
	my($projectfile_F)  = 
         "/home/spiderman2/implementation/staff.pop5/staff.pop5_F.projectfile";
	my($projectfile_R)  =
         "/home/spiderman2/implementation/staff.pop5/staff.pop5_R.projectfile";
	my($short_vector_F) =
         "/home/spiderman2/implementation/staff.pop5/staff.pop5_F.short_vector";
	my($short_vector_R) =
         "/home/spiderman2/implementation/staff.pop5/staff.pop5_R.short_vector";

	mysystem( $iter,
		  "train -V $short_vector_F $opts -p $projectfile_F",
		  "train -V $short_vector_R $opts -p $projectfile_R",
		  "$outdir/staff.pop5.$outtag" );
    }

    ### Commented this out at the suggestion of David Ho. (the data is no good)
    #if( intersection_exists( \@keywords, "ALL", "ALL_pop5", "staff.pop5.2" ))#8
    #{
    #   my($datadir) = 
    #	    "/home/spiderman2/implementation/staff.pop5.2/P0698-3-A04-ver1";
    #
    #	mysystem( $iter,
    #		  "train $opts" .
    #		    " -C /home/blackstone16/gena/abi/CO*/P0698-3-A04final.txt" .
    #		    " -d $datadir",
    #		  "$outdir/staff.pop5.2.$outtag" );
    #}

    if( intersection_exists( \@keywords, "ALL", "ALL_pop5", "tigr.pop5" ) ) # 12
    {
	my($projectfile_F)  =
	    "/home/spiderman4/tigr.pop5/tigr.pop5.F.projectfile";      
	my($projectfile_R)  = 
	    "/home/spiderman4/tigr.pop5/tigr.pop5.R.projectfile";
	my($short_vector_F) = 
	    "/home/spiderman4/tigr.pop5/pHOS2_F.short_vector";
	my($short_vector_R) =
	    "/home/spiderman4/tigr.pop5/pHOS2_R.short_vector";

	mysystem( $iter,
		  "train -V $short_vector_F $opts -p $projectfile_F",
		  # NOTICE: "test" was appended to file name
		  "train -V $short_vector_R $opts -p $projectfile_R",
		  "$outdir/tigr.pop5.$outtag" );
    }

    if( intersection_exists( \@keywords, "ALL","ALL_pop5","whitehead.pop5"))# 14
    {
	my($consensus) = "/home/blackstone16/gena/abi/CO*/L9*txt";
	my($datadir)   = "/home/spiderman2/implementation/whitehead.pop5/L906";

	mysystem( $iter,
		  "train $opts -C $consensus -d $datadir",
		  "$outdir/whitehead.pop5.$outtag" );
    }

    # added 8/21/2001
    if( intersection_exists( \@keywords, "ALL", "ALL_377_dt","Project_377_dt")) 
    {
	mysystem( $iter, "train " .
		  " -p /home/spiderman2/implementation/baylor/Project_377_dt",
		  "$outdir/Project_377_dt.$outtag" );
    }
    
    if( intersection_exists( \@keywords, "ALL", "ALL_377_dt", "tim-data" )) 
    {
	my($consendir) = "/home/blackstone16/gena/abi/CONSENSUSes";
	my($spi4data)  = "/home/spiderman4";
	
	mysystem( $iter, 
              "train " .
		  " -C $consendir/114N19final.fasta" .
		  " -d $spi4data/tim-data/projects/114N19/chromat_dir", 
              "train " .
		  " -C $consendir/201F1final.fasta" .
		  " -d $spi4data/tim-data/projects/201F1/chromat_dir", 
              "train " .
		  " -C $consendir/314E14final.fasta" .
		  " -d $spi4data/tim-data/projects/316E14/chromat_dir", 
	      "$outdir/tim-data.$outtag" );
    }

    {
	my($consendir) = "/home/blackstone16/gena/abi/CONSENSUSes";
	my($spi1data)  = "/home/spiderman1/implementation";

	if( intersection_exists( \@keywords, "ALL", "ALL_377_dt", "438" )) 
	{
	    mysystem( $iter, "train " .
		      " -C $consendir/438final.txt" .
		      " -d $spi1data/377/438/438x_583", 
		      "$outdir/438.$outtag" );
	}

	if( intersection_exists( \@keywords, "ALL", "ALL_377_dt", "AL" )) 
	{
	    mysystem( $iter, "train " .
		      " -C $consendir/ALfinal.txt" .
		      " -d $spi1data/377/AL.chromat/ALx_958", 
		      "$outdir/AL.$outtag" );
	}

	if( intersection_exists( \@keywords, "ALL", "ALL_377_dt", "B27" )) 
	{
	    mysystem( $iter, "train " .
		      " -C $consendir/B27final.txt" .
		      " -d $spi1data/377/B27", 
		      "$outdir/B27.$outtag" );
	}

    }

    #if( $iter==0 )
    if( 0 )
    {
	print(
	 "Do you agree with the above set of commands to be executed? (y/n)");
	$ans = <STDIN>;
	chomp( $ans );
	if( $ans ne 'y' )
        {
	    print( "try again\n" );
	    exit;
        }
    }
}# Loop over "scripts"
#===============================================================================


# &mysystem( 0, "pwd", "ls -a -l", "g.dat" ); # debug test 

sub keywords_are_legal
{
    my($at,$bt) = @_;
    my(@a,@b);
    @a = @$at;
    @b = @$bt;
    #print( "a=", @a, "\n" );
    #print( "b=", @b, "\n" );
    my($ok) = 1;
    foreach $a_elem (@a) 
    {  
	if( !&element_of( $a_elem, @b ) ) 
        {
	    print( " ERROR: ", $a_elem, " is not a legal keyword\n" );
	    $ok = 0;
        }
    }
    return $ok;
}


sub intersection_exists # (@@)   \@a, @b   # my attempt at prototypes
{
    my($at,@bt) = @_;
    my(@a,@b);
    @a = @$at;
    @b = @bt;
    #print( "a=", @a, "\n" );
    #print( "b=", @b, "\n" );
    foreach $a_elem (@a) {  if( &element_of( $a_elem, @b ) ) { return 1; } }
    return 0;
}


sub element_of
{
    my($val) = shift(@_);
    my(@array) = @_;
    foreach $elem (@array)  {  if( $elem eq $val ) { return 1; }  }
    return 0;
}


sub mysystem
{
    $isize = @_;
    if( $isize<3 ) { die( "isize<3\n" ) }
    $flag     = shift(@_);
    $out_file = pop(@_);
    $isize = @_;
    if( !(-e $out_file) )		# Does file exist?
    {
	for( $i=0; $i<$isize; $i++ )
        {
	    if( $i==0 ) { $cmd = $_[$i] . " >  " . $out_file; }
	    else        { $cmd = $_[$i] . " >> " . $out_file; }
	    if( !$flag )  { print( "     ", "$cmd\n" ); }
	    else          { system( $cmd ); }
        }
    }
    else
    {
	for( $i=0; $i<$isize; $i++ )
        {
	    if( $i==0 ) { $cmd = $_[$i] . " >  " . $out_file; }
	    else        { $cmd = $_[$i] . " >> " . $out_file; }
	    if( !$flag )  { print( "DONE:", "$cmd\n" ); }
        }
    }   
}

