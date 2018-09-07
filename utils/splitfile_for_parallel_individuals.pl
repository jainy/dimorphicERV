#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  May 2017
# email  :  jainythomas1@gmail.com
# Pupose :  split a file to multiple files and make the list in another file
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Copy;
use File::Path qw(make_path remove_tree);
use Data::Dumper;

my $version = "4.0";
my $scriptname = "splitfile_jt_v4.0.pl";
my $changelog = "
#   - v1.0 = 21 May 2017
#	- v2.0 = 17 June 2017
#  				added the option sort the files based on the individual name
#	- v3.0 = 29 August 2017
#				split based on the number of individuals is an option
#	- v4.0 = 7 September 2018
#				deprecated the split by line or bytes option, retaining only split by individual option
#
#
#
#
#
\n";
my $usage = "\nUsage [$version]: 
    perl $scriptname -f <file that needs to be split> -s <sort the file or not, give yes or no> -n <number of individuals> [-o name of outputfile with the list of files][-pr name of the splitfiles] [-v] [-c] [-h] 
	
	MANDATORY ARGUMENT:	
 
    -f,--file  (string) file 
    -n, --number of individuals (STRING) split by the number of individuals	
    -s,--presort (STRING) sort the file based on the individual name, also prints the sorted files  
    
    OPTIONAL ARGUMENTS:
    
    -pr,--prefix (STRING) name of the split files(aa,ab,ac will be added to the prefix) available only for lines/bytes splitting  
    -p, --path   (STRING)  path where the splitbyindividuals folder is created                
    
    -c,--chlog   (BOOL)   Print updates
    -v,--v       (BOOL)   Print version if only option
    -h,--help>   (BOOL)   Print this usage\n\n";
   
#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$line,$bytes,$prefix,$presort,$path,$indinumber,$out,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            'l=s' => \$line,
            'b=s' => \$bytes,
           'pr=s' => \$prefix,
            's=s' => \$presort,
            'o=s' => \$out,
            'n=s' => \$indinumber,
            'p=s' => \$path,	
            'c'   => \$chlog, 
            'h'   => \$help,
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $file)&&  (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $file) ||  ($help));
$prefix = "$file.sorted.split" if ((! $prefix) && ($presort));
$out = "$prefix.sorted.list.txt" if ((! $out) && ($presort)) ;
$prefix = "$file.split" if (! $prefix);
$out = "$file.split.list.txt" if (! $out);
$path = "./splitbyindividuals" if (! $path);
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my %sortindividual;
my $sortedfile;
my @listofsplitfiles;
if ($presort eq "yes") {
	print STDERR "the files is sorting now based on the file name";
	$sortedfile = "$file.sorted.txt";
	&load_hash();
	&sorthashbyvalue();
	print STDERR "Starting to split based on individuals....\n";
	&split_byindivi($indinumber) if ($indinumber);
}
if ($presort eq "no") {
		
	if ($indinumber) {
		print STDERR "splitting based on individuals....\n";
		&load_hash();
		&split_byindivi($indinumber);
	}
}
print STDERR "DONE..\n";
exit;
#-------------------------------------------------------------------------------
#----------------------------------- SUB ---------------------------------------
#-------------------------------------------------------------------------------


sub split_byindivi {
	
	my ($indinu) = shift;
	my $currentindividual;
	my %keyindividual;
	my $numofindi;
	my @elementstoprint;
	my $count =0;
	my $lastindivi;
	
	make_path ("$path") if ($path);
	foreach my $name (sort { $sortindividual{$a} cmp  $sortindividual{$b}} keys %sortindividual) {
		$currentindividual = $sortindividual{$name};
		
		push (@elementstoprint,$name);
		#print $fp "$name\n";
		$keyindividual{$currentindividual} = 1;
		$numofindi = keys (%keyindividual);
		
		if ($numofindi > "$indinu") {
			$count++;
			#print "$count\n";
			$lastindivi = pop (@elementstoprint);
			%keyindividual= ();
			&print_array($count,@elementstoprint);
			@elementstoprint =();
			push (@elementstoprint,$lastindivi);
		} else {
			
			#print Dumper %keyindividual;
			next;
		}	
	}
	$count++;
	#print "final count is $count\n";
	&print_array($count,@elementstoprint);
	&print_filename(@listofsplitfiles);
}

sub load_hash {
	open (my $fh,"<",$file) or die ("cannot open file $file to read $!\n");
	while (<$fh>) {
		chomp (my $data = $_);
		my @col = split (/\s+/,$data);
		#my @firstcol = split (/\./,$col[0]);
		my $individual = $col[0];
		$sortindividual{$data} = $individual;
	}
	return (%sortindividual);
}

sub sorthashbyvalue {
	open (my $hp,">","$file.sorted.txt") || die ("failed to open file to write (sub sorthashbyvalue) sorted file $!\n");
	foreach my $name (sort { $sortindividual{$a} cmp  $sortindividual{$b}} keys %sortindividual) {
		print $hp "$name\n";
	}
	close $hp;
}

sub print_array {
	my ($num,@array) = @_;
	open (my $ah,">","$path/$file.$num.individuals.sorted.txt") || die ("cannot open file (sub print_array) to write $!\n");
	push (@listofsplitfiles,"$file.$num.individuals.sorted.txt");
		foreach my $dataline (@array) {
			print $ah "$dataline\n";
		}
	close ($ah);
}

sub print_filename {
	my @list = @_;
	my @sortedsplitfiles = sort (@list);
	#print "the array contains @sortedsplitfiles\n";	
	open (my $fhout,">",$out) || die ("cannot open file $out to write (sub print_filename) $!\n");
		foreach my $filename (@sortedsplitfiles) {
			print $fhout "$filename\n";
		}
	close $fhout;
}