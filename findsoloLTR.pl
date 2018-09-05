#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  July 2017
# email  :  jainythomas1@gmail.com
# Pupose :  to predict the genotype of the TE based on its location and read coverage
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Data::Dumper;
#use MIME::Lite;
use DBI;


my $version = "3.2";
my $scriptname = "verifydeletions.pl";
my $changelog = "
#   - v1.0 = 30 July 2017 
#	- v2.0 = 10 August 2017
#				fixed a bug in the pipeline when the read depth is zero for the flanking region
#	- v3.1 = included email, changed the ratio for calling heterozygote to equal or less than 50,changed how diffratio is calcuated	
#	- v3.2 = 15 May 2018
#				included mappability scores calculation,using read depth estimation of reads with more than 99.99% probability of mappabablity (>=30 score)
#				
#
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <table> -f <file with ltr cordinates> -bl <location of bamfiles> [-m yes to find mappabilty] [-u Username] [-pd password][-db mysql database][-mt mysql table][-p <path of the outputdirectory>][-o <output file>] [-v] [-c] [-h] 
	
    MANDATORY ARGUMENT:	
    -t,--table (string) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file  (string) file containing accesion information (output from the script get_coverage_coordinates.pl script on bedtools)
    -bl,--bamlocation (STRING) location of bam files
    	  
    OPTIONAL ARGUMENTS:
    -m, --mappability  (STRING)  yes if need to find mappability; no if not needed
    -mt,--table 		(STRING) mysql table e.g.	hg19wgEncodeCrgMapabilityAlign100mer_index
    												wgEncodeCrgMapabilityAlign100merhg38_lo_index
    -p,--path         	(STRING) output directory name (path)
                            	 Default = <current working directory>
    -db --mysqldbinfo 	(STRING) ex. command: DBI:mysql:database=jainys_db;host=localhost;port=22;     
    -u --user 			(STRING) Username for mysql databasemy e.g	jainy;
    -pd,--password 		(STRING) password for mysql database e.g. wysql123
    -p,--path   (STRING) output directory name (path)
                         Default = <current working directory>
    -o,--output (STRING) output file
    -i,--igv    (BOOL)   get IGV files for the regions with 250 bp flank
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$path,$rawout,$table,$verbose,$mappability,$bamlocation,$mysqltable,$mysqldb,$user,$password,$igv,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            'p=s' => \$path,
            'o=s' => \$rawout,
            't=s' => \$table,
            'm'   => \$mappability,
            'bl=s'=> \$bamlocation,
            'mt=s'=> \$mysqltable,
            'db=s'=> \$mysqldb,
            'u=s' => \$user,
            'pd=s'=> \$password,
            'i'	  => \$igv,
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table) && (! $file) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $table) || (! $file) ||  ($help));
my $cwd = getcwd();
$path = $cwd if (!$path) ;
$rawout = "$path/$file.readdepth.output.txt" if (!$rawout);
my $logfile = "$path/$file.log";
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#my $bamlocation = "/kbod2/WGS_DATA/SGDP_bams_public";#vader server

my %bamfile = ();
my $bamid;
my $elegenomeloc;
my $fivepgenomeloc;
my $threepgenomeloc;
my $readdepthpath;
my %allreaddepth;
my @ratios = ();
my $eledepth;
my $fivepdepth;
my $threepdepth;
my $uniqueid;
my $individual;
my $IGVloc;
my @logvalues;
my %cordinatescore;
my $dbh;

%bamfile = load_file ($table);

if ($mappability) {
	
	my $dsn = "DBI:mysql:database=$mysqldb;host=localhost;port=22";
	$dbh = DBI->connect($dsn, $user, $password,
						{'RaiseError' => 1});

}
open (my $fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\s+/,$line);
		$individual = $col[0];
		my $chr = $col[1] ;
		my $elestart = $col[2];
		my $eleend = $col[3];
		my $length = ($eleend- $elestart);
		my $fivepstart = $col[2] - 250;
		my $fivepend = $col[2] - 3;
		my $threepend = $col[3] + 250;
		my $threepstart = $col[3] + 3;
		$IGVloc = $chr.":".$fivepstart."-".$threepend;
		$elegenomeloc = $chr.":".$elestart."-".$eleend;
		my $avgscore = &store_mapscores($elegenomeloc) if ($mappability);
		my $gmscore = $length."=>".$avgscore if ($mappability);
		$fivepgenomeloc = $chr.":".$fivepstart."-".$fivepend;
		$threepgenomeloc = $chr.":".$threepstart."-".$threepend;
		my $tetype = $col[4];
		&find_bamid();
		$uniqueid = $individual.".".$tetype.".".$elegenomeloc;
		#get IGV
		&get_IGV if ($igv);
		# get read depth
		@ratios = ();
		make_path  ("$path/ReadDepth/$tetype");
		$readdepthpath = "$path/ReadDepth/$tetype";
		$eledepth = &avg_readdepth($elegenomeloc);
		push (@ratios,{'element' => $eledepth});
		$fivepdepth = &avg_readdepth($fivepgenomeloc);
		push (@ratios,{'5pflank' => $fivepdepth});
		$threepdepth = &avg_readdepth($threepgenomeloc);
		push (@ratios,{'3pflank' => $threepdepth});
		push (@ratios,{'gmscore' => $gmscore});
		$allreaddepth{$uniqueid} = [@ratios];
		
		
	}

&print_hash;
if (@logvalues) {
	print STDERR "Please check logfile for values for which no read depth was calculated...\n";
	&print_array($logfile,@logvalues);
}
#&email;
exit;
#-----------------------------------------------------------------------------
#----------------------------------- SUB -------------------------------------
#-----------------------------------------------------------------------------
sub load_file {
	my ($file1) = @_;
	my %sbamfile;
	open (my $th, "<", $file1) or confess "\n ERROR (main): could not open to read $file1 $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @namebam = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
			my $name = $namebam[0];
			$sbamfile{$name} = $namebam[1]; #loading to the hash
		}
	return (%sbamfile);
	close $th;	
}
sub find_bamid {
	if (exists ($bamfile { $individual } ))  {
		$bamid = $bamfile{$individual};
		print STDERR " the file now analysing is $bamid \n" if ($verbose);
	}
	else {
		die "There is someting wrong with the name, cant locate the bamid \n";
	}
}
sub avg_readdepth {
	my ($genomeloc) = @_;
	my $avg;
	my $total = 0;
	my $count = 0;
	system("samtools depth -a -r $genomeloc -q 20 -Q 30 $bamlocation/$bamid.38.sorted.bam > $readdepthpath/$individual.$genomeloc.readdepth.txt") == 0 or die ("failed extract readdepth via samtools\n"); 
	my $readdepthfile = "$readdepthpath/$individual.$genomeloc.readdepth.txt";
	my $size = (-s $readdepthfile);
	if ($size == 0) {
		push (@logvalues,$readdepthfile);
		$avg = 0;
		return ($avg);
	} else {
		open (my $rh, "<", "$readdepthfile") or confess "\n ERROR (main): could not open to read $readdepthfile $!\n";
			while(<$rh>) {
				chomp(my $line = $_);
				$count++ if  (!/^\s+?$/);#count non-blank lines
				my @columns = split(/\s+/,$line);
				$total += $columns[2];
			}
			$avg = $total/$count;
		return ($avg);
		close ($rh);
	} 
}
sub get_IGV {
	make_path  ("$path/IGV/$individual");
	system("samtools view -b -o $path/IGV/$individual/$uniqueid.bam $bamlocation/$bamid.38.sorted.bam $IGVloc") == 0 or die ("unable to run command on $uniqueid \n");
	print STDERR " 	Extracting 	$uniqueid using samtools done\n" if ($verbose);
	system ("samtools index -b $path/IGV/$individual/$uniqueid.bam ") == 0 or die ("unable to create index file of $uniqueid.bam \n");
	print STDERR " indexing the bamfile done \n" if ($verbose) ;
}

sub print_hash {
	open (my $fhout,">","$rawout") || die ("cannot open file $rawout to write $!\n");
	if ($mappability) {
		print $fhout "name\t5'Read_depth\telement_Read_depth\t3'Read_depth\tpredicted_genotype\t%Read_depth\tgmap_score_avg.weightedavg\n" ;
	} else {
		print $fhout "name\t5'Read_depth\telement_Read_depth\t3'Read_depth\tpredicted_genotype\t%Read_depth\n" ;
	}
	foreach my $name (sort keys %allreaddepth) {
		#print "$name\n";
		my $fivepratio;
		my $threepratio;
		my $elementratio;
		my $avgflankratio;
		my $predigenoty;
		my $fiftypflaratio;
		my $tenflaratio;
		my $diffratio;
		my $gmapscore;
		
		foreach my $hash_ref (@{$allreaddepth{$name}}) {
			#print "$hash_ref\n";
			foreach my $type (keys %{$hash_ref}) {
				#print " $type => ${$hash_ref}{$type}\n\n";
				$fivepratio = ${$hash_ref}{$type} if ($type eq '5pflank');
				$threepratio = ${$hash_ref}{$type} if ($type eq '3pflank');
				$elementratio = ${$hash_ref}{$type} if ($type eq 'element');
				if ($mappability) {
					$gmapscore = ${$hash_ref}{$type} if ($type eq 'gmscore') ;
				}
			}
		}
		if (($fivepratio == 0) && ($threepratio == 0) && ($elementratio == 0)) {
			$predigenoty = "2";
			$diffratio = "undef";
			if ($mappability) {
				print $fhout "$name\t$fivepratio\t$elementratio\t$threepratio\t$predigenoty\t$diffratio\t$gmapscore\n";
			} else {
				print $fhout "$name\t$fivepratio\t$elementratio\t$threepratio\t$predigenoty\t$diffratio\n";	
			}
		} else {
			$avgflankratio = ($fivepratio + $threepratio)/2;
			$fiftypflaratio = ($avgflankratio * 0.5);#it .40 before
			$tenflaratio = ($avgflankratio * 0.10);
		 
			unless ($avgflankratio == 0) {
				#$diffratio = (($avgflankratio - $elementratio)/$avgflankratio)*100;
				$diffratio = ($elementratio/$avgflankratio)*100;
			} else {
				$diffratio = "undef";
			}
			if ($elementratio > $fiftypflaratio) {
				$predigenoty = "2";
				#$diffratio = ($elementratio - $fiftypflaratio);
			} elsif (($elementratio <= $fiftypflaratio) && ($elementratio >= $tenflaratio)) {
				$predigenoty = "1";
				#$diffratio = ($elementratio - $twentyflaratio);
			} elsif ($elementratio < $tenflaratio) {
				$predigenoty = "0";
				#$diffratio = ($twentyflaratio - $elementratio);
			}
				$fivepratio = sprintf("%.2f",$fivepratio);
				$threepratio = sprintf("%.2f",$threepratio);
				$elementratio = sprintf("%.2f",$elementratio);
				$diffratio = sprintf("%.2f",$diffratio) if ($diffratio ne "undef");
				if ($mappability) {
					print $fhout "$name\t$fivepratio\t$elementratio\t$threepratio\t$predigenoty\t$diffratio\t$gmapscore\n" ;
				} else {
					print $fhout "$name\t$fivepratio\t$elementratio\t$threepratio\t$predigenoty\t$diffratio\n" ;
				}
		
		}
	}
close ($fhout);
}


sub print_array {
	my ($outfile, @array)= @_;
	open (my $ah,">","$outfile") || die ("cannot open file $outfile to write $!\n");
		foreach my $dataline (@array) {
			print $ah "$dataline\n";
		}
	close ($ah);
}
sub store_mapscores {
	my ($gcordinate) = @_;
	#print STDERR "$gcordinate\n";
	my $mappingscore;
	if (exists ($cordinatescore{$gcordinate})) {
		my $value = $cordinatescore{$gcordinate};
		$mappingscore = $value;
		#push (@mappingscores,$value);
	} else {
		my ($chro,$start5,$end5) = split /[:,-]/, $gcordinate;
		#print STDERR "$chro,$start5,$end5\n";
		$mappingscore = &extractgmscores_mysql($chro,$start5,$end5);
		#push (@mappingscores,$mappingscore);
	}
	return $mappingscore;
}
sub extractgmscores_mysql {
	my ($rchr,$rstart,$rend) = @_;
	#my $length = $rend - $rstart;
	my $gposition = $rchr.":".$rstart."-".$rend;
	#print "$gposition\n";
	my @startdata =();
	my @enddata=();
	my $length_raw;
	#my $len_avg;
	my $avg;
	my $weightedavg;
	my $twoscore;
	#retrieve closest numbers from the table using a command 
	my $sth = $dbh->prepare("(SELECT * FROM $mysqltable WHERE chromo='$rchr' and start <= '$rstart' order by start desc limit 1) 
								union 
							(SELECT * from $mysqltable where chromo='$rchr' and end >= '$rend' order by start asc limit 1)") ;
	$sth->execute();
	my $line = 0;
	while (my $row = $sth->fetchrow_hashref()) {
		my $startd = $row->{'start'};
		push (@startdata,$startd);
		my $endd = $row->{'end'};
		push (@enddata,$endd);
		
		print "Found a row: chr = $row->{'chromo'},$row->{'start'},$row->{'end'},score = $row->{'score'}\n" if ($verbose);
		print "-"x50,"\n"if ($verbose);
		$line++;
	}
	$sth->finish();
	my $tstart = $startdata[0];
	my $tend = $enddata[0] if ($line == 1);
	$tend = $enddata[1] if ($line == 2);
	#retrieve whole raw from mysql table using perlDBI  
	$sth = $dbh->prepare("SELECT * FROM $mysqltable WHERE chromo='$rchr' AND start >= '$tstart' AND end <= '$tend'");
	$sth->execute();
	my $count = 0;
	my $sum = 0;
	my $weightedsum = 0;
	
	while (my $row = $sth->fetchrow_hashref()) {
		my $score = $row->{'score'};
		my $rawstart = $row->{'start'};
		my $rawend = $row->{'end'};
		$rawstart = $rstart if ($rawstart == $tstart);
		$rawend = $rend if ($rawend == $tend);
		my $length_raw = ($rawend - $rawstart);
		my $weighscore = $score*$length_raw;
		print "Found a row: chr = $row->{'chromo'},$row->{'start'},$row->{'end'},score = $row->{'score'}\n" if ($verbose);
		print "*"x50,"\n"if ($verbose);
		$count++;
		$sum += $score;
		$weightedsum += $weighscore;
	}
	$sth->finish();
	if ($count >= 1) {
		$avg = $sum/$count ;
		$avg = sprintf("%.3f",$avg);
		$weightedavg = $weightedsum/$count;
		$weightedavg = sprintf("%.3f",$weightedavg);
		#$len_avg = $length."=>".$avg;
		#$len_avg = "$tstart-$tend => $avg";
		$twoscore = $avg.",".$weightedavg;
		$cordinatescore{$gposition} = $avg.",".$weightedavg;
		
	}
	return ($twoscore);	
}
sub email {
	my $to = 'jainyt@genetics.utah.edu';
	my $cc = 'jainyt@genetics.utah.edu';
	my $from = 'jainy@yoda.genetics.utah.edu';
	my $subject = 'findsoloLTR output ';
	my $message = "Results are ready for $path/$file";

	my $msg = MIME::Lite->new(From     => $from,
							  To       => $to,
							  Cc       => $cc,
							  Subject  => $subject,
							  Data     => $message);
	$msg->send;
	print "Email Sent Successfully\n";
}
