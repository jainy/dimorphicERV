#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  Jan 2017
# email  :  jainythomas1@gmail.com
# Purpose :  helps to identify if a proviral sequence exists in place of soloLTR in populations only doable if a bam file avaiable with with paired end reads and if the location is known
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Bio::SearchIO; 
use Bio::SeqIO;
use Bio::DB::Fasta;	
use MIME::Lite;
use Data::Dumper;
use File::Copy;
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(firstidx);
use DBI;

my $version = "18.0_part1";
my $scriptname = "findprovirus_1.pl";
my $changelog = "
#   - v1.0 = 16 Jan 2017 
#	- v2.0 = 24 Jan 2017
#				modified the script so that so that path of the output can be given at the commandline
#				using the picardtools version 2.8.1 and java 1.8
#	- v2.1 = 1 	Feb 2017
#				included concatenate function for contigs and singlets
#	- v2.5 = 1 	Feb 2017
#				the output will be directed to a folder with region within each ltr
#				included the verbose option
#   - v2.6 = 13 Feb 2017
#				can be selected in the script only IGV run or extracted reads needs to be done or both	
#   - v3.0 = 8 	May 2017
#				perform blast with envelope and gag sequences supplied in the command line
#	- v4.0 = 10	May 2017
#				included subroutines
#	- v4.1 = 16 May 2017
#				will not write the reads file in picard tools script
#				changed name to findprovirus.v4.1_forHERVK.pl from Run_samtools_picard_cap3_selection_blast.v4.0_coverage.pl
#	- v5.3 = 14 June 2017
#				deprecated the region option
#               trying to extract only discordant reads using picard tools so that time for extraction can be reduced#
#				found that it does not reduce the time at least 3 hours go throught the file evenif the reads that needs to be collected are less
#				rewrote the script so extraction of reads done only one time for an individual
#   - v5.4 = 20 June 2017 
#				made the array unique for when collecting reads from all bam
#	- v5.5 = 26 June 2017
#				fixed one bug that will the job for zero sized files
#   - v5.6 = 11 July 2017
#				fixed the bug that the script will die if there are no discordant reads identified for an individual
#   - v5.7 = 13 July 2017
#				introduced the quality features while extracting reads, will extract reads only that are above the mapping quality of 20 or above			
#	- v6.0 = 25 August 2017
#				To find number of mates of discordant reads that have hit to env or gag sequences
#				bug: found discordant reads were extracting without taking flanking. So need rerun the with the flanking
#	- v7.0 = 29 August 2017	
#				introduced a prediction report based on the length of the alignment to the genomic soloLTR locus and number of discordant reads that have hit to TE sequence
#
#	- v8.0 = 6 September 2017
#				read depth analysis is done as part of the pipeline and provide read depth info in the prediction report
#	- v9.0 = 26 September 2017
#				to do assembly with out the discordant reads and blast it genomic extract
#	- v10.0 = 10 October 2017
#				Trying SPAde instead of cap3 assembly
#				giving whole internal sequence for blast will help in catching stuff if there is internal truncation in the copy
# 				deprecating TE blast to assembled sequences with discordant pairs- only blast to mates of discordant reads and mapped reads to genomic sequence kept
#  - v12.4 = 26 October 2017
#				successfully implemented SPADE for assembling the mapped reads. (options for restart and continue available)
#				successfully implemented SPAde for assembling all reads including discordant reads. Split discordant reads into 1 and 2 and find teh corresponding mate from the unpaired end file and join them to the c
#				corresponding read sequence
#  - v12.5 = 26 October 2017
#				SPAdes misserably failed for some locations and did not produce any assembled sequence. In some locations, where is not enough coverage sits in recent TE it also failed not able to identify the insert distance.
#				no scaffold.fasta produced in such case.So going to bring back CAP3 where SPAdes fails.
#				
#  - v13.0 = 31 October 2017
#				SPAdes is very slow. So going to change the assembler to Minia.	
#	-v14.0 = 8 	November 2017
#				Minia cannot run when there are multiple minia running. So going to split the pipeline to two parts. One part will run CAP3 and will give the output. Those regions cannot be solved by Cap3 will be ran by SPAde for better genotyping 
#	-v14.5 = 12 February 2018
#				incorporated genome mappability score, for that extracting the reference position of the discordant reads using bamtobed
#				for read depth and trying to normalise the read depth (by adding the read depth of all LTRs of an individual) and calculate read depth dividing by the normalised value		
#				introduced the mysql access feature for getting the genome mappability scores
#	-v15.0 = 22 February 2018
#				not time wise not good when calculate the scores for each multiple read, instead going to find boundaries of overlapping region that reads mapps to and get scores for such regions
#	-v15.1 = 26 February 2018
#				 most of the time extracting same region from the mysql, if possible modify the code so that sql extraction is done once for a loci and the scores stored in table and distributed accordingly from the table
#	-v15.2 = 26 February 2018			
#				created a new table with index for mysql table for faster extraction; fixed illegal division by zero error so modified the script sub normalise read depth
#	-v15.4 = 1 	March 2018
#				 was able to sort of aoH because I changed the way it was loaded into the hash, fixed the bug where previous values were printed for the sql scores, no bug when extracting overlapped boundaries
#
#	-v16.0 = May 15 2018
#				changed from avg flank depth to avg element depth
#	-v17.0 = June 3 2018
#				introduced calculation of average of all scores, now the prediction depends only on the no of discordant reads, and average read, ability to form soloLTRs are only reported
#				To do:
#				instead of calculating try to get this from the output by formating the output to give query coverage
#	-v18.0 = June 13 2018
#				changed mappability to optional step, add the provision for providing user, password in the commandline

\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <BAM ID table> -f <file with ltr cordinates> -bl <location of bamfiles> [-p <path of the outputdirectory>][-g <path of the genome>][-m mapscores][-te <TEseq>][-u Username] [-pd password][-db mysql database][-mt mysql table][-b] [-i] [-e] [-rd] [-x] [-v] [-c] [-h] [-s]
	
    MANDATORY ARGUMENT:	
    -t,--table 	      (STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file         (STRING) file containing cordinates of solo LTR (sampleID, chr, start, end (of solo LTR),unique identifier,length,strand ) 
    -bl,--bamlocation (STRING) location of bam files
    -b,--both   		(BOOL)   run extraction and assembly with sliced bam file for viewing IGV
    -st,--seqtkpro	  (STRING) path to the seqtk
    -pc,--picard	  (STRING) path to the picardtools
    -cp,--cap3		  (STRING) path to the cap3 assembler
    -bp,--blast		  (STRING) path to the blast 
    -bd,--bedtools	  (STRING) path to the bedtools 
    
    OPTIONAL ARGUMENTS:  
    -p,--path   		(STRING) output directory name (path)
                         	 	 Default = <current working directory>
    -te,--teseq 		(STRING) consensus internal sequence of the provirus
    -g,	--genome		(STRING) path of the genome
    -rd,--readepth  	(BOOL)   read depth analysis
    -x, --extract   	(BOOL)	 extract genomic sequence               
    -i,--igv    		(BOOL)   only IGV has to be run
    -e,--reads  		(BOOL)   run picard tools and cap3 for extracting and assembling reads 
    -m --mapscores 		(BOOL) 	 Need to calculate mappability scores	
    -db --mysqldbinfo 	(STRING) database=jainys_db, (needs to prepare before running script)   
    -u --user 			(STRING) Username for mysql database,(needs to prepare before running script) 
    -pd,--password 		(STRING) password for mysql database,(needs to prepare before running script) 
    -mt,--mysqltable    (STRING) mysql table containing mappability scores (needs to prepare before running script)  
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help   		(BOOL)   Print this usage\n\n";
   
   
#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$table,$path,$teseq,$extract,$out,$igv,$reads,$readdepth,$both,$GENOME,$mapscores,$mysqltable,$mysqldb,$user,$password,$verbose,$bamlocation,$seqtkpro,$picardpro,$BLASTpro,$CAP3pro,$bedtoolspro,$help,$v,$chlog);
GetOptions ('t=s' => \$table,
            'f=s' => \$file,            
            'p=s' => \$path,
            'te=s'=> \$teseq,
            'x'   => \$extract,
            'g=s' => \$GENOME,
            'o=s' => \$out,
            'rd'  => \$readdepth,
            'i'   => \$igv,
            'e'   => \$reads,
            'b'   => \$both,
            'm'   => \$mapscores,	
            'mt=s'=> \$mysqltable,
            'db=s'=> \$mysqldb,
            'u=s' => \$user,
            'pd=s'=> \$password,
            'bl=s'=> \$bamlocation,
            'st=s'=> \$seqtkpro,
            'pc=s'=> \$picardpro,
            'cp=s'=> \$CAP3pro,
            'bp=s'=> \$BLASTpro,
            'bd=s'=> \$bedtoolspro,
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table) && (! $file) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $file) || !($table) || ($help));
my $cwd = getcwd();
$path = $cwd if (!$path);
$out = "$path/$file.prediction_alleles.txt" if (! $out);

#check if the following tools are installed and the right path is given
#my $seqtkpro = "/home/jainy/software/seqtk";
#my $picardpro = "/home/jainy/software/picard-2.9.2";
#my $CAP3pro = "/home/jainy/software/CAP3";
#my $BLASTpro = "/home/jainy/software/ncbi-blast-2.6.0+/bin";
#my $bamlocation = "/kbod2/WGS_DATA/SGDP_bams_public";
#my $bedtoolspro = "/home/jainy/software/bedtools2/bin";
#only if mappability scores are calculated, my sql table needs to be constructed and indexed prior to running script
#my $mysqltable = "wgEncodeCrgMapabilityAlign100merhg38_lo_index";


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my %bamfile = ();
my $bamid ;
my $dispath;
my $uniqueid;
my $individual;
my $genomeloc;
my $ltr;
my @alldisreads =();
my %hashindividual =();
my %discordantmatelist =();
my $discordmatefile;
my $disrenamedfile;
my @nodisreads;
my @nodisreadswithERV;
my $nudiscortefile;
my $dismatereadblastout;
my $extractgenomicseqout;
my @dbIDs;
my $db;
my $gchr;
my $gstart;
my $gend;
my $gstrand;
my $gexstart;
my $gexend;
my $renamedfile;
my $glength;
my %genomicseqblast;
my @predictiondetails;
my @noreadsregion;
my $gextractloc;
my $discornu;
my $elegenomeloc;
my $fivepgenomeloc;
my $threepgenomeloc;
my @avgflankdepthall;
my @logvalues;
my $readdepthpath;
my %allreaddepth;
my $nodisreadidlist = "$file.ids_with_no_discordantreads.txt";
my $noreadsidlist = "$file.ids_with_no_reads.txt";#most likely Y regions in females
my $nodisreadswithERVs = "$file.ids_with_no_discordantreadsERVs.txt";
my $renamedmappedreadfile;
my %bedfile;
my %avgrdindvi;
my %rd_uid;
my @person;

%bamfile = &load_hash ($table);
	
open (my $fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\t/,$line);
		my $count = @col;
		print "the number of columns in each line of the $file is $count, \n" if ($verbose);
		$individual = $col[0];
		push (@person,$individual) if (! @person);
		$ltr = $col[4];
		my ($chr,$start,$end,$strand) = ($col[1],$col[2],$col[3],$col[6]);
		my $elestart = $start; #need to change depending on cordinates initially used
		my $eleend = $end;#need to change depending on on cordinates initially used
		my $fivepstart = $elestart - 250;#250 before
		my $fivepend = $elestart - 3;#3 before
		my $threepend = $eleend + 250;
		my $threepstart = $eleend + 3;
		$elegenomeloc = $chr.":".$elestart."-".$eleend;
		$fivepgenomeloc = $chr.":".$fivepstart."-".$fivepend;
		$threepgenomeloc = $chr.":".$threepstart."-".$threepend;
		#since these coordinates represent the element coordinates, need to get 100 bp flanking for extracting all the discordant reads
		$start = $start - 100;
		$end = $end + 100;
		$genomeloc = $chr.":".$start."-".$end;
		&find_bamid();
		$uniqueid = $individual.".".$ltr.".".$genomeloc;
		$genomicseqblast{$uniqueid} = [@predictiondetails];
		if ($readdepth) {
				if ($person[0] ne $individual) {
				@avgflankdepthall = ();
				$person[0] = $individual;
			}
			make_path  ("$path/ReadDepth/$ltr");
			$readdepthpath = "$path/ReadDepth/$ltr";
			my $eledepth = &avg_readdepth($elegenomeloc);
			my $fivepdepth = &avg_readdepth($fivepgenomeloc);#  may be no need to calculate these flanking read depth since I am not using it
			my $threepdepth = &avg_readdepth($threepgenomeloc);
			my $avgflankdepth = ($fivepdepth + $threepdepth)/2;
			#going to take the average of read depth of only solo LTRs
			push (@avgflankdepthall,{$genomeloc => $eledepth});#storing avg read depth at  a locus for an individual
			$allreaddepth{$individual} = [@avgflankdepthall];#storing avg read depth for all locus for an individual
			$rd_uid{$uniqueid} = $eledepth;
		}
		make_path  ("$path/IGV/$individual")  if (($igv) || ($both));
		make_path  ("$path/ExtractedReads/$individual") if (($reads) || ($both));
		if (($igv) || ($both)) {	
			system("samtools view -b -o $path/IGV/$individual/$uniqueid.bam $bamlocation/$bamid.38.sorted.bam $genomeloc") == 0 or die ("unable to run command on $uniqueid \n");
			print STDERR " 	Extracting 	$uniqueid using samtools done\n" if ($verbose);
			system ("samtools index -b $path/IGV/$individual/$uniqueid.bam ") == 0 or die ("unable to create index file of $uniqueid.bam \n");
			print STDERR " indexing the bamfile done \n" if ($verbose) ;
		}	
			##Extracting reads for the assembling the reads
		if (($reads) || ($both)) {	
			##Extracting discordant reads
			make_path  ("$path/Discordantreads/$individual");
			#unless (-e "$path/Discordantreads/$individual/$uniqueid.dismapped.reads.fasta") {
				system ("samtools view -b -q 30 -F 3854 $bamlocation/$bamid.38.sorted.bam $genomeloc > $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.bam") == 0 or die ("unable to extract readsby F 3854 flag $!");
				system ("samtools bam2fq $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.bam | $seqtkpro/seqtk seq -A -q30 > $path/Discordantreads/$individual/$uniqueid.dismapped.reads.fasta") == 0 or die ("unable to convert discordant bam file  to fasta $uniqueid \n");
				
			#}
			# identifying the reference coordinates of discordant reads
			system ("$bedtoolspro/bedtools bamtobed -i $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.bam > $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.reads.bed") == 0 or die ("unable to convert bamtobed $uniqueid \n");
			#make a subroutine for loading the info into a hash-uniqueid-individual-readID==>reference position
			&read_location("$path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.reads.bed");

			#identifying mates of discordant reads
			$dispath = "$path/Discordantreads/$individual";
			my $mappedisreads = "$uniqueid.dismapped.reads.fasta";
			%discordantmatelist = &load_readIDs ($dispath,$mappedisreads);
			my $totaldisreads = keys %discordantmatelist;
			push (@{$genomicseqblast{$uniqueid}},{'totaldisreads' => $totaldisreads});
			#print discordant read mates to a file
			make_path  ("$path/Discordantreads/$file/dismateIDLists");	
			$discordmatefile = "$path/Discordantreads/$file/dismateIDLists/$uniqueid.discordantmatesreadIDlist.txt";
			&print_hash(%discordantmatelist);
			#loading the all the discordant reads that needs to be extracted to a hash
			&collectreads_indi();
		}
	}
#print Dumper %discordantmatelist, "\n";	
#print Dumper %hashindividual, "\n";
#print Dumper %bedfile, "\n";					
close $fh;
make_path  ("$path/Allreadsforbam/$file");
my $indiallreads = "$path/Allreadsforbam/$file"; 
&printreads_indi();
&extractreads_bam();
make_path ("$path/Discordantreads/$file/discordantmatesonly");
my $dismateseqfilepath = "$path/Discordantreads/$file/discordantmatesonly"; 
my $dismateIDlistpath = "$path/Discordantreads/$file/dismateIDLists";
my %readswithte;
&extractdiscordmates();
#print Dumper %allreaddepth;
&normalise_readepth();
#print Dumper %avgrdindvi;
#print Dumper %rd_uid;
&calculate_readepth(%avgrdindvi);
if ($extract) {   
	# index the genome and connect to the fasta file
	my $reindex;
	my $indexfile = "$GENOME.index";
		if (-e $indexfile) {
			$reindex = 0;
			print "\t Genome previously indexed - Skipping indexing...\n\t (if you want to reindex the genome, delete $indexfile)\n";
		} else {
			$reindex = 1;
			print "\t  Genome not indexed - Indexing $GENOME...\n";
		}
	$db = Bio::DB::Fasta->new( $GENOME, -reindex=>$reindex) or die "\t    ERROR - Failed to create Bio::DB::Fasta object from $GENOME $!\n";
	#create list of the ID of the genome file  
	@dbIDs = $db->get_all_ids();
}
#connecting to mysql database at yoda
my $dsn;
my $dbh;
if ($mapscores) {

	$dsn = "DBI:mysql:database=$mysqldb;host=localhost;port=22";
#my $user = "jainy";
#my $password = "wysql123";
	$dbh = DBI->connect($dsn, $user, $password,
						{'RaiseError' => 1});
}
my %overlapscore;
						
open ($fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\t/,$line);
		my $count = @col;
		print "now analysing $line \n";
		$individual = $col[0];
		$ltr = $col[4];
		($gchr,$gstart,$gend,$gstrand) =($col[1],$col[2],$col[3],$col[6]);
		#here cordinates represents the coordinates of the element so need to take 100 bp flanking to catch all the discordant reads mapping to that position
		$gstart = $gstart - 100;
		$gend = $gend + 100;
		$genomeloc = $gchr.":".$gstart."-".$gend;
		$gexstart = $gstart + 50;
		$gexend = $gend - 50;
		$glength = $gexend - $gexstart;
		$gextractloc = $gchr.":".$gexstart."-".$gexend;# element cordinates with 100 bp flanking
		&find_bamid();
		$uniqueid = $individual.".".$ltr.".".$genomeloc;
		make_path  ("$path/ExtractedReads/$ltr") if (($reads) || ($both));
		make_path  ("$path/MappedReads/$ltr") if (($reads) || ($both));
		my $epath = "$path/ExtractedReads/$ltr";
		my $mpath = "$path/MappedReads/$ltr";
		my $dmatefile; 
		if (($reads) || ($both)) {	
			#extracting only the mapped reads
			system("samtools view -b -q 30 -F 4 $bamlocation/$bamid.38.sorted.bam $genomeloc > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.bam ") == 0 or die ("unable to extract bam at $uniqueid \n");
			system ("samtools bam2fq $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.bam > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fastq") == 0 or die ("unable to convert to fasta $uniqueid \n");
			copy("$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fastq", "$path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.fastq") or die "Copy failed fastq $uniqueid:$!"; 
			system ("$seqtkpro/seqtk seq -A -q30 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fastq > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fasta") == 0 or die ("unable to run seqtk on fastq $uniqueid \n");
			if (-e "$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
				unless (-z "$path/Discordantreads/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
					#copy the corresponding the discordant mate reads, 
					copy("$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta", "$path/ExtractedReads/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta") or die "Copy failed:$!"; 
					#concatenate mapped reads and discordant mates for CAP3 assembly
					system ("cat $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fasta $path/ExtractedReads/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta > $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta" ) == 0 or die ("unable to concatenate dismatereads on $uniqueid\n");
					my $fisizeAR = -s "$path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta";#concatenated file
					if ($fisizeAR > 0) {		   
						system ("$CAP3pro/cap3 $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta > $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.asemb.fasta ") == 0 or 
						die ("unable to assemble all reads for $uniqueid \n");
						print STDERR " Assembling fasta done \n" if ($verbose);
						#Rename file 
						my $assembledfile = "$uniqueid.allreads.concatenated.fasta.cap.contigs";
						&renameseq_filename ($assembledfile,$epath);
					} 
				}
			} else {
				warn "No discordant reads were identified $uniqueid\n";
				system ("$CAP3pro/cap3 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fasta > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fasta.assb.fasta") == 0 or 
 				die ("unable to assemble mapped reads for $uniqueid \n");
				my $assembledfile = "$uniqueid.onlymappedreadIDs.fasta.cap.contigs";
				&renameseq_filename ($assembledfile,$epath);
				push (@nodisreads,$uniqueid);
			}	
		}
		if ($teseq) {
			#preparing for blasting the discordant mates to env/gag
			make_path  ("$path/QuantifyDisc/$ltr");
			if (-e "$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
				copy("$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta", "$path/QuantifyDisc/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta") or die "Copy failed:$!"; 			
				my $dismatereads = "$path/QuantifyDisc/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta";
				$dismatereadblastout = "$path/QuantifyDisc/$ltr/$uniqueid.discordantmatesreadIDlist.blast.out";
				my $filesize = -s $dismatereads;
				if ($filesize > 0) {
					print STDERR "performing BLAST on individual discordant reads..\n" if ($verbose);
					unless (-e "$dismatereads.nhr") {
						system ("$BLASTpro/makeblastdb -in $dismatereads -dbtype nucl") == 0 or die ("unable to makeblastdb on $dismatereads \n");
					}
					#system ("$BLASTpro/blastn -db $dismatereads -query $teseq -evalue 0.0001 -out $dismatereadblastout") == 0 or die ("unable to to perform blast $uniqueid \n");
					system ("$BLASTpro/blastn -db $dismatereads -query $teseq -evalue 0.0001 -outfmt 6 -out $dismatereadblastout.tabular.out") == 0 or die ("unable to to perform tabular blast $uniqueid \n");#fortableoutput
					%readswithte = &parse_blast("$dismatereadblastout.tabular.out");
					#print Dumper %readswithte;
					$discornu = keys (%{$readswithte{$uniqueid}});
					#print "number of discordant reads with TE $discornu\n";
					print "list of discordant reads position of reads that have mate in TE\n";
					my ($arrarefposi) = &compare_twoHOH(\%readswithte,\%bedfile);
					#print Dumper %$posirefs;
					print Dumper @$arrarefposi if ($verbose);
					#my $overlaps = &extract_overlap_bound(%$posirefs);
					my $overlaps = &extract_overlap_bound_aoH(@$arrarefposi);
					print "extracted overlapping regions of the discordant reads\n"if ($verbose);
					print Dumper @$overlaps if ($verbose);
					print "extracting scores from genome mappability scores\n" if (($verbose) && ($mapscores));
					my $mapscoref = &extractgmscores($overlaps,$gchr) if ($mapscores);
					print "extracted scores from genome mappability scores\n" if (($verbose) && ($mapscores));
					print Dumper  @$mapscoref, "\n" if (($verbose) && ($mapscores));
					#need to output the name of the reads, identify the position of those reads, check the mappability score and output them 
					push (@{$genomicseqblast{$uniqueid}},{'nu_of_discreadTE' => $discornu});
					push (@{$genomicseqblast{$uniqueid}},{'mapscore_discread' => [@$mapscoref]}) if ($mapscores);
				} else {
					$discornu = 0;
					push (@{$genomicseqblast{$uniqueid}},{'nu_of_discreadTE' => $discornu});
					push (@nodisreadswithERV,$uniqueid);
				}
			} 
		}
		if ($extract) {
			&extract_genomicseq();
			my $paligned;
			my $pmappaligned;
			print STDERR "performing BLAST on the allpairedend reads assembly with genomic seq......\n" if ($verbose);
			make_path  ("$path/Genomeblast/$ltr");
			$renamedfile = "$path/ExtractedReads/$ltr/Renamedcontigs/$uniqueid.rename.fasta";			
			my $gblastout = "$path/Genomeblast/$ltr/$uniqueid.gblast.out";
			if (-e $renamedfile) {
				my $sizeoffil = (-s $renamedfile);
				if ($sizeoffil > 0) {
					unless (-e "$renamedfile.nhr") {
						system ("$BLASTpro/makeblastdb -in $renamedfile -dbtype nucl") == 0 or die ("unable to makeblastdb on $renamedfile \n");
					}
					#system ("$BLASTpro/blastn -db $renamedfile -query $extractgenomicseqout -evalue 0.0001 -out $gblastout") == 0 or die ("unable to to perform gblast $uniqueid \n");
					system ("$BLASTpro/blastn -db $renamedfile -query $extractgenomicseqout -evalue 0.0001 -outfmt 6 -out $gblastout.tabular.out") == 0 or die ("unable to to perform gtabular blast $uniqueid \n");#fortableoutput
					$paligned = &parse_blast_foralen("$gblastout.tabular.out");
					push (@{$genomicseqblast{$uniqueid}},{'perallreadsaligned' => $paligned});
				} else {
					$paligned = 0;
					push (@{$genomicseqblast{$uniqueid}},{'perallreadsaligned' => $paligned});
				}
			}		
		}	
	}
#print Dumper %genomicseqblast;
&print_predictionAoH(%genomicseqblast) if ($extract);
print_array($nodisreadidlist,@nodisreads);
print_array($noreadsidlist,@noreadsregion);
print_array($nodisreadswithERVs,@nodisreadswithERV);
&email();
close $fh;
# Disconnect from the database.
$dbh->disconnect();
exit;

#-------------------------------------------------------------------------------
#----------------------------------- SUB ---------------------------------------
#-------------------------------------------------------------------------------

sub load_hash {
	my ($file) = @_;
	my %sbamfile;
	open (my $th, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @namebam = split(/\t/,$data); # splitting the data based on tab and storing into the arrray
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
		die "there is something wrong with the bamids check inputfiles $!\n";
	}
}
sub renameseq_filename {
	my ($contigfile,$fpath) = @_;
	make_path  ("$fpath/Renamedcontigs");
	#copy ("$contigfile","$path/$contigfile.mod.fasta") or die "Copy failed to $path: $!";
	open (my $bhout, ">","$fpath/Renamedcontigs/$uniqueid.rename.fasta") or die "\n ERROR (main): could not open to read $fpath/Renamedcontigs/$uniqueid.rename.fasta $!\n";
	open (my $bh, "<", "$fpath/$contigfile") or confess "\n ERROR (main): could not open to read $contigfile $!\n";
		while(my $dataline = <$bh>) {
			chomp($dataline);
			#print STDERR "$line\n";
				if ($dataline =~ m/^\>\w+\d+/) {
					$dataline =~ s/^\>(\w+\d+)/\>$uniqueid\.$1/;
					print $bhout "$dataline\n";
				}
				else {
					print $bhout "$dataline\n";
				}
		}
	close $bh;
	close $bhout;
}
sub load_readIDs {
	my ($dpath,$mappedid) = @_;
	open (my $idh,"<","$dpath/$mappedid") || die ("unable to open $mappedid $! \n ");
	my %idlist;
	while (<$idh>) {
		 my $head = $_;
		 chomp $head;
		if ($head =~ /^\>(.*)\/(\d)/) {
			$head = "$1\/2" if ($2 == 1);
			$head = "$1\/1" if ($2 == 2);
			$idlist{$head} =1;
		} elsif ($head =~ /^\>(.*)/) {
			$head = "$1\/0";
		 	$idlist{$head} =1;
		} else {
			next;
		}
	}
	return (%idlist);
	close $idh;
}
sub print_hash {
	my %hashtoprint = @_;
	open (my $hp,">","$discordmatefile") || die ("failed to open file to write discordant mate list $!\n");
	foreach my $mate (sort keys %hashtoprint) {
		print $hp "$mate\n";
	}
	close $hp;
}
sub collectreads_indi {#collect all reads for an individual
	@alldisreads =();
	foreach my $read (sort keys %discordantmatelist) {
		$read = substr $read, 0,-2;#remove /1or2 from the file
		push (@{$hashindividual{$individual}},$read);
		(@{$hashindividual{$individual}}) = uniq (@{$hashindividual{$individual}});#for making the array unique
	}
	return (%hashindividual);
}
sub printreads_indi {
	foreach my $indi (sort keys %hashindividual) {
		open (my $ih, ">","$indiallreads/$indi.allreadIDs.txt" ) or die ("cannot write file $indi.allreadIDs.txt $!\n");
		foreach my $reads (@{$hashindividual{$indi}}) {
			print $ih "$reads\n";
		}
		close $ih;
	}
}
sub extractreads_bam {
	my @indifiles = `ls $path/Allreadsforbam/$file`;
	my $bam_id;
	foreach my $indivi (@indifiles) {
		my @allreadsfile_name = split (/\./,$indivi);
		my $allreadindi = $allreadsfile_name[0];
		if (exists ($bamfile{$allreadindi})) {
			$bam_id = $bamfile{$allreadindi};
		} else {
			die "there is something wrong with the bamids check inputfiles (sub extractreads_bam) $!\n";
		}
		#extract reads using picard tools
		unless (-e "$path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam") {
			system ("java -jar $picardpro/picard.jar FilterSamReads INPUT=$bamlocation/$bam_id.38.sorted.bam FILTER=includeReadList READ_LIST_FILE=$path/Allreadsforbam/$file/$allreadindi.allreadIDs.txt WRITE_READS_FILES=false OUTPUT=$path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam") == 0 or die ("unable to run picard tools in $file on $allreadindi \n");
			system ("samtools bam2fq $path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam | $seqtkpro/seqtk seq -A -q30 > $path/Allreadsforbam/$file/$allreadindi.allreadIDs.fasta") == 0 or die ("unable to convert to fasta $allreadindi \n");

		}
	}
}
sub extractdiscordmates {
	my @dismates = `ls $dismateIDlistpath`;
	foreach my $matefile (@dismates) {
		chomp $matefile;
		my @matefilename = split (/\./,$matefile);
		my $indivifilename = $matefilename[0];
		my %mates = ();
		open (my $mh, "<", "$dismateIDlistpath/$matefile") or confess "\n ERROR (main): could not open to read $matefile $!\n";
		while (my $dataline = <$mh>) { 
			chomp $dataline; 
			$mates{$dataline} = 1; 
		}
		close ($mh);
		#print Dumper %mates, "\n";
		if (-e "$path/Allreadsforbam/$file/$indivifilename.allreadIDs.fasta") {
			my $readio_obj = Bio::SeqIO->new(-file 	 => "$path/Allreadsforbam/$file/$indivifilename.allreadIDs.fasta", 
											 -format => 'fasta') 
										 or die "\t    ERROR - Failed to create SeqIO FH object from $indivifilename.allreadIDs.fasta $!\n";  
			my $outreadio_obj = Bio::SeqIO->new(-file   => ">$dismateseqfilepath/$matefile.fasta",
												-format => 'fasta') 
											 or die "\t    ERROR - Failed to create SeqIO FH object from $dismateseqfilepath/$matefile.fasta $!\n";  
			while (my $seq = $readio_obj->next_seq() ){
				my $header = $seq->display_id;
				if (exists $mates{$header}) {
					$outreadio_obj->write_seq($seq);
				} else {
					next;
				}
			}
		}
	}
}
sub check_file {
	my ($filename,$summaryfile) = @_;
	open (my $oh, ">>", $summaryfile ) || die ("cannot open file $summaryfile to append\n");
	my $size = -s $filename;
	my @names = split (/\//,$filename);
	my $finalname = $names[-1];
	print $oh "$finalname\t$size\n";
	close $oh;
}
sub parse_blast {
	#parse blast tabular output unique hit with highest blast score
	my $blast = shift;
	my %d = ();
	open (my $bl, "<", $blast) or confess "\n ERROR (sub parse_blast): could not open to read $blast $!\n";
	LINE: while (my $line = <$bl>) {
		chomp $line;
		next LINE if ($line !~ /\w/);
		my @col = split(/\t/,$line); # query_id	subject_id	pct_identity	aln_length	n_of_mismatches	gap_openings	q_start	q_end	s_start	s_end	e_value	bit_score
		my $read = $col[1];
		my $score = $col[11];
		$d{$uniqueid}{$read} = $score;
		if (exists ($d{$uniqueid}{$read})) {
		 	my $currentscore = $d{$uniqueid}{$read};
		 	if ($currentscore > $score){
		 		next;
		 	} else {
		 		$d{$uniqueid}{$read} = $score;
		 	}
		}
	}
	#my $noofdisread = keys (%d);
	#return ($noofdisread);
	return (%d);
	close ($bl);
}
sub parse_blast_foralen {
	#parse blast tabular output unique hit with highest alignment length
	my $blast2 = shift;
	my %s = ();
	open (my $sl, "<", $blast2) or confess "\n ERROR (sub parse_blast_foralen): could not open to read $blast2 $!\n";
	while (my $line = <$sl>) {
		chomp $line;
		next LINE if ($line !~ /\w/);
		my @col = split(/\s+/,$line); # query_id	subject_id	pct_identity	aln_length	n_of_mismatches	gap_openings	q_start	q_end	s_start	s_end	e_value	bit_score
		my $read = $col[1];
		
		my $pidentity = $col[2];
		my $score = $col[11];
		
		$s{$line} = $score;
		if (exists ($s{$line})) {
		 	my $currentscore = $s{$line};
		 	if ($currentscore > $score){
		 		next;
		 	} else {
		 		$s{$line} = $score;
		 	}
		}
	}
	my ($line1,$bigscore) = &sorthashbyvalue(%s);
	my @colu = split(/\t/,$line1);
	my $alength = $colu[3];
	my $percaligned = (($alength/$glength)*100);#instead of calculating try to get this from the output by formating the output to give query coverage
	$percaligned = sprintf ("%.2f",$percaligned);
	#print STDERR "the max perc aligned is $percaligned\n";
	return ($percaligned);
	close $sl; 
}
sub sorthashbyvalue {
	my %contighash = @_;
	my $contig_name;
	my $readnumber;
	foreach my $name (sort { $contighash{$b} <=>  $contighash{$a}} keys %contighash) {
		$contig_name = $name;
		$readnumber = $contighash{$name};
		last;
	}
	return ($contig_name,$readnumber);
}
sub print_keyvalue {
	my ($filepath,%hashtoprint) = @_;
	open (my $kp,">","$filepath") || die ("failed to open file to write $filepath $!\n");
	foreach my $mate (sort keys %hashtoprint) {
		print $kp "$mate\t$hashtoprint{$mate}\n";
	}
	close $kp;
}
sub extract_genomicseq {
	make_path  ("$path/ExtractGenomicsequences");
	$extractgenomicseqout = "$path/ExtractGenomicsequences/$gextractloc.extract.seq.fa";
	unless (-e $extractgenomicseqout) {
		my $gout = Bio::SeqIO->newFh(-format => 'Fasta',
									-file   => ">$extractgenomicseqout") 
									or 	die "\t    ERROR - Failed to create SeqIO FH object from $extractgenomicseqout $!\n";
		my $log = "$extractgenomicseqout.log";
		open (LOG, ">>$log") or die "\t    ERROR - can not create log file $log $!\n";
		my $newId = join('_',$gchr,$gexstart,$gexend);	
			if ("@dbIDs" =~ m/(\S*)($gchr)(\S*)/) {
				my $subSeq = $db->seq($gchr,$gexstart,$gexend);	# extract target sequence
				my $seqobj = Bio::Seq->new( -display_id => $newId,
												   -seq => $subSeq); # create object with target
				print $gout $seqobj;	# print it out (in fasta format)
			} else {
				print LOG "$gchr was not found in $GENOME\n";
			}
	}
	close LOG;
}
sub print_predictionAoH {
	my (%genomelocation) = @_;
	open (my $fhout,">",$out) || die ("cannot open file $out to write $!\n");
	print $fhout "#element\t#total_discordant\t#discordantmates_maps_ERV\t#percent_of_assembly_alignedto_refseq\t#read_depthratio\t#genotype_prediction\t#variable_allele\t#mapability_scores(cordinates=>scores)\n" ;
	foreach my $element (sort keys %genomelocation) {
 		my $nudiscreadTE;
		my $perallreadsalign;
		my $perreadfladepth;
		my $totdismates;
		my @scoreslist;
		my $AVGavgscore;
		foreach my $hash_ref (@{$genomelocation{$element}}) {
			foreach my $type (keys %{$hash_ref}) {
				#print STDERR "$element $type => ${$hash_ref}{$type}\n";
				$nudiscreadTE = ${$hash_ref}{$type} if ($type eq 'nu_of_discreadTE');
				$perallreadsalign = ${$hash_ref}{$type} if ($type eq 'perallreadsaligned');
				$perreadfladepth = ${$hash_ref}{$type} if ($type eq 'read_depthratio');
				$totdismates = ${$hash_ref}{$type} if ($type eq 'totaldisreads');
				@scoreslist = @{${$hash_ref}{$type}} if ($type eq 'mapscore_discread');#"$tstart-$tend => $avg";
				my $SUMscores;
				my $numelem = @scoreslist;
				foreach my $SCORE (@scoreslist) {
					my ($cordType, $AVGscore) = split/=>/,$SCORE;
					$SUMscores += $AVGscore;
				}
				$AVGavgscore = $SUMscores/$numelem unless ($numelem == 0);
				
				#$permappedreadsalign = ${$hash_ref}{$type} if ($type eq 'mappedreadsaligned');
			}
		}
		if ((($nudiscreadTE >= 1) && ($nudiscreadTE < 4)) && ($perallreadsalign >= 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSP\t(P?)\t$AVGavgscore\t@scoreslist\n" ;
		}  elsif ((($nudiscreadTE >= 1) && ($nudiscreadTE < 4)) && ($perallreadsalign <= 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSP\t(P?)\t$AVGavgscore\t@scoreslist\n" ;
		} elsif ((($nudiscreadTE >= 4) && ($nudiscreadTE < 20)) &&  ($perallreadsalign < 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSP\t(S?)\t$AVGavgscore\t@scoreslist\n" ;
		} elsif ((($nudiscreadTE >= 4) && ($nudiscreadTE < 20)) &&  ($perallreadsalign >= 95)) { 
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSP\tSP\t$AVGavgscore\t@scoreslist\n" ;
		} elsif (($nudiscreadTE >= 20) && ($perallreadsalign < 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tPP\tPP\t$AVGavgscore\t@scoreslist\n" ;
		} elsif (($nudiscreadTE >= 20) && ($perallreadsalign > 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSP\t(SP)\t$AVGavgscore\t@scoreslist\n" ;
		} elsif (($nudiscreadTE == 0) && ($perallreadsalign >= 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSS\tSS\tNoavgscore\tno scores\n" ;
		} elsif (($nudiscreadTE == 0) &&  ($perallreadsalign == 0)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tNA\\tno readstNoavgscore\tno scores\n" ;
		} elsif (($nudiscreadTE == 0) && ($perallreadsalign < 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSS\tnested TE?\tNoavgscore\tno scores\n" ;#polymorphic TE inside
		} elsif (($nudiscreadTE == 0) && ($perallreadsalign >= 95)) {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tSS\tSS\tNoavgscore\tno scores\n" ;
		} else {
			print $fhout "$element\t$totdismates\t$nudiscreadTE\t$perallreadsalign\t$perreadfladepth\tND\tND\tNoavgscore\t@scoreslist\n";
		}
	}
	close $fhout;
}
sub print_array {
	my ($filetoprint,@array)= @_;
	open (my $ah,">","$path/$filetoprint") || die ("cannot open file (sub print_array) $filetoprint to write $!\n");
		foreach my $dataline (@array) {
			print $ah "$dataline\n";
		}
	close ($ah);
}
sub avg_readdepth {
	my ($genomelocatn) = @_;
	my $avg;
	my $total = 0;
	my $count = 0;
	system("samtools depth -a -r $genomelocatn -q 20 -Q 20 $bamlocation/$bamid.38.sorted.bam > $readdepthpath/$individual.$genomelocatn.readdepth.txt") == 0 or die ("failed extract readdepth via samtools\n"); 
	my $readdepthfile = "$readdepthpath/$individual.$genomelocatn.readdepth.txt";
	my $size = (-s $readdepthfile);
	if ($size == 0) {
		push (@logvalues,$readdepthfile);
		$avg = 0;
		#return ($avg);
	} else {
		open (my $rh, "<", "$readdepthfile") or confess "\n ERROR (main): could not open to read $readdepthfile $!\n";
			while(<$rh>) {
				chomp(my $line = $_);
				$count++ if  (!/^\s+?$/);#count non-blank lines
				my @columns = split(/\s+/,$line);
				$total += $columns[2];
			}
			$avg = $total/$count;
		close ($rh);
		
	} 
	return ($avg);
}
sub read_location {
	my ($bdfile) = @_;
	print STDERR "the file input is $bdfile\n" if ($verbose);
	open (my $dh, "<", $bdfile) or confess "\n ERROR (sub read_location): could not open to read $bdfile $!\n";
		while (my $data = <$dh>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @namebd = split(/\s+/,$data); # splitting the data based on space and storing into the arrray
			my $rdname = $namebd[3];#uniqueid-individual-readID==>reference position
			my $refposn = $namebd[0].":".$namebd[1]."-".$namebd[2];
			#print STDERR " read name is $rdname for $uniqueid and $individual and refposition is $refpos\n" if ($verbose);
			$bedfile{$uniqueid}{$rdname} = $refposn; #loading to the hash
		}
	return (%bedfile);
	#print Dumper %bedfile;
	close $dh;	
}
#need to output the name of the reads, identify the position of those reads, check the mappability score and output them 
sub compare_twoHOH {
	my ($readhash,$bedhash) = @_;
	my %readhash = %$readhash;
	my $head;
	my @referencepostions = ();#directly loaded into arrays as it is done for each unique id
	foreach my $uid (sort keys (%readhash)){
		#print STDERR "1st key in readhash is $uid\n";
		foreach my $readname (sort keys $readhash{$uid}) {
			if ($readname =~ /^(.*)\/(\d)/) {
				$head = "$1\/2" if ($2 == 1);
				$head = "$1\/1" if ($2 == 2);
				#$idlist{$head} =1;
			} 
			if (exists (${$bedhash}{$uid}{$head})) {
				#print STDERR " $head exists in bedhash\n";
				my $refposition = ${$bedhash}{$uid}{$head};
				#print STDERR " $refpos is the ref position of $head\n";
				my ($rchr,$rstart,$rend) = split (/[:-]/,$refposition);
				push (@referencepostions,{'start' => $rstart,'end' => $rend });
			} else {
				next;
			}
		}
	}
	my @sortedrefe = sort {$a->{start} <=>$b->{start}} @referencepostions;
	return (\@sortedrefe);
}

sub extract_overlap_bound_aoH {#used when cordinates are sorted
	my @rposi = @_;
	my %overlaps;
	my $min1 = 0;
	my $max1 = 0;
	my $region = 1;
	my @AoH = ();
	my $startcord;
	my $endcord;
	foreach my $genloc (@rposi) {
		#print "$genloc=\n";
		foreach my $key (keys %{$genloc}) {
			#print "$key => $genloc->{$key} \n";
			$startcord = $genloc->{$key} if ($key eq 'start');
			$endcord = $genloc->{$key} if ($key eq 'end');
		}	
		
			#print "$startcord = $endcord\n";
			if (($min1 == 0) && ($max1 == 0)) {
				$min1 = $startcord;
				$max1 = $endcord;
				push (@AoH, {$min1 => $max1});
				print "min and max loaded are $min1 = $max1\n";
				#next;
				#now min1 is the first value of the table751
			} else {
				for my $i ( 0 .. $#AoH ) {
					for my $astart ( sort (keys %{ $AoH[$i] } )) {
						my $aend=$AoH[$i]{$astart};
						#print  "in the array $start = $end ";
						if (($startcord >= $astart  ) && ($endcord >= $aend  ) && ($startcord < $aend  )) {
							$min1 = $astart;
							$max1 = $endcord;
							splice(@AoH,$i,1, {$min1 => $max1});
						} elsif ( (($startcord < $astart) && ($endcord < $astart)) || (($startcord > $aend) && ($endcord > $aend)) ) {
							my ($p);
							foreach my $hasref (@AoH) {
								#print "$p\n";
								foreach my $akey (keys %$hasref) {
									my $value = $hasref->{$akey};
									if (($startcord >=  $akey) && ($startcord <= $value)) {
										$p = 1;	
										last;
									} else {
										next;
									}
								}
							}
							push (@AoH, {$startcord => $endcord}) if (!defined $p);
						} 
					}
				}
			}
	}
		#print Dumper @AoH;	
	return (\@AoH);
}
sub extractgmscores {
	my ($cordinates,$rchr) = @_;
	my @mappingscores = ();
	
	foreach my $overlapp (@$cordinates) {
		foreach my $rstart (keys %$overlapp) {
			
			my @startdata =();
			my @enddata=();
			my $rend = $overlapp->{$rstart};
			my $givencordinate = $rchr.":".$rstart."-".$rend;
			if (exists ($overlapscore{$givencordinate})) {
				my $value = $overlapscore{$givencordinate};
				push (@mappingscores,$value);
			} else {
				my $len_avg;
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
				my $length = $tend - $tstart;
				while (my $row = $sth->fetchrow_hashref()) {
					my $score = $row->{'score'};
					print "Found a row: chr = $row->{'chromo'},$row->{'start'},$row->{'end'},score = $row->{'score'}\n" if ($verbose);
					print "*"x50,"\n"if ($verbose);
					$count++;
					$sum += $score;
				}
				$sth->finish();
				if ($count >= 1) {
					my $avg = $sum/$count ;
					$avg = sprintf("%.3f",$avg);
					#$len_avg = "$length => $avg";
					$len_avg = "$tstart-$tend => $avg";
					$overlapscore{$givencordinate} = $len_avg;
				}
				push (@mappingscores,$len_avg);
			}
		}
	}
	@mappingscores = uniq (@mappingscores);
	return (\@mappingscores);
}


sub normalise_readepth {
	foreach my $sampleindi (sort keys %allreaddepth) {
		print STDERR "-"x50,"\n"if ($verbose);
		print STDERR "$sampleindi\n"if ($verbose);
		my $sumrdindi = 0;
		my $totalloci = 0;
		foreach my $hash_ref (@{$allreaddepth{$sampleindi}}) {
			#print STDERR "$hash_ref\n";
			foreach my $pos (keys %{$hash_ref}) {#I need to exclude chrY from calculation or will have to add the sex info
				my $rd = ${$hash_ref}{$pos};
				print STDERR "$pos ==> $rd\n"if ($verbose);
				
				if ($pos =~ /^(chrY)\:(.*)/i) {
					next;
				} else {
					$sumrdindi += $rd ;
					$totalloci++ ;
				}
			}
			
		}
			my $indvidualavgrd = ($sumrdindi/$totalloci);
			$avgrdindvi{$sampleindi} = $indvidualavgrd;
	}
	return (%avgrdindvi);
}
sub calculate_readepth {
	my (%r) = @_;
	my $depthdiff;
	my $avrd_in;
	foreach my $unigenoloc (keys %rd_uid) {
		print STDERR "*"x50,"\n"if ($verbose);
		print STDERR "$unigenoloc\n"if ($verbose);
		my @parts = split (/\./,$unigenoloc);
		my $rdelement = $rd_uid{$unigenoloc};
		my $person = $parts[0];
		print STDERR "$person => $rdelement\n" if ($verbose);
		if (exists ($r{$person})) {
			$avrd_in = $r{$person};
			print STDERR "$person => $avrd_in\n"if ($verbose);	
		} else {
			warn ("avg readdepth for $person could not be identified\n");
		}
		if (! $rdelement == 0) {
			$depthdiff = ($rdelement/$avrd_in);
			
			$depthdiff = sprintf("%.2f",$depthdiff);
			print STDERR "$person => $depthdiff\n" if ($verbose);	
			push (@{$genomicseqblast{$unigenoloc}},{'read_depthratio' => $depthdiff});
		} else {
			$depthdiff = "undef";
			push (@{$genomicseqblast{$unigenoloc}},{'read_depthratio' => $depthdiff});
		}	
			
	}
		
}
sub email {
	my $to = 'jainythomas1@gmail.com';
	my $cc = 'jainyt@genetics.utah.edu';
	my $from = 'jainy@yoda.genetics.utah.edu';
	my $subject = 'script for env gag presence';
	my $message = "Results are ready for $path/$file";

	my $msg = MIME::Lite->new(From     => $from,
							  To       => $to,
							  Cc       => $cc,
							  Subject  => $subject,
							  Data     => $message);
	$msg->send;
	print "Email Sent Successfully\n";
}
