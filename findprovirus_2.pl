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
#use Cwd;
use Bio::SearchIO; 
use Bio::SeqIO;
use Bio::DB::Fasta;
use MIME::Lite;
use Data::Dumper;
use File::Copy;
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(firstidx);

my $version = "14.0_part2";
my $scriptname = "findprovirus.v14.0_part2forHERVW.pl";
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
#	- v14.0 = 10 November 2017
# 				split the script into two. First script performs an initial analysis. Second part of the script runs on the output of the first one and tries a bunch of different assemblers to see if soloLTR allele could be reconstructed
#				
#				
#			

\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <BAM ID table> -f <prediction output file from the first run> -p <path of the outputdirectory> -bl <location of bamfiles>[-v] [-c] [-h] [-s]
	
    MANDATORY ARGUMENT:	
    -t,--table 		  (STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file   	  (STRING) output containing genotype prediction
    -p,--path   	  (STRING) output directory name (path where extracted reads, mapped reads, extracted genomic sequences are found)	  
    -bl,--bamloc	  (STRING) location of bam files
    -st,--seqtkpro	  (STRING) path to the seqtk
    -pc,--picard	  (STRING) path to the picardtools
    -cp,--cap3		  (STRING) path to the cap3 assembler
    -bp,--blast		  (STRING) path to the blast 
    -bd,--bedtools	  (STRING) path to the bedtools
    -bu,--bamutils	  (STRING) path of bamutils
    -sp,--spade		  (STRING) path of spade
    -mp,--minia		  (STRING) path of minia assembler
    
    OPTIONAL ARGUMENTS:  
    -c,--chlog  	(BOOL)   Print updates
    -v,--v      	(BOOL)   Print version if only option
    -s,--verbose	(BOOL)   The script will talk to you
    -h,--help>  	(BOOL)   Print this usage\n\n";
   
#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($pfile,$table,$path,$out,$bamlocation,$seqtkpro,$BLASTpro,$CAP3pro,$bedtoolspro,$bamUtilpro,$SPAdepro,$miniapro,$verbose,$help,$v,$chlog);
GetOptions ('t=s' => \$table,
            'f=s' => \$pfile,            
            'p=s' => \$path,
            'o=s' => \$out,
            'bl=s'=> \$bamlocation,
            'st=s'=> \$seqtkpro,
            'cp=s'=> \$CAP3pro,
            'bp=s'=> \$BLASTpro,
            'bd=s'=> \$bedtoolspro,
            'bu=s'=> \$bamUtilpro,
            'sp=s'=> \$SPAdepro,
            'mp=s'=> \$miniapro,
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table) && (! $pfile) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $pfile) || !($table) || !($path) ||($help));
#my $cwd = getcwd();
#$path = $cwd if (!$path);
$out = "$path/Refined.prediction_alleles.txt" if (! $out);

#The user has to add the path
#my $seqtkpro = "/home/jainy/software/seqtk";
#my $picardpro = "/home/jainy/software/picard-2.9.2";
#my $CAP3pro = "/home/jainy/software/CAP3";
#my $BLASTpro = "/home/jainy/software/ncbi-blast-2.6.0+/bin";
#my $bamlocation = "/kbod2/WGS_DATA/SGDP_bams_public";
#my $SPAdepro = "/home/jainy/software/SPAdes-3.11.1-Linux/bin";
#my $bamUtilpro = "/home/jainy/software/bamUtil_1.0.13/bamUtil/bin";
#my $bamUtilpro = "/home/jainy/software/bamUtil";
#my $bedtoolspro = "/home/jainy/software/bedtools2/bin";
#my $miniapro = "/home/jainy/software/minia-v2.0.7-Source/build/bin";
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
my $extractgenomicseqout;
my $gchr;
my $gstart;
my $gend;
my $gstrand;
my $gexstart;
my $gexend;
my $renamedfile;
my $glength;
my @predictiondetails;
my $gextractloc;
my $renamedmappedreadfile;
my %predictiondata;
my $dmatefile;
my @secondassembly; 

%bamfile = &load_hash ($table);
open (my $ph, "<", $pfile) or confess "\n ERROR (main): could not open to read $pfile $!\n";
	while(<$ph>) {
		chomp (my $line = $_);
		my @col = split(/\t/,$line);
		$uniqueid = $col[0];
		$predictiondata{$line} = [@secondassembly];#loading into hash
		my @filename = split (/\./,$col[0]);
		$individual = $filename[0];
		$ltr = $filename[1];
		$genomeloc = $filename[2];
		($gchr,$gstart,$gend) = split (/[:-]/,$filename[2]);
		$gexstart = $gstart + 50;
		$gexend = $gend - 50;
		$glength = $gexend - $gexstart;
		$gextractloc = $gchr.":".$gexstart."-".$gexend;
		my $genoprediction = $col[6];
		if ($genoprediction =~ /^\(.\?\)/) {
			&find_bamid();
			my $epath = "$path/ExtractedReads/$ltr";
			my $mpath = "$path/MappedReads/$ltr";
			#sort the bam file on name
			system ("samtools sort -n -o $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.bam $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.bam") == 0 or die ("unable to name sort bam file $uniqueid \n");
			#Extract mapped reads from the bam
			system ("$bamUtilpro/bam bam2FastQ --in $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.bam --firstOut $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDsBU.namesorted.R1.fastq --secondOut $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDsBU.namesorted.R2.fastq --unpairedOut $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq --noReverseComp") == 0 or die ("unable to extract reads using bamUtil $uniqueid \n");
			unlink ("$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDsBU.namesorted.R1.fastq","$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDsBU.namesorted.R2.fastq") or warn "Could not unlink : $!";
			#extracting the R1 and R2 from bedtools as there is a bug with bamutil for R1 and R2
			system ("$bedtoolspro/bedtools bamtofastq -i $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.bam -fq $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq -fq2 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq") == 0 or die ("unable to extract reads using bamUtil $uniqueid \n");
			
			system ("$seqtkpro/seqtk seq -A -q20 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta") == 0 or die ("unable to run seqtk on R1 $uniqueid \n");
			system ("$seqtkpro/seqtk seq -A -q20 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta") == 0 or die ("unable to run seqtk on R2 $uniqueid \n");
			system ("$seqtkpro/seqtk seq -A -q20 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta") == 0 or die ("unable to run seqtk on Up $uniqueid \n");
			
			#concatenate reads for minia assembly
			system ("cat $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta > $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.concatenated.fasta" ) == 0 or die ("unable to concatenate R1,R2 and Upfasta on $uniqueid\n");
			
			copy("$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq", "$path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq") or die "Copy failed R1 fastq $uniqueid:$!"; 
			copy("$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq", "$path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq") or die "Copy failed R2 fastq $uniqueid:$!"; 
			copy("$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq", "$path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq") or die "Copy failed Up fastq $uniqueid:$!"; 
			
			copy("$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fastq", "$path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.fastq") or die "Copy failed mapped fastq $uniqueid:$!"; 
			
			my $fisizR1 = -s "$path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq";
			my $fisizR2 = -s "$path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq";
			if (($fisizR1 > 0) || ($fisizR2 > 0))  {
				#Assemble only the mapped reads		   
				system ("$SPAdepro/spades.py -1 $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq -2 $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq -s $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq --careful -o $path/MappedReads/$ltr/$uniqueid.mappedSPAdeout") == 0 or 
				system ("$SPAdepro/dipspades.py -1 $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fastq -2 $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fastq -s $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fastq -o $path/MappedReads/$ltr/$uniqueid.mappeddipSPAdeout") == 0 or 
				system ("$miniapro/minia -in $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.fastq -kmer-size 55 -abundance-min 3 -out $path/ExtractedReads/$ltr/$uniqueid.mappedreads.concatenated.fasta_k55_ma3") == 0 or 
				#system ("$CAP3pro/cap3 $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.concatenated.fasta > $path/MappedReads/$ltr/$uniqueid.onlymappedreadIDs.concatenated.assembld.fasta") == 0 or 
				die ("unable to assemble mapped reads for $uniqueid \n");
				print STDERR " Assembling fasta done \n" if ($verbose);
					#Rename scaffolds.fasta/capcontigs 
				if (-e "$path/MappedReads/$ltr/$uniqueid.mappedSPAdeout/scaffolds.fasta" ) {
					copy("$path/MappedReads/$ltr/$uniqueid.mappedSPAdeout/scaffolds.fasta", "$path/MappedReads/$ltr/$uniqueid.mapped.scaffolds.fasta") or die "Copy failed scaffolds.fasta $uniqueid:$!";
				
					print STDERR " 	processing of	$uniqueid is done\n" if ($verbose);	
					#rename the fasta sequence with its filename
				
					my $assembledfile = "$uniqueid.mapped.scaffolds.fasta";
					&renameseq_filename ($assembledfile,$mpath);
				}  elsif (-e "$path/MappedReads/$ltr/$uniqueid.mappeddipSPAdeout/dipspades/consensus_contigs.fasta") {# if scaffold.fasta was not created with SPAdes run
					copy("$path/MappedReads/$ltr/$uniqueid.mappeddipSPAdeout/dipspades/consensus_contigs.fasta", "$path/MappedReads/$ltr/$uniqueid.consensus_contigs.fasta") or die "Copy failed consensus_contigs.fasta $uniqueid:$!";
					my $assembledfile = "$uniqueid.consensus_contigs.fasta";
					&renameseq_filename ($assembledfile,$mpath);
				} elsif (-e "$path/MappedReads/$ltr/$uniqueid.mappedreads.concatenated.fasta_k55_ma3.contigs.fa") {
					my $assembledfile = "$uniqueid.mappedreads.concatenated.fasta_k55_ma3.contigs.fa";
					&renameseq_filename ($assembledfile,$mpath);
				} 
			} 
			if (-e "$path/ExtractedReads/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
				unless (-z "$path/ExtractedReads/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
					$dmatefile = "$uniqueid.discordantmatesreadIDlist.txt.fasta";
					&splitdismatetopairs($epath,$dmatefile);
					#concatenate mapped and discordant reads that is split based on read pair info 
					system ("cat $epath/$dmatefile.r1 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R1.fasta > $path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D1R1.fasta") == 0 or die ("unable to concatenate discordantsplitreads $uniqueid \n");
					system ("cat $epath/$dmatefile.r2 $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.R2.fasta > $path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D2R2.fasta") == 0 or die ("unable to concatenate discordantsplitreads $uniqueid \n");
					#concatenate mapped reads and discordant mates for minia assembly
					system ("cat $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.concatenated.fasta $path/ExtractedReads/$ltr/$uniqueid.discordantmatesreadIDlist.txt.fasta > $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta" ) == 0 or die ("unable to concatenate dismatereads on $uniqueid\n");
					my $fisizeDR1 = -s "$path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D1R1.fasta";#concatenated file
					my $fisizeDR2 = -s "$path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D2R2.fasta";#concatenated file
					if (($fisizeDR1 > 0) || ($fisizeDR2 > 0)) {		   
						#Assemble the reads to contigs
						system ("$SPAdepro/spades.py -1 $path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D1R1.fasta -2 $path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D2R2.fasta -s $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.DUp.fasta --careful --only-assembler -o $path/ExtractedReads/$ltr/$uniqueid.allreadsSPAdeout") == 0 or 
						system ("$SPAdepro/dipspades.py -1 $path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D1R1.fasta -2 $path/ExtractedReads/$ltr/$uniqueid.allreads.namesorted.D2R2.fasta -s $path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.DUp.fasta --only-assembler --expect-rearrangements -o $path/ExtractedReads/$ltr/$uniqueid.allreadsdipSPAdeout") == 0 or 
						system ("$miniapro/minia -in $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta -kmer-size 55 -abundance-min 3 -out $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta_k55_ma3") == 0 or 
						#system ("$CAP3pro/cap3 $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta > $path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.asemb.fasta ") == 0 or 
						die ("unable to assemble allreads for $uniqueid \n");
						print STDERR " Assembling fasta done \n" if ($verbose);
						if (-e "$path/ExtractedReads/$ltr/$uniqueid.allreadsSPAdeout/scaffolds.fasta") {
							#Rename scaffolds.fasta 
							copy("$path/ExtractedReads/$ltr/$uniqueid.allreadsSPAdeout/scaffolds.fasta", "$path/ExtractedReads/$ltr/$uniqueid.allreads.scaffolds.fasta") or die "Copy failed scaffolds.fasta $uniqueid:$!";
							#rename the fasta sequence with its filename
							#my $allpath = "$path/ExtractedReads/$ltr";
							my $assembledfile = "$uniqueid.allreads.scaffolds.fasta";
							&renameseq_filename ($assembledfile,$epath);
						} elsif (-e "$path/ExtractedReads/$ltr/$uniqueid.allreadsdipSPAdeout/dipspades/consensus_contigs.fasta") { 
							copy("$path/ExtractedReads/$ltr/$uniqueid.allreadsdipSPAdeout/dipspades/consensus_contigs.fasta", "$path/ExtractedReads/$ltr/$uniqueid.allreads.consensus_contigs.fasta") or die "Copy failed consensus_contigs.fasta $uniqueid:$!";
							my $assembledfile = "$uniqueid.allreads.consensus_contigs.fasta";
							&renameseq_filename ($assembledfile,$epath);
						} elsif (-e "$path/ExtractedReads/$ltr/$uniqueid.allreads.concatenated.fasta_k55_ma3.contigs.fa") { 
							my $assembledfile = "$uniqueid.allreads.concatenated.fasta_k55_ma3.contigs.fa";
							&renameseq_filename ($assembledfile,$epath);
						}
					} 
				}
			} else {
				warn "No discordant reads were identified $uniqueid\n";
				if (-e "$path/MappedReads/$ltr/Renamed_Assembledseq/$uniqueid.rename.fasta") {
					copy("$path/MappedReads/$ltr/Renamed_Assembledseq/$uniqueid.rename.fasta","$path/ExtractedReads/$ltr/Renamed_Assembledseq/$uniqueid.rename.fasta") or die "Copy failed renamed mapped reads.fasta $uniqueid:$!";
				} 
				
			}	
			my $paligned;
			my $pmappaligned;
			print STDERR "performing BLAST on the mapped eads assembly with genomic seq......\n";
			make_path  ("$path/Genomeblast_S2/$ltr");
			#$summarygenomicseqbl = "$path/Genomeblast/Summarygenomicseqblast.txt";
			$renamedmappedreadfile = "$path/MappedReads/$ltr/Renamed_Assembledseq/$uniqueid.rename.fasta";
			$extractgenomicseqout = "$path/ExtractGenomicsequences/$gextractloc.extract.seq.fa";
			my $gmappedblastout = "$path/Genomeblast_S2/$ltr/$uniqueid.gmappedblast.out";
			if (-e $renamedmappedreadfile) {
				my $sizefil = (-s $renamedmappedreadfile);
				if ($sizefil > 0) {
					unless (-e "$renamedmappedreadfile.nhr") {
						system ("$BLASTpro/makeblastdb -in $renamedmappedreadfile -dbtype nucl") == 0 or die ("unable to makeblastdb on $renamedmappedreadfile \n");
					}
					#system ("$BLASTpro/blastn -db $renamedmappedreadfile -query $extractgenomicseqout -evalue 0.0001 -out $gmappedblastout") == 0 or die ("unable to to perform gmappblast $uniqueid \n");
					system ("$BLASTpro/blastn -db $renamedmappedreadfile -query $extractgenomicseqout -evalue 0.0001 -outfmt 6 -out $gmappedblastout.tabular.out") == 0 or die ("unable to to perform gmapptabular blast $uniqueid \n");#fortableoutput
					$pmappaligned = &parse_blast_foralen("$gmappedblastout.tabular.out");
					push (@{$predictiondata{$line}},{'mappedreadsaligned' => $pmappaligned});
				} else {
					$pmappaligned = 0;
					push (@{$predictiondata{$line}},{'mappedreadsaligned' => $pmappaligned});
				}
			}		
			$renamedfile = "$path/ExtractedReads/$ltr/Renamed_Assembledseq/$uniqueid.rename.fasta";			
			my $gblastout = "$path/Genomeblast_S2/$ltr/$uniqueid.gblast.out";
			if (-e $renamedfile) {
				my $sizeoffil = (-s $renamedfile);
				if ($sizeoffil > 0) {
					unless (-e "$renamedfile.nhr") {
						system ("$BLASTpro/makeblastdb -in $renamedfile -dbtype nucl") == 0 or die ("unable to makeblastdb on $renamedfile \n");
					}
					#system ("$BLASTpro/blastn -db $renamedfile -query $extractgenomicseqout -evalue 0.0001 -out $gblastout") == 0 or die ("unable to to perform gblast $uniqueid \n");
					system ("$BLASTpro/blastn -db $renamedfile -query $extractgenomicseqout -evalue 0.0001 -outfmt 6 -out $gblastout.tabular.out") == 0 or die ("unable to to perform gtabular blast $uniqueid \n");#fortableoutput
					$paligned = &parse_blast_foralen("$gblastout.tabular.out");
					push (@{$predictiondata{$line}},{'perallreadsaligned' => $paligned});
				} else {
					$paligned = 0;
					push (@{$predictiondata{$line}},{'perallreadsaligned' => $paligned});
				}
			}		
		} else {
			next;
		}	
	}
#print Dumper %predictiondata;
&print_predictionAoH(%predictiondata);
&email();
close $ph;
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
		die " Not able to identify the bamid of the individual. Please check the input files $! \n" ;
	}
}
sub renameseq_filename {
	my ($contigfile,$fpath) = @_;
	make_path  ("$fpath/Renamed_Assembledseq");
	#copy ("$contigfile","$path/$contigfile.mod.fasta") or die "Copy failed to $path: $!";
	open (my $bhout, ">","$fpath/Renamed_Assembledseq/$uniqueid.rename.fasta") or die "\n ERROR (main): could not open to read $fpath/Renamed_Assembledseq/$uniqueid.rename.fasta $!\n";
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
		$d{$read} = $score;
		if (exists ($d{$read})) {
		 	my $currentscore = $d{$read};
		 	if ($currentscore > $score){
		 		next;
		 	} else {
		 		$d{$read} = $score;
		 	}
		}
	}
	my $noofdisread = keys (%d);
	return ($noofdisread);
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
		my $upperlimit = $glength +10;
		my $lowerlimit = $glength - 10;
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
	my $percaligned = (($alength/$glength)*100);
	$percaligned = sprintf ("%.2f",$percaligned);
	print STDERR "the max perc aligned is $percaligned\n";
	return ($percaligned); 
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
sub print_predictionAoH {
	my (%genomelocation) = @_;
	open (my $fhout,">",$out) || die ("cannot open file $out to write $!\n");
	print $fhout "#element\t#total_discordant\t#discordantmates_maps_ERV\t#percent_of_assembly_alignedto_refseq\t#read_depthratio\t#genotype_prediction\t#variable_allele\t#mapability_scores(cordinates=>scores)\t#perce_mappedreads_aligntorefseq_otherassembly\t#perce_allreads_align_otherassembly\t#predictionof'S'allel\n" ;
	foreach my $element (sort keys %genomelocation) {
 		
		my $perallreadsalign;
		my $permappedreadsalign;
 		#print "$element\t @{$genomelocation{$element}}\n";
# 		my $numelements = @{$genomelocation{$element}};
		foreach my $hash_ref (@{$genomelocation{$element}}) {
			foreach my $type (keys %{$hash_ref}) {
				#print STDERR "$element $type => ${$hash_ref}{$type}\n";
				$perallreadsalign = ${$hash_ref}{$type} if ($type eq 'perallreadsaligned');
				$permappedreadsalign = ${$hash_ref}{$type} if ($type eq 'mappedreadsaligned');
			}
		}
		if (($permappedreadsalign >= 95) && ($perallreadsalign >= 95)) {# changed 'or' to 'and'
			print $fhout "$element\t$permappedreadsalign\t$perallreadsalign\tS\n" ;
		}  elsif (($permappedreadsalign < 95) && ($perallreadsalign < 95)) {
			print $fhout "$element\t$permappedreadsalign\t$perallreadsalign\n" ;
		}  elsif (($permappedreadsalign == 0) && ($perallreadsalign == 0)) {
			print $fhout "$element\t$permappedreadsalign\t$perallreadsalign\tNA\tno reads\n" ;
		}  else {
			print $fhout "$element\t$permappedreadsalign\t$perallreadsalign\tND\n";
		}
	}
}
sub splitdismatetopairs {
	my ($dmpath,$dmfile) = @_;
	my %dmlistR1 = ();
	my %dmlistR2 = ();
	my %upaddlist = ();
	my @seqArrayR1 = ();
	my @seqArrayR2 = ();
	my @seqArrayUp = ();
	my $readio_obj = Bio::SeqIO->new(-file 	 => "$dmpath/$dmfile", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $dmfile $!\n";  
	while (my $seq = $readio_obj->next_seq() ){#loading the reads in the discordant reads into a hash and also into separate arrays
		my $header = $seq->display_id;
		if ($header =~ /^(.*)\/(\d)/) {
			if ($2 == 1) {
				$dmlistR1{$1}=1;
				push (@seqArrayR1,$seq);
			} elsif ($2 == 2) {
				$dmlistR2{$1}=1;
				push (@seqArrayR2,$seq);
			}
		} 
	}
	my $DUpfile = "$uniqueid.onlymappedreadIDs.namesorted.DUp.fasta";
	my $Upreadio_obj = Bio::SeqIO->new(-file 	 => "$path/ExtractedReads/$ltr/$uniqueid.onlymappedreadIDs.namesorted.Up.fasta", 
									 -format => 'fasta') 
								 or die "\t    ERROR - Failed to create SeqIO FH object from $dmfile $!\n";  
	while (my $seq = $Upreadio_obj->next_seq() ){
		my $header = $seq->display_id;
		if (exists ($dmlistR1{$header})) {#check if the read in unpaired file is present in the discordant read1 list has remove it discordant list
			delete ($dmlistR1{$header});
			$seq->display_id("$header/2");#modify display_id add /1 or /2
			push (@seqArrayR2,$seq);
		} elsif (exists ($dmlistR2{$header})) {
			delete ($dmlistR2{$header});
			$seq->display_id("$header/1");
			push (@seqArrayR1,$seq);
		} else {
			push (@seqArrayUp, $seq);
			next;
		}
	}
	
	 
	if (%dmlistR1) {# puttting the rest of the rest of the sequences in the discordant mate file to unpaired sequences
		foreach my $read (keys %dmlistR1){
			my $ind = firstidx { $_->display_id eq "$read/1" } @seqArrayR1;
			
			my $upseq = splice (@seqArrayR1,$ind,1);
			
			$read = substr $read, 0,-2;
			$upseq->display_id("$read");
			push (@seqArrayUp, $upseq);
		}
	} 
	
	if (%dmlistR2) {
		foreach my $read (keys %dmlistR2){
			my $ind = firstidx { $_->display_id eq "$read/2" } @seqArrayR2;
			
			my $upseq = splice (@seqArrayR2,$ind,1);
			
			$read = substr $read, 0,-2;
			$upseq->display_id("$read");
			
			push (@seqArrayUp, $upseq);
		}
	
	}
	
	
	&sortfasta_header($dmpath,"$dmfile.r1",@seqArrayR1);
	&sortfasta_header($dmpath,"$dmfile.r2",@seqArrayR2);
	&sortfasta_header($dmpath,"$DUpfile",@seqArrayUp);
	
}
sub sortfasta_header {
	my ($spath,$filetosort,@seqarray) = @_;
		
	my $sortedReadio_obj = Bio::SeqIO->new(-file   => ">$spath/$filetosort",
										      -format => 'fasta') 
									         or die "\t    ERROR - Failed to create SeqIO FH object from $spath/$filetosort $!\n";  
	@seqarray = sort { ($a->display_id cmp $b->display_id) } @seqarray;
	foreach my $seque (@seqarray) {
		$sortedReadio_obj->write_seq($seque);	
	}
}
sub email {
	my $to = 'jainythomas1@gmail.com';
	my $cc = 'jainyt@genetics.utah.edu';
	my $from = 'jainy@yoda.genetics.utah.edu';
	my $subject = 'script for env gag presence';
	my $message = "Results are ready for $pfile";

	my $msg = MIME::Lite->new(From     => $from,
							  To       => $to,
							  Cc       => $cc,
							  Subject  => $subject,
							  Data     => $message);
	$msg->send;
	print "Email Sent Successfully\n";
}
