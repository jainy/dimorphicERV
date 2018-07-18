# dimorphicERV

DimorphicERV is the integrated name of two pipelines called findprovirus and findsoloLTR. 
These pipelines identify proviral deletions resulting from LTR recombination events
using whole genome resequencing paired-end data aligned to the reference genome. 

Both pipeline requires the start and end cordinates of solo LTR or proviral copies of the HERV family.


findprovirus pipeline
------------------------------------------------------------------------------------------------------------
This pipeline helps to identify solo-LTR to provirus variants.

Inorder to obtain the mappability scores of HERV regions, an indexed mysql table of the mappability scores for the genome has to be generated. The script uses the data from the table to calculate the mappability scores of the regions where discordant reads are mapped.

Usage: perl findprovirus_1.pl -t BAM ID table -f file with ltr cordinates -bl location of bamfiles -b [-p path of the outputdirectory][-g path of the genome][-te TEseq][-m ][-u Username] [-pd password][-db mysql database][-mt mysql table] [-i] [-e] [-x] [-v] [-c] [-h] [-s]
	
    MANDATORY ARGUMENT:	
	-t,--table 	      	(STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file         	(STRING) file containing cordinates of solo LTR (sampleID, chr, start, end (of solo LTR),unique identifier,length,strand ) 
    -bl,--bamlocation 	(STRING) location of bam files  
    -b,--both   		(BOOL)   run extraction and assembly with sliced bam file for viewing IGV
    
    OPTIONAL ARGUMENTS:  
    -p,--path   		(STRING) output directory name (path)
                         	 	 Default = current working directory
    -te,--teseq 		(STRING) consensus internal sequence of the provirus 
    -g,--genome   		(STRING) path of the genome
    -rd,--readepth  	(BOOL)   read depth analysis
    -x, --extract   	(BOOL)	 extract genomic sequence               
    -i,--igv    		(BOOL)   only IGV has to be run
    -e,--reads  		(BOOL)   extraction and assembling of reads 
    -m --mapscores 		(BOOL) 	 Need to calculate mappability scores	
    -db --mysqldbinfo 	(STRING) database=jainys_db, (needs to prepare before running script)   
    -u --user    		(STRING) Username for mysql database,(needs to prepare before running script) 
    -pd,--password 		(STRING) password for mysql database,(needs to prepare before running script) 
    -mt,--mysqltable	(STRING) mysql table containing mappability scores (needs to prepare before running script)  
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help    		(BOOL)   Print this usage


Output:The predictions are in the *.prediction_alleles.txt 

	Col 1: HERV name, location and individual name
	Col 2: Total number of discordant reads identified
	Col 3: Number of discordants reads whose mates have significant homology with the internal sequence of provirus
	Col 4: Percent of the de novo assembled contig aligned to reference solo LTR allele (best hit when using BLAST)
	Col 5: Ration of average read depth at the solo LTR to the average of read depths of all solo LTRs
	Col 6: Predicted genotype (S for solo LTR, P for Provirus)
	Col 7: Allele for which the prediction is not well-supported
	Col 8: Average mappability score
	Col 9: Cordinates of non-overlapping regions where discordant reads are mapped and the respective average mappability scores


   

if you want to try an alternate assembler, the following script is recommended to run on the output from first (*.prediction_alleles.txt).

Usage: perl findprovirus_2.pl -t BAM ID table -f prediction output file from the first run -p path of the outputdirectory -bl location of bamfiles[-v] [-c] [-h] [-s]
	
    MANDATORY ARGUMENT:	
    -t,--table   	(STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file   	(STRING) output from findprovirus_1 script (*.prediction_alleles.txt)
    -p,--path   	(STRING) output directory name (path where extracted reads, mapped reads, extracted genomic sequences are found from the previous run)	  
    -bl,--bamloc	(STRING) location of bam files
    
    OPTIONAL ARGUMENTS:  
    -c,--chlog  	(BOOL)   Print updates
    -v,--v      	(BOOL)   Print version if only option
    -s,--verbose	(BOOL)   The script will talk to you
    -h,--help  	(BOOL)   Print this usage

Output: The output is Refine.prediction_alleles.txt. Reports if able to assemble solo LTR allele using an alternate assembler.



findsoloLTR pipeline
------------------------------------------------------------------------------------------------------------
This pipeline identifies provirus to solo LTR variants.

Usage:perl findsoloLTR.pl -t table -f file with ltr cordinates -bl location of bamfiles [-m] [-u Username] [-pd password][-db mysql database][-mt mysql table][-p path of the outputdirectory][-o output file] [-v] [-c] [-h] 
	
    MANDATORY ARGUMENT:
    -t,--table  		(STRING) file with first column needs to be the sampleIDs, second column BAMIDs
    -f,--file   		(STRING) file containing cordinates of solo LTR (sampleID, chr, start, end (of solo LTR),unique identifier,length,strand ) 
    -bl,--bamlocation 	(STRING) location of bam files
      	  
    OPTIONAL ARGUMENTS:
    -m, --mappability  	(BOOL)  if need to find mappability 
    -mt,--table 		(STRING) mysql table e.g.	hg19wgEncodeCrgMapabilityAlign100mer_index/wgEncodeCrgMapabilityAlign100merhg38_lo_index
    -p,--path         	(STRING) output directory name (path) Default = current working directory
    -db, --mysqldbinfo	(STRING) ex. jainys_db
    -u, --user  		(STRING) Username for mysql databasemy e.g	jainy;
    -pd,--password		(STRING) password for mysql database e.g. wysql123
    -o,--output  		(STRING) output file
    -i,--igv    		(BOOL)   get IGV files for the regions with 250 bp flank
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help   		(BOOL)   Print this usage

Output: The output is *.readdepth.output.txt. 

	Col 1: HERV name, location and individual name
	Col 2: Average Read depth at 5' 250 bp
	Col 3: Average Read depth across the HERV
	Col 4: Average Read depth at 3' 250 bp
	Col 5: Predicted genotype (2 (two provirus alleles), 1 (1 provirus and 1 solo LTR), 0 (two solo LTR alleles))
	Col 6: Percent Read depth  (Read depth at the HERV/(Average of Read depth at flanks))*100
	Col 7: Mappability score at the provirus region, length of the provirus => average score and weighted score where score is mulitplied by its corresponding length
	

Requirements
------------------------------------------------------------------------------------------------------------
These are required for running findprovirus_1 and findsoloLTR scripts. These softwares has to be installed and the location has to be added to the script

samtools (http://www.htslib.org/)
seqtk (https://github.com/lh3/seqtk)
picardtools (https://github.com/broadinstitute/picard)
cap3(http://seq.cs.iastate.edu/cap3.html)
blast(ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
bedtools(https://github.com/arq5x/bedtools2/releases/download/)

for findprovirus_2
SPAdes(http://cab.spbu.ru/software/spades/)
bamUtil(https://github.com/statgen/bamUtil)


Preparing the input files
------------------------------------------------------------------------------------------------------------


1.To generate the catalogue of solo LTR and proviral elements of a particular family the tool
One Code to Find Them All (Bailly-Bechet et al. 2014) is used on the repeatmasker output.
 
2. The output of One Code to Find Them All is parsed to identify the boundaries of each copy and provide a unique name to each copy using rename_mergedLTRelements.pl (see utils folder).

Typical Usage: perl rename_mergedLTRelements.pl -f file that needs to renamed -ltr name of ltr=length of LTR -int name of internalsequence -ilen length of internal sequence -rn name that you would like give [-v version] [-c change log] [-h help]

	MANDATORY ARGUMENT:
    -f,--file          	(STRING) file
    -ltr,--ltrname     	(STRING) Name of the ltr that needs to be rejoined=length of the ltr that can be classified as soloLTR (~5-10 bp length less than consensus length) 
								e.g. -ltr MER66C=550 -ltr MER66B=481 -ltr MER66D=479
	-int,--intname  	(STRING) Name of the internal erv sequence that needs to be rejoined
	-ilen,--lenint		(STRING) length of the total internal sequence that can be classified as complete (~5-10 bp length less than consensus length)
    -rn,--rename		(STRING) the name that you would like to give to the element (eg. HERVH or HERVW or HERV17)
    
    OPTIONAL ARGUMENTS:
    -o,--out    		(STRING) name of the output file
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help    		(BOOL)   Print this usage
    
Output: A bed format file will be generated cordinates and with a unique name stating whether its a solo LTR or a 2 LTR provirus

3. Use the command 'grep' "soloLTR" and "2LTR" to create two text files containing solo LTRs and provirus
4. To add the sample IDs that need to be tested as the first column, makelist.pl provided in the util folder is used.
5. To speed up the script, findprovirus_1.pl, the output file from makelist.pl is split to multiple files using  called 'splitfile_for_parallel.pl' (provided in the util folder).  Using 'parallel' (https://www.gnu.org/software/parallel/parallel_tutorial.html) multiple jobs are launched.  

Questions
------------------------------------------------------------------------------------------------------------
Please contact Jainy Thomas (jainyt@genetics.utah.edu) for questions or support.
