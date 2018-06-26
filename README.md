# dimorphicERV

DimorphicERV is the integrated name of two pipelines called findprovirus and findsoloLTR. 
These pipelines are used to identify proviral deletions resulting from LTR recombination events
from whole genome resequencing paired-end data aligned to the reference genome. 

Both pipeline requires the start and end cordinates of solo LTR or proviral copies of the HERV family.

To generate the catalogue of solo LTR and proviral elements of a particular family the tool
One Code to Find them All (Bailly-Bechet et al. 2014) is used. The script called rename_mergedLTRelementsdotpl is used to identify the boundaries of element and provide a unique name to each copy.

Typical Usage: perl rename_mergedLTRelements.pl -f file that needs to renamed -ltr name of ltr=length of LTR -int name of internalsequence -ilen length of internal sequence -rn name that you would like give [-v version] [-c change log] [-h help]

	MANDATORY ARGUMENT:
    -f,--file          	(STRING) file
    -ltr,--ltrname     	(STRING) Name of the ltr that needs to be rejoined=length of the ltr that can be classified as soloLTR (~5-10 bp length less than consensus length) 
								e.g. -ltr MER66C=550 -ltr MER66B=481 -ltr MER66D=479
	-int,--intname  	(STRING) Name of the internal erv sequence that needs to be rejoined
	-ilen,--lenint		(STRING) length of the total internal sequence that can be classified as complete (~5-10 bp length less than consensus length)
    -rn,--rename		(STRING) the name that you would like to give to the element (eg. HERVH or HERVW or HERV17)
    
    OPTIONAL ARGUMENTS:
    -o,--out			(STRING) name of the output file
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help			(BOOL)   Print this usage
    
Output: A bed format file will be generated cordinates and with a unique name stating whether its a solo LTR or a 2 LTR provirus
Create a text file containing only solo LTRs by using grep "soloLTR" and create a text file containing only 2 LTR provirus by using grep "2LTR".
Create a input file for the next two scripts by adding the sample IDs to be tested as the first column.


findprovirus pipeline
------------------------------------------------------------------------------------------------------------
Find solo-LTR to provirus variants using findprovirus pipeline.

Need to generate a mysql (mariadb) database and an indexed table if need to obtain the mappability scores.
Typical Usage: perl findprovirus_1.pl -t BAM ID table -f file with ltr cordinates -bl location of bamfiles -b [-p path of the outputdirectory][-g path of the genome][-m mapscores][-te TEseq][-u Username] [-pd password][-db mysql database][-mt mysql table] [-i] [-e] [-x] [-v] [-c] [-h] [-s]
	
    MANDATORY ARGUMENT:	
    -t,--table 	      	(STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file         	(STRING) file containing cordinates of solo LTR (sampleID, chr, start, end (of solo LTR),unique identifier,length,strand ) 
    -bl,--bamlocation 	(STRING) location of bam files  
    -b,--both   		(BOOL)   run extraction and assembly with sliced bam file for viewing IGV
    
    OPTIONAL ARGUMENTS:  
    -p,--path   		(STRING) output directory name (path)
                         	 	 Default = current working directory
    -te,--teseq 		(STRING) consensus internal sequence of the provirus 
    -g,	--genome		(STRING) path of the genome
    -rd,--readepth  	(BOOL)   read depth analysis
    -x, --extract   	(BOOL)	 extract genomic sequence               
    -i,--igv    		(BOOL)   only IGV has to be run
    -e,--reads  		(BOOL)   extraction and assembling of reads 
    -m --mapscores 		(BOOL) 	 Need to calculate mappability scores	
    -db --mysqldbinfo 	(STRING) database=jainys_db, (needs to prepare before running script)   
    -u --user 			(STRING) Username for mysql database,(needs to prepare before running script) 
    -pd,--password 		(STRING) password for mysql database,(needs to prepare before running script) 
    -mt,--mysqltable	(STRING) mysql table containing mappability scores (needs to prepare before running script)  
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help			(BOOL)   Print this usage\n\n"


Output:The predictions are in the *.prediction_alleles.txt 
Note: parallel can be used to speed up this step. The input file can be split to multiple files. Then parallel is used when running the above command.     

if you want to try an alternate assembler the following script is recommended to run.


Typical Usage: perl findprovirus_2.pl -t BAM ID table -f prediction output file from the first run -p path of the outputdirectory -bl location of bamfiles[-v] [-c] [-h] [-s]
	
    MANDATORY ARGUMENT:	
    -t,--table 		(STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file   	(STRING) output containing genotype prediction
    -p,--path   	(STRING) output directory name (path where extracted reads, mapped reads, extracted genomic sequences are found)	  
    -bl,--bamloc	(STRING) location of bam files
    
    OPTIONAL ARGUMENTS:  
    -c,--chlog  	(BOOL)   Print updates
    -v,--v      	(BOOL)   Print version if only option
    -s,--verbose	(BOOL)   The script will talk to you
    -h,--help		(BOOL)   Print this usage

Output: The output is in Refine.prediction_alleles.txt. Reports if able to assemble solo LTR allele using an alternate assembler.

findsoloLTR pipeline
------------------------------------------------------------------------------------------------------------
Find provirus to solo LTR variants using findsoloLTR pipeline.

Typical Usage:perl findsoloLTR.pl -t table -f file with ltr cordinates -bl location of bamfiles [-m yes to find mappabilty] [-u Username] [-pd password][-db mysql database][-mt mysql table][-p path of the outputdirectory][-o output file] [-v] [-c] [-h] 
	
    MANDATORY ARGUMENT:	
    -t,--table			(STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file			(STRING) file containing accesion information (output from the script get_coverage_coordinates.pl script on bedtools)
    -bl,--bamlocation 	(STRING) location of bam files
      	  
    OPTIONAL ARGUMENTS:
    -m, --mappability  	(STRING)  yes if need to find mappability 
    -mt,--table 		(STRING) mysql table e.g.	hg19wgEncodeCrgMapabilityAlign100mer_index/wgEncodeCrgMapabilityAlign100merhg38_lo_index
    -p,--path         	(STRING) output directory name (path) Default = current working directory
    -db, --mysqldbinfo 	(STRING) ex. command: DBI:mysql:database=jainys_db;host=localhost;port=22;     
    -u, --user 			(STRING) Username for mysql databasemy e.g	jainy;
    -pd,--password 		(STRING) password for mysql database e.g. wysql123
    -o,--output 		(STRING) output file
    -i,--igv    		(BOOL)   get IGV files for the regions with 250 bp flank
    -c,--chlog  		(BOOL)   Print updates
    -v,--v      		(BOOL)   Print version if only option
    -s,--verbose		(BOOL)   The script will talk to you
    -h,--help  			(BOOL)   Print this usage

Output: The output is *.readdepth.output.txt.  

