# rnaSEQ_pre-process_data_small_machines
Complete pipeline to pre-process RNAseq paired end set of files from scratch working in machines without high computational capabilities. Currently working in Ubuntu 23.04 Lunar Lobster, 16 cores, 32Gb RAM, 1TB SSD. <br/>

This script combines python and bash tools with the aim to get a final csv with counts of reads against a reference genome. In this example I have precompiled with STAR the genome for five different plant species I work more frequently.

# Prerequisites
There are some specific programs not instaled by default but here you have minor instructions you can use to include without further problems with your sudo privileges in some specific cases.

1) You must have installed fastQC: pip install sequana-fastqc
2) You must have cutadapt installed: sudo apt-get install cutadapt
3) You must have installed sortmerna (more information here): https://github.com/sortmerna/sortmerna
4) You must have installed trimmomatic (http://archive.ubuntu.com/ubuntu/pool/universe/t/trimmomatic/trimmomatic_0.39+dfsg-1_all.deb)
5) You must have star installed (more information here): https://github.com/alexdobin/STAR
6) You must have subread installed in order to use featureCounts more information here): https://github.com/ShiLab-Bioinformatics/subread
merge-paired-reads.sh and unmerge-paired-reads.sh currently in folder with data (linked to /usr/share/sormerna/scripts)<br />
To generate the database indexes for sortmerna (directly copy the fasta files to the rRNA_databases folder and then generate the indexes)
$ indexdb_rna --ref ./silva-bac-16s-id90.fasta,./silva-bac-16s-id90-db:./silva-bac-23s-id98.fasta,./silva-bac-23s- id98-db:./silva-arc-16s-id95.fasta,./silva-arc-16s-id95-db:./silva-arc-23s-id98.fasta,./silva-arc-23s-id98- db:./silva-euk-18s-id95.fasta,./silva-euk-18s-id95-db:./silva-euk-28s-id98.fasta,./silva-euk-28s-id98-db:. /rfam-5s-database-id98.fasta,./rfam-5s-database-id98-db:./rfam-5.8s-database-id98.fasta,./rfam-5.8s-database-id98-db: <br />
If the indexes are not generated correctly, nothing will work correctly.
7) You must have multiQC installed: pip install multiqc
8) IMPORTANT: all STAR annotation genomes must be included in "/usr/local/bin/STAR_GENOMES/" with the specific genome name. More details in the code

These options are recommended ones, but not mandatory. You have alternatives for some of the programs mentioned above.

# How to use
Of course you have many possibilities. The simpliest, copy the .py file in the folder containing all the paired end files and execute in the terminal. You will be asked about the genome to use and how many processors can you employ.<br/>

![Initial_options](https://github.com/dlazarogi/rnaSEQ_pre-process_data_small_machines/assets/77571239/488e5d05-0c1c-493b-9d45-1a888333c4cd)

# Workflow
This pipe line is designed to take the raw data from the sequencer and follow the logical steps established to get to perform an RNAseq analysis with the data
for different species. Important, the files must be compressed as they come from the sequencer.<br />

* The first step is to quality check the fastqc scripts and move the results to a folder
* The second step is to join the pairs of files into one to perform a sequence deletion that can match with contaminating ribosomal sequences with SortmeRNA.
* The third step is to pass the files sequentially through SortmeRNA and get a new pair of paired end files.
* The fourth step is to separate the paired end files back into two separate files each with its own name
* The fifth step is to perform a new fastQC to know the quality of the parsed readings and enter them in a QC/SortmeRNA directory.
* The sixth step is removing adapters with trimmomatic. Install to /usr/local/bin (basically copies the entire contents of the folder in that directory).
* The seventh step is to do new fastQC for the trimmomatic output and enter the results in the QC/trimmomatic directory.
* The eighth step is the annotation of the STAR sequences against the reference genome. If you have it installed you don't need to have ../STAR/source/
* The ninth step is to generate .bai files. We use samtools.
* The tenth step is to count the number of reads that will align overlap with the ones we have mapped. Using feature accounts. Finally you get some information when everything is finished:<br/>
**Preprocess of data finished successfully<br/>
Please review the QC data for multiple files in the following directories:<br/>
./QC/raw/multiqc_report.html -----------> multiQC for raw data<br/>
./QC/SortmeRNA/multiqc_report.html -----> multiQC for reads from raw data remaining after SortmeRNA<br/>
./QC/Trimmomatic/multiqc_report.html ---> multiQC for SortmeRNA data remaining after Trimmomatic<br/>
./SortmeRNA ----------------------------> log files of processes and fq files with rejected reads to check<br/>
./Trimmomatic --------------------------> Summary of reads for each group after processed with Trimmomatic<br/>
./STAR_results/ ------------------------> Include all .bam and .bai files to review<br/>
./STAR_results/log ---------------------> Log progress for each condition processed and annotated<br/>
featureCounts_results.txt.summary is interesting to check in order to see how many information was finaly incuded.<br/>
featureCounts_results.csv is the file to be used in R**.
* The eleventh step would already be the DE analysis in my case I use own scripts in R.

# Notes
The aim is to use an easy code to automatize the pre-processing using the tools quoted before (fastQC,sortmerna,trimmomatic,STAR,and subread, multiQC) using your own pre-compiled STAR genomes included in the corresponding folder. I prefer to have them in "/usr/local/bin/STAR_GENOMES/", but you can adapt the code accordingly.<br/>
Good point is that as the steps are executed in paired end files, you can check the progress despite you did not complete the pre-processing. At least you can supervise the evolution of the process if you are impatient. Besides, if the process crashes for some reason you can control where it failed, so can move the last paired files from the "original_files" folder to the upper level and re-launch, keeping all the other processed files correctly stored. Of course you will have to detete the intermediate files that are stored as I do not move the files to destination until the process is complete. <br/>
Sometimes the SWAP is increasingly filled from time to time. There is a reminder during the process to indicate at least once a day to empty the  SWAP memory.
Some good practices is to keep in mind the path length. The script take this into account, and if the paht plus the file names are longer than 255 characters, the process will exit.
