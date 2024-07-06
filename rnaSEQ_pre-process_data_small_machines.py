#!/usr/bin/python3

# Prerequisites
# The name of the files including the full path to them must not exceed.
# all 255 characters long, or there will be problems.
#1) You must have cutadapt installed (sudo in general for Ubuntu, but pip also ok): sudo apt-get install cutadapt or pip install cutadapt
#2) You must have star installed: sudo apt-get install star
#3) You must have subread installed in order to use featureCounts: sudo apt-get install subread
# 4) You must have multiqc installed to be able to do the evaluation of all processed files: pip install multiqc
# 4) merge-paired-reads.sh and unmerge-paired-reads.sh currently in folder with data (linked to /usr/share/sormerna/scripts in V2)
# To generate the database indexes for sortmerna (directly copy the fasta files to the rRNA_databases folder and then generate the indexes)
# $ indexdb_rna --ref ./silva-bac-16s-id90.fasta,./silva-bac-16s-id90-db:./silva-bac-23s-id98.fasta,./silva-bac-23s- id98-db:./silva-arc-16s-id95.fasta,./silva-arc-16s-id95-db:./silva-arc-23s-id98.fasta,./silva-arc-23s-id98- db:./silva-euk-18s-id95.fasta,./silva-euk-18s-id95-db:./silva-euk-28s-id98.fasta,./silva-euk-28s-id98-db:. /rfam-5s-database-id98.fasta,./rfam-5s-database-id98-db:./rfam-5.8s-database-id98.fasta,./rfam-5.8s-database-id98-db:
# If the indexes are not generated correctly, nothing will work correctly.
# This pipe line is designed to take the raw data from the sequencer and follow the
# logical steps established to get to perform an RNAseq analysis with the data
# for different species. Important, the files must be uncompressed, there is a moment that
# it is necessary to have them and I prefer to do it from scratch.
# The first step is to quality check the fastqc scripts and move the results to a folder
# The second step is to join the pairs of files into one to perform a sequence deletion
# that can match with contaminating ribosomal sequences with SortmeRNA.
# The third step is to pass the files sequentially through SortmeRNA and get a new pair of paired end files.
# The fourth step is to separate the paired end files back into two separate files each with its own name
# The fifth step is to perform a new fastQC to know the quality of the parsed readings and enter them in a
# QC/SortmeRNA directory.
# The sixth step is removing adapters with trimmomatic. Install to /usr/local/bin (basically copies the entire contents of the folder
# in that directory).
# The seventh step is to do new fastQC for the trimmomatic output and enter the results in the
# QC/trimmomatic directory.
# The eighth step is the annotation of the STAR sequences against the reference genome. If you have it installed you don't need to have ../STAR/source/
# The ninth step is to generate .bai files. We use samtools.
# The tenth step is to count the number of reads that will align overlap with the ones we have mapped. Using feature accounts.
# The eleventh step would already be the DE analysis in R.


import signal
import sys
import glob
import os
import subprocess
import time
import shutil
import pathlib
import shlex
##########################Force process cancelation############################

def custom_ctrl_c_handler(sig, frame):
    print("\nUps. You pressed Ctrl+C. Maybe something went wrong for you.")
    # Add any custom cleanup or additional actions here
    sys.exit(0)  # Exit the program gracefully

# Set the custom Ctrl+C handler
signal.signal(signal.SIGINT, custom_ctrl_c_handler)

###############################################################################


###############################check dependencies##############################
# List of libraries to check and install if missing
libraries_to_check = ["glob", "subprocess", "time","time","shutil","pathlib","shlex","multiqc"]

# Function to check and install libraries
def check_and_install_libraries(libraries):
    for lib in libraries:
        try:
            importlib.import_module(lib)  # Try to import the library
            #print(f"{lib} is already installed.")
        except ImportError:
            print(f"{lib} is not installed. Installing...")
            try:
                subprocess.check_call(["pip", "install", lib])  # Use pip to install the library
                print(f"{lib} has been successfully installed.")
            except Exception as e:
                print(f"Error installing {lib}: {e}")


################################################################################


#############################check path adequation##############################
# Specify the folder path
folder_path = os.getcwd()

# Iterate through files in the folder
for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)
    
    # Check if the total length of the file path is greater than 255 characters
    if len(file_path) > 255:
        print(f"File path {file_path} exceeds 255 characters.")
    else:
        pass
##########################################################################################
print('\033[1;37;44m _____  _   _                                      _            _ _            \033[0;0m\n'+
'\033[1;37;44m|  __ \| \ | |   /\                               (_)          | (_)           \033[0;0m\n'  +  
'\033[1;37;44m| |__) |  \| |  /  \ ______ ___  ___  __ _   _ __  _ _ __   ___| |_ _ __   ___ \033[0;0m\n'+
'\033[1;37;44m|  _  /| . ` | / /\ \______/ __|/ _ \/ _` | | \'_ \| | \'_ \ / _ \ | | \'_ \ / _ \\\033[0;0m\n'+
'\033[1;37;44m| | \ \| |\  |/ ____ \     \__ \  __/ (_| | | |_) | | |_) |  __/ | | | | |  __/\033[0;0m\n'+
'\033[1;37;44m|_|  \_\_| \_|_/    \_\    |___/\___|\__, | | .__/|_| .__/ \___|_|_|_| |_|\___|\033[0;0m\n'+
'\033[1;37;44m                                        | | | |     | |                        \033[0;0m\n'+
'\033[1;37;44m                                        |_| |_|     |_|                        \033[0;0m')

print("\n \nPipeline under development by David Lázaro Gimeno")
print("\ncontact: david{dot}lazaro{dot}gimeno@gmail.com, david{dot}lazaro{dot}gimeno@umu.se")
print("\nPlease make sure your files are compressed as .fq.gz, otherwise they process will fail")

print("\n \nHow many threads can you use?")
threads = input()

print("""
1) Arabidopsis thaliana  2) Picea abies v1.0
3) Populus tremula x tremuloides 4) Populus tremula (Potra) v2.2
5) Populus trichocarpa""")

genome=input("Select one genome of reference: ")

if genome =="1":
    genome="/usr/local/bin/STAR_GENOMES/Tair10/" 
    feature="/usr/local/bin/gff/Tair10_3.gtf"
elif genome =="2":
    genome="/usr/local/bin/STAR_GENOMES/Picea_annotation_STAR_original/" 
    feature="/usr/local/bin/gff/Pabies01-gene.gtf"
elif genome =="3":
    genome="/usr/local/bin/STAR_GENOMES/T89/" 
    feature="/usr/local/bin/gff/Potrx01-genome.gtf"
elif genome =="4":
    genome="/usr/local/bin/STAR_GENOMES/Potra/" 
    feature="/usr/local/bin/gff/Potra02_genescopy.gff" # Replaced ID by gene_id in all document. Worked
elif genome =="5":
    genome="/usr/local/bin/STAR_GENOMES/Potri/" 
    feature="/usr/local/bin/gff/Ptrichocarpa_210_v3.0.gtf"
else:
    print("\n\nYour option is not valid. \nThe script will stop. \nTry again\n")
    exit(1)

#sortme_rna path
sort ="/usr/share/sortmerna/rRNA_databases/"

print(f"\nOK, we're going to use {threads} threads... analyzing files in batch process")

day = time.strftime("%y/%m/%d")
hour = time.strftime("%H:%M:%S")
print ('\n\nUncompressing original files at \033[0;94;48m' + hour + '\033[0;0m on \033[0;94;48m' + day + '\033[0;0m')
cwd = os.getcwd()

#Uncompress.gz files
lista =  glob.glob('*.gz')
drlist1 = sorted (lista)
os.makedirs('./original_files', exist_ok = True)
while len(drlist1) > 0:
        compr1 = drlist1[0]
        compr2 = drlist1[1]
        p = subprocess.Popen('gunzip -f -k '+ compr1,
                            shell=True, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        out, err = p.communicate()
        shutil.move(os.path.join (cwd,compr1), os.path.join (cwd + '/original_files', compr1)) 
        p = subprocess.Popen('gunzip -f -k '+ compr2,
                            shell=True, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        out, err = p.communicate()
        shutil.move(os.path.join (cwd,compr2), os.path.join (cwd + '/original_files', compr2))



        day = time.strftime("%y/%m/%d")
        hour = time.strftime("%H:%M:%S")
        print ('\nFiles uncompressed and moved to original_files folder at \033[0;94;48m' + hour + '\033[0;0m on \033[0;94;48m' + day + '\033[0;0m')

        day = time.strftime("%y/%m/%d")
        hour = time.strftime("%H:%M:%S")
        #Quality control of raw sequences with .fq extension
        print('\nfastQC for raw data\n' + "The process starts at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day + '\033[0;0m')

        print ('Files to be processed:')
        raw =  glob.glob('*.fq')
        drlist = sorted (raw)
        while len(drlist) > 0:
            str1 = drlist[0]
            print(str1)
            del drlist[0:1]
        p = subprocess.Popen('fastqc -t {threads} *.fq',
                             shell=True, 
                             stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE)
        out, err = p.communicate()

        os.makedirs('./QC/raw', exist_ok = True)
        for ext in ["zip","html"]:
             for file in glob.glob('./*.{}'.format(ext)):
                 shutil.move(os.path.join (cwd,file), os.path.join (cwd + '/QC/raw', file)) #This way files will be overwrite if already exist

        day = time.strftime("%y/%m/%d")
        hour = time.strftime("%H:%M:%S")
        print("The raw lectures were QC completed at \033[0;94;48m" + hour + " on " + day + '\033[0;0m')

        #Work with couple of files
        lista =  glob.glob('*.fq')

        #Order the list so it can be processed by on a logical sequence
        #merge by paris and move he raw files to a folder and leave the required ones available
        drlist = sorted (lista)

        #merge-paire-reads.sh
        while len(drlist) > 0:
            str1 = drlist[0]
            #print (str1)
            str1_0 = os.path.splitext(str1)[0] #Name without extension
            str2 = drlist[1]
            #print (str2)
            str2_0 = os.path.splitext(str2)[0]
            print('\nThe files \033[0;91;48m' + str1 + '\033[0;0;48m and \033[0;91;48m' + str2  + '\033[0;0;48m to be merged')
            day = time.strftime("%y/%m/%d")
            hour = time.strftime("%H:%M:%S")
            print("The process starts at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day + '\033[0;0m')
            p = subprocess.Popen ('bash /usr/share/sortmerna/scripts/merge-paired-reads.sh '+ (str1) +' '+ (str2) + ' '+ (str1_0+str2_0)+'.fq',
                shell=True,
                executable='/bin/bash',
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE)
            out, err = p.communicate()
            os.remove(str1)
            os.remove(str2)

            print('\nThe files \033[0;92;48m' + str1 + '\033[0;0m and \033[0;92;48m'+ str2 + '\033[0;0m have ben removed. Original files are safe')
            day = time.strftime("%y/%m/%d")
            hour = time.strftime("%H:%M:%S")
            print("The merged files were completed at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day +'\033[0;0m')

            del drlist[0:2]

        #sortmeRNA
        lista = glob.glob('*.fq')
        drlist = sorted (lista)

        while len(drlist) > 0:
           str1 = drlist[0]
           str1_0 = os.path.splitext(str1)[0] #Name without extension
           print ('\nThe file \033[0;92;48m' + str1 + '\033[0;0m is being processed by sortmeRNA')
           day = time.strftime("%y/%m/%d")
           hour = time.strftime("%H:%M:%S")
           print("The process starts at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day + '\033[0;0m')
           dir_short = ('/usr/share/sortmerna/rRNA_databases/')
           lista_sortme = glob.glob(dir_short+'*.fasta')
        #print (lista_sortme). I need to improve this part
           cmd = (
           f'sortmerna '
           f'--ref {sort}rfam-5.8s-database-id98.fasta,{sort}rfam-5.8s-database-id98-db:'
           f'{sort}silva-bac-16s-id90.fasta,{sort}silva-bac-16s-id90-db:'
             f'{sort}rfam-5s-database-id98.fasta,{sort}rfam-5s-database-id98-db:'
             f'{sort}silva-bac-23s-id98.fasta,{sort}silva-bac-23s-id98-db:'
             f'{sort}silva-arc-16s-id95.fasta,{sort}silva-arc-16s-id95-db:'
             f'{sort}silva-euk-18s-id95.fasta,{sort}silva-euk-18s-id95-db:'
             f'{sort}silva-arc-23s-id98.fasta,{sort}silva-arc-23s-id98-db:'
             f'{sort}silva-euk-28s-id98.fasta,{sort}silva-euk-28s-id98-db '
             f'--reads {str1} --paired_in -a {{threads}} --log --fastx '
             f'--aligned {str1_0}_rRNA --other {str1_0}_subset_sortme_rna_p'
         )

           p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
           out, err = p.communicate()
           #del drlist[0:1]
        # Move results to the sortmeRNA folder
           os.makedirs('./SortmeRNA', exist_ok = True)
           for ext in [".log","_rRNA.fq"]:
             for file in glob.glob('./*{}'.format(ext)):
                 shutil.move(os.path.join (cwd,file), os.path.join (cwd + '/SortmeRNA', file)) #file, './SortmeRNA'
           day = time.strftime("%y/%m/%d")
           hour = time.strftime("%H:%M:%S")
           print ('\nResults moved to folder SortmeRNA\nNow is time to recover the two files\n \n\033[0;91;48mPlease remember to empty sap memory as it is overloaded after some time\nsudo swapoff -a and sudo swapon -a in a different terminal twice per day\nduring the pre processing is good enough.\033[0;0m')

           print("\nFiles by SortmeRNA were completed at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day + '\033[0;0m')
   
           #Define the two names to split
           half=int(len(str1_0)/2)
           first_file = (str1_0[0:half])# First file
           print('\033[0;92;48m'+first_file)
           second_file = (str1_0[half:]) # Second file
           print(second_file + '\033[0;0m')
           for ext in ['_rna.fq']:
               p = subprocess.Popen ('bash /usr/share/sortmerna/scripts/unmerge-paired-reads.sh ' + str1_0 + '_subset_sortme_rna_p.fq ' + first_file + '_subset_sortme_rna_p.fq ' + second_file + '_subset_sortme_rna_p.fq',
                    shell=True,
                    executable='/bin/bash',
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
               out, err = p.communicate()
               os.remove(str1_0 + '_subset_sortme_rna_p.fq') # COMBINED PARSED FILE
               print('\nRemoved paired files\n\nResults are new imput for fastQC')
               #fastQC for SortmeRNA and MV results to QC
               day =time.strftime("%y/%m/%d")
               hour = time.strftime("%H:%M:%S")
               print("\nfastQC analysis starts at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day + '\033[0;0m')
               p = subprocess.Popen('fastqc -t {threads} ' + first_file + '_subset_sortme_rna_p.fq ' +second_file + '_subset_sortme_rna_p.fq',
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
               out, err = p.communicate()
               os.makedirs('./QC/SortmeRNA', exist_ok = True)
           for ext in ["zip","html"]:
               for file in glob.glob('./*.{}'.format(ext)):
                   shutil.move(os.path.join (cwd,file), os.path.join (cwd + '/QC/SortmeRNA', file)) #This way files will be overwrite if already exist
           os.remove(first_file+second_file+'.fq')
           print('Results moved to QC/SortmeRNA')



        #Generate trimmomatic file executable on terminal. Check if I can replace -threads 6 by -threads {threads}
        #Teitur indicated to modify TruSeq3-PE.fa:2:30:10:2:KeepBothReads by TruSeq3-PE.fa:2:30:10:2:True
           day = time.strftime("%y/%m/%d")
           hour =time.strftime("%H:%M:%S")
           print("\nNow is time to trim the sequences with Trimmomatic \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day+ '\033[0;0m')
           file= open ('trim.sh','w')
           file.write ('#!/bin/bash\njava -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads '+repr(threads)+' -phred33 -summary summary.txt '+
           first_file + '_subset_sortme_rna_p.fq ' +second_file + '_subset_sortme_rna_p.fq '+
           first_file + '_paired.fq.gz ' + first_file + '_unpaired.fq.gz ' +
           second_file + '_paired.fq.gz ' + second_file + '_unpaired.fq.gz ' +
           'ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36' )# recuerda cambiar PE por SE cuando trabajes con single end
           file.close()

           for ext in ['.sh']:
                p = subprocess.Popen ('bash ./trim.sh',
                    shell=True,
                    executable='/bin/bash',
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
                out, err = p.communicate()
                os.remove('trim.sh') #It is better to delete it and generate a new one afterwards. Evaluate to generate a hidden file
           day = time.strftime("%y/%m/%d")
           hour = time.strftime("%H:%M:%S")
           summary_trim = "summary.txt"
           os.rename(cwd + '/' + summary_trim,cwd + '/' + first_file+second_file+'_'+summary_trim)
           os.makedirs('./Trimmomatic', exist_ok = True)
           for ext in ["txt"]:
               for file in glob.glob('./*.{}'.format(ext)):
                   shutil.move(os.path.join (cwd,file), os.path.join (cwd + '/Trimmomatic', file)) #This way files will be overwrite if already exist
           ren1 = first_file + '_unpaired.fq.gz'
           ren2 =  second_file + '_unpaired.fq.gz'
           shutil.move(cwd + '/' + ren1 , cwd + '/Trimmomatic/' + ren1)
           shutil.move(cwd + '/' + ren2 , cwd +'/Trimmomatic/' + ren2)
           print('Results moved to Trimmomatic')
           print("\nTrimmomatic process finished at \033[0;94;48m"+ hour + "\033[0;0m  on \033[0;94;48m" + day + '\033[0;0m ')
        #fastQC for trimmomatic output results, but before I have to move the results of unpaired.fq.gz
           day = time.strftime("%y/%m/%d")
           hour = time.strftime("%H:%M:%S")
           print("\nfastQC analysis starts at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day + '\033[0;0m')
           p = subprocess.Popen('fastqc -t {threads} ' + first_file + '_paired.fq.gz ' +second_file + '_paired.fq.gz',
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
           out, err = p.communicate()
           os.makedirs('./QC/Trimmomatic', exist_ok = True)
           for ext in ["zip","html"]:
               for file in glob.glob('./*.{}'.format(ext)):
                   shutil.move(os.path.join (cwd,file), os.path.join (cwd + '/QC/Trimmomatic', file)) #This way files will be overwrite if already exist
           print('Results moved to QC/Trimmomatic')

           os.remove(first_file + '_subset_sortme_rna_p.fq') # Input files for Trimmomatic no more necessary. Maybe can be interesting to keep them, but many big files. Evaluate
           os.remove(second_file + '_subset_sortme_rna_p.fq') # Input files for Trimmomatic no more necessary. Maybe can be interesting to keep them, but many big files. Evaluate

        #STAR
           day = time.strftime("%y/%m/%d")
           hour = time.strftime("%H:%M:%S")
           print("\nSequence annotation with STAR starts at \033[0;94;48m"+ hour + "\033[0;0m on \033[0;94;48m" + day + '\033[0;0m')
           p = subprocess.Popen('STAR --runThreadN 4 --genomeDir '+ repr(genome) + '--readFilesIn ./' +
           first_file +'_paired.fq.gz'+ ',./'+ second_file + '_paired.fq.gz ' + '--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate',
                     shell=True, 
                     stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)
           out, err = p.communicate()
           os.remove (first_file + '_paired.fq.gz')
           os.remove (second_file + '_paired.fq.gz')
           ren1 = 'SJ.out.tab'
           ren2 = 'Aligned.sortedByCoord.out.bam'
           ren3 = 'Log.final.out'
           ren4 = 'Log.out'
           ren5 = 'Log.progress.out'
           os.rename(cwd + '/' + ren1,cwd + '/' + first_file+second_file+'_'+ ren1)
           os.rename(cwd + '/' + ren2,cwd + '/' + first_file+second_file+'_'+ ren2)
           os.rename(cwd + '/' + ren3,cwd + '/' + first_file+second_file+'_'+ ren3)
           os.rename(cwd + '/' + ren4,cwd + '/' + first_file+second_file+'_'+ ren4)
           os.rename(cwd + '/' + ren5,cwd + '/' + first_file+second_file+'_'+ ren5)
           os.makedirs('./STAR_results/Log', exist_ok = True)
           for ext in ["out", "tab"]:
            for file in glob.glob('./*.{}'.format(ext)):
                shutil.move(os.path.join (cwd,file), os.path.join (cwd + '/STAR_results/Log/', file)) #This way files will be overwrite if already exist
           del drlist[0:1]


   
    # end of general loop
        del drlist1[0:2]

#Convert files from .BAM to .BAI and save
print("\nIt's time to index the bam files")
lista = glob.glob('*.bam')
drlist = sorted (lista)
print (drlist)
while len(drlist) > 0:
    str1 =  drlist[0]
    day = time.strftime("%y/%m/%d")
    hour = time.strftime("%H:%M:%S")
    str1 = drlist[0]
    print('\nThe file \033[0;93;48m'+  str1  + '\033[0;0m starts at \033[0;94;48m' + hour + '\033[0;0m on \033[0;94;48m' + day +'\033[0;0m')
    p = subprocess.Popen('samtools index ' + str1,
        shell=True, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE)
    out, err = p.communicate()
    print('\nFile processed')
    del drlist[0:1]

#Move .bam and .bai files to STAR_results
for ext in ["bam","bai"]:
    for file in glob.glob('./*.{}'.format(ext)):
        shutil.move(os.path.join (cwd,file), os.path.join (cwd + '/STAR_results/', file)) #This way files will be overwriten if already exist


#Feature Counts with all .bam files generated
day = time.strftime("%y/%m/%d")
hour = time.strftime("%H:%M:%S")
print('Feature counts starts at \033[0;94;48m' + hour + '\033[0;0m on \033[0;94;48m' + day + '\033[0;0m')
p = subprocess.Popen('featureCounts -a '+ repr(feature)+' -o feautreCounts_results.csv ./STAR_results/*.bam -T '+ repr(threads) + '--verbose',
    shell=True, 
    stdout=subprocess.PIPE, 
    stderr=subprocess.PIPE)
out, err = p.communicate()


#MultiQC for all the steps where I have performed fastQC
print("\nNow is time to complete the results\nmultiQC for raw data")

file= open ('multiQC.sh','w')
file.write ('#!/bin/bash\nmultiqc -f -d -o ./QC/raw ./QC/raw' )
file.close()

for ext in ['.sh']:
        p = subprocess.Popen ('bash ./multiQC.sh',
            shell=True,
            executable='/bin/bash',
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        os.remove('multiQC.sh')

print("\nmultiQC for SortmeRNA data")
file= open ('multiQC.sh','w')
file.write ('#!/bin/bash\nmultiqc -f -d  -o ./QC/SortmeRNA ./QC/SortmeRNA' )
file.close()

for ext in ['.sh']:
        p = subprocess.Popen ('bash ./multiQC.sh',
            shell=True,
            executable='/bin/bash',
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        os.remove('multiQC.sh')

print("\nmultiQC for Trimmomatic data")
file= open ('multiQC.sh','w')
file.write ('#!/bin/bash\nmultiqc -f -d  -o ./QC/Trimmomatic ./QC/Trimmomatic' )
file.close()

for ext in ['.sh']:
        p = subprocess.Popen ('bash ./multiQC.sh',
            shell=True,
            executable='/bin/bash',
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        os.remove('multiQC.sh')
# Move the original files to the main directory
source = './original_files/'
dest1 = './'
files = os.listdir(source)
for f in files:
        shutil.move(source+f, dest1)

os.removedirs('./original_files')


day = time.strftime("%y/%m/%d")
hour = time.strftime("%H:%M:%S")
print("Cool! you reached the end at \033[0;96;48m" + hour + " on " + day + '\033[0;0m')

print ('\033[0;96;48m\nPreprocess of data finished successfully\n\nPlease review the QC data for multiple files in the following directories:')
print('./QC/raw/multiqc_report.html -----------> multiQC for raw data')
print('./QC/SortmeRNA/multiqc_report.html -----> multiQC for reads from raw data remaining after SortmeRNA ')
print('./QC/Trimmomatic/multiqc_report.html ---> multiQC for SortmeRNA data remaining after Trimmomatic')
print('./SortmeRNA ----------------------------> log files of processes and fq files with rejected reads to check')
print('./Trimmomatic --------------------------> Summary of reads for each group after processed with Trimmomatic')
print('./STAR_results/ ------------------------> Include all .bam and .bai files to review with')
print('./STAR_results/log ---------------------> Log progress for each condition processed and annotated')
print('\nfeatureCounts_results.txt.summary is interesting to check in order to see how many information was finaly incuded.')
print('featureCounts_results.csv is the file to be used in R\033[0;0m')

#review Rscript to perform properly the data analysis afterwards
#subprocess.call (["/usr/bin/Rscript", "--vanilla", "./Script_Differential_Expression_1.R"])


