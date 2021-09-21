#!/bin/bash
#SBATCH --job-name=homework_1                           # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=2:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/vjp98982/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=vjp98982@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="home/vjp98982/homework_1"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi


# curl the gff file
curl -s ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz | gunzip -c > ${OUTDIR}/ecoli_MG1655.gff


# cut out the third column (specifying the nucleotide type), pipe to grep
# each line containing "CDS" is piped to wc
# wc counts the number of lines output from grep and stores the total in results.txt
cut -f3 ${OUTDIR}/ecoli_MG1655.gff | grep "CDS" | wc -l > ${OUTDIR}/results.txt
