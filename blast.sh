#!/bin/bash
#SBATCH --job-name=testBLAST		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=2:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/cbergman/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=cbergman@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/work/gene8940/cbergman/"                       # replace cbergman in the following line with your myid
DATADIR="/work/gene8940/instructor_data"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load blast module
module load BLAST+/2.9.0-gompi-2019b

# run blast against local copy of NCBI nucleotide database
blastn -num_threads 4 -query $DATADIR/sample.fasta -db /db/ncbiblast/nt/06042020/nt -out $OUTDIR/sample.fa.blastn.${SLURM_JOB_ID}.tsv -outfmt 6 -max_target_seqs=2
