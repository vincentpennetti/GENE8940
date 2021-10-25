#!/bin/bash
#SBATCH --job-name=homework_4                           # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=6:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/vjp98982/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=vjp98982@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/work/gene8940/vjp98982/homework_4"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load necessary modules for file processing
module load SRA-Toolkit/2.9.6-1-centos_linux64
module load BWA/0.7.17-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load BCFtools/1.10.2-GCC-8.3.0

# curl the refseq E coli MG1655 genome sequence
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > ${OUTDIR}/GCF_000005845.2_ASM584v2.fna

# dowload and extract paired end illumina reads from ncbi for e coli C600
prefetch -O ${OUTDIR} SRR8082143
fastq-dump --split-files --gzip ${OUTDIR}/SRR8082143.sra -O ${OUTDIR}

# construct a BWA index for the E. coli MG1655 refseq reference genome
bwa index ${OUTDIR}/GCF_000005845.2_ASM584v2.fna

# map reads to reference with bwa and output as sorted BAM file
bwa mem -t 6 ${OUTDIR}/GCF_000005845.2_ASM584v2.fna ${OUTDIR}/SRR8082143_1.fastq.gz ${OUTDIR}/SRR8082143_2.fastq.gz | samtools view - -O BAM | samtools sort - > ${OUTDIR}/SRR8082143_pipe.sorted.bam

# generate index for sorted .bam file using samtools index
samtools index ${OUTDIR}/SRR8082143_pipe.sorted.bam

#####################
# call variants
# i) reads with mapping quality >60
bcftools mpileup -Oz --threads 6 --min-MQ 60 -f ${OUTDIR}/GCF_000005845.2_ASM584v2.fna ${OUTDIR}/SRR8082143_pipe.sorted.bam > ${OUTDIR}/SRR8082143.sorted.mpileup.vcf.gz

# call the variants
bcftools call -Oz -mv --threads 6 --ploidy 1 ${OUTDIR}/SRR8082143.sorted.mpileup.vcf.gz > ${OUTDIR}/SRR8082143.sorted.mpileup.call.vcf.gz

# exclude calls with a quality score < 40 and a depth of <10 reads
bcftools filter -Oz -i 'QUAL<40 && DP<10' ${OUTDIR}/SRR8082143.sorted.mpileup.call.vcf.gz >  ${OUTDIR}/SRR8082143.sorted.mpileup.call.filter.vcf.gz
#########################

# the above as a single line pipe, avoids wasting compute time on compression and uncompression
#bcftools mpileup -Oz --threads 6 --min-MQ 60 -f ${OUTDIR}/GCF_000005845.2_ASM584v2.fna ${OUTDIR}/SRR8082143_pipe.sorted.bam | \
#bcftools call -Oz -mv --threads 6 --ploidy 1  ${OUTDIR}/SRR8082143.sorted.mpileup.vcf.gz | \
#bcftools filter -Oz -e 'QUAL<40 || DP<10'  ${OUTDIR}/SRR8082143.sorted.mpileup.call.vcf.gz >  ${OUTDIR}/SRR8082143.sorted.mpileup.call.filter.vcf.gz

# create IGV readable index file
bcftools index ${OUTDIR}/SRR8082143.sorted.mpileup.call.filter.vcf.gz
