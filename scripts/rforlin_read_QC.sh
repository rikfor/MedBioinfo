#!/bin/bash
echo "script start: download and initial sequencing read quality control"
date

echo "get run_accessions for correct samples"
sqlite3 -batch -noheader -csv /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "select run_accession from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='rforlin'" > rforlin_run_accessions.txt

echo "loading sra-tool"
module load sra-tools/2.11.0

echo "making directory to store fastq-files"
mkdir ../data/sra_fastq/

echo "download sample fastq-files"
cat rforlin_run_accessions.txt | srun --cpus-per-task=1 --time=00:30:00 xargs fastq-dump --readids --gzip \
--outdir ../data/sra_fastq/ --disable-multithreading --split-e

echo "load seqkit and check stats for each fastq-file, and compare it to the previous read counts"
module load seqkit

seqkit stats ../data/sra_fastq/* --threads 4 > ../analyses/seqkit_stats.txt

sqlite3 -batch -csv /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "select run_accession, read_count, base_count from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='rforlin';" > ../analyses/author_fastq_readstats.txt

echo "check if duplicates in samples (we want to work with deduplicated files)"
seqkit rmdup -s -i ../data/sra_fastq/*.fastq.gz | wc -l > ../analyses/duplicated_file_counts.txt

echo "check if reads have already been trimmed"

#Check for read 1:
seqkit grep -r -p "^AGATCGGAAGAGCACACGTCTGAACTCCAGTCA|AGATCGGAAGAGCACACGTCTGAACTCCAGTCA$" ../data/sra_fastq/*_1.fastq.gz | wc -l > adaptors_matching_R1.txt

#Check for read 2:
seqkit grep -r -p "^AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT$" ../data/sra_fastq/*_2.fastq.gz | wc -l > adaptors_matching_R2.txt

echo "--------------------------------------------------------"

echo "Running FastQC"
mkdir ../analyses/fastqc
module load fastqc

srun --cpus-per-task=4 --time=00:30:00 xargs -I{} -a rforlin_run_accessions.txt fastqc --outdir ../analyses/fastqc/ \
--threads 2 --noextract ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 


echo "Looking at the .html-files from FastQC - I think the reads have been trimmed from adaptor sequences since we get not % of Adaptor sequences ('Adaptor content')"
echo "Per base quality score is also really good - no low quality scores here."
echo "-------------------------------------------------------"
echo "Merging paired end reads"

module load flash2
mkdir ../data/merged_pairs

#srun --cpus-per-task=2 flash2 --threads=2 -z --output-directory=../merged_pairs/ --output-prefix=ERR6913149.flash \
#../data/sra_fastq/ERR6913149_1.fastq.gz ../data/sra_fastq/ERR6913149_2.fastq.gz 2>&1 | tee -a rforlin_flash2.log

echo "Proportion of reads that was merged: ~80% of my reads were merged successfully"
echo "Using seqkit stat: the range goes from 35-292 read-lengths"
echo "Checking out the histogram files, we see that the majority of reads are between 100-150 base pairs, then there is a steep decline"

echo "Now running flash on all sequencing runs"
srun --cpus-per-task=2 --time=00:30:00 xargs -a rforlin_run_accessions.txt -n 1 -I{} flash2 --threads=2 -z --output-directory=../data/merged_pairs/ \
--output-prefix={}.flash ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a rforlin_flash2.log 


#Comparing base pairs in unmerged vs merged reads - we seem to have losr some information (~15-20% av bases)
#seqkit stats -a ERR69131*.fastq.gz
#stats -a  ERR69131*.flash.extendedFrags.fastq.gz


echo "Use read mapping to check for PhiX contamination"
mkdir ../data/reference_seqs
efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna

module load bowtie2
mkdir ../data/bowtie2_DBs
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_DBs/PhiX_bowtie2_DB

srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/PhiX_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz  -S ../analyses/bowtie/rforlin_merged2PhiX.sam \
--threads 8 --no-unal 2>&1 | tee ../analyses/bowtie/rforlin_bowtie_merged2PhiX.log

echo "No hits observed!"

echo "Now, we create a DB and align with SARS-COV-2 instead"
efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fna
srun bowtie2-build -f ../data/reference_seqs/SC2_NC_045512.fna ../data/bowtie2_DBs/SC2_bowtie2_DB

srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/SC2_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz  -S ../analyses/bowtie/rforlin_mergedSC2.sam \
--threads 8 --no-unal 2>&1 | tee ../analyses/bowtie/rforlin_bowtie_mergedSC2.log

echo "0.01% alignment rate"

echo "--------------------------------------------"
echo "multiQC-time!"
module load multiqc

srun multiqc --force --title "rforlin sample sub-set" ../data/merged_pairs/ ../analyses/fastqc/ ./rforlin_flash2.log ../analyses/bowtie/

echo "ahhh enjoy that bird's eye view now"

date
echo "script end."
