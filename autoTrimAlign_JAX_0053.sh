#!/bin/sh

SAMPLE_LIST=~/BioinfoPrograms/AutoTrimAlign/samplelist_JAX_0053
MONITOR_DIR=~/BioinfoPrograms/AutoTrimAlign/JAX_0053/
TRIMMED_FASTQS=~/BioinfoPrograms/AutoTrimAlign/TrimmedFastQs/JAX_0053/
SAM_FILES=~/BioinfoPrograms/AutoTrimAlign/SAM_Files/JAX_0053/
FASTQC_REPORTS=~/BioinfoPrograms/autoTrimAlign/FastQC_Reports/JAX_0053/
GATK=~/BioinfoPrograms/My_Tools/gatk_tools/GenomeAnalysisTK.jar
FASTQC=~/BioinfoPrograms/My_Tools/FastQC/fastqc
TRIM_GALORE=/usr/bin/trim_galore
BWA=/usr/bin/bwa
SAMTOOLS=/usr/bin/samtools
XHMM=~/BioinfoPrograms/My_Tools/xhmm_tools/build/execs/xhmm
REF_FASTA=~/BioinfoPrograms/My_Tools/Reference_Genomes/hs37d5/hs37d5.fa

old_samples=$(cat ${SAMPLE_LIST})
echo "$old_samples" > old_samples_JAX_0053.txt
sort old_samples_JAX_0053.txt > old_samples_sorted_JAX_0053.txt
rm old_samples_JAX_0053.txt
current_samples=$(ls -1 ${MONITOR_DIR})
echo "$current_samples" > current_samples_JAX_0053.txt
sort current_samples_JAX_0053.txt > current_samples_sorted_JAX_0053.txt
rm current_samples_JAX_0053.txt
comm -13 old_samples_sorted_JAX_0053.txt current_samples_sorted_JAX_0053.txt > new_samples_JAX_0053.txt
rm old_samples_sorted_JAX_0053.txt
rm current_samples_sorted_JAX_0053.txt
new_samples=$(cat 'new_samples_JAX_0053.txt')
rm new_samples_JAX_0053.txt

for sample in $new_samples; do
        mkdir ${TRIMMED_FASTQS}${sample}
        cd ${MONITOR_DIR}${sample}
        if [ "$?" = "0" ]; then
                for fastqc in $(ls *_R1_*.fastq.gz | rev | cut -c 10- | rev | uniq); do
                        trim_stem=` echo ${fastqc} | rev | cut -c 8- | rev`
                        trim_num=` echo ${fastqc} | cut -d_ -f8`

                        ${FASTQC} -o ${FASTQC_REPORTS} ${trim_stem}_R1_${trim_num}.fastq.gz ${trim_stem}_R2_${trim_num}.fastq.gz
                        ${TRIM_GALORE} -o ${TRIMMED_FASTQS}${sample}/ --paired ${trim_stem}_R1_${trim_num}.fastq.gz ${trim_stem}_R2_${trim_num}.fastq.gz
                done
        
                cd ${TRIMMED_FASTQS}${sample}
                if [ "$?" = "0" ]; then         

                        for fastq in $(ls *_R1_*_val_1.fq.gz | rev | cut -c 13- | rev | uniq); do
                                stem=` echo ${fastq} | rev | cut -c 8- | rev`
                                num=` echo ${fastq} | cut -d_ -f8`

                        	for read in $num; do
                                	f1=${stem}_R1_${num}_val_1.fq.gz
                                	f2=${stem}_R2_${num}_val_2.fq.gz
                                	key=` zcat ${f1} | head -1 | cut -d: -f1-4`
                                	sm=` echo ${stem} | cut -d_ -f1`
                                	asy=` echo ${f1} | cut -d_ -f3`
                                	bar=` echo ${f1} | cut -d_ -f4`
                                	pl=` echo ${f1} | cut -d_ -f5`
                                	ln=` echo ${f1} | cut -d_ -f6`
                                	index=` zcat ${f1} | head -1 | cut -d: -f10`
                                	id=${key}:${index}
                                	cwd=`pwd`

                                	# Local BWA test run
					${BWA} mem -R "@RG\tID:${id}\tPL:Illumina\tLB:${sm}\tSM:${sm}\tCN:CGRL" -t 4 ${REF_FASTA} ${cwd}/${f1} ${cwd}/${f2} > ${SAM_FILES}${sm}_${asy}_${bar}_${pl}_${ln}_${num}.sam
                        	done
			done
                else
                        echo "Trimmed Sample Directory not found!" 1>&2
                        exit 1
                fi
        else
                echo "Raw FASTQ Sample Directory not found!" 1>&2
                exit 1
        fi 
done

echo "$current_samples" > ${SAMPLE_LIST};

