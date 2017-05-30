#!/bin/sh

# AutoTrimAlign
# 2/21/17
# This is a simple shell script designed to be run in conjunction with a cronjob for automatic processing 
# of fastq files in regular intervals (i.e. daily at midnight)

# Create variable for list of samples
SAMPLE_LIST=~/filepath/sample_list.txt
# Create variable with path of directory to be monitored
MONITOR_DIR=~/filepath/
# Create variable for trimGalore output directory 
TRIMMED_FASTQS=~/filepath/
# Create variable for .sam file output directory
SAM_FILES=~/filepath/
# Create variable for FASTQC output directory
FASTQC_REPORTS=~/filepath/
# Create variable for GATK  .jar file location
GATK=~/filepath/GenomeAnalysisTK.jar
# Create variable for FASTQC binary file location
FASTQC=~/filepath/fastqc
# Create variable for TrimGalore binary file location
TRIM_GALORE=/filepath/trim_galore
# Create variable for BWA aligner binary file location
BWA=/filepath/bwa
# Create variable for Reference Genome file location
REF_FASTA=~/filepath/hs37d5.fa

# Read list of old samples to variable
old_samples=$(cat ${SAMPLE_LIST})
# Create temporary file for list of old samples
echo "$old_samples" > old_samples.txt
# Sort old samples and create temporary file for sorted old samples
sort old_samples.txt > old_samples_sorted.txt
# Remove temporary file for old samples
rm old_samples.txt
# Create variable by listing all files in monitored directory
current_samples=$(ls -1 ${MONITOR_DIR})
# Create temporary file from list of current samples in monitored directory
echo "$current_samples" > current_samples.txt
# Sort list of current samples and create temporary file of sorted current samples
sort current_samples.txt > current_samples_sorted.txt
# Remove temporary file for current samples
rm current_samples.txt
# Create temporary file for new samples by comparing temporary files for sorted old samples versus sorted current samples
comm -13 old_samples_sorted.txt current_samples_sorted.txt > new_samples.txt
# Remove temporary file for sorted old samples
rm old_samples_sorted.txt
# Remove temporary file for current samples
rm current_samples_sorted.txt
# Create variable from reading temporary file for new samples
new_samples=$(cat 'new_samples.txt')
# Remove temporary file for new samples
rm new_samples.txt

# Loop through all new samples
for sample in $new_samples; do
	# Create a subdirectory for TrimGalore output for each sample
        mkdir ${TRIMMED_FASTQS}${sample}
	# Change directory to sample directory
        cd ${MONITOR_DIR}${sample}
	# Error check to ensure in proper directory
        if [ "$?" = "0" ]; then
		# Loop through file stems for only R1 of each set of paired end reads for sample to avoid redundancy
                for fastqc in $(ls *_R1_*.fastq.gz | rev | cut -c 10- | rev | uniq); do
                        # Save file stem from the read as variable
			trim_stem=` echo ${fastqc} | rev | cut -c 8- | rev`
			# Save read number for reads > 100 bp as variable
                        trim_num=` echo ${fastqc} | cut -d_ -f8`
			# Run FastQC on paired end reads and output to FastQC directory
                        ${FASTQC} -o ${FASTQC_REPORTS} ${trim_stem}_R1_${trim_num}.fastq.gz ${trim_stem}_R2_${trim_num}.fastq.gz
			# Run TrimGalore on paired end reads and output trimmed reads to directory
                        ${TRIM_GALORE} -o ${TRIMMED_FASTQS}${sample}/ --paired ${trim_stem}_R1_${trim_num}.fastq.gz ${trim_stem}_R2_${trim_num}.fastq.gz
                done
        
		# Change directory to TrimGalore output directory with trimmed reads for sample
                cd ${TRIMMED_FASTQS}${sample}
                # Error trapping to ensure in proper directory
		if [ "$?" = "0" ]; then         
			
			# Loop through trimmed fastqs for only R1 of each set of paired end reads for sample to avoid redundancy
                        for fastq in $(ls *_R1_*_val_1.fq.gz | rev | cut -c 13- | rev | uniq); do
                                # Save file stem from the read as variable
				stem=` echo ${fastq} | rev | cut -c 8- | rev`
				# Save read number for reads > 100 bp as variable
                                num=` echo ${fastq} | cut -d_ -f8`

                        	# Loop through different read numbers for reads > 100 bp
				for read in $num; do
					# Create variables for each mate pair of the read
                                	f1=${stem}_R1_${num}_val_1.fq.gz
                                	f2=${stem}_R2_${num}_val_2.fq.gz
                                	# Create a variable for the current directory
                                	cwd=`pwd`

                                	# Run BWA mem alignment specifying 4 threads and output .sam file to directory
					${BWA} mem -t 4 ${REF_FASTA} ${cwd}/${f1} ${cwd}/${f2} > ${SAM_FILES}${stem}.sam
                        	done
			done
                else
			# Error Trapping Exit Message
                        echo "Trimmed Sample Directory not found!" 1>&2
                        exit 1
                fi
        else
                # Error Trapping Exit Message
		echo "Raw FASTQ Sample Directory not found!" 1>&2
                exit 1
        fi 
done

# Repopulate sample list with newly processed files to avoid re-alignment
echo "$current_samples" > ${SAMPLE_LIST};

