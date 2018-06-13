################################################################################
# Title: crossResultsCNV
# Author: Nicholas Bense
# Date: 6/26/17
# Description: Python script to cross reference and annotate raw results from the XHMM, Codex, and Clamms germline 
# CNV callers to detect matching events that are discovered in the same gene regions. This is intended to allow for comparison
# of the respective quality metric scores used in determining the likelihood of an event being true positive or false positive 
# as well as the correlation of phenotype data in a germ-line CNV calling pipeline.
################################################################################

# Import packages
from collections import defaultdict
import subprocess
import re

# Declare dictionaries
phenotype_dict = {}
xhmm_dict = defaultdict(list)
codex_dict = defaultdict(list)
clamms_dict = defaultdict(list)

# Create filehandles
xhmm_file = 'PartC_XHMM_Results.xcnv'
temp_xhmm_file = 'temp.' + xhmm_file
codex_file = 'PartC_CODEX_Results.txt'
temp_codex_file = 'temp.' + codex_file
clamms_file = 'PartC_CLAMMS_Results.txt'
temp_clamms_file = 'temp.' + clamms_file
master_bed_file = 'master.bed'
phenotype_file = 'Phenotypes_PartC.txt'
output_file = 'PartC_Crossed_Results_Final.txt'

# Open output file for writing
output = open(output_file,'w')

# Read in phenotype file
phenotype_lines = open(phenotype_file).read().splitlines()
for line in phenotype_lines:
    fields = re.split(r'\t+', line)
    # Populate dictionary with sample ID and phenotype
    phenotype_dict[fields[0]] = fields[1]
   
# Preprocess XHMM results file for desired fields and annotate with master bed file and bedtools intersect   
subprocess.call(["sed '1d' " + str(xhmm_file) + "| sed 's/^.\{15\}//g' | sed 's/\-N//g' | sed 's/\-/    /' | sed 's/\:/    /' | awk 'BEGIN {OFS=\"\t\"} {print $3, $4, $5, $1, $2, $12, $6}' | sed 's/^/chr/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g'> temp." + str(xhmm_file)], shell=True)
xhmm_formatted = subprocess.Popen("bedtools intersect -wa -wb -a " + str(master_bed_file) + ' -b ' + str(temp_xhmm_file) + " | uniq | awk 'BEGIN {OFS=\"\t\"} {print $1,$6,$7,$8,$4,$9,$10,$11}' | sort | uniq | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g'", shell=True, stdout=subprocess.PIPE)   
    
# Read in processed XHMM results file and populate dictionary with desired fields    
xhmm_lines = xhmm_formatted.stdout.readlines()
for line in xhmm_lines:
    fields = re.split(r'\t+', line)
    interval = fields[0] + ":" + fields[1] + "-" + fields[2]
    # Create dictionary key by concatenating SampleID, Gene, and Type of Event (DEL/DUP)
    xhmm_key = str(fields[3] + fields[4] + fields[5])
    xhmm_data = [fields[3], fields[4], fields[5], interval, fields[6], fields[7].rstrip(), phenotype_dict[fields[3]]]
    xhmm_dict[xhmm_key].append(xhmm_data)  

# Preprocess Codex results file for desired fields and annotate with master bed file and bedtools intersect
subprocess.call(["sed '1d' " + str(codex_file) + "| sed 's/\-N//g' | sed 's/dup/DUP/' | sed 's/del/DEL/' | awk '!($7 ~ \"-\")' | awk 'BEGIN {OFS=\"\t\"} {print $2, $5, $6, $1, $4, $11, $12, $7}' | sed 's/^/chr/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' > temp." + str(codex_file)], shell=True)
codex_formatted = subprocess.Popen("bedtools intersect -wa -wb -a " + str(master_bed_file) + ' -b ' + str(temp_codex_file) + " | uniq | awk 'BEGIN {OFS=\"\t\"} {print $1,$6,$7,$8,$4,$9,$10,$11,$12}' | sort | uniq | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g'", shell=True, stdout=subprocess.PIPE)

# Read in processed Codex results file and populate dictionary with desired fields  
codex_lines = codex_formatted.stdout.readlines()
for line in codex_lines:
    fields = re.split(r'\t+', line)
    interval = fields[0] + ":" + fields[1] + "-" + fields[2]
    # Create dictionary key by concatenating SampleID, Gene, and Type of Event (DEL/DUP)
    codex_key = str(fields[3] + fields[4] + fields[5])
    codex_data = [fields[3], fields[4], fields[5], interval, fields[6], fields[7], fields[8].rstrip(), phenotype_dict[fields[3]]]
    codex_dict[codex_key].append(codex_data)

# Preprocess Clamms results file for desired fields and annotate with master bed file and bedtools intersect
subprocess.call(["sed 's/\-N//g' " + str(clamms_file) + "| awk 'BEGIN {OFS=\"\t\"} {print $1, $2, $3, $5, $6, $9, $10}' | sed 's/^/chr/g' | sed 's/^chrX/chr23/g' | sed 's/^chrY/chr24/g' > temp." + str(clamms_file)], shell=True)
clamms_formatted = subprocess.Popen("bedtools intersect -wa -wb -a " + str(master_bed_file) + ' -b ' + str(temp_clamms_file) + " | uniq | awk 'BEGIN {OFS=\"\t\"} {print $1,$6,$7,$8,$4,$9,$10,$11}' | sort | uniq | sed 's/^chr23/chrX/g' | sed 's/^chr24/chrY/g'", shell=True, stdout=subprocess.PIPE)

# Create dictionary key by concatenating SampleID, Gene, and Type of Event (DEL/DUP)
clamms_lines = clamms_formatted.stdout.readlines()
for line in clamms_lines:
    fields = re.split(r'\t+', line)
    interval = fields[0] + ":" + fields[1] + "-" + fields[2]
    length_kb = (float(fields[2]) - float(fields[1]))/1000
    # Create dictionary key by concatenating SampleID, Gene, and Type of Event (DEL/DUP)
    clamms_key = str(fields[3] + fields[4] + fields[5])
    clamms_data = [fields[3], fields[4], fields[5], interval, fields[6], fields[7].rstrip(), length_kb, phenotype_dict[fields[3]]]
    clamms_dict[clamms_key].append(clamms_data)  

# Cross reference the dictionaries and print the results
for key in xhmm_dict: 
    if key in codex_dict and key in clamms_dict:
        output.write(str(xhmm_dict[key][0][0])+ "\t" + str(xhmm_dict[key][0][6]) + "\t" + str(xhmm_dict[key][0][1]) + "\t" + str(xhmm_dict[key][0][2]) + 
                    "\t" + str(xhmm_dict[key][0][3]) + "\t" + str(xhmm_dict[key][0][4])+ "\t" + str(xhmm_dict[key][0][5]) + "\t" + str(codex_dict[key][0][3]) +
                    "\t" + str(codex_dict[key][0][4]) + "\t" + str(codex_dict[key][0][5])+ "\t" + str(codex_dict[key][0][6]) + "\t" + str(clamms_dict[key][0][3]) +
                    "\t" + str(clamms_dict[key][0][4]) + "\t" + str(clamms_dict[key][0][5])+ "\t" + str(clamms_dict[key][0][6]) +"\n")
    elif key in codex_dict and key not in clamms_dict:
        output.write(str(xhmm_dict[key][0][0])+ "\t" + str(xhmm_dict[key][0][6]) + "\t" + str(xhmm_dict[key][0][1]) + "\t" + str(xhmm_dict[key][0][2]) + 
                    "\t" + str(xhmm_dict[key][0][3]) + "\t" + str(xhmm_dict[key][0][4]) + "\t" + str(xhmm_dict[key][0][5]) + "\t" + str(codex_dict[key][0][3]) +
                    "\t" + str(codex_dict[key][0][4]) + "\t" + str(codex_dict[key][0][5]) + "\t" + str(codex_dict[key][0][6]) + "\t" + "NA" +
                    "\t" + "0" + "\t" + "0" + "\t" + "NA" + "\n")
    elif key not in codex_dict and key in clamms_dict:
        output.write(str(xhmm_dict[key][0][0])+ "\t" + str(xhmm_dict[key][0][6]) + "\t" + str(xhmm_dict[key][0][1]) + "\t" + str(xhmm_dict[key][0][2]) + 
                    "\t" + str(xhmm_dict[key][0][3]) + "\t" + str(xhmm_dict[key][0][4]) + "\t" + str(xhmm_dict[key][0][5]) + "\t" + "NA" +
                    "\t" + "0" + "\t" + "0" + "\t" + "NA" + "\t" + str(clamms_dict[key][0][3]) + "\t" + str(clamms_dict[key][0][4]) + "\t" + str(clamms_dict[key][0][5])+ 
                    "\t" + str(clamms_dict[key][0][6]) +"\n")
    else:
        output.write(str(xhmm_dict[key][0][0])+ "\t" + str(xhmm_dict[key][0][6]) + "\t" + str(xhmm_dict[key][0][1]) + "\t" + str(xhmm_dict[key][0][2]) + 
                    "\t" + str(xhmm_dict[key][0][3]) + "\t" + str(xhmm_dict[key][0][4]) + "\t" + str(xhmm_dict[key][0][5]) + "\t" + "NA" +
                    "\t" + "0" + "\t" + "0" + "\t" + "NA" + "\t" + "NA" +"\t" + "0" + "\t" + "0" + "\t" + "NA" +"\n")
for key in codex_dict: 
    if key not in xhmm_dict and key in clamms_dict:
        output.write(str(codex_dict[key][0][0])+ "\t" + str(codex_dict[key][0][7]) + "\t" + str(codex_dict[key][0][1]) + "\t" + str(codex_dict[key][0][2]) + 
                    "\t" + "NA" + "\t" + "0" + "\t" + "NA" + "\t" + str(codex_dict[key][0][3]) + "\t" + str(codex_dict[key][0][4]) + "\t" + str(codex_dict[key][0][5]) + 
                    "\t" + str(codex_dict[key][0][6]) + "\t" + str(clamms_dict[key][0][3]) + "\t" + str(clamms_dict[key][0][4]) + "\t" + str(clamms_dict[key][0][5])+ 
                    "\t" + str(clamms_dict[key][0][6]) +"\n")
    if key not in xhmm_dict and key not in clamms_dict:
        output.write(str(codex_dict[key][0][0])+ "\t" + str(codex_dict[key][0][7]) + "\t" + str(codex_dict[key][0][1]) + "\t" + str(codex_dict[key][0][2]) + 
                    "\t" + "NA" + "\t" + "0" + "\t" + "NA" + "\t" + str(codex_dict[key][0][3]) + "\t" + str(codex_dict[key][0][4]) + "\t" + str(codex_dict[key][0][5]) + 
                    "\t" + str(codex_dict[key][0][6]) + "\t" + "NA" + "\t" + "0" + "\t" + "0" + "\t" + "NA" +"\n")
for key in clamms_dict:
    if key not in xhmm_dict and key not in codex_dict:
        output.write(str(clamms_dict[key][0][0])+ "\t" + str(clamms_dict[key][0][7]) + "\t" + str(clamms_dict[key][0][1]) + "\t" + str(clamms_dict[key][0][2]) + 
                    "\t" + "NA" + "\t" + "0" + "\t" + "NA" + "\t" + "NA" + "\t" + "0" + "\t" + "0" + "\t" + "NA" + "\t" + str(clamms_dict[key][0][3]) + 
                    "\t" + str(clamms_dict[key][0][4]) + "\t" + str(clamms_dict[key][0][5]) + "\t" + str(clamms_dict[key][0][6]) +"\n")

# Delete temporary files
subprocess.call(["rm " + temp_xhmm_file], shell=True)
subprocess.call(["rm " + temp_codex_file], shell=True)
subprocess.call(["rm " + temp_clamms_file], shell=True)   
