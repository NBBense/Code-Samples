'''
Created on Jan 27, 2017

@author: bensen
'''
# Import PyVCF
import vcf

# Allow for cross platform compatibility
import os
os.path.join('usr', 'bin', 'spam')

# Declare dictionaries
g1k_dict = {}
esp_dict = {}
exac_dict = {}

# Accept user input for 1000 Genomes Project .vcf file
g1k_input_file = input("\nPlease type in the file path/name to the 1000 Genomes Project (1KG) .vcf file you would like to parse and press 'Enter': ")

# Code in development to incorporate tabix indexing
#g1k_chrNum = input("\nPlease type in the chromosome number of the 1000 Genomes Project .vcf file you would like to parse and press 'Enter': ")

# Accepts user input for Exome Sequencing Project .vcf file
esp_input_file = input("\nPlease type in the file path/name to the Exome Sequencing Project (ESP) .vcf file you would like to parse and press 'Enter': ")

# Code in development to incorporate tabix indexing
#esp_chrNum = input("\nPlease type in the chromosome number of the Exome Sequencing Project (ESP) .vcf file you would like to parse and press 'Enter': ")

# Accepts user input for Exome Aggregation Consortium .vcf file
exac_input_file = input("\nPlease type in the file path/name to the Exome Aggregation Consortium (ExAC) .vcf file you would like to parse and press 'Enter': ")
#exac_chrNum = input("\nPlease type in the chromosome number of the Exome Aggregation Consortium (ExAC) .vcf file you would like to parse and press 'Enter': ")

# Accepts user input to specify output file
output_file = input("\nPlease type in the file path/name of your desired output file and press 'Enter': ")
output = open(output_file, 'w')

# Creates a header for field referencing
output.write("##CHROM" + "\t" + "POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "1KG_AC" + "\t" + "1KG_AN" + "\t" + 
             "1KG_EAS_AF" + "\t" + "1KG_EUR_AF" + "\t" + "1KG_AFR_AF" + "\t" + "1KG_AMR_AF" + "\t" +
             "1KG_SAS_AF" + "\t" + "1KG_TOTAL_AF" + "\t" + "ESP_EA_AC" + "\t" + "ESP_AA_AC" + "\t" + "ESP_TOTAL_AC" + "\t" +
             "ESP_EA_AN" + "\t" + "ESP_AA_AN" + "\t" + "ESP_TOTAL_AN" + "\t" + "ESP_EA_AF" + "\t" + "ESP_AA_AF" + "\t" + 
             "ESP_TOTAL_AF" + "\t" + "ExAC_AC_AFR" + "\t" + "ExAC_AC_AMR" +  "\t" + "ExAC_AC_EAS" + "\t" + "ExAC_AC_FIN" + "\t" +
             "ExAC_AC_NFE" + "\t" + "ExAC_AC_OTH" + "\t" + "ExAC_AC_SAS" + "\t" + "ExAC_AC" + "\t" + "ExAC_AN_AFR" + "\t" +
             "ExAC_AN_AMR" + "\t" + "ExAC_AN_EAS" + "\t" + "ExAC_AN_FIN" + "\t" + "ExAC_AN_NFE" + "\t" + "ExAC_AN_OTH" + "\t" +
             "ExAC_AN_SAS" + "\t" + "ExAC_AN" + "\t" + "ExAC_AF_AFR" + "\t" + "ExAC_AF_AMR" + "\t" + "ExAC_AF_EAS" + "\t" + 
             "ExAC_AF_FIN" + "\t" + "ExAC_AF_NFE" + "\t" + "ExAC_AF_OTH" + "\t" + "ExAC_AF_SAS" + "\t" + "ExAC_AF" + "\n")

# Opens an instance of PyVCF reader to parse the 1000 Genome file 
g1k_input_reader = vcf.Reader(open(g1k_input_file, 'r'))

# Parse each 1000 Genome variant record loaded to the reader from the file
for g1k_record in g1k_input_reader:
    
# Select only single allele variant records to parse
    if len(g1k_record.ALT) == 1:
        
# Recursive code in development for multiple allele variants
        #g1k_altMultiples = len(g1k_record.ALT)
        #g1k_counter = 0
        #while counter < length:
        
# Retrieve and format the desired fields from the variant record and populate the dictionary with the position as key and record as value
        g1k_data = [g1k_record.CHROM, g1k_record.POS, g1k_record.REF, g1k_record.ALT[0], g1k_record.INFO['AC'][0], 
                      g1k_record.INFO['AN'], "%.9f" % float(g1k_record.INFO['EAS_AF'][0]), "%.9f" % float(g1k_record.INFO['EUR_AF'][0]), "%.9f" % float(g1k_record.INFO['AFR_AF'][0]), 
                      "%.9f" % float(g1k_record.INFO['AMR_AF'][0]), "%.9f" % float(g1k_record.INFO['SAS_AF'][0]), "%.9f" % float(g1k_record.INFO['AF'][0])]
        g1k_dict[g1k_record.POS] = g1k_data
        
# Recursive code in development for multiple allele variants
        #g1k_multDict[g1k_counter] = g1k_dict

# Opens an instance of PyVCF reader to parse the Exome Sequencing Project file 
esp_input_reader = vcf.Reader(open(esp_input_file, 'r'))

# Parse each Exome Sequencing Project variant record loaded to the reader from the file
for esp_record in esp_input_reader:

# Select only single allele variant records to parse    
    if len(esp_record.ALT) == 1:

# Recursive code in development for multiple allele variants
        #esp_length = len(esp_record.INFO['EA_AC'])
        #esp_counter = 0
        #while esp_counter != esp_length: 

# Retrieve and format the desired fields from the variant record and populate the dictionary with the position as key and record as value
        esp_data = [esp_record.CHROM, esp_record.POS, esp_record.REF, esp_record.ALT[0], esp_record.INFO['EA_AC'][0], 
                    esp_record.INFO['AA_AC'][0], esp_record.INFO['TAC'][0], int(esp_record.INFO['EA_AC'][0]) + int(esp_record.INFO['EA_AC'][1]),
                    int(esp_record.INFO['AA_AC'][0]) + int(esp_record.INFO['AA_AC'][1]), int(esp_record.INFO['TAC'][0]) + int(esp_record.INFO['TAC'][1]), 
                    "%.9f" % (float(esp_record.INFO['MAF'][0])/100), "%.9f" % (float(esp_record.INFO['MAF'][1])/100), "%.9f" % (float(esp_record.INFO['MAF'][2])/100)]
        esp_dict[esp_record.POS] = esp_data
        
# Recursive code in development for multiple allele variants
        #esp_multDict[esp_counter] = esp_dict

# Code to incorporate ExAC files once format is corrected/adjusted for
exac_input_reader = vcf.Reader(open(exac_input_file, 'r'))
for exac_record in exac_input_reader:  
    if len(exac_record.ALT) == 1:
        if int(exac_record.INFO['AN_AFR'][0]) > 0:
            AF_AFR = "%.9f" % (float(exac_record.INFO['AC_AFR'][0])/float(exac_record.INFO['AN_AFR'][0]))
        else:
            AF_AFR = 0    
        if int(exac_record.INFO['AN_AMR'][0]) > 0: 
            AF_AMR = "%.9f" % (float(exac_record.INFO['AC_AMR'][0])/float(exac_record.INFO['AN_AMR'][0]))
        else:
            AF_AMR = 0    
        if int(exac_record.INFO['AN_EAS'][0]) > 0: 
            AF_EAS = "%.9f" % (float(exac_record.INFO['AC_EAS'][0])/float(exac_record.INFO['AN_EAS'][0]))
        else:
            AF_EAS = 0    
        if int(exac_record.INFO['AN_FIN'][0]) > 0: 
            AF_FIN = "%.9f" % (float(exac_record.INFO['AC_FIN'][0])/float(exac_record.INFO['AN_FIN'][0]))
        else:
            AF_FIN = 0    
        if int(exac_record.INFO['AN_NFE'][0]) > 0: 
            AF_NFE = "%.9f" % (float(exac_record.INFO['AC_NFE'][0])/float(exac_record.INFO['AN_NFE'][0]))
        else:
            AF_NFE = 0    
        if int(exac_record.INFO['AN_OTH'][0]) > 0: 
            AF_OTH = "%.9f" % (float(exac_record.INFO['AC_OTH'][0])/float(exac_record.INFO['AN_OTH'][0]))
        else:
            AF_OTH = 0    
        if int(exac_record.INFO['AN_SAS'][0]) > 0: 
            AF_SAS = "%.9f" % (float(exac_record.INFO['AC_SAS'][0])/float(exac_record.INFO['AN_SAS'][0]))
        else:
            AF_SAS = 0
        
        exac_data = [exac_record.CHROM, exac_record.POS, exac_record.REF, exac_record.ALT[0], exac_record.INFO['AC_AFR'][0],
                exac_record.INFO['AC_AMR'][0], exac_record.INFO['AC_EAS'][0], exac_record.INFO['AC_FIN'][0], exac_record.INFO['AC_NFE'][0],
                exac_record.INFO['AC_OTH'][0], exac_record.INFO['AC_SAS'][0], exac_record.INFO['AC'][0], exac_record.INFO['AN_AFR'][0],
                exac_record.INFO['AN_AMR'][0], exac_record.INFO['AN_EAS'][0], exac_record.INFO['AN_FIN'][0], exac_record.INFO['AN_NFE'][0],
                exac_record.INFO['AN_OTH'][0], exac_record.INFO['AN_SAS'][0], exac_record.INFO['AN'][0], 
                AF_AFR, AF_AMR, AF_EAS, AF_FIN, AF_NFE, AF_OTH, AF_SAS, "%.9f" % float(exac_record.INFO['AF'][0]) ]
    exac_dict[exac_record.POS] = exac_data
    
# Loop through the variants in the 1000 Genomes Project
for key in g1k_dict:
    
# Print out the desired fields if the variant is a match between the 1000 Genomes Project, Exome Sequencing Project, and Exome Aggregation Consortium
    if key in esp_dict and key in exac_dict:
        output.write(str(g1k_dict[key][0]) + "\t" + str(g1k_dict[key][1]) + "\t" + str(g1k_dict[key][2]) + "\t" + str(g1k_dict[key][3]) + 
            "\t" + str(g1k_dict[key][4]) + "\t" + str(g1k_dict[key][5]) + "\t" + str(g1k_dict[key][6]) + "\t" + str(g1k_dict[key][7]) + 
            "\t" + str(g1k_dict[key][8]) + "\t" + str(g1k_dict[key][9]) + "\t" + str(g1k_dict[key][10]) + "\t" + str(g1k_dict[key][11]) +  
            "\t" + str(esp_dict[key][4]) + "\t" + str(esp_dict[key][5]) + "\t" + str(esp_dict[key][6]) + "\t" + str(esp_dict[key][7]) + 
            "\t" + str(esp_dict[key][8]) + "\t" + str(esp_dict[key][9]) + "\t" + str(esp_dict[key][10]) + "\t" + str(esp_dict[key][11]) + 
            "\t" + str(esp_dict[key][12]) + "\t" + str(exac_dict[key][4]) + "\t" + str(exac_dict[key][5]) + "\t" + str(exac_dict[key][6]) +
            "\t" + str(exac_dict[key][7]) + "\t" + str(exac_dict[key][8]) + "\t" + str(exac_dict[key][9]) + "\t" + str(exac_dict[key][10]) + 
            "\t" + str(exac_dict[key][11]) + "\t" + str(exac_dict[key][12]) + "\t" + str(exac_dict[key][13]) + "\t" + str(exac_dict[key][14]) + 
            "\t" + str(exac_dict[key][15]) + "\t" + str(exac_dict[key][16]) + "\t" + str(exac_dict[key][17]) + "\t" + str(exac_dict[key][18]) + 
            "\t" + str(exac_dict[key][19]) + "\t" + str(exac_dict[key][20]) + "\t" + str(exac_dict[key][21]) + 
            "\t" + str(exac_dict[key][22]) + "\t" + str(exac_dict[key][23]) + "\t" + str(exac_dict[key][24]) + 
            "\t" + str(exac_dict[key][25]) + "\t" + str(exac_dict[key][26]) + "\t" + str(exac_dict[key][27]) + "\n")
# Print out the desired fields if the variant is a match between the 1000 Genomes Project and Exome Aggregation Consortium
    elif key not in esp_dict and key in exac_dict:
        output.write(str(g1k_dict[key][0]) + "\t" + str(g1k_dict[key][1]) + "\t" + str(g1k_dict[key][2]) + "\t" + str(g1k_dict[key][3]) + 
            "\t" + str(g1k_dict[key][4]) + "\t" + str(g1k_dict[key][5]) + "\t" + str(g1k_dict[key][6]) + "\t" + str(g1k_dict[key][7]) + 
            "\t" + str(g1k_dict[key][8]) + "\t" + str(g1k_dict[key][9]) + "\t" + str(g1k_dict[key][10]) + "\t" + str(g1k_dict[key][11]) + 
            "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
            "\t" + str(exac_dict[key][4]) + "\t" + str(exac_dict[key][5]) + "\t" + str(exac_dict[key][6]) + "\t" + str(exac_dict[key][7]) + 
            "\t" + str(exac_dict[key][8]) + "\t" + str(exac_dict[key][9]) + "\t" + str(exac_dict[key][10]) + "\t" + str(exac_dict[key][11]) + 
            "\t" + str(exac_dict[key][12]) + "\t" + str(exac_dict[key][13]) + "\t" + str(exac_dict[key][14]) + "\t" + str(exac_dict[key][15]) + 
            "\t" + str(exac_dict[key][16]) + "\t" + str(exac_dict[key][17]) + "\t" + str(exac_dict[key][18]) + "\t" + str(exac_dict[key][19]) + 
            "\t" + str(exac_dict[key][20]) + "\t" + str(exac_dict[key][21]) + "\t" + str(exac_dict[key][22]) + "\t" + str(exac_dict[key][23]) + 
            "\t" + str(exac_dict[key][24]) + "\t" + str(exac_dict[key][25]) + "\t" + str(exac_dict[key][26]) + "\t" + str(exac_dict[key][27]) + "\n")    

# Recursive code in development for multiple allele variants
#    if len(g1k_dict[key][3]) = 1: 
  
# Print out the desired fields if the variant is a match between the 1000 Genomes Project and Exome Sequencing Project but not in the Exome Aggregation Consortium  
    elif key in esp_dict and key not in exac_dict:
        output.write(str(g1k_dict[key][0]) + "\t" + str(g1k_dict[key][1]) + "\t" + str(g1k_dict[key][2]) + "\t" + str(g1k_dict[key][3]) + 
                        "\t" + str(g1k_dict[key][4]) + "\t" + str(g1k_dict[key][5]) + "\t" + str(g1k_dict[key][6]) + "\t" + str(g1k_dict[key][7]) + 
                        "\t" + str(g1k_dict[key][8]) + "\t" + str(g1k_dict[key][9]) + "\t" + str(g1k_dict[key][10]) + "\t" + str(g1k_dict[key][11]) + 
                        "\t" + str(esp_dict[key][4]) + "\t" + str(esp_dict[key][5]) + "\t" + str(esp_dict[key][6]) + "\t" + str(esp_dict[key][7]) + 
                        "\t" + str(esp_dict[key][8]) + "\t" + str(esp_dict[key][9]) + "\t" + str(esp_dict[key][10]) + "\t" + str(esp_dict[key][11]) + 
                        "\t" + str(esp_dict[key][12]) + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
                        "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
                        "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\n")
    
# Print out the desired fields if the variant is only found in the 1000 Genomes Project   
    else:
        output.write(str(g1k_dict[key][0]) + "\t" + str(g1k_dict[key][1]) + "\t" + str(g1k_dict[key][2]) + "\t" + str(g1k_dict[key][3]) + 
                        "\t" + str(g1k_dict[key][4]) + "\t" + str(g1k_dict[key][5]) + "\t" + str(g1k_dict[key][6]) + "\t" + str(g1k_dict[key][7]) + 
                        "\t" + str(g1k_dict[key][8]) + "\t" + str(g1k_dict[key][9]) + "\t" + str(g1k_dict[key][10]) + "\t" + str(g1k_dict[key][11]) + 
                        "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
                        "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
                        "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
                        "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "\n")

# Loop through the variants in the Exome Sequencing Project
for key in esp_dict:
    
# Print out the desired fields if the variant is not in the 1000 Genomes project but a match between the Exome Sequencing Project and Exome Aggregation Consortium
    if key not in g1k_dict and key in exac_dict:
        output.write(str(esp_dict[key][0]) + "\t" + str(esp_dict[key][1]) + "\t" + str(esp_dict[key][2]) + "\t" + str(esp_dict[key][3]) + 
            "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
            "\t" + str(esp_dict[key][4]) + "\t" + str(esp_dict[key][5]) + "\t" + str(esp_dict[key][6]) + "\t" + str(esp_dict[key][7]) + 
            "\t" + str(esp_dict[key][8]) + "\t" + str(esp_dict[key][9]) + "\t" + str(esp_dict[key][10]) + "\t" + str(esp_dict[key][11]) + 
            "\t" + str(esp_dict[key][12]) + "\t" + str(exac_dict[key][4]) + "\t" + str(exac_dict[key][5]) + "\t" + str(exac_dict[key][6]) +
            "\t" + str(exac_dict[key][7]) + "\t" + str(exac_dict[key][8]) + "\t" + str(exac_dict[key][9]) + "\t" + str(exac_dict[key][10]) + 
            "\t" + str(exac_dict[key][11]) + "\t" + str(exac_dict[key][12]) + "\t" + str(exac_dict[key][13]) + "\t" + str(exac_dict[key][14]) + 
            "\t" + str(exac_dict[key][15]) + "\t" + str(exac_dict[key][16]) + "\t" + str(exac_dict[key][17]) + "\t" + str(exac_dict[key][18]) + 
            "\t" + str(exac_dict[key][19]) + "\t" + str(exac_dict[key][20]) + "\t" + str(exac_dict[key][21]) + "\t" + str(exac_dict[key][22]) + 
            "\t" + str(exac_dict[key][23]) + "\t" + str(exac_dict[key][24]) + "\t" + str(exac_dict[key][25]) + "\t" + str(exac_dict[key][26]) + 
            "\t" + str(exac_dict[key][27]) + "\n")

# Print the variants that are in the Exome Sequencing Project but not matches to the 1000 Genomes Project or Exome Aggregation Consortium    
    elif key not in g1k_dict and key not in exac_dict:
        output.write(str(esp_dict[key][0]) + "\t" + str(esp_dict[key][1]) + "\t" + str(esp_dict[key][2]) + "\t" + str(esp_dict[key][3]) + 
            "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
            "\t" + str(esp_dict[key][4]) + "\t" + str(esp_dict[key][5]) + "\t" + str(esp_dict[key][6]) + "\t" + str(esp_dict[key][7]) + 
            "\t" + str(esp_dict[key][8]) + "\t" + str(esp_dict[key][9]) + "\t" + str(esp_dict[key][10]) + "\t" + str(esp_dict[key][11]) + 
            "\t" + str(esp_dict[key][12]) + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
            "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
            "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "\n")

# Print out the desired fields if the variant is only found in the Exome Aggregation Consortium
for key in exac_dict:
    if key not in g1k_dict and key not in esp_dict:
        output.write(str(exac_dict[key][0]) + "\t" + str(exac_dict[key][1]) + "\t" + str(exac_dict[key][2]) + "\t" + str(exac_dict[key][3]) + 
            "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
            "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + 
            "\t"+ str(exac_dict[key][4]) + "\t" + str(exac_dict[key][5]) + "\t" + str(exac_dict[key][6]) + "\t" + str(exac_dict[key][7]) + 
            "\t" + str(exac_dict[key][8]) + "\t" + str(exac_dict[key][9]) + "\t" + str(exac_dict[key][10]) + "\t" + str(exac_dict[key][11]) + 
            "\t" + str(exac_dict[key][12]) + "\t" + str(exac_dict[key][13]) + "\t" + str(exac_dict[key][14]) + "\t" + str(exac_dict[key][15]) + 
            "\t" + str(exac_dict[key][16]) + "\t" + str(exac_dict[key][17]) + "\t" + str(exac_dict[key][18]) + "\t" + str(exac_dict[key][19]) + 
            "\t" + str(exac_dict[key][20]) + "\t" + str(exac_dict[key][21]) + "\t" + str(exac_dict[key][22]) + "\t" + str(exac_dict[key][23]) + 
            "\t" + str(exac_dict[key][24]) + "\t" + str(exac_dict[key][25]) + "\t" + str(exac_dict[key][26]) + "\t" + str(exac_dict[key][27]) + "\n")

# Maybe implement default dict for nested dictionaries? 
#dict.setdefault(key, default=None)    

    
    