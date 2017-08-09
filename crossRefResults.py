'''
Created on Jun 19, 2017

@author: bensen
'''
from collections import defaultdict
import re
phenotype_dict = {}
xhmm_dict = {}
d = defaultdict(list)
xhmm_file = 'XHMM_PARTC_76.txt'
phenotype_file = 'PARTC_Phenotypes.txt'

outfilename = 'Crossed_Results.txt'
outfile = open(outfilename,'w')

phenotype_lines = open(phenotype_file).read().splitlines()
for line in phenotype_lines:
    fields = re.split(r'\t+', line)
    phenotype_dict[fields[0]] = fields[1]
    
xhmm_lines = open(xhmm_file).read().splitlines()
for line in xhmm_lines:
    fields = re.split(r'\t+', line)
    xhmm_data = [fields[0], fields[1], fields[2], fields[3], fields[4], fields[7], fields[9]]
    print (xhmm_data)
    xhmm_key = fields[0]+fields[1]+fields[2]
    print (xhmm_key)
    d[xhmm_key].append(xhmm_data)

codex_lines = open(codex_file).read().splitlines()
for line in codex_lines:
    fields = re.split(r'\t+', line)
    codex_data = [fields[0], fields[1], fields[2], fields[3], fields[4], fields[7], fields[9]]
    codex_key = fields[0]+fields[1]+fields[2]
    codex_dict[codex_key] = codex_data    
    
first=dict(a=1,b=3,c=5)
second=dict(d=1,b=5,f=5)
third=dict(f=1,q=3,a=5)
import pandas as pd
df=pd.DataFrame.from_dict(dict(first=first,second=second,third=third),)
print(df)