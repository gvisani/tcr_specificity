
import os
import pandas as pd
import json

df = pd.read_csv('ternary_templates_v2.tsv', sep='\t')

adict = {}
for i, row in df.iterrows():

    if row['mhc_class'] == 2:
        continue
    
    mhc_allele = row['mhc_allele']

    # put allele in NetMHCPan format
    if '*' in mhc_allele:
        mhc_allele = 'HLA-' + mhc_allele.replace('*', '')
    elif mhc_allele.startswith('H2'):
        mhc_allele = 'H-2-' + mhc_allele[3:]
    else:
        continue
    
    pep_length = len(row['pep_seq'])

    if mhc_allele not in adict:
        adict[mhc_allele] = set()
    
    adict[mhc_allele].add(pep_length)

for k, v in adict.items():
    adict[k] = sorted(list(v))

with open('mhc_to_pep_length.json', 'w') as f:
    json.dump(adict, f, indent=4)
