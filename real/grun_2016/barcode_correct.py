"""
pre-correct CEL-Seq barcodes prior to running kallisto
input: a fastq file where each read contains a cell barcode
input: a text file with a list of valid barcodes extdata/barcodes.txt
output: the same fastq file where each "N" in a barcode is replaced with the closest match from the valid barcode list
"""

# assumes the barcode is in the first 8 characters of a line

import re

bc_list = open("extdata/barcodes.txt","r").readlines()
bc_list = tuple(i.replace("\n","") for i in bc_list)
#bc_list = ('CAGTCTCG','TAGTCTCG','ACAGTATC')

#bc = 'NAGTCTCG'
#match_types = {}

def barcode_correct(bc,bc_list):
    """
    Given a barcode and a list of valid barcodes, checks if the barcode matches a valid barcode by replacing any "N" characters. If so, the corrected barcode is returned. Otherwise, it returns None
    """
    if 'N' not in bc:
        return None
    else:
        pattern = re.compile(bc.replace('N','[ACGT]'))
        matches = (pattern.match(x) for x in bc_list)
        matches = tuple(x.string for x in matches if x)
        if len(matches)==1: return matches[0]
        else: return None


ifile_name = "data/original/debug/test1.fastq"
ofile_name = "data/original/debug/test1_corrected.fastq"

counter = 0
num_chg = 0
with open(ifile_name,'r') as ifile, open(ofile_name,'w') as ofile:
    for line in ifile:
        counter += 1
        if counter % 4 != 2:
            ofile.write(line)
        else:
            bc_correct = barcode_correct(line[:8], bc_list)
            if bc_correct:
                num_chg += 1
                ofile.write(bc_correct+line[8:])
            else:
                ofile.write(line)
