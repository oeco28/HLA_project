#!/usr/bin/env python

"""find_unclassified_ONTBarcodes.py: Assign unclassified ONT reads to appropriate barcode while allowing for mismatches"""

__author__      = "Omar E. Cornejo"
__copyright__   = "Copyright 2024, The HLA/Microbiome project"
__credits__     = ["Omar Cornejo", "Obed Garcia", "Sam Bogan"]
__license__     = "GPL"
__version__     = "1.0"
__maintainer__  = "Omar E Cornejo"
__email__       = "omcornej@ucsc.edu"
__status__      = "Beta"

import sys
from fuzzysearch import find_near_matches
from Bio import SeqIO
import numpy as np
import argparse

parser=argparse.ArgumentParser(
    description='''The script will identify sequences in the unclassified pile that match the provided barcode''', epilog='''The intended usage is python find_unclassified_ONTBarcodes.py barcode <string> max_dist <integer> query_file <string> output_file <string>''')
    
#parser.add_argument('', type=int, default=242, help='Barcode sequence')
#parser.add_argument('bar', nargs='*', default=242, help='maximum edit distance')
#parser.add_argument('', type=int, default=242, help='input file')
#parser.add_argument('', type=int, default=242, help='output file')
#args=parser.parse_args()


barcode = sys.argv[1]
my_dist = sys.argv[2]
myquery_file = sys.argv[3]
output_file = sys.argv[4]
output_index = print(output_file,".index",sep="")

print("Analyzing barcode:", barcode)


def process_dna(dna_sequences):
  # Iterate and perform operations on DNA sequences in the array
  for sequence in dna_sequences:
    print(f"DNA Sequence: {sequence}")

# Replace "your_fastq_file.fastq" with the path to your actual file
dna_sequences = []  # Create an empty list to store sequences
with open(myquery_file) as fastq_file:
  for record in SeqIO.parse(fastq_file, "fastq"):
    dna_sequence = str(record.seq).upper()  # Extract and convert to uppercase
    dna_sequences.append(dna_sequence)


myfindings = [] # Create an empty list to store the results of the search

def barcode_search(query):
        # Performs the search on each sequence
        return myfindings.append(find_near_matches(barcode, query, max_l_dist=my_dist))

processed_search = [barcode_search(number) for number in dna_sequences]


res = [i for i in range(len(myfindings)) if myfindings[i] != []];
clean_findings = [ele for ele in myfindings if ele != []];

# Because clean_findings is a complex output and I don't have the time to work it out and separate elements into single lists, I will save each list independently and then you have to run a simple paste command in linux to create a final (useful) outcome
# This is the command:
# paste -d "/t" "output_file" "index_file" > "Final_barcode_out.tab"

with open(output_file, 'w') as f:
    for line in clean_findings:
        f.write(f"{line}\n")

with open(output_index, 'w') as f:
    for line in res:
        f.write(f"{line}\n")

