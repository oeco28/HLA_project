#!/usr/bin/env python

"""findBarcodes.py: Assign unclassified ONT reads to appropriate barcode while allowing for mismatches"""

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
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
warnings.simplefilter("ignore")


#parser=argparse.ArgumentParser(
#    description='''The script will identify sequences in the unclassified pile that match the provided barcode''',
#    epilog='''usage: $python find_unclassified_ONTBarcodes.py barcode <string> max_dist <integer> query_file <string> output_file <string>''')


parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-b", "--barcode", default="", type=str, help="Barcode sequence")
parser.add_argument("-d", "--max_distance", default=1, type=int, help="Maximum edit distance from barcode to mathcing sequence in query")
parser.add_argument("-i", "--input_file", default="", type=str, help="path and name of input file")
parser.add_argument("-o", "--output_file", default="", type=str, help="path and name of output file")
args = parser.parse_args()


# Set up parameters
barcode = args.barcode
my_dist = int(args.max_distance)
myquery_file = args.input_file
output_file_p = args.output_file
myoutput_file = output_file_p + ".out"
output_index = output_file_p + ".index"
output_fastq = output_file_p + ".fastq"




#Inform the user which barcode they are using
print("Analyzing barcode:", barcode)

# Function to read fastq sequence

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


#Create a list to save findings 
myfindings = [] # Create an empty list to store the results of the search

# Find barcode in reads and save into list
def barcode_search(query):
        # Performs the search on each sequence
	return myfindings.append(find_near_matches(barcode, query, max_l_dist = my_dist))
 
processed_search = [barcode_search(number) for number in dna_sequences]


# A cleaning step to make sure we keep only non-empty elements of the list (index and output of matching algorithm
res = [i for i in range(len(myfindings)) if myfindings[i] != []];
clean_findings = [ele for ele in myfindings if ele != []];



def extract_reads(fastq_file, read_indices, output_file):
  """
  Extracts reads based on indices from a FASTQ file and saves them to a new file.

  Args:
      fastq_file (str): Path to the input FASTQ file.
      read_indices (list): List of integer indices of reads to extract (1-based indexing).
      output_file (str): Path to the output FASTQ file.
  """

  with open(fastq_file) as fastq_in, open(output_file, "w") as fastq_out:
    # Validate indices (assuming 1-based indexing)
    num_reads = sum(1 for _ in SeqIO.parse(fastq_in, "fastq"))
    for index in read_indices:
      if index < 1 or index > num_reads:
        print(f"WARNING: Index {index} is out of range. Skipping...")
        continue

    # Reset fastq_in to beginning for proper iteration
    fastq_in.seek(0)

    # Extract and write specific reads
    read_counter = 0  # Track current read position
    for record in SeqIO.parse(fastq_in, "fastq"):
      read_counter += 1
      if read_counter in read_indices:
        SeqIO.write(record, fastq_out, "fastq")

if __name__ == "__main__":
  # Replace with your file paths and desired indices (1-based)
  fastq_file = "your_fastq_file.fastq"
  read_indices = [10, 303, 501]
  output_file = "extracted_reads.fastq"

  extract_reads(myquery_file, res, output_fastq)
  print("Extraction of reads completed. Output saved to:", output_fastq)


# Save explanatory output to files
# The output_file will contain the record of what reads contain a barcode, what are the coordinates of the match and the edit distance between barcode and read
# The output_index will contain the an index of the barcode (based on the order in the original fastq) that matches a barcode

with open(myoutput_file, 'w') as f:
    for line in clean_findings:
        f.write(f"{line}\n")

with open(output_index, 'w') as f:
    for line in res:
        f.write(f"{line}\n")

print("Process completed. Detailed output saved to:", myoutput_file, ", and index of reads saved to:",output_index)

# Because clean_findings is a complex output and I don't have the time to work it out and separate elements into single lists, I will save each list independently and then you have to run a simple paste command in linux to create a final (useful) outcome
# This is the command:
# paste -d "/t" "output_file" "index_file" > "Final_barcode_out.tab"
