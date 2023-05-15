import argparse
import os
import pandas as pd

# Define command-line arguments
parser = argparse.ArgumentParser(description='Convert a GTF file to GFF format')
parser.add_argument('input_file', metavar='input.gtf', help='input GTF file')
args = parser.parse_args()

# Determine the output file name and path
output_file = os.path.splitext(args.input_file)[0] + '.gff'

# Read in the GTF file using pandas
df = pd.read_csv(args.input_file, sep='\t', comment='#', header=None)

# Rename columns to match GFF format
df = df.rename(columns={0: 'seqid', 1: 'source', 2: 'type', 3: 'start', 4: 'end',
                        5: 'score', 6: 'strand', 7: 'phase', 8: 'attributes'})

# Map "transcript" entries to "mRNA"
df.loc[df['type'] == 'transcript', 'type'] = 'mRNA'

# Remove any other entries in the "type" column that are not part of the GFF format
valid_types = set(['gene', 'mRNA', 'exon', 'CDS', 'start_codon', 'stop_codon',
                   'five_prime_UTR', 'three_prime_UTR', 'noncoding_exon'])
df = df[df['type'].isin(valid_types)]

# Write out the GFF file
df.to_csv(output_file, sep='\t', index=False, header=False)

