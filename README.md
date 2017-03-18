# Miscellaneous Scripts

### 1. contig_stats.pl
This script will calculate and output a number of qualitative length statistics for any de novo transcriptome assembly.

Example usage: $ perl contig_stats.pl input.fasta

### 2. count_fasta.py
This script will count the number of each nucleotide (A,C,G,T,N), add them together, give the length of the sequence, and GC%. This scripts uses BioPython and can take a multifasta file. Specify your input file name in the script.

Example usage: $ python count_fasta.py

### 3. extract_features.py
This script can be used to extract CDS or other features from a genbank file. Each feature will be output to a separate fasta file.

Example usage: $ python extract_features.py filename.gbk

### 4. extract_seqs_from_fasta.py
This script extracts fasta sequences provided in comma-separated list from a provided input fasta file.

Example usage: $ perl extract_seqs_from_fasta.pl IDlist.csv input.fasta output.fasta

### 5. histogram_fasta.pl
This script will calculate histogram intervals for contig lengths and GC content for de novo transcriptome assemblies.

Example usage: $ perl count_fasta.pl input.fasta

### 6. mean_distance_from_alignments.py
This script will calculate pairwise mean patristic and Jukes-Cantor genetic distances for a directory of sequence alignments and outputs a table .

Example usage: $ python mean_distance_from_alignment.py in_dir aln_format

### 7. File format conversion
  fasta2phylip.pl, nexus2fasta.pl, phylip2fasta.pl
  
Example usage: $ perl fasta2phylip.pl in.fa out.phy
