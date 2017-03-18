#! /usr/bin/python

from Bio import SeqIO
import os
import sys

filename = sys.argv[1]

filebase = os.path.splitext(filename)[0]

handle = open(filename, 'r')

for record in SeqIO.parse(handle,'genbank'):
    seq = str(record.seq)
    for feature in record.features:
        if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA':
            if 'gene' in feature.qualifiers:
                geneName = feature.qualifiers['gene'][0]
            elif 'product' in feature.qualifiers:
                geneName = feature.qualifiers['product'][0]
            else:
                print('ERROR when parsing feature:')
                print(feature.qualifiers)
                quit()
            
            geneFile = open(filebase + '_' + record.id + '_' + geneName + '.fa','w')
            geneFile.write('>')
            geneFile.write(os.path.basename(filebase) + '_' + geneName + '\n')
            
            geneFile.write(seq[feature.location.start.position: feature.location.end.position])
            geneFile.write('\n\n')
