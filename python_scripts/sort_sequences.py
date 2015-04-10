#! /bin/usr/env python
# D.J. Bennett
# 19/03/2015
"""
Rename sequences to species names.
Bring BLAST results and primer pairs together.
"""

# PARAMETERS
selection = 0  # index of the outgroup seq
exclude = ['LM405059.1']  # list of Accessions to ignore

# PACKAGES
import os
from Bio import SeqIO


# FUNCTION
def getFasta(sequence, name=None):
    '''Take SeqRecord and return fasta with name: "species_accession"'''
    if not name:
        # get sp name
        name = sequence.annotations['organism']
        name = name.replace(' ', '_')
    # add the accession
    name += '__' + sequence.id
    sequence.id = name
    sequence.description = ''
    return(sequence.format('fasta'))


# DIRS
output_dir = '1_sequences'
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
input_dir = 'data'
blast_dir = os.path.join(input_dir, 'BLAST_results')

# PROCESS
# dict of infiles and outfiles
files = {'p1': {'ingroup': os.path.join(blast_dir,
                                        'p1_930F_Cyclophyllidea.gb'),
                'outgroup': os.path.join(blast_dir,
                                         'p1_cestode6_Caryophyllidea.gb'),
                'sample': os.path.join(input_dir, 'p1.fasta'),
                'outfile': os.path.join(output_dir, 'p1.fasta')},
         'p2': {'ingroup': os.path.join(blast_dir,
                                        'p2_600F_Cyclophyllidea.gb'),
                'outgroup': os.path.join(blast_dir,
                                         'p2_600F_Caryophyllidea.gb'),
                'sample': os.path.join(input_dir, 'p2.fasta'),
                'outfile': os.path.join(output_dir, 'p2.fasta')},
         'p3': {'ingroup': os.path.join(blast_dir,
                                        'p3_wormA_Cyclophyllidea.gb'),
                'outgroup': os.path.join(blast_dir,
                                         'p3_wormA_Caryophyllidea.gb'),
                'sample': os.path.join(input_dir, 'p3.fasta'),
                'outfile': os.path.join(output_dir, 'p3.fasta')}}
for key in files.keys():
    # read in GB ingroup sequences
    sequences_gb = []
    with open(files[key]['ingroup'], 'rb') as f:
        for record in SeqIO.parse(f, "gb"):
            sequences_gb.append(record)
    # convert to fasta
    sequences_fasta = []
    for s in sequences_gb:
        if s.id not in exclude:
            sequences_fasta.append(getFasta(s))
    # read in GB outgroup sequences
    outgroup_sequences = []
    with open(files[key]['outgroup'], 'rb') as f:
        for record in SeqIO.parse(f, "gb"):
            outgroup_sequences.append(record)
    # choose one to add to sequences_gb
    outgroup = outgroup_sequences[selection]
    # rename
    outgroup = getFasta(outgroup, name='outgroup')
    # add to sequences_fasta
    sequences_fasta.append(outgroup)
    # read in primer pair samples
    samples = []
    with open(files[key]['sample'], 'rb') as f:
        for record in SeqIO.parse(f, "fasta"):
            samples.append(record)
    # rename as 'sample' and add to sequences_fasta
    for sequence in samples:
        sequences_fasta.append(getFasta(sequence, name='sample'))
    # write out as fasta
    with open(files[key]['outfile'], 'wb') as f:
        for s in sequences_fasta:
            f.write("{0}\n".format(s))
