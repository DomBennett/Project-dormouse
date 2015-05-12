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
from parts import parts  # import dictionary of different sequence parts


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


def getIngroupSeqs(infile):
    # read in GB ingroup sequences
    sequences = []
    with open(infile, 'rb') as f:
        for record in SeqIO.parse(f, "gb"):
            sequences.append(record)
    # convert to fasta
    sequences_fasta = []
    for s in sequences:
        if s.id not in exclude:
            sequences_fasta.append(getFasta(s))
    return(sequences_fasta)


def getOutgroupSeq(infile):
    if os.path.isfile(infile):
        sequences = []
        with open(infile, 'rb') as f:
            for record in SeqIO.parse(f, "gb"):
                sequences.append(record)
        # choose one to add to sequences_gb
        outgroup = sequences[selection]
        # rename
        outgroup = getFasta(outgroup, name='outgroup')
        # add to sequences_fasta
        return(outgroup)
    else:
        return(None)

if __name__ == '__main__':
    # DIRS
    output_dir = '1_sequences'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    input_dir = 'data'
    blast_dir = os.path.join(input_dir, 'BLAST_results')
    sample_dir = os.path.join(input_dir, 'samples')

    # PROCESS
    for key in parts.keys():
        print key
        print os.path.join(blast_dir, parts[key]['ingroup'])
        if not os.path.isfile(os.path.join(blast_dir, parts[key]['ingroup'])):
            next
        sequences = getIngroupSeqs(os.path.join(blast_dir,
                                                parts[key]['ingroup']))
        # read in GB outgroup sequences
        outgroup = getOutgroupSeq(os.path.join(blast_dir,
                                               parts[key]['outgroup']))
        if outgroup:
            sequences.append(outgroup)
        # read in primer pair samples
        samples = []
        with open(os.path.join(sample_dir, parts[key]['sample']), 'rb') as f:
            for record in SeqIO.parse(f, "fasta"):
                samples.append(record)
        # rename as 'sample' and add to sequences_fasta
        for sequence in samples:
            sequences.append(getFasta(sequence, name='sample'))
        # write out as fasta
        with open(os.path.join(output_dir, parts[key]['outfile']), 'wb') as f:
            for s in sequences:
                f.write("{0}\n".format(s))
