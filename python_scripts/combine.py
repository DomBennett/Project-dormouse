#! /bin/usr/env python
# D.J. Bennett
# 28/03/2015
"""
Combine alignments from each pair into single supermatrix
"""

# PACKAGES
import os
import re
from Bio import SeqIO
from Bio import AlignIO

# PARAMETERS
sample_ids_f = ['sample__930f', 'sample__600f', 'sample__wormA']
sample_ids_r = ['sample__cestode6', 'sample__A27', 'sample__1270R']


# FUNCTIONS
def readSequences(infile):
    if not os.path.isfile(infile):
        raise IOError('No {0}, are you sure you`ve run previous stage?'.
                      format(infile))
    sequences = []
    with open(infile, 'rb') as f:
        for record in SeqIO.parse(f, "fasta"):
            # rename sample r and f ids to share same ID
            if record.id in sample_ids_f:
                record.id = 'sample__f'
                record.description = ''
            if record.id in sample_ids_r:
                record.id = 'sample__r'
                record.description = ''
            sequences.append(record)
    return(sequences)


def trimSequences(sequences):
    starts = []
    ends = []
    for s in sequences:
        starts.append(re.search('[atcgATCG]', str(s.seq)).start())
        ends.append(len(s) - re.search('[atcgATCG]', str(s.seq)[::-1]).end())
    res = []
    for s in sequences:
        res.append(s[max(starts):min(ends)])
    return(res)

# DIRS
inout_dir = '2_alignments'

# READ
p1_seqs = readSequences(os.path.join(inout_dir, 'p1_alignment.fasta'))
p2_seqs = readSequences(os.path.join(inout_dir, 'p2_alignment.fasta'))
p3_seqs = readSequences(os.path.join(inout_dir, 'p3_alignment.fasta'))

# TRIM
p1_seqs = trimSequences(p1_seqs)
p2_seqs = trimSequences(p2_seqs)
p3_seqs = trimSequences(p3_seqs)

# MATCH IDS INTO SINGLE DICTIONARY
sequences = {}
for s in p1_seqs:
    sp, seqid = s.id.split('__')
    if seqid in sequences.keys():
        sequences[seqid]['p1'].append(s)
    else:
        sequences[seqid] = {'p1': [s], 'p2': [], 'p3': []}

for s in p2_seqs:
    sp, seqid = s.id.split('__')
    if seqid in sequences.keys():
        sequences[seqid]['p2'].append(s)
    else:
        sequences[seqid] = {'p1': [], 'p2': [s], 'p3': []}

for s in p3_seqs:
    sp, seqid = s.id.split('__')
    if seqid in sequences.keys():
        sequences[seqid]['p3'].append(s)
    else:
        sequences[seqid] = {'p1': [], 'p2': [], 'p3': [s]}

# CONSTRUCT SUPERMATRIX
supermatrix = []
for key in sequences:
    if sequences[key]['p1'] and sequences[key]['p2'] and sequences[key]['p3']:
        # stick together, always use the first element
        s1 = sequences[key]['p1'][0]
        s2 = sequences[key]['p2'][0]
        s3 = sequences[key]['p3'][0]
        s = s1 + s2 + s3
        # rename outgroup
        sp, _ = s.id.split('__')
        if 'outgroup' == sp:
            s.id = 'outgroup'
            s.description = ''
        supermatrix.append(s)

# OUTPUT
alignment = AlignIO.MultipleSeqAlignment(supermatrix)
outfile = os.path.join(inout_dir, 'supermatrix.phy')
with open(outfile, "w") as f:
    # write out using PhylipWriter in order to extend id_width
    AlignIO.PhylipIO.PhylipWriter(f).write_alignment(alignment, id_width=40)
