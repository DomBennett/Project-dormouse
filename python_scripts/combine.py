#! /bin/usr/env python
# D.J. Bennett
# 28/03/2015
"""
Combine alignments from each pair into single supermatrix
"""

# PACKAGES
import os
import re
import sys
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from parameters import parts as all_parts
from parameters import part_names
from parameters import min_parts
from parameters import byid
from parameters import normalise_name


# PARAMETER PARSE
if min_parts > len(part_names):
    sys.exit('min_parts is greater than part_names')
# create parts dict from all_parts and part_names
parts = {}
for part in part_names:
    parts[part] = all_parts[part]


# FUNCTIONS
def readSequences(infile):
    if not os.path.isfile(infile):
        raise IOError('No {0}, are you sure you`ve run previous stage?'.
                      format(infile))
    sequences = []
    with open(infile, 'rb') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append(record)
    if normalise_name:
        for s in sequences:
            sp, seqid = s.id.split('__')
            sp = re.split('_', sp)[0:2]
            sp = '_'.join(sp)
            s.id = sp + '__' + seqid
    return(sequences)


def getPartitions(parts, part_names):
    """Return partition.txt"""
    begin = 1
    end = 0
    text = ''
    cgene = part_names[0].split('_')[0]
    ngene = 1
    for part in part_names:
        gene = part.split('_')[0]
        print gene
        print begin
        print end
        if cgene == gene:
            end += parts[part]['length']
        else:
            cgene = gene
            text += 'DNA, gene{0} = {1}-{2}\n'.\
                format(ngene, begin, end)
            ngene += 1
            begin = end + 1
            end += parts[part]['length']
    if ngene > 1:
        text += 'DNA, gene{0} = {1}-{2}\n'.\
            format(ngene, begin, end)
    return(text)


def trimSequences(sequences):
    """Return sequences without tail-ends"""
    starts = []
    ends = []
    for s in sequences:
        starts.append(re.search('[atcgATCG]', str(s.seq)).start())
        ends.append(len(s) - re.search('[atcgATCG]', str(s.seq)[::-1]).end())
    res = []
    for s in sequences:
        res.append(s[max(starts):min(ends)])
    return(res)


def getSeqDict(parts, part_names):
    """Return dictionary of seq IDs with sequences for each part"""
    seqdict = {}
    for key in parts.keys():
        # read
        alignment_file = key + '.fasta'
        alignment_file = os.path.join(input_dir, alignment_file)
        print alignment_file
        if not os.path.isfile(alignment_file):
            print 'not file'
            part_names.pop(part_names.index(key))
            continue
        seqs = readSequences(alignment_file)
        # trim
        seqs = trimSequences(seqs)
        # add length of trimmed seq to parts dict
        parts[key]['length'] = len(seqs[0])
        for s in seqs:
            sp, seqid = s.id.split('__')
            if sp == 'sample':
                seqid = 'sample'
            if byid:
                group_factor = seqid
            else:
                group_factor = sp
            if group_factor in seqdict.keys():
                if key in seqdict[group_factor].keys():
                    seqdict[group_factor][key].append(s)
                else:
                    seqdict[group_factor][key] = [s]
            else:
                seqdict[group_factor] = {key: [s]}
    return(seqdict, part_names)


def getSupermatrix(seqdict, parts, part_names):
    """Return large alignment of each sequence part combined"""
    supermatrix = []
    ngaps = 0
    for key in seqdict:
        temp_names = seqdict[key].keys()
        if len(temp_names) >= min_parts:
            # stick together, use random seq
            s = None
            seqids = []
            for part_name in part_names:
                if part_name in seqdict[key].keys():
                    ps = random.sample(seqdict[key][part_name], 1)[0]
                    sp, psid = ps.id.split('__')
                    if psid not in seqids:
                        seqids.append(psid)
                else:
                    # create seqrecord of '-'
                    ps = SeqRecord(Seq('-' * parts[part_name]['length']))
                    ps.id = s
                    ngaps += parts[part_name]['length']
                if s:
                    s += ps
                else:
                    s = ps
            # rename outgroup + sample
            sp = re.sub("_sp.*$", "", sp)
            sp = re.sub("_cf.*$", "", sp)
            if len(seqids) == 1:
                s.id = sp + '__' + seqids[0]
            else:
                s.id = sp + '__' + "_".join(seqids)
            s.description = ""
            if key != 'sample':
                if 'outgroup' == sp:
                    s.id = 'outgroup'
                    s.description = ''
            else:
                s.id = 'sample'
            supermatrix.append(s)
    return(supermatrix, ngaps)

if __name__ == '__main__':
    # DIRS
    input_dir = '2_alignments'
    output_dir = '3_supermatrix'
    # MATCH IDS INTO SINGLE DICTIONARY
    print part_names
    seqdict, part_names = getSeqDict(parts, part_names)
    # CONSTRUCT SUPERMATRIX
    supermatrix, ngaps = getSupermatrix(seqdict, parts, part_names)
    print len(supermatrix)
    print ngaps
    ngaps_psp = float(ngaps)/len(supermatrix)
    # GET PARTITIONS
    partition_text = getPartitions(parts, part_names)
    # OUTPUT
    alignment = AlignIO.MultipleSeqAlignment(supermatrix)
    print('Supermatix of [{0}] length and [{1}] species generated with [{2}] \
    gaps per species'.format(alignment.get_alignment_length(), len(alignment),
                             ngaps_psp))
    outfile = os.path.join(input_dir, 'supermatrix.phy')
    with open(outfile, "w") as f:
        # write out using PhylipWriter in order to extend id_width
        AlignIO.PhylipIO.PhylipWriter(f).write_alignment(alignment,
                                                         id_width=100)
    # OUTPUT PARITIONS
    if partition_text:
        outfile = os.path.join(input_dir, 'partitions.txt')
        with open(outfile, 'w') as file:
            file.write(partition_text)
