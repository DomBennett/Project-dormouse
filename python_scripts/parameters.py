#! /bin/usr/env python
# D.J. Bennett
# 04/07/2015
"""
Parameter settings
"""

normalise_name = False  # normalise the name to Genus_species
byid = False  # combine by seqid or species name
min_parts = 3  # minimum number of parts in a sequence
# Use this to select which genes to include, form: GENE_PART
part_names = ['18S_p1', '18S_p2', '18S_p3', '28S_p1', '28S_p2']
#part_names = ['28S_p1', '28S_p2']

# Dictionary of the file locations of each gene section
parts = {'18S_p1': {'ingroup': '18S_p1_930F_Cyclophyllidea.gb',
                    'outgroup': '18S_p1_cestode6_Caryophyllidea.gb',
                    'sample': '18S_p1.fasta',
                    'outfile': '18S_p1.fasta'},
         '18S_p2': {'ingroup': '18S_p2_600F_Cyclophyllidea.gb',
                    'outgroup': '18S_p2_600F_Caryophyllidea.gb',
                    'sample': '18S_p2.fasta',
                    'outfile': '18S_p2.fasta'},
         '18S_p3': {'ingroup': '18S_p3_wormA_Cyclophyllidea.gb',
                    'outgroup': '18S_p3_wormA_Caryophyllidea.gb',
                    'sample': '18S_p3.fasta',
                    'outfile': '18S_p3.fasta'},
         '28S_p1': {'ingroup': '28S_p1_Cyclophyllidea.gb',
                    'outgroup': '28S_p1_Caryophyllidea.gb',
                    'sample': '28S_p1.fasta',
                    'outfile': '28S_p1.fasta'},
         '28S_p2': {'ingroup': '28S_p2_Cyclophyllidea.gb',
                    'outgroup': '28S_p2_Caryophyllidea.gb',
                    'sample': '28S_p2.fasta',
                    'outfile': '28S_p2.fasta'}}
