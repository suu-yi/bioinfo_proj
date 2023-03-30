#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 17:10:58 2023
@author: citra

############### UPDATE USAGE

Description:
    This script is used to improve the annotations of the reference genome 
    created by roary, after fix_pan_core.py has updated it with PA/Uniprot tags 
    from gff files.
    Uses the first hit from blast results to update unlabelled prokka ids from 
    a reference genome and makes a list of the update ids.

Input:
    diamond blast result
    reference genome
    
Output:
    updated reference genome
    list of replaced ids and their replacements

Usage:
    python moblast_annot.py
    


"""

import argparse

parser = argparse.ArgumentParser(description='Replace prokka IDs with first hit blast results')
parser.add_argument('-b','--blast', type=str, metavar='', help='Blast format output file from diamond blastx')
parser.add_argument('-i', '--input', type=str, metavar='', help='Input reference genome in fasta format')
parser.add_argument('-o', '--output', type=str, metavar='', default='new_reference.fa', help='Name for your updated reference')

args = parser.parse_args()


###############################################################################

names = []
prokka2 = ''

with open(args.blast, 'r') as blast:
    for line in blast:
        ids = line.replace('\\t', '\t').split()[:3:2]
        prokka = ids[0]
        if prokka != prokka2:
            names.append(ids)
        prokka2 = prokka 

out = open(args.output, 'w')

count = 0

with open(args.input, 'r') as panin:
    for line in panin:
        if line.startswith('>'):
            found = False
            #print(line)
            old = line.strip('>').split()[0]
            #print(old)
            for entry in names:
                if old == entry[0]:
                    print(line.rstrip().replace(old, entry[1]), file=out)
                    print(entry[1])
                    count = count + 1
                    found = True
                    break    
            if not found:
                print(line.rstrip(), file = out)
        else:
            print(line.rstrip(), file =out)
            
out.close()

with open('new_blast_tags.tsv', 'w') as listout:
    for entry in names:
        print( '\t'.join(entry), file=listout)
        
