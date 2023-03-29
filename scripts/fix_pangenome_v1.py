#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 09:11:18 2023

@author: su

Run:
    python fix_pangenome.py
    or 
    python fix_pangenome.py -i <pangenome.filename>

Input:
    - gff files fed into roary
    - pangenome reference from roary
    
Output:
    - pangenome with updated headers
    - list of seqs not updated
    - list of updated
    - fasta of all unknown seqs that can be used in blastx or blastn
    - log file of stats


"""

# argparse and file management

import os
import argparse

parser = argparse.ArgumentParser(description='Replace prokka IDs with locustags in the prokka .gff files')
parser.add_argument('-i','--input', type=str, default='pan_genome_reference.fa', metavar='', help='Pangenome reference fasta file')

args = parser.parse_args()

gff_count = 0
PAtags = 0
unitags = 0

file_list = []
taglist = []

cwd_path = os.getcwd()
for file in os.listdir(cwd_path):
    f_name, f_ext = os.path.splitext(file)
    if f_ext == '.gff':
        gff_count = gff_count + 1
        file_list.append(file)

        gff_in = open(file,'r')

        version = gff_in.readline() # pass version line
        seq_region = gff_in.readline() # pass seq region line

        gene_details = []
        PAlist = []
        unilist = []
        notag = []

        for line in gff_in:
            if line.startswith("##"):
                break
            else:
                description = line.rstrip().split("\t")[8].split(";")
                # print(description)
                gene_details.append(description)

        for entry in gene_details:
            pair = []
            if "PA" in entry[1]:
                tag = entry[1].split(':')[4]
                ID = entry[0].strip('ID=')
                pair.append(ID)
                pair.append(tag)
                PAlist.append(pair)
            elif len(entry) >= 5:
                if "UniProtKB" in entry[4]:
                    # print(entry)
                    tag = entry[4].split(':')[4]
                    ID = entry[0].strip('ID=')
                    pair.append(ID)
                    pair.append(tag)
                    unilist.append(pair)
            else:    
                notag.append(entry)
        
        taglist = taglist + PAlist + unilist
        PAtags = PAtags + len(PAlist)
        unitags = unitags + len(unilist)
        
########### import pangenome fasta file data (WIP)

pangen = open(args.input, 'r')

names = [] # list of headers in pangenome
seq_frag = []
seq_list = [] # list of complete seqs

# seq_frag is used since the sequence can be in multiple lines
# these sections of a sequence will need to be concatenated

for line in pangen:
    if line.startswith('>'):
        nam = line.strip('>\n').split(' ') # remove '>' and newline
        names.append(nam) #list of list?
        if seq_frag:    #find the list of seq fragments, if not empty
            sequence = "".join(seq_frag)     #join seq fragments
            seq_list.append(sequence) #add sequences to the list
        seq_frag = [] # empties the list for new sequences
    else:       #find more exiting sequences
        seq = line.rstrip()
        seq_frag.append(seq)
if seq_frag:    #find the remaining seq fragments in the file
    sequence = "".join(seq_frag)
    seq_list.append(sequence.upper())   #add completed sequence into the list
    
########## match ID in pangenome and replace with new PA tag

# add info to list of list????? use in loop for each gff file

out = open(f'new_{args.input}', 'w')
unk_out = open(f'notags_seqs_{args.input[:-3]}.fa', 'w')
unk_list = open(f'notags_list_{args.input[:-3]}.txt', 'w')
updated = open(f'updated_list_{args.input[:-3]}.txt', 'w')

# header
updated.write('Prokka ID' + '\t' + 'Locustag' + '\t' + 'gene_name\n')

count = 0
not_grp = 0  # in gff but not group_ for name
unk = 0
# def add_info(gene_details):
for num in range(len(names)): # index pangen seqs
    found = False
    gen_id = names[num][0] # prokka id
    for tag in taglist: # prokka id and tag
        if gen_id == tag[0]:
            count = count + 1
            found = True
            print('>' + tag[1] + '\t' + names[num][1], file=out) # >PA/Unitag gene/group
            print(seq_list[num], file=out) # corresponding seq
            print(names[num][0] + '\t' + tag[1] + '\t' + names[num][1], file=updated )
            break
    if not found:
        print('\t'.join(names[num]), file=unk_list)
        print('>' + '\t'.join(names[num]), file=out)
        print(seq_list[num], file=out)
        print('>' + '\t'.join(names[num]), file=unk_out)
        print(seq_list[num], file=unk_out)
        unk = unk + 1
        if not names[num][1].startswith('group'):
            not_grp = not_grp + 1

############ log file

with open(f'output_{args.input[:-3]}.log', 'w') as logs:
    print(f'>>> Using {gff_count} gff files', file=logs)
    for n in file_list:
        print(n, file=logs)
    print(f'\n>>> Pangenome reference: {args.input}\nFound {len(names)} genes', file=logs)
    print(f'Updated {count} with locustag or UniProtKB accession', file=logs)
    print(f'{not_grp} have gene names\n{unk-not_grp} remaining have prokka ID and group_????', file=logs)
    print(f'>>>Output files\nnew_{args.input}\nnotags_seqs_{args.input[:-3]}.fa\nnotags_list_{args.input[:-3]}.txt\nupdated_list_{args.input[:-3]}.txt\noutput_{args.input[:-3]}.log')

out.close()
unk_list.close()
unk_out.close()
updated.close()

