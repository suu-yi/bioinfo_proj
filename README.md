# RNA-seq analysis pipeline for *Pseudomonas aeruginosa* in clinical samples

This is an attempt to build a RNA-seq analysis pipeline for P. aeruginosa using pangenome reference.


```mermaid
graph TD
        id3[NCBI sequence data\n.] --> id4[pan/core reference genome\n.]
        id4 --> id5[annotation]
  id0((raw reads)) --> id1((decontaminate))
        id1 --> id2((human depleted reads))
        id2 & id5 --> id6(mapping/alignment\n.) --> id7(normalization\n.) --> id8(visualisation\n.)

classDef default fill:#bbf,stroke:#f66,stroke-width:2px,color:#fff;
classDef ref fill:#f96;
classDef de fill:#A2D9CE;
class id3,id4,id5 ref;
class id6,id7,id8 de;
```

## Introduction
This is a work in progress... only private view for now

## Pan/core genome reference

### Overview

#### File structure

```bash

00_data
# nucleotide and amino acid sequences of 21 PA strains
# nuc from ncbi, aa from pseudomonas.com

01_prokka
# result files from prokka, gff files are used for roary
# only files from prokka seem to work with roary??

02_roary
# roary results from prokka gff files
# reference pangenome

03_clean_ref
# coregenome needs to be extracted
# replaced prokka ids with PAtags from gff files
# blast remaining unlabelled seqs
# replace remaining with first hit

```

### Data collection

21 strains of *Pseudomonas aeruginosa* (PA) organisms from KEGG GENOME Database were used. The genome data were collected from GenBank 
through the link provided by the KEGG record, and their corresponding protein sequence data were downloaded from [pseudomonas.com](http://pseudomonas.com). 

[Record of PA strains](00_data/PA_strains.csv)

Retrieving files with URLs:

```bash
# download genome assembly from ncbi
wget -i PA_genome_list.txt

# get protein sequence files from pseudomonas.com
wget -i PA_AA_list.txt
```

### Creating a pangenome

```bash
# SOME PREAMBLE ABOUT RENAMING FILES:

# rename AA protein files with loop
# first remove the Pseudomonas_aeruginosa_ and keep the strain names for easy ref
ls Pseudomonas* | while read line; do mv ${line} ${line#P*nosa_} ; done
# then add the corresponding genome assembly file name (genbank assembly ID) done manually , no loops for this

# make a list of all genbank genome ids
ls GCA* > genome_list.txt
# extract just the id
sed s/(.*.[12])_.*_.*.fna/1/ genome_list.txt > genome_list.txt
# add the strain names to ids manually, use underscore to join
# rename AA protein files by adding the genbank ids
cat strain_id.txt | while read line; do AA=$(echo ${line:16}.faa); mv ${AA} ${line}.faa ; done

# simplify the .faa headers for prokka
sed 's/ref|.*|//' Pseudomonas_aeruginosa_PAO1_107.faa
```

#### Prokaryotic genome annotation

- **Prokka**
    - version 1.14.6
    - installed with conda version 4.12.0

```bash
# running prokka ### looped
for f in 00_data/*cleaned.faa ; do prokka --kingdom Bacteria --outdir 01_prokka/$(basename ${f} _cleaned.faa) --genus Pseudomonas 
--notrna --proteins ${f} --cpus 8 --prefix $(basename ${f} _cleaned.faa) ${f:0:15}*genomic.fna ; done

```

#### Creating the pangenome

- **Roary**
    - version 3.13.0
- module installed on lunar aurora

```bash

# run roary: 8 threads, verbose, use mafft
roary -p 8 -v -e --mafft *.gff
```

### Creating reference genomes

#### Extracting the core genome

The file `gene_presence_absence.csv` contains information for each gene in the pangenome and which strains contains it. This is 
used to extract the list of gene names in the core genome, which is then used with the python script `fix_pan_core.py` to create 
the reannotated core genome.

```bash
 # extracting core genome
awk 'BEGIN {FS="\""}; {if ($8==21) print $2}' gene_presence_absence.csv > list_coregene
```

#### Reference genome reannotation

- python scripts:
    - fix_pan_core.py
    - moblast_annot.py
        - (combine the `awk` process into `fix_pan_core.py` ?)
        - (combine both py scripts?)

#### - Reannotate with PAtags

The reference pangenome file from roary is in a fasta format, with each sequence header containing a unique prokka assigned ID and 
a gene name. Most of the gene names are also named as group_???? by prokka, and since these do not provide much information about 
the sequences, these will need to be reannotated. The .gff files from prokka, which were used in roary, also contain the PA tags, 
and these files will be used with the  `fix_pan_core.py` script to update the prokka IDs with PA tags that can be searched on 
pseudomonas.com

There are some sequences which would not have PA tags from the gff files, these sequences will be saved in another file which can 
be put into diamond later.

```bash
# for lunarc aurora cluster: try biopython
module add GCC/7.3.0-2.30 OpenMPI/3.1.1
module add Biopython/1.73-Python-3.6.6

# extracting core genome
awk 'BEGIN {FS="\""}; {if ($8==21) print $2}' gene_presence_absence.csv > list_coregene

# run script to update prokka IDs and create coregenome from list_coregene
python fix_pan_core.py

# if using fix_pangenome_v1.py and create coregenome with seqtk
module add seqtk/1.2
seqtk subseq pan_genome_reference_simpleheader.fa list_coregene > coregenome.fa

```

#### - More annotations

The file of sequences that were not updated are searched through all 21 of the previously downloaded `.faa` protein files from 
[pseudomonas.com](http://pseudomonas.com) with `diamond blastx` 

- diamond blast
    - version 2.1.4
    - installed with conda version 4.12.0 on lunarc aurora
    
    ```bash
    # make diamond database from all PA.faa files
    cat *.faa > all_PA.faa
    diamond makedb --in all_PA.faa --db PA_db --log
    
    # run diamond blastx ()
    diamond blastx -d PA_db.dmnd -q ../notags_seqs_pan_genome_reference.fa -f 6 -o notag_seqs.out
    ```
    

The diamond blast results output file is used to find the sequences that matched and the first hit for each sequence with identity 
\> 70% and  e-value < 0.001.

Using `moblast_annot.py` , the prokka IDs in the reference genome are updated with the new PA tags found by diamond blast.

```bash
# run moblast_annot.py
python moblast_annot.py
```

### Exploring the reference genome
- Using the script from roary [link](https://github.com/sanger-pathogens/Roary/blob/21ffb84504fd55d256eca90a47e3f2f5a9012c5c/contrib/roary_plots/roary_plots.py) to explore the pangenome with plots.

```bash
python roary_plots.py tree.file gene_presence_absence.csv --format pdf --labels
```


### Reproducibility

- **GNU Wget:** version 1.14 built on linux-gnu
- **Prokka:** version 1.14.6
- **Roary:** version 3.13.0
- **DIAMOND**: version 2.1.4
- **seqtk:** version 1.2
- **GNU Awk:** version 4.0.2
- **Python:** version 3.6.6

## Human reads removal from sequenced sample data

### Synopsis

- The sample reads listed in the table were depleted of human reads using a combination of two different methods of detecting human reads: taxonomy classification method with the software Kraken2 (version 2.1.1) and alignment method software bowtie2 (version 2.4.4).
- The two-step method is used to ensure all human reads are removed from the samples. Different methods of detecting human reads in microbial sequencing datasets have been tested by [Bush et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7478626/)
- Kraken2 software was used for the first step in detecting human reads. Using the .kraken output, the sequence ID for the reads that were not assigned by kraken2 as ‘Homo sapien’ are saved as a list and used with seqtk (version1.2) subseq command to extract non-human reads from the sample reads files.
- The subsequent reads are then mapped to the human genome GRCh37 from NCBI using bowtie2. SAM flags are interpreted using the [Picard utility](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7478626/) in the resulting SAM file output from bowtie2. SAMtools (version 1.15.1) is used to find reads flagged as unmapped and both reads unmapped for paired-end reads. These are then extracted into gzip compressed FASTQ files, completing the second step of removing human reads.
- A second kraken report is made for the final cleaned product.

### FastQC
Initial quality of the sequences
MultiQC for all thge reports from fastQC


### Kraken2

- version 2.1.1
- Database:
    - The database used for kraken2 was a pre-made index (dated 12/9/2022) from the standard collection containing archae, bacteria, plasmid and human data.
    - Index downloaded from: [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)
- Running Kraken2
    - Paired-end reads and single-end reads are processed separately
    
    ```bash
    # kraken paired-end reads
    for f in ../00_fastq/*1_001.fastq.gz;
    do base_nm=$(basename ${f} 1_001.fastq.gz);
    
    kraken2 --use-names \
            --threads 64 \
            --db ../database/ \
            --gzip-compressed \
            --report ${base_nm}.report.kraken \
            --output ${base_nm}.kraken \
            --paired ${f} ${f%%1_001.fastq.gz}2_001.fastq.gz;
    echo "...Finished run for ${base_nm} at $(date)";
    done
    
    # kraken single-end reads
    for f in ../00_fastq/*.gz;
    do base_nm=$(basename ${f} .fastq.gz);
    
    kraken2 --use-names \
            --threads 64 \
            --db ../database/ \
            --gzip-compressed \
            --report ${base_nm}.report.kraken \
            --output ${base_nm}.kraken \
            ${f};
    echo "Finished run for ${base_nm} at $(date)";
    done
    ```
    

### Seqtk

- version 1.2
- using kraken output files, extract list of sequence ids NOT assigned as 'Homo sapiens'
- use seqtk subseq command with list of non-human ids to extract reads from fastq

```bash
# extract seq id into list
for file in ../kraken_results/*1_001.kraken;
do bs_name=$(basename ${file} .kraken);
grep -v 'Homo sapiens' ${file} | cut -f2 > list_cleaned_${bs_name};
rm ${file};

# seqtk subseq the listed ids from original fastq into new cleaned_fastq files
seqtk subseq ../20*/${bs_name}.fastq.gz list_cleaned_${bs_name} | gzip > cleaned_${bs_name}.fastq.gz;
# optional: remove large files to save space
rm ../20*/${bs_name}.fastq.gz;
echo "${bs_name} ...done at $(date)";
done
```

### Bowtie2

- version 2.4.4
- Bowtie pre-made index *H. sapiens,* GRch37
- Running Bowtie2
    
    ```bash
    # bowtie2 paired-end
    for f in ../kraken/cleaning/*1_001.fastq.gz;
    do bs_name=$(basename ${f} 1_001.fastq.gz);
    bowtie2 -t -p 16 -x GRCh37/GRCh37 \
            -1 ${f} \
            -2 ${f%%1_001.fastq.gz}2_001.fastq.gz \
            -S ${bs_name}.sam;
    echo ">>>>> ${bs_name} ...done at $(date)";
    done
    
    # bowtie2 single-end reads
    for f in ../kraken/cleaning/*1_001.fastq.gz;
    do bs_name=$(basename ${f} 1_001.fastq.gz);
    bowtie2 -t -p 32 -x ../bowtie2/GRCh37/GRCh37 \
            -U ${f} \
            -S ${bs_name}.sam;
    rm ${f};
    echo "${bs_name} ...done at $(date)";
    done
    ```
    

### SAMtools

- version 1.15.1
- SAMtools is utilised for the second step of removing human reads, filtering for reads unmapped to human genome using SAM output from bowtie2 and creating the final FASTQ files
- Running SAMtools

```bash
# paired-end
for f in ../bowtie2/cleaned*;
do bs_name=$(basename ${f} .sam);
# convert SAM to BAM
samtools view -b ${f} -o ${bs_name}.bam;
# sort in order of names
samtools sort -n ${bs_name}.bam -o ${bs_name}.bam;
samtools fastq -f 0x4 \
        -1 ${bs_name}1.fastq.gz \
        -2 ${bs_name}2.fastq.gz \
        -0 /dev/null \
        -s /dev/null -n ${bs_name}.bam;
echo " ${bs_name} ...done at $(date)";
done

# single-end reads
for f in ../bowtie2/cleaned*;
do bs_name=$(basename ${f} .sam);
# convert SAM to BAM
samtools view -b ${f} -o ${bs_name}.bam;
# optional: rm large files
rm ${f};
samtools fastq -f 0x4 \
        -0 ${bs_name}1.fastq.gz \
        -s /dev/null -n ${bs_name}.bam;
echo " ${bs_name} ...done at $(date)";
done
```

### FastQC

- version v0.12.1
- installed with conda

```bash
# noextract will not extract zip files
fastqc --noextract -o fastqc seq_1.fastq.gz seq_2.fastq.gz 
```

### kallisto 

- version 0.46.0
- module is also available on lunarc

### Create the index

Using the reference pan/core-genome, create the database index used by kallisto

### Mapping reads to the reference

The `kallisto quant` command maps the sample reads to the reference and outputs the results into their 
corresponding directories.

- **Single end reads**
    
    (Different command for SE reads)
    
- **Paired end**
```bash

echo ">>>>>Starting step1: build kallisto index at $(date)"
kallisto index --make-unique -i PA_pan.idx newtag_pan_genome_reference.fa

echo "Done at $(date)"

echo ">>>>>Starting step2: kallisto quant at $(date)"

for f in ../DEMULTIPLEX_2020_19_R1_cleaned/*R1.fastq.gz;
do bs_name=$(basename ${f} R1.fastq.gz);
kallisto quant -i PA_pan.idx \
        -t 16 \
        -o ${bs_name} \
        ${f} \
        ${f%%R1.fastq.gz}R2.fastq.gz;
echo "${bs_name} ...done at $(date)";
done

echo "COMPLETED at $(date)"
```
