# pipeline notes

```mermaid
graph TD
	id3[NCBI sequence data\n.] --> id4[pan/core reference genome\n.]
	id4 --> id5[annotation]
  id0((raw reads)) --> id1((decontaminate))
	id1 --> id2((human depleted reads))
	id2 & id5 --> id6(mapping/alignment\n.) --> id7(normalization\n.) --> id8(quantification\n.)

classDef default fill:#bbf,stroke:#f66,stroke-width:2px,color:#fff;
classDef ref fill:#f96;
classDef de fill:#A2D9CE;
class id3,id4,id5 ref;
class id6,id7,id8 de;
```

```mermaid
graph LR
id[reference PA data\n.] --> id9[prokka\n.] --> id10[roary\n.] --> id6[pan/core\ngenome\n .] --> id11[annoated reference]
id0((clinical samples\n.)) --> id1((kraken2)) --> id2((\nseqtk\n\n .))
id2 --> id3((bowtie2))
id3 --> id4((samtools))
id4 --> id5((\nfastqc\n\n.))
id11 & id5 --> id7(kallisto\n.) --> id8(deseq2?\n .)

classDef default fill:#bbf,stroke:#f66,stroke-width:2px,color:#fff;
classDef ref fill:#f96;
classDef de fill:#A2D9CE;
class id,id6,id9,id10,id11 ref;
class id7,id8 de;
```

### bioinfo methods?

- comparing methods(?)how/why certain tools are chosen?:
    - decontamination process
    - pan/core genome softwares
        - why pan/core genome
        - compare with using PA reference PAO1 genome?
    - alignment softwares

version control/reproducibility: github

***how do these apply to PA?***

why not use established databases such as BACTOME?

- further(?):
    - reference-free *de novo* transcriptome?
    - metatranscriptomics?