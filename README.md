# MetaGeneMark-2

MetaGeneMark-2: Improved Gene Prediction in Metagenomes

Karl Gemayel, Alexandre Lomsadze* and Mark Borodovsky*

Georgia Institute of Technology, Atlanta, Georgia, USA

*joint last authors

Reference: PAPER LINK


## Overview
MetaGeneMark-2 is an unsupervised metagenomic gene finder. It improves on MetaGeneMark by adding models for better gene start prediction, as well as automatic selection of genetic code (4 or 11). The models for gene start prediction are based in part on the work done for GeneMarkS-2; they include Shine-Dalgarno RBS, non-Shine-Dalgarno (or non-canonical) RBS, and bacterial and archaeal promoter models (for use in cases of leaderless transcription).

## Installing MetaGeneMark-2 locally
Running MetaGeneMark-2 using automatic genetic code detection is done through the `run_mgm.pl` script found in `src`.

     cd src;
     pf_makefile=Makefile.macos    # NOTE: change based on operating system
     make -f $pf_makefile


This generates a binary =gmhmmp2=.

## Running MetaGeneMark-2
Running MetaGeneMark-2 with automatic genetic code detection is done using `run_mgm.pl`. The following files should be in the same directory: `run_mgm.pl`, `gmhmmp2`, `mgm2_11.mod`, `mgm2_4.mod`. MetaGeneMark-2 can then be run (from anywhere) using:

    $path_to_binary/run_mgm.pl --seq [name]  --out [name]

    Required options:
         --seq  [name]            nucleotide sequence of metagenome in FASTA format.
         --out  [name]            output file with coordinates of predicted protein coding genes.

    Output options:

           --nt  [name]           output file with nucleotide sequences of predicted genes in FASTA format.
           --aa  [name]           output file with protein sequences of predicted genes in FASTA format.
           --format  [gtf]        format of output file with gene coordinates: gtf or gff3.
           --clean                delete temporay files

    Other parameters:
          --verbose
