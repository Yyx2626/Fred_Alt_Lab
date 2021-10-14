# GRO_SEQ PIPELINE

Author: Nia Kyritsis, Jianqiao Hu, Zhou Du, et al. @ Boston Children's Hospital / Harvard Medical School


### prerequisite

* python (3.7.4)
* bowtie2 (2.4.2)
* samtools (0.1.19)
* homer (4.9.1)
* bedtools (2.29.0)
* gfold
* bx-python
* RSeQC
  * bam2wig.py (3.0.1) by Liguo Wang
* UCSC binary utilities (http://hgdownload.soe.ucsc.edu/admin/exe/)
  * wigToBigWig (v4)


This pipeline is a series of smaller scripts/packages to analyze GRO-seq libraries. It uses the program homer to identify transcripts, and the RSeQC package to convert to bw files for viewing in IGV. There are several versions of this pipeline, changed over the years for several members. The -o parameter will the file names of all output files with the relevant extensions. The pipeline makes an /alignment folder in the directory the command is run in, and all output files are placed here.

For single-end or paired-end to be treated as single-end (if given two files, R2 is reverse complemented and it is treated like single end data):
```bash
python3 GROSeqPL_SE_R2opt.py -f1 read1.fq.gz -f2 read2.fq.gz -o outputnames -cg /path_to_bowtie2files --chromInfo /path_to_chromosomesizes.txt 
```

For single end, no -f1 or -f2, -f1 becomes -f:
```bash
python3 GROSeqPL_SE_R2opt.py -f read1.fq.gz -o outputnames -cg /path_to_bowtie2files --chromInfo /path_to_chromosomesizes.txt
```

For joining R1 and R2 normal paired-end:
```bash
python3 GROSeqPL_joined.py -f1 read1.fq.gz -f2 read2.fq.gz -o outputnames -cg /path_to_bowtie2files --chromInfo /path_to_chromosomesizes.txt
```
