# ChIP-SEQ PIPELINE

Author: Nia Kyritsis, Zhou Du, et al. @ Boston Children's Hospital / Harvard Medical School


### prerequisite

* bowtie2 (2.4.2)
* samtools (0.1.19)
* macs2 (2.2.5)
* UCSC binary utilities (http://hgdownload.soe.ucsc.edu/admin/exe/)
  * bedGraphToBigWig (v4)


Pipeline used to analyze ChIP-seq libraries; runs very similarly to the GRO-seq pipeline, except that it uses MACS2 to peak characterization and analysis. MACS2 analysis can also take an optional control library, which is sometimes also called the “input” library without IP by lab members. If someone does have a control library, this should be the first library run through the pipeline, and the control bam file is then given to the following commands.

Running with no control or running the control first:
```bash
python3 Paired_Analysis.py -f1 read1.fq.gz -f2 read2.fq.gz -d analysis -o name_of_output_files -g /path_to_bowtie2files --chromInfo /path_to_chromsizes.txt
```

With control bam file:
```bash
python3 Paired_Analysis.py -f1 read1.fq.gz -f2 read2.fq.gz -d analysis -c control_file.bam -o name_of_output_files -g /path_to_bowtie2files --chromInfo /path_to_chromsizes.txt
```

