import os, sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='ChIP Analysis script. Calls Bowtie2 and Macs2 to generate BW files of the peaks')
    parser.add_argument("-f1", type = str, required = True, help = "Read file 1" )
    parser.add_argument("-f2", type = str, required = False, help = "Read file 2" )
    parser.add_argument("-d", dest = "firstdirectory", type = str, required = True, help = "1st directory: outputs will be placed here" )
    parser.add_argument("-o", dest = "output", type = str, required = True, help = "Output file name for the R1 file" )
    parser.add_argument("-g", dest = "genome", type = str, required = True, help = "Genome used in the ChIP analysis" )
    parser.add_argument("-m", dest = "module", default="both" ,type = str, required = False, help = "Picks which modules to run {alignment, bigwig, both}. DEFAULT: both" )
    parser.add_argument("--chromInfo", dest = "chromInfo", type = str, required = True, help = "chromInfo file" )
    parser.add_argument("-c", dest = "control", type = str, required = False, help = "bam file for control sample")

    args = parser.parse_args()
    return args

def read_concat(args):
    index_1 = args.f1.find("L001")
    cat_R1 = args.f1.replace(args.f1[index_1:],"L00*_R1*")
    index_2 = args.f2.find("L001")
    cat_R2 = args.f2.replace(args.f2[index_2:],"L00*_R2*")
    args.f1 = args.f1.replace("L001_","")
    args.f2 = args.f2.replace("L001_","")

    os.system("cat {0}> {1}".format(cat_R1,args.f1))
    os.system("cat {0}> {1}".format(cat_R2,args.f2))
    os.system("rm -f {0}".format(cat_R1[:-4]))
    return

def fq_to_bam(args):
    if args.f2:
        os.system("bowtie2 -x {4} -1 {0} -2 {1} -p 6 --non-deterministic -S {2}/{3}.sam".format(args.f1,args.f2,args.firstdirectory,args.output,args.genome))
    else:
        os.system("bowtie2 -x {3} -f -U {0} -p 6 --non-deterministic -S {1}/{2}.sam".format(args.f1,args.firstdirectory,args.output,args.genome))
    os.system("samtools view -bhS {0}/{1}.sam -o {0}/{1}.bam".format(args.firstdirectory, args.output))
    os.system("samtools sort {0}/{1}.bam {0}/{1}.sort".format(args.firstdirectory, args.output))
    os.system("mv {0}/{1}.sort.bam {0}/{1}.bam".format(args.firstdirectory, args.output))
    os.system("samtools index {0}/{1}.bam".format(args.firstdirectory, args.output))
    #os.system("rm -f {0}/{1}.sam".format(args.firstdirectory, args.output))


def wig_convert(args):
    if args.genome == "mm9" or args.genome == "hg19":
        os.system("python3 bam2bw.py -b {1}/{0}.bam -g {2} -s F".format(args.output,args.firstdirectory,args.chromInfo))
    else:
        if args.control:
            os.system("macs2 callpeak -t {0}/{2}.bam -c {3} -n {1}/{2} --nomodel --keep-dup all --extsize 51 --nolambda -B --SPMR -g mm --verbose 0 --broad".format(args.firstdirectory, args.firstdirectory, args.output, args.control))
        else:
            #os.system("python3 bam2bw.py -b {1}/{0}.bam -g {2} -s F".format(args.output,args.firstdirectory,args.chromInfo))
            os.system("macs2 callpeak -t {0}/{2}.bam -n {1}/{2} --nomodel --keep-dup all --extsize 51 --nolambda -B --SPMR -g mm --verbose 0 --broad".format(args.firstdirectory, args.firstdirectory, args.output))
        #os.system("bedtools sort {0}/{1}_treat_pileup.bdg {0}/{1}_treat_pileup.bdg".format(args.firstdirectory, args.output)) #####
        os.system("bedGraphToBigWig {0}/{1}_treat_pileup.bdg {2} {0}/{1}_treat_pileup.bw".format(args.firstdirectory, args.output, args.chromInfo))
        os.system("samtools view -b -F 0x10 {0}/{1}.bam -o {0}/{1}.pos.bam".format(args.firstdirectory, args.output))
        os.system("samtools view -b -f 0x10 {0}/{1}.bam -o {0}/{1}.neg.bam".format(args.firstdirectory, args.output))
        os.system("macs2 callpeak -t {0}/{2}.pos.bam -n {1}/{2}.pos --nomode --keep-dup all --extsize 51 --nolambda -B --SPMR -g mm --verbose 0 --broad".format(args.firstdirectory, args.firstdirectory, args.output))
        os.system("macs2 callpeak -t {0}/{2}.neg.bam -n {1}/{2}.neg --nomode --keep-dup all --extsize 51 --nolambda -B --SPMR -g mm --verbose 0 --broad".format(args.firstdirectory, args.firstdirectory, args.output))
        #os.system("bedtools sort {0}/{1}.pos_treat_pileup.bdg {0}/{1}.pos_treat_pileup.bdg".format(args.firstdirectory, args.output))
        #os.system("bedtools sort {0}/{1}.neg_treat_pileup.bdg {0}/{1}.neg_treat_pileup.bdg".format(args.firstdirectory, args.output))
        os.system("bedGraphToBigWig {0}/{1}.pos_treat_pileup.bdg {2} {0}/{1}.treat_pileup_pos.bw".format(args.firstdirectory, args.output, args.chromInfo))
        os.system("bedGraphToBigWig {0}/{1}.neg_treat_pileup.bdg {2} {0}/{1}.treat_pileup_neg.bw".format(args.firstdirectory, args.output, args.chromInfo))

def GROSeqPL(args):
    os.system("mkdir -p {0}".format(args.firstdirectory))
    #if (os.path.isfile(args.f1) and "L001" in args.f1):
    #    read_concat(args)
    if args.module == "both":
        fq_to_bam(args)
        wig_convert(args)
    if args.module == "alignment":
        fq_to_bam(args)
    if args.module == "bigwig":
        wig_convert(args)
    os.system("rm -f {0}/{1}.sam".format(args.firstdirectory, args.output))

def main():
    args = parse_args()
    GROSeqPL(args)

main()
