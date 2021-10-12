import sys, os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='GRO-Seq Pipeline')
    parser.add_argument("-f", dest = "file", type = str, required = False, help = "Read input file" )
    parser.add_argument("-f1", dest = "file_1", type = str, required = False, help = "Read input file 1 of pair" )
    parser.add_argument("-f2", dest = "file_2", type = str, required = False, help = "Read input file 2 of pair" )
    parser.add_argument("-o", dest = "output", type = str, required = True, help = "Output file name" )
    parser.add_argument("-g", dest = "genome", type = str, required = False, help = "Genome file to align to (do not include .fa)" )
    parser.add_argument("-cg", dest = "custom_genome", type = str, required = False, help = "Custom Genome file to align to (includes .fa)" )
    parser.add_argument("--chromInfo", dest = "chromInfo", type = str, required = True, help = "chromInfo file" )
    args = parser.parse_args()
    return args

def read_concat(args):
    index_1 = args.file_1.find("L001")
    cat_R1 = args.file_1.replace(args.file_1[index_1:],"L00*_R1*")
    index_2 = args.file_2.find("L001")
    cat_R2 = args.file_2.replace(args.file_2[index_2:],"L00*_R2*")
    args.file_1 = args.file_1.replace("L001_","")
    args.file_2 = args.file_2.replace("L001_","")

    os.system("cat {0}> {1}".format(cat_R1,args.file_1))
    os.system("cat {0}> {1}".format(cat_R2,args.file_2))
    os.system("rm -f {0}".format(cat_R1[:-4]))
    return

def reverse_r2(args):

    #os.system("gunzip {0}".format(args.file_2))
    #unzipped2 = args.file_2[:-3]
    rev2file = "rev_" + args.file_2[:-3]
    os.system("seqtk seq -r {0} > {1}".format(args.file_2, rev2file))
    new_r2 = rev2file + ".gz"
    os.system("gzip {0}".format(rev2file))
    args.file_2 = new_r2
    return

def fq_to_bam(args):
    if args.custom_genome:
        # path, folder = os.path.split(args.custom_genome)
        args.genome = args.custom_genome

    if args.file_1:
        os.system("bowtie2 -x {2} -U {0},{3} -p 6 --non-deterministic -S alignment/{1}.sam".format(args.file_1, args.output, args.genome, args.file_2))
    else:
        os.system("bowtie2 -x {2} -U {0} -p 6 --non-deterministic -S alignment/{1}.sam".format(args.file, args.output, args.genome))

    os.system("samtools view -bhS alignment/{0}.sam -o alignment/{0}.bam".format(args.output))
    os.system("samtools sort alignment/{0}.bam alignment/{0}.sort".format(args.output))
    os.system("mv alignment/{0}.sort.bam alignment/{0}.bam".format(args.output))
    os.system("samtools index alignment/{0}.bam alignment/{0}.bam.bai".format(args.output))
    #os.system("samtools stats {0}.bam > {1}.stats.txt".format(args.output, args.output))
    os.system("rm -rf alignment/{0}.sam".format(args.output))
    return

def wig_convert(args):
    os.system("python3 bam2bw.py -b alignment/{0}.bam -g {1} -s T".format(args.output,args.chromInfo))
    os.system("samtools view -b -F 0x10 alignment/{0}.bam -o alignment/{0}.pos.bam".format(args.output))
    os.system("samtools view -b -f 0x10 alignment/{0}.bam -o alignment/{0}.neg.bam".format(args.output))
    return

def homer(args):
    if args.custom_genome:
        os.system("makeTagDirectory conv_rpkm/{0} alignment/{0}.bam -genome {1}.fa".format(args.output,args.custom_genome))
    else:
        os.system("makeTagDirectory conv_rpkm/{0} alignment/{0}.bam -genome {1}".format(args.output,args.genome))
    # os.system("/usr/local/homer/bin/findPeaks conv_rpkm/{0} -style groseq -o auto -uniqmap /usr/local/homer/data/uniqmap/mm9-uniqmap".format(args.output))
    os.system("findPeaks conv_rpkm/{0} -style groseq -o auto".format(args.output))
    return

def convergent(args):
    os.system("grep -v \"#\" conv_rpkm/{0}/transcripts.txt | awk '{{OFS=\"\t\"}}{{print $2, $3, $4, $1, $6, $5}}' > conv_rpkm/{0}.transcript.bed".format(args.output))
    os.system("bedtools sort -i conv_rpkm/{0}.transcript.bed > see".format(args.output))
    os.system("mv see conv_rpkm/{0}.transcript.bed".format(args.output))
    os.system("grep -v \"+\" conv_rpkm/{0}.transcript.bed > conv_rpkm/{0}.transcripts.minus.bed".format(args.output))
    os.system("grep \"+\" conv_rpkm/{0}.transcript.bed > conv_rpkm/{0}.transcripts.plus.bed".format(args.output))
    os.system("bedtools multiinter -i conv_rpkm/{0}.transcripts.plus.bed conv_rpkm/{0}.transcripts.minus.bed |grep \"1,2\" |awk \'{{OFS=\"\\t\"}}{{ if ($3-$2>100){{print $1, $2, $3}} }}\' > conv_rpkm/{0}.transcripts.convergent.bed".format(args.output))
    return

def gmean_cal(args):
    os.system("awk '{{OFS=\"\\t\"}}{{printf \"%s\\t%s\\t%s\\t%s_%s_%s\\t.\\t.\\n\", $1, $2, $3, $1, $2, $3}}' conv_rpkm/{0}.transcripts.convergent.bed > conv_rpkm/{0}.transcripts.convergent.6.bed".format(args.output))
    os.system("samtools view alignment/{0}.pos.bam | gfold count -annf BED -ann conv_rpkm/{0}.transcripts.convergent.6.bed -tag stdin -o conv_rpkm/{0}.pos.cnt".format(args.output))
    os.system("samtools view alignment/{0}.neg.bam | gfold count -annf BED -ann conv_rpkm/{0}.transcripts.convergent.6.bed -tag stdin -o conv_rpkm/{0}.neg.cnt".format(args.output))
    # os.system("samtools view alignment/{0}.pos.bam | gfold count -annf BED -ann mm9.transcript.longest.GROseq.bed -tag stdin -o conv_rpkm/{0}.pos3.cnt".format(args.output))
    # os.system("samtools view alignment/{0}.neg.bam | gfold count -annf BED -ann mm9.transcript.REV.longest.GROseq.bed -tag stdin -o conv_rpkm/{0}.neg3.cnt".format(args.output))
    os.system("python COVT.gmean_cal.py conv_rpkm/{0} > conv_rpkm/{0}.gmean.cnt".format(args.output))
    os.system("sort -gr -k 4 conv_rpkm/{0}.gmean.cnt > conv_rpkm/{0}.gmean.cnt.sort".format(args.output))
    os.system("python Add.col.py conv_rpkm/{0}.gmean.cnt.sort > conv_rpkm/{0}.gmean.bed".format(args.output))
    return

def GROSeqPL(args):
    os.system("mkdir -p alignment conv_rpkm")
    if (os.path.isfile(args.file_1) and "L001" in args.file_1):
        read_concat(args)
    if args.file_2:
        reverse_r2(args)
    fq_to_bam(args)
    wig_convert(args)
    homer(args)
    convergent(args)
    gmean_cal(args)

def main():
    args = parse_args()
    GROSeqPL(args)

main()
