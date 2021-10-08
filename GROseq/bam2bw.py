"""
Convert bam file to bigwig file
"""

import sys, os
import argparse
import glob

def parse_args():
    parser = argparse.ArgumentParser(description='Convert BAM file to BIGWIG files')
    parser.add_argument("-b", dest = "bamfile", type = str, required = True,
                              help = "bam file" )
    parser.add_argument("-g", dest = "genome", type = str, default = "mm9", required = True,
                              help = "specify the genome: mm9, hg19 or abs path (mm9)" )
    parser.add_argument("-s", "--strandspecific", dest = "strandspecific", 
                              type = str, default = "F", choices={'F', 'T'},
                              help = "Divide into strand plus and minus bigwig files (F)" )
    parser.add_argument("-n", "--normalization", dest = "normalization", 
                              type = str, default = "T", choices={'F', 'T'},
                              help = "Do normalization of bigwig to coverage \
                              of 10 million 100nt reads (T)" )
    #parser.add_argument("-q", "MAP_QUAL", dest = "MAP_QUAL", type = int, default = 30,
    #                          help = "Minimum mapping quality for an alignment to be called 'uniquely mapped'. default=30" )
    args = parser.parse_args()
    return args

def parseFILE(args):
    if args.genome == 'mm9':
        genomefile = '/path/to/chromInfo_mm9.txt'
    elif args.genome == 'hg19':
        genomefile = '/path/to/chromInfo_hg19.txt'
    else:
        genomefile = args.genome
    normdef = ''
    if args.normalization == 'T':
        normdef = ' -t 1000000000 '

    if args.strandspecific == 'F':
        prefix = args.bamfile.replace('.bam', '')
        os.system('bam2wig.py -i %s -o %s -s %s %s' % (args.bamfile, prefix, genomefile, normdef))
        os.system('rm -rf %s.wig' % prefix)
    else:
        prefix_pos = args.bamfile.replace('.bam', '.pos')
        os.system('samtools view -F 0x10 -hb -o %s.bam %s' % (prefix_pos, args.bamfile))
        os.system('samtools sort %s.bam -o %s.sort > %s.sort.bam' % (prefix_pos, prefix_pos, prefix_pos))
        os.system('mv %s.sort.bam %s.bam' % (prefix_pos, prefix_pos))
        os.system('samtools index %s.bam' % prefix_pos)
        os.system('bam2wig.py -i %s.bam -o %s -s %s %s' % (prefix_pos, prefix_pos, genomefile, normdef))

        # Need to add minus symbol in minus bw file
        prefix_neg = args.bamfile.replace('.bam', '.neg')
        os.system('samtools view -f 0x10 -hb -o %s.bam %s' % (prefix_neg, args.bamfile))
        os.system('samtools sort %s.bam -o %s.sort > %s.sort.bam' % (prefix_neg, prefix_neg, prefix_neg))
        os.system('mv %s.sort.bam %s.bam' % (prefix_neg, prefix_neg))
        os.system('samtools index %s.bam' % prefix_neg)
        os.system('bam2wig.py -i %s.bam -o %s -s %s %s' % (prefix_neg, prefix_neg, genomefile, normdef))
        fout = open('%s.2.wig' % prefix_neg, 'w')
        for line in open('%s.wig' % prefix_neg):
            if line.startswith('variableStep'):
                fout.write(line)
            else:
                fout.write('%s\t-%s\n' % (line.split()[0], line.split()[1]))
        fout.close()
        os.system('wigToBigWig -clip %s.2.wig %s %s.bw' % (prefix_neg, genomefile, prefix_neg))

        #os.system('rm -rf %s.bam* %s.bam*' % (prefix_pos, prefix_neg))
        os.system('rm -rf %s.wig %s.wig' % (prefix_pos, prefix_neg))
        os.system('rm -rf %s.2.wig' % (prefix_neg))

def main():
    args = parse_args()
    parseFILE(args)

main()
