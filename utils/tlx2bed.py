"""
Convert TLX file to BED/BEDGRAPH file
"""

import sys, os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Convert TLX file to BED or BEDGRAPH format')
    parser.add_argument("-f", dest = "tlxfile", type = str, required = True,
                              help = "TLX file" )
    parser.add_argument("-o", dest = "outprefix", type = str,
                              help = "Output file prefix. Default: TLX file prefix" )
    parser.add_argument("-t", dest = "type", type = str, default = "both", choices={'bed', 'bedgraph', 'both'},
                              help = "Output files: BED, BEDGRAPH, BOTH. Default: both" )
    parser.add_argument("-g", dest = "genome", type = str, default = "mouse", choices={'mouse', 'human'},
                              help = "specify the genome: mouse or human" )
    parser.add_argument("--v3", action='store_true',
                              help = "The tlx file is generated by HTGTS pipeline V3. Default: False" )
    args = parser.parse_args()
    return args

def filenames(args):
    # convert the windows line ending to unix line ending
    os.system("perl -pi -e 's/\r\n|\n|\r/\n/g' %s" % args.tlxfile)
    # parse output filenames 
    if args.outprefix:
        output_prefix = args.outprefix
    else:
        output_prefix = args.tlxfile.replace(".tlx", "")
    return output_prefix

def parseFILE(args):
    # Parse filenames
    output_prefix = filenames(args)
    bedtmp = output_prefix + '.tmp.bed'
    bedname = output_prefix + '.bed'
    bedgraphname = output_prefix + '.bedgraph'
    bedgraphname_pos = output_prefix + '.pos.bedgraph'
    bedgraphname_neg = output_prefix + '.neg.bedgraph'

    # convert to sorted BED format
    strandconvt = {"1": "+", "-1": "-"}
    bedfile = open(bedtmp, 'w')
    for line in open(args.tlxfile):
        if line.startswith('Qname'): continue
        l = line.strip().split('\t')
        if len(l) < 5: continue
        if l[2] == 'Adapter': continue
        if args.v3:
            try:
                bedfile.write("%s\t%d\t%s\t%s\t.\t%s\n" % (l[2], int(l[3])-1, l[3], l[0], strandconvt[l[4]]))
            except:
                pass
        else:
            bedfile.write("%s\t%d\t%s\t%s\t.\t%s\n" % (l[1], int(l[2])-1, l[2], l[0], strandconvt[l[3]]))
    bedfile.close()
    os.system('bedtools sort -i %s > %s' % (bedtmp, bedname))

    # convert combined bedgraph file
    if args.genome == 'mouse':
        genomefile = '/path/to/chromInfo_mm9.txt'
    elif args.genome == 'human':
        genomefile = '/path/to/chromInfo_hg19.txt'
    fout = open("%s" % bedgraphname, 'w')
    fout.write('track type=bedGraph name="%s" visibility=full color=70,0,199\n' % bedgraphname)
    fout.close()
    os.system('bedtools genomecov -bg -i %s -g %s >> %s' % (bedname, genomefile, bedgraphname))
    os.system('gzip %s' % (bedgraphname))
    # for positive strand
    fout = open("%s" % bedgraphname_pos, 'w')
    fout.write('track type=bedGraph name="%s" visibility=full color=199,0,30\n' % bedgraphname_pos)
    fout.close()
    os.system('bedtools genomecov -strand + -bg -i %s -g %s >> %s' % (bedname, genomefile, bedgraphname_pos))
    os.system('gzip %s' % (bedgraphname_pos))
    # for the minus strand
    fout = open("%s" % bedgraphname_neg, 'w')
    fout.write('track type=bedGraph name="%s" visibility=full color=0,30,200\n' % bedgraphname_neg)
    fout.close()
    os.system('bedtools genomecov -strand - -scale -1 -bg -i %s -g %s >> %s' % (bedname, genomefile, bedgraphname_neg))
    os.system('gzip %s' % (bedgraphname_neg))

    # delete files
    os.system('rm -rf %s' % (bedtmp))
    if args.type == 'bed':
        os.system('rm -rf %s.gz %s.gz %s.gz' % (bedgraphname, bedgraphname_pos, bedgraphname_neg))
    if args.type == 'bedgraph':
        os.system('rm -rf %s' % bedname)

def main():
    args = parse_args()
    parseFILE(args)

main()

