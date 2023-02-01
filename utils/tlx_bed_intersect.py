import sys
import os

###python3 tlx_bed_intersect.py XFC1_AltN118_new.tlx NC_000078_anno_mm9_+1_imgt.bed extract_comp

def main():

    tlx = sys.argv[1]
    bed = sys.argv[2]
    operation = sys.argv[3] ###extract, extract_comp, count

    if operation == "count":
        tlxbed = tlx[:-3] + "bed"
        os.system("/storage/researchers/yiwen/tlx2BED.pl {0} {1}".format(tlx, tlxbed))
        outputbed = tlx[:-4] + "_count.bed"
        os.system("bedtools intersect -c -a {0} -b {1} > {2}".format(bed, tlxbed, outputbed))
    
    elif operation == "extract":
        tlxbed = tlx[:-3] + "bed"
        os.system("/storage/researchers/yiwen/tlx2BED.pl {0} {1}".format(tlx, tlxbed))
        outputbed = tlx[:-4] + "_intersect.bed"
        os.system("bedtools intersect -u -a {0} -b {1} > {2}".format(tlxbed, bed, outputbed))
        outputtlx = tlx[:-4] + "_intersect.tlx"
        os.system("/storage/researchers/yiwen/pullTLXFromBED.pl {0} {1} {2}".format(tlx, outputbed, outputtlx))

    elif operation == "extract_comp":
        tlxbed = tlx[:-3] + "bed"
        os.system("/storage/researchers/yiwen/tlx2BED.pl {0} {1}".format(tlx, tlxbed))
        outputbed = tlx[:-4] + "_intersect.bed"
        os.system("bedtools intersect -v -a {0} -b {1} > {2}".format(tlxbed, bed, outputbed))
        outputtlx = tlx[:-4] + "_intersect.tlx"
        os.system("/storage/researchers/yiwen/pullTLXFromBED.pl {0} {1} {2}".format(tlx, outputbed, outputtlx))

main()
