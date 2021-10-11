import sys, os
import pandas as pd
import argparse

def main():
	parser = argparse.ArgumentParser(description='''
Downsampling normalization for tlx file
''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', metavar='input.tlx', help='input tlx file', required=True)
	parser.add_argument('-n', '--normalize', metavar='INT', help='normalize to ? total junctions', default=1e6)
	parser.add_argument('-s', '--seed', metavar='INT', help='random seed (for reproducibility)', default=1234567)
	parser.add_argument('-o', '--output', metavar='output.tlx', help='output tlx file', default='(STDOUT)')
	
	args = parser.parse_args()
	if args.output == '(STDOUT)':
		args.output = sys.stdout
	
	records = pd.read_csv(args.input, sep='\t')
	#if sys.argv[2] == 'GC':
	#    records = records[records['V_MUTATION']>0]
	#elif sys.argv[2] == 'nonGC':
	#    records = records[records['V_MUTATION']==0]
	records = records.sample(int(args.normalize), random_state=int(args.seed))
	records.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
	main()

