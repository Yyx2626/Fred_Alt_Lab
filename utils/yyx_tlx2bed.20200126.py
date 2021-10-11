
#### Usage: cat <input.tlx> | python this.py <which_part>
#### Options:
####    <which_part> can be: junction (default) | bait | prey
#### Output: STDOUT   bed format 6 columns

import sys


which_part = 'junction'
if len(sys.argv) > 1:
	which_part = sys.argv[1]
if which_part not in ['junction', 'bait', 'prey']:
	print('Error: cannot recognize which_part="{}", should be one of "junction"|"bait"|"prey"'.format(which_part), file=sys.stderr)
	sys.exit(-1)


def fields_to_field2idx(fields):
	ans = {}
	for i, x in enumerate(fields):
		ans[x] = i
	return ans


is_headline = True
fields = []
field2idx = {}
rowidx = -1
for line in sys.stdin:
	rowidx += 1
	line = line.rstrip('\r\n')
	F = line.split('\t')
	if is_headline:
		is_headline = False
		fields = F
		field2idx = fields_to_field2idx(fields)
		continue
	else:
		name = F[field2idx['Qname']]
		score = str(rowidx)
		strand = F[field2idx['Strand']]
		if strand != '':
			strand = '+' if ( int(strand) > 0 ) else '-'
		chr = F[field2idx['Rname']]
		if which_part == 'junction':
			pos = F[field2idx['Junction']]
			if pos == '':
				continue
			pos = int(pos)
			start = str(pos-1)
			end = str(pos)
		elif which_part == 'prey':
			pos = F[field2idx['Rstart']]
			if pos == '':
				continue
			pos = int(pos)
			start = str(pos-1)
			end = F[field2idx['Rend']]
		elif which_part == 'bait':
			chr = F[field2idx['B_Rname']]
			pos = F[field2idx['B_Rstart']]
			if pos == '':
				continue
			pos = int(pos)
			start = str(pos-1)
			end = F[field2idx['B_Rend']]
			strand = '+' if ( int(F[field2idx['Strand']]) > 0 ) else '-'
		print('\t'.join([chr, start, end, name, score, strand]))
