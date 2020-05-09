
#### Usage: cat <input.tlx> | python3 this.py
#### Input: STDIN   <input.tlx>
####     should have at least these columns:
####       B_Qstart, B_Qend, Qstart, Qend, Seq
#### Output: STDOUT   append 5 columns
####     pre  bait  mid  prey  post


import sys


is_headline = True
fields = []
field2idx = {}
for line in sys.stdin:
	line = line.rstrip('\r\n')
	F = line.split('\t')
	if is_headline:
		is_headline = False
		fields = F
		for i,x in enumerate(F):
			field2idx[x] = i
		for now_colname in ['B_Qstart', 'B_Qend', 'Qstart', 'Qend', 'Seq']:
			if now_colname not in field2idx:
				print('Error: cannot find colname "{}"'.format(now_colname), file=sys.stderr)
				sys.exit(-1)
		output_fields = fields
		output_fields.extend(['pre', 'bait', 'mid', 'prey', 'post'])
		print('\t'.join(output_fields))
		continue
	else:
		seq = F[field2idx['Seq']]
		Bstart  = int(F[field2idx['B_Qstart']])
		Bend    = int(F[field2idx['B_Qend']])
		Pstart  = int(F[field2idx['Qstart']])
		Pend    = int(F[field2idx['Qend']])
		bait = seq[(Bstart-1):Bend]
		prey = seq[(Pstart-1):Pend]
#		if Bstart < Pstart:
		pre = seq[:(Bstart-1)]
		mid = seq[Bend:(Pstart-1)]
		post = seq[Pend:]
		if Pstart < Bstart:
			print('Warning: prey start < bait start ...', file=sys.stderr)
		F.extend([pre, bait, mid, prey, post])
		print('\t'.join(F))
		
