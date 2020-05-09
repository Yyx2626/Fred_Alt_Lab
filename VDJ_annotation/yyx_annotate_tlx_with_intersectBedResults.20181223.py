
#### Usage: python3 this.py <input.tlx> <anno1.bed> [anno2.bed] ...

import sys


input_tlx_filename = sys.argv[1]


### read input.tlx into memory
is_headline = True
original_NC = -1
storage = []
rowidx = -1
with open(input_tlx_filename, 'r') as fin:
	for line in fin:
		rowidx += 1
		line = line.rstrip('\r\n')
		F = line.split('\t')
		if is_headline:
			is_headline = False
			original_NC = len(F)
		if len(storage) != rowidx:
			print('Error: something goes wrong', file=sys.stderr)
			sys.exit(-1)
		storage.append(F)
original_NR = len(storage)
#print(storage)   # debug


### parse each annotation file
anno_idx = 0
for anno_filename in sys.argv[2:]:
	anno_idx += 1
	print('Now parse anno{} file {}'.format(anno_idx, anno_filename), file=sys.stderr)
	storage_idx = 0
	storage[storage_idx].append('anno{}_overlap_bps'.format(anno_idx))
	storage[storage_idx].append('anno{}_overlap_features'.format(anno_idx))
	with open(anno_filename, 'r') as fin:
		for line in fin:
			line = line.rstrip('\r\n')
			F = line.split('\t')
			name = F[3]
			rowidx = int(F[4])
			overlap_bp = F[len(F)-1]
			overlap_feature = F[9]
			
			# for each input tlx row, I will append 2 columns for each annotation file: overlap_feature_bps, overlap_feature_names
			# so anno_idx=1 -> output original_NC+2*1 columns; anno_idx=2 -> output original_NC+2*2 columns
			while storage_idx < rowidx:
				storage_idx += 1
				if storage_idx >= original_NR:
					print('Warning: anno{} file {} has more rows than input tlx file, I will skip them'.format(anno_idx, anno_filename), file=sys.stderr)
					break
				storage[storage_idx].append('0')
				storage[storage_idx].append('.')
#				print(len(storage[storage_idx]), file=sys.stderr)   # debug
			if storage_idx >= original_NR:
				break
			if name != storage[storage_idx][0]:
				print('Warning: rowidx={} name {} in anno{} not match storage_idx={} name {} in input tlx ?!'.format(rowidx, name, anno_idx, storage_idx, storage[storage_idx][0]), file=sys.stderr)
			if overlap_bp != '0':
				if storage[storage_idx][original_NC+2*anno_idx-2] == '0':
					storage[storage_idx][original_NC+2*anno_idx-2] = overlap_bp
					storage[storage_idx][original_NC+2*anno_idx-1] = overlap_feature
				else:
					storage[storage_idx][original_NC+2*anno_idx-2] += ',' + overlap_bp
					storage[storage_idx][original_NC+2*anno_idx-1] += ',' + overlap_feature


### output to STDOUT
for record in storage:
	print('\t'.join(record))
	