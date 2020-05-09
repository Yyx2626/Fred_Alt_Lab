
#### Usage: cat <input.tlx> | python3 this.py <D.fa> [score_cutoff (default:5)]
#### Input: STDIN   <input.tlx>
####     should have at least these columns:
####       B_Qstart, B_Qend, Qstart, Qend, Seq
#### Output: STDOUT   append 5+2 columns
####     pre  bait  mid  prey  post   mid_D_score  mid_D_annotate


import sys
#import Bio.pairwise2
import Bio.SeqIO
#import Bio.Seq
#import Bio.Alphabet
import re


## 2018-12-23 turn back to use continuous_LCS_DP (LCS=longest continuous substring, not allow gap) instead of Bio.pairwise2.align.globalms

def continuous_LCS_DP(query, subjt, caseSensitive=False):
	import numpy as np
	if not caseSensitive:
		query = query.upper()
		subjt = subjt.upper()
		
	len1 = len(query)
	len2 = len(subjt)
	max_LCS_len = 0
	max_len_i, max_len_j = -1, -1
	DPmat = np.zeros((len1+1, len2+1), dtype='int')
	for i in range(1, len1+1):
		for j in range(1, len2+1):
			if query[i-1] == subjt[j-1]:
				DPmat[i,j] = DPmat[i-1,j-1] + 1
				if DPmat[i,j] > max_LCS_len:
					max_LCS_len = DPmat[i,j]
					max_len_i, max_len_j = i, j
	
	return max_LCS_len, max_len_i, max_len_j
	# max_match_shift = j - i


trailing_branket_pattern = re.compile('\(.*$')


D_ref_fa_filename = sys.argv[1]

score_cutoff = 5
if len(sys.argv) > 2:
	score_cutoff = int(sys.argv[2])


### read in D reference fasta file
D_ref_hash = {}
D_rC_hash = {}   # reverseComplement
for seq_record in Bio.SeqIO.parse(D_ref_fa_filename, 'fasta'):
	now_ID = seq_record.id
	now_ID = now_ID.replace('lcl|', '')
	D_ref_hash[now_ID] = seq_record.seq
	D_rC_hash[now_ID] = seq_record.seq.reverse_complement()


def best_align(query, ref_hash):
	max_score = -1
	max_score_hits = []
	for k, v in ref_hash.items():
#		alignments = Bio.pairwise2.align.globalms(v, query, 1, 0, -.5, -.1, penalize_end_gaps=[True, False], one_alignment_only=True)
#		score = -2
#		if len(alignments) > 0 and len(alignments[0]) > 2:
#			score = alignments[0][2]
		## 2018-12-23 turn back to use continuous_LCS_DP (LCS=longest continuous substring, not allow gap) instead of Bio.pairwise2.align.globalms
		alignment = continuous_LCS_DP(query, str(v))
		score = alignment[0]
		if score > max_score:
			max_score = score
			max_score_hits = [k]
		elif abs(score-max_score) < 1e-5:
			max_score_hits.append(k)
	return max_score, max_score_hits

Dalign_storage = {}
def Dalign(query):
	if query in Dalign_storage:
		return Dalign_storage[query]
	D_score, D_hit = best_align(query, D_ref_hash)
	D_rC_score, D_rC_hit = best_align(query, D_rC_hash)
	max_score = max([D_score, D_rC_score])
	max_score_hits = []
	if abs(D_score - max_score) < 1e-5:
		max_score_hits.extend(sorted(D_hit))
	if abs(D_rC_score - max_score) < 1e-5:
		max_score_hits.extend([x+' (rC)' for x in sorted(D_rC_hit)])
	Dalign_storage[query] = [max_score, max_score_hits]
	return max_score, max_score_hits


is_headline = True
fields = []
field2idx = {}
rowidx = 0
for line in sys.stdin:
	rowidx += 1
	if rowidx % 10000 == 0:
		print('Now process {}-th read ...'.format(rowidx), file=sys.stderr)
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
		output_fields.extend(['pre', 'bait', 'mid', 'prey', 'post', 'mid_D_score', 'mid_D_annotate'])
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
		
		## alignment
#		score, mid_annotate = Dalign(Bio.Seq.Seq(mid, Bio.Alphabet.IUPAC.unambiguous_dna))
		## 2018-12-23 turn back to use continuous_LCS_DP (LCS=longest continuous substring, not allow gap) instead of Bio.pairwise2.align.globalms
		score, mid_annotate = Dalign(mid)
		if score < score_cutoff:
			mid_annotate = '-'
			score = '-'
		else:
			score = str(score)
		
		F.extend([pre, bait, mid, prey, post, score, ', '.join(mid_annotate)])
		print('\t'.join(F))
		
