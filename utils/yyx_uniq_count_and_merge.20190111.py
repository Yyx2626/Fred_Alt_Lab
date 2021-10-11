#!/usr/bin/env python3
##### Usage: python this.py <empty_fill> <should_split_output> <output_prefix> <file1> <file2> [file3 ...]
#####   <file1> can be 'filename' or 'filename:key_col_idx' or 'filename:col1-col2,col3,...' or 'filename:col:fieldname_prefix'
#####   (should have headline, columns separated by '\t', key_col_idx is 1-based)
##### Output:
#####   <output_prefix>.log  =  STDERR
#####   <output_prefix>.tsv  (separated by '\t')
#####      shared_field(s) , count_in_file_1, append_fields_in_file_1, count_in_file_2, append_fields_in_file_2, ...
#####   <output_prefix>.1=???.tsv, <output_prefix>.1x2.tsv, <output_prefix>.1x2x3.tsv ...  if <should_split_output> = True


import sys
import re

empty_line_pattern = re.compile('^\s*$')
all_number_pattern = re.compile('^[0-9]+$')
number_to_number_pattern = re.compile('^([0-9]+)-([0-9]+)$')

def yyx_parse_number_arg(str):
	ans_vec = str.split(';')
	for i in range(len(ans_vec)):
		now_str = ans_vec[i]
		now_str_vec = now_str.split(',')
		now_vec = []
		for now_str in now_str_vec:
			s1 = all_number_pattern.search(now_str)
			s2 = number_to_number_pattern.search(now_str)
			if s1:
				now_vec.append(int(now_str))
			elif s2:
				now_vec.extend(range(int(s2.group(1)), int(s2.group(2))+1))
#			elif now_str in retain:
#				now_vec.append(now_str)
			else:
				print("Warning: cannot parse '" + now_str + "', which will be skipped")
		ans_vec[i] = now_vec
	return ans_vec


key_list = []
key2count = {}
# key2count[key][file_idx] = count_in_file
key2otherContents = {}
# key2otherContents[key][file_idx] = other_fields
fileIdx2otherNC = {}
output_fieldnames = []

def read_and_store_input(filename, key_col_idxes, file_idx, fieldname_prefix, empty_fill):
	with open(filename, 'r') as fin:
		is_headline = True
		fieldnames = []
		NC = 0
		rowidx = 0
		for line in fin:
			rowidx += 1
			line = line.rstrip()
			if empty_line_pattern.search(line):
				continue
			F = line.split('\t')
			if is_headline:
				is_headline = False
				fieldnames = F
				NC = len(fieldnames)
				if file_idx == 1:
					Fk = (F[i] for i in key_col_idxes)
					output_fieldnames.extend(Fk)
				output_fieldnames.append('count_in_file_{}'.format(file_idx))
				Fo = (F[i] for i in range(NC) if i not in key_col_idxes)
				other_fieldnames = [fieldname_prefix + x for x in Fo]
				output_fieldnames.extend(other_fieldnames)
				fileIdx2otherNC[file_idx] = len(other_fieldnames)
			else:
				if len(F) < NC:
					print('Warning: not enough columns (NC={}) for Line {} in file {} : "{}"'.format(NC, rowidx, filename, line), file=sys.stderr)
					for i in range(len(F), NC):
						F.append(empty_fill)
				for i in range(len(F)):
					if F[i] == '':
						F[i] = empty_fill
				Fk = (F[i] for i in key_col_idxes)
				key = '\t'.join(Fk)
				if key not in key2count:
					key_list.append(key)
					key2count[key] = {}
					key2otherContents[key] = {}
				if file_idx not in key2count[key]:
					key2count[key][file_idx] = 0
				key2count[key][file_idx] += 1
				if file_idx not in key2otherContents[key]:
					Fo = (F[i] if i < len(F) else '' for i in range(NC) if i not in key_col_idxes)
					key2otherContents[key][file_idx] = '\t'.join(Fo)

def output(empty_fill, should_split_output, output_prefix):
	global output_fieldnames, key_list, file_num, key2count, key2otherContents, fileIdx2otherNC, file_labels
	
	def bool_vec_to_x_str(bool_vec):
		is_first = True
		ans = ''
		for i in range(len(bool_vec)):
			if bool_vec[i] > 0:
				if is_first:
					is_first = False
				else:
					ans += 'x'
				ans += str(i+1)
		if ans == '':
			ans = '0'
		return ans
	
	def bool_str_to_x_str(bool_str):
		return bool_vec_to_x_str(list(map(int, bool_str.split('\t'))))
		
	count_hash = {}
	bool_hash = {}
	for key in key_list:
		count_vec = []
		bool_vec = []
		for file_idx in range(1, file_num+1):
			if file_idx in key2count[key]:
				count_vec.append(key2count[key][file_idx])
				bool_vec.append(1)
			else:
				count_vec.append(0)
				bool_vec.append(0)
		# for summary
		count_str = '\t'.join([str(x) for x in count_vec])
		bool_str = '\t'.join([str(x) for x in bool_vec])
		if count_str not in count_hash:
			count_hash[count_str] = 0
		count_hash[count_str] += 1
		if bool_str not in bool_hash:
			bool_hash[bool_str] = 0
		bool_hash[bool_str] += 1
	summary_str = 'Summary:\n'
	summary_str += '\n'
	summary_str += 'Bool\n'
	summary_str += 'row_count\t'
	summary_str += '\t'.join(['{}={}'.format(x, file_labels[x-1]) for x in range(1, file_num+1)])
	summary_str += '\n'
	for kv in sorted(bool_hash.items(), key=lambda kv: kv[1]):
		summary_str += '{}\t{}\n'.format(kv[1], kv[0])
	summary_str += '\n'
	summary_str += 'Frequence\n'
	summary_str += 'row_count\t'
	summary_str += '\t'.join(['{}={}'.format(x, file_labels[x-1]) for x in range(1, file_num+1)])
	summary_str += '\n'
	for kv in sorted(count_hash.items(), key=lambda kv: kv[1]):
		summary_str += '{}\t{}\n'.format(kv[1], kv[0])
	
	log_filename = output_prefix + '.log'
	with open(log_filename, 'w') as flog:
		print(summary_str, file=flog)
	print(summary_str, file=sys.stderr)
	
	out_file_handles = {}
	with open(output_prefix + '.tsv', 'w') as fout:
		print('\t'.join(output_fieldnames), file=fout)
		if should_split_output:
			for k in bool_hash.keys():
				x_str = bool_str_to_x_str(k)
				if 'x' not in x_str:
					x_str = x_str + '=' + file_labels[int(x_str)-1]
				out_file_handles[k] = open(output_prefix + '.' + x_str + '.tsv', 'w')
				print('\t'.join(output_fieldnames), file=out_file_handles[k])
		
		try:
			for key in key_list:
				row = [key]
				bool_vec = []
				for file_idx in range(1, file_num+1):
					if file_idx in key2count[key]:
						row.append(str(key2count[key][file_idx]))
						row.append(key2otherContents[key][file_idx])
						bool_vec.append(1)
					else:
						row.append('0')
						row.append('\t'.join([empty_fill] * fileIdx2otherNC[file_idx]))
						bool_vec.append(0)
				print('\t'.join(row), file=fout)
				if should_split_output:
					bool_str = '\t'.join([str(x) for x in bool_vec])
#					k = bool_vec_to_x_str(bool_vec)
					print('\t'.join(row), file=out_file_handles[bool_str])
		finally:
			if should_split_output:
				for k in bool_hash.keys():
					out_file_handles[k].close()


def str2bool(x):
	try:
		x = int(x)
		return x > 0
	except ValueError:
		return x.lower() in ('true', 't', 'yes', 'y')


if len(sys.argv) <= 4:
	print('No file is input. Exit.', file=sys.stderr)
	sys.exit(-1)
	
empty_fill, should_split_output, output_prefix = sys.argv[1:4]
should_split_output = str2bool(should_split_output)

file_idx = 0
prefix_slash_pattern = re.compile('^.*/')
post_dot_pattern = re.compile('[.].*?$')

file_labels = []
for now_param in sys.argv[4:]:
	file_idx += 1
	P = now_param.split(':')
	filename = P[0]
	key_col_idxes = 1
	if len(P) > 1:
		key_col_idxes = [x-1 for x in yyx_parse_number_arg(P[1])[0]]
	fieldname_prefix = ''
	if len(P) > 2:
		fieldname_prefix = P[2]
		tmp = prefix_slash_pattern.sub('', fieldname_prefix)
		file_labels.append(tmp)
	else:
		tmp = prefix_slash_pattern.sub('', filename)
		tmp = post_dot_pattern.sub('', tmp)
		file_labels.append(tmp)
	
	print('# Now process {}-th file {} with field_prefix={} ...'.format(file_idx, filename, fieldname_prefix), file=sys.stderr)
	read_and_store_input(filename, key_col_idxes, file_idx, fieldname_prefix, empty_fill)
	
file_num = file_idx
print('# Now output ...', file=sys.stderr)
output(empty_fill, should_split_output, output_prefix)
print('# Successfully exit.', file=sys.stderr)




