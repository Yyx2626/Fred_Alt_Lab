#!/usr/bin/env python3

import sys, os, os.path, subprocess, time
import argparse


### flushing print, reference: https://mail.python.org/pipermail/python-list/2015-November/698426.html
def _print(*args, **kwargs):
	file = kwargs.get('file', sys.stdout)
	print(*args, **kwargs)
	file.flush()


def check_elapsed_time(start_time):
	end_time = time.time()
	elapsed_time = end_time - start_time
	day = int(elapsed_time / (3600*24))
	hour = int(elapsed_time % (3600*24) / 3600)
	min = int(elapsed_time % 3600 / 60)
	sec = elapsed_time % 60
	elapsed_time = ''
	if day>0 : elapsed_time += '{}day '.format(day)
	if hour>0: elapsed_time += '{}h'.format(hour)
	if min>0 : elapsed_time += '{}min'.format(min)
	if sec>0 or elapsed_time == '': elapsed_time += '{:.2f}s'.format(sec)
	_print('[PYTHON-TIME] ' + elapsed_time, file=sys.stderr)


def check_is_empty_dir(dirname):
	if os.path.isdir(dirname):
		for nowDirName, subdirList, fileList in os.walk(dirname):
			if nowDirName == '.':
				return len(subfileList) == 0
			else:
				continue
	else:
		_print('Warning: Not a dir is attemped to be checked', file=sys.stderr)
		return None
		

def exist_file_or_dir(filenames, prompt_str, mode='any'):
	if mode not in ('any', 'all'):
		_print('Error: mode should be either "any" or "all", in exist_file_or_dir()', file=sys.stderr)
		return None
	is_mode_all = False
	if mode == 'any':
		is_mode_all = True
		
	if isinstance(filenames, str):
		filenames = [filenames]
	not_None_count = 0
	for filename in filenames:
		if filename is None:
			continue
		not_None_count += 1
		if os.path.isdir(filename):
			if not check_is_empty_dir(filename):
				_print('[CHECK-EXIST] Dir ' + filename + ' has already existed, and not empty. ' + prompt_str, file=sys.stderr)
				if not is_mode_all:
					return True
			else:
				if is_mode_all:
					return False
		elif os.path.isfile(filename) and os.path.getsize(filename) >= 100:
			# os.path.getsize(x) may also be os.stat(x).st_size
			_print('[CHECK-EXIST] File ' + filename + ' has already existed. ' + prompt_str, file=sys.stderr)
			if not is_mode_all:
				return True
		else:
			if is_mode_all:
				return False
	if not_None_count > 0:
		return is_mode_all
	return False
	

def check_final_file_then_remove_intermediate_file(final_filenames, intermediate_filenames, mode='all'):
	if mode not in ('any', 'all'):
		_print('Error: mode should be either "any" or "all", in check_final_file_then_remove_intermediate_file()', file=sys.stderr)
		return
		
	if isinstance(intermediate_filenames, str):
		intermediate_filenames = [intermediate_filenames]
	if exist_file_or_dir(final_filenames, 'So remove intermediate files...', mode=mode):
		for filename in intermediate_filenames:
			if filename is None:
				continue
			if os.path.exists(filename):
				_print('[PYTHON-REMOVE] ' + filename, file=sys.stderr)
				os.remove(filename)

def check_file_then_exec_command(filenames, command, should_time=False, error_stop=False, not_run=False):
	start_time = time.time()
	
	_print('[SUBPROCESS-CALL] ' + command, file=sys.stderr)
	if exist_file_or_dir(filenames, 'Skip this above command...', mode='any'):
		return
	
	if not not_run:
#		returnValue = os.system('/bin/bash -c ' + command)
		returnValue = subprocess.call(['/bin/bash', '-c', command])
		if returnValue != 0:
			if error_stop:
				_print('Error: when exec last command, return value = {}'.format(returnValue), file=sys.stderr)
				sys.exit(returnValue)
	
	if should_time:
		check_elapsed_time(start_time)


def stop_if_file_not_exist(filenames, mode='any'):
	if mode not in ('any', 'all'):
		_print('Error: mode should be either "any" or "all", in stop_if_file_not_exist()', file=sys.stderr)
		return None
	is_mode_all = False
	if mode == 'any':
		is_mode_all = True
		
	if isinstance(filenames, str):
		filenames = [filenames]
	checkFileNumber = 0
	missingFileNumber = 0
	for filename in filenames:
		if filename is None:
			continue
		checkFileNumber += 1
		if not os.path.isfile(filename):
			# os.path.getsize(x) may also be os.stat(x).st_size
			_print('[CHECK-EXIST] File ' + filename + ' does not exist.', file=sys.stderr)
			missingFileNumber += 1
		else:
			_print('[CHECK-EXIST] File ' + filename + ' exists.  Good.', file=sys.stderr)
	if missingFileNumber > 0:
		if not is_mode_all:
			_print('[STOP-NOT-EXIST] Error: requested {} file(s) is missing. Terminate!'.format(missingFileNumber), file=sys.stderr)
			sys.exit(missingFileNumber)
		elif missingFileNumber == checkFileNumber:
			_print('[STOP-NOT-EXIST] Error: requested file(s) is missing. Terminate!', file=sys.stderr)
			sys.exit(missingFileNumber)



def trim_bed(input_bedfile, output_bedfile, remain_length, towards_end='5', allow_elongate_length=True):
	print('Now trim ' + input_bedfile + ' to length={}'.format(remain_length) + ', output to ' + output_bedfile, file=sys.stderr)
	with open(input_bedfile, 'r') as fin:
		with open(output_bedfile, 'w') as fout:
			for line in fin:
				F = line.rstrip().split('\t')
				F[1] = int(F[1])
				F[2] = int(F[2])
				if remain_length > 0 and allow_elongate_length or F[2] - F[1] > remain_length:
					if F[5] == '-':
						if towards_end == '5':
							F[1] = F[2] - remain_length
						else:
							F[2] = F[1] + remain_length
					else:
						if towards_end == '5':
							F[2] = F[1] + remain_length
						else:
							F[1] = F[2] - remain_length
				F[1] = str(F[1])
				F[2] = str(F[2])
				print('\t'.join(F), file=fout)


def bed_to_genomecov_bdg(bedfile, output_bdgfile, genomeChrSizeFile, other_options=''):
	command = 'bedtools genomecov ' + other_options + ' -bg -i ' + bedfile + ' -g ' + genomeChrSizeFile + ' >' + output_bdgfile
	check_file_then_exec_command([output_bdgfile], command)


def bdg_extract_multiple(input_file, output_bdgfile, multiply_factor, chr='', start='', end=''):
	command = 'perl yyx_bdg_extract_multiply.20200120.pl ' + input_file + ' ' + output_bdgfile + ' ' + str(multiply_factor)
	if chr != '':
		command = command + ' ' + chr
		if start != '':
			command = command + ' ' + start
			if end != '':
				command = command + ' ' + end
	check_file_then_exec_command([output_bdgfile], command)


def bam_to_bw(bamfile, output_prefix, genomeChrSizeFile, trim, normalize_to, total_number, should_negative=False):
	bam2bedfile = output_prefix + '.bamtobed.bed'
	trim_bedfile = output_prefix + '.trim_{}.bed'.format(trim)
	genomecov_bdgfile = output_prefix + '.genomecov.bdg'
	
	## bam -> bed
	command = 'bedtools bamtobed -i ' + bamfile + ' | bedtools sort -i - >' + bam2bedfile
	check_file_then_exec_command([bam2bedfile], command)
	bedfile = bam2bedfile
	if exist_file_or_dir(trim_bedfile, 'Skip trim bed file.'):
		pass
	else:
		trim_bed(bam2bedfile, trim_bedfile, trim)
	bedfile = trim_bedfile
	
	## bed -> bdg
	bed_to_genomecov_bdg(bedfile, genomecov_bdgfile, genomeChrSizeFile)
	bdgfile = genomecov_bdgfile
	bwfile = output_prefix + '.trim_{}.genomecov.bw'.format(trim)
	if normalize_to > 0:
		if total_number <= 0:
			raise ZeroDivisionError('total_number should be greater than zero')
		norm_bdgfile = output_prefix + '.trim_{}.norm_div_{}_mult_{}.bdg'.format(trim, total_number, normalize_to)
		bwfile = output_prefix + '.trim_{}.norm_div_{}_mult_{}.bw'.format(trim, total_number, normalize_to)
		bdg_extract_multiple(genomecov_bdgfile, norm_bdgfile, normalize_to/total_number * (-1 if should_negative else 1) )
		bdgfile = norm_bdgfile
	
	## bdg -> bw
	command = 'bedGraphToBigWig ' + bdgfile + ' ' + genomeChrSizeFile + ' ' + bwfile
	check_file_then_exec_command([bwfile], command)

def get_bam_total(bamfile, output_prefix):
	flagstat_file = bamfile.replace(".bam", ".flagstat")
	if not os.path.exists(flagstat_file):
		flagstat_file = output_prefix + '.flagstat'
		command = 'samtools flagstat ' + bamfile + ' >' + flagstat_file
		check_file_then_exec_command([flagstat_file], command)
	with open(flagstat_file, 'r') as fin:
		for line in fin:
			F = line.split(' ')
			if F[3] == 'mapped':
				return int(F[0])

def split_bam_strand(bamfile, output_pos_bamfile, output_neg_bamfile):
	command = 'samtools view -b -f 0x10 ' + bamfile + ' >' + output_neg_bamfile
	check_file_then_exec_command([output_neg_bamfile], command)
	command = 'samtools view -b -F 0x10 ' + bamfile + ' >' + output_pos_bamfile
	check_file_then_exec_command([output_pos_bamfile], command)


def parse_args(should_show_content=True):
	parser = argparse.ArgumentParser(description='Convert BAM file to BigWig files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', metavar='genome.chromsize', required=True, help='genome chrom sizes file')
	parser.add_argument('-i', '--input', metavar='INPUT', required=True, help='input bam file')
	parser.add_argument('-o', '--output', metavar='OUTPUT', required=True, help='output prefix')
	parser.add_argument('-q', '--mapQ', metavar='INT', type=int, default=0, help='only keep reads with mapQ >= ?')
	parser.add_argument('-f', '--require-flag', metavar='INT', type=int, default=0, help='reads required flag ?')
	parser.add_argument('-F', '--filter-flag', metavar='INT', type=int, default=0, help='reads filtering flag ?')
	parser.add_argument('-s', '--strand', metavar='True|False|Both', default='Both', help='should divide strand plus and minus?')
	parser.add_argument('-n', '--normalize', metavar='INT', type=int, default=0, help='should scale to ? , 0 for NO')
	parser.add_argument('-T', '--total', metavar='INT', type=int, default=0, help='total read number for normalization , 0 for using mapped R1 in .flagstat')
	parser.add_argument('-t', '--trim', metavar='INT', type=int, default=0, help='trim read to remain 5\' ?-bp , 0 for NO')
	
	args = parser.parse_args()
	
	if should_show_content:
		for key in sorted(args.__dict__.keys()):
			print('{} = {}'.format(key, args.__dict__.get(key)), file=sys.stderr)
	
	return args
	
	
def main():
	### reference: Yyx_system_command_functions.20160607.pl, tlx2bed_v3.py
	start_time = time.time()
	_print('[PYTHON-START] ' + time.ctime(), file=sys.stderr)
	
	print('Now parse command-line arguments ...', file=sys.stderr)
	args = parse_args()
	
	bamfile = args.input
	output_prefix = args.output
	
	pos_bamfile = output_prefix + '.pos.bam'
	neg_bamfile = output_prefix + '.neg.bam'
	
	## deal with option mapQ, require_flag, filter_flag (bam read filtering)
	if args.mapQ + args.require_flag + args.filter_flag > 0:
		output_prefix = output_prefix + '.q{}_f{}_F{}'.format(args.mapQ, args.require_flag, args.filter_flag)
		filtered_bamfile = output_prefix + '.bam'
		command = 'samtools view -b -q {} -f {} -F {} {} >{}'.format(args.mapQ, args.require_flag, args.filter_flag, bamfile, filtered_bamfile)
		check_file_then_exec_command([filtered_bamfile], command)
		bamfile = filtered_bamfile
		pos_bamfile = output_prefix + '.pos.bam'
		neg_bamfile = output_prefix + '.neg.bam'
	
	## deal with option total == 0
	total_number = args.total
	if total_number <= 0:
		total_number = get_bam_total(bamfile, output_prefix)
		print('total is set to {}, according to flagstat'.format(total_number), file=sys.stderr)
	
	## deal with option strand
	if args.strand == 'False' or args.strand == 'Both':
		bam_to_bw(bamfile, output_prefix, args.genome, args.trim, args.normalize, total_number)
	if args.strand == 'True' or args.strand == 'Both':
		split_bam_strand(bamfile, pos_bamfile, neg_bamfile)
		bam_to_bw(pos_bamfile, output_prefix + '.pos', args.genome, args.trim, args.normalize, total_number)
		bam_to_bw(neg_bamfile, output_prefix + '.neg', args.genome, args.trim, args.normalize, total_number, True)

	_print('[PYTHON-END] ' + time.ctime(), file=sys.stderr)
	check_elapsed_time(start_time)


if __name__ == '__main__':
	main()
