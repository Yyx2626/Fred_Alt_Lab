
#### Usage: python3 this.py <VDJ.bed> <D.fa> <scripts_dir> <input.tlx> <output_prefix>


#### 2020-02-11, we think the default score threshold = 5 of D annotation in VDJ recombination (yyx_annotate_tlx_midD_LCS.20190116.py) may be too stringent.
####   We decide to change it to 3, and then implement another script to iteratively calculate D usage in VDJ


import sys, os, os.path, subprocess, time


### flushing print, reference: https://mail.python.org/pipermail/python-list/2015-November/698426.html
def _print(*args, **kwargs):
	file = kwargs.get('file', sys.stdout)
	print(*args, **kwargs)
	file.flush()


### reference: Yyx_system_command_functions.20160607.pl, tlx2bed_v3.py
start_time = time.time()
_print('[PYTHON-START] ' + time.ctime(), file=sys.stderr)


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
	if mode == 'all':
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
		
			
#import string
#def generateRandomString(length, charSet=string.ascii_uppercase+string.ascii_lowercase):
#	return ''.join(random.choice(charSet) for _ in range(length))


def list_assign(inputList, idx, value):
	while len(inputList) <= idx:
		inputList.append(None)
	
	inputList[idx] = value
	return inputList







bed_filename, D_ref_filename, scripts_root_dirname, input_tlx_filename, output_prefix = sys.argv[1:]

date_postfix = time.strftime('%Y%m%d')

### reference: work_flow.annotate_GSE82126_tlx.20190116.sh, work_flow.annotate_Alt225_tlx.20190116.sh

output_filenames = []
parts = ('junction', 'prey', 'bait')

is_tlx_gzipped = False
if input_tlx_filename.endswith('.gz'):
	is_tlx_gzipped = True

list_assign(output_filenames, 4, '{}.sequence_segmented.{}.tlx'.format(output_prefix, date_postfix))
list_assign(output_filenames, 6, '{}.intersectBed_annotated.{}.tlx'.format(output_prefix, date_postfix))
list_assign(output_filenames, 8, '{}.annotate_tlx_midD_LCS.{}.tlx'.format(output_prefix, date_postfix) )
list_assign(output_filenames, 10, '{}.HTGTS_annotate_merged.{}'.format(output_prefix, date_postfix) )
list_assign(output_filenames, 11, output_filenames[10] + '.tsv' )
list_assign(output_filenames, 12, output_filenames[10] + '.log' )
list_assign(output_filenames, 15, '{}.HTGTS_VDJ_annotated.{}.tsv'.format(output_prefix, date_postfix) )
#used_scripts_filenames = [
#		'yyx_tlx2bed.20181223.py', 'yyx_sequence_segment_tlx.20181221.py', 'yyx_annotate_tlx_with_intersectBedResults.20181223.py',
#		'yyx_annotate_tlx_midD_LCS.20190116.py', 'yyx_uniq_count_and_merge.20190111.py'
#	]
if not exist_file_or_dir(output_filenames[15], 'So skip all HTGTS annotation part...', mode='any'):
	stop_if_file_not_exist(input_tlx_filename)
	
	stop_if_file_not_exist('{}/Yyx_check_col_num.pl'.format(scripts_root_dirname))
	command = '{}/Yyx_check_col_num.pl '.format(scripts_root_dirname) + input_tlx_filename
	if is_tlx_gzipped:
		command = '{}/Yyx_check_col_num.pl <(zcat {} | head -n2)'.format(scripts_root_dirname, input_tlx_filename)
	_print('[SUBPROCESS-CHECK-OUTPUT] ' + command, file=sys.stderr)
	outBytes = subprocess.check_output(['/bin/bash', '-c', command])
	outStr = outBytes.decode('utf-8')
#	_print('[DEBUG] ' + outStr, file=sys.stderr)
	headline_colnum_str = outStr.split('\n')[1]
	_print('  Original column number = ' + headline_colnum_str, file=sys.stderr)
	original_colnum = int(headline_colnum_str)

	if not exist_file_or_dir(output_filenames[6], 'So skip intersectBed part...', mode='any'):
		stop_if_file_not_exist('{}/yyx_tlx2bed.20200126.py'.format(scripts_root_dirname))
		for i in range(len(parts)):
			part = parts[i]
			tmp_bed_filename = '{}.{}.{}.tmp.bed'.format(output_prefix, part, date_postfix)
			out_bed_filename = '{}.{}.{}.bed'.format(output_prefix, part, date_postfix)
			list_assign(output_filenames, i, out_bed_filename)
			command = ('z' if is_tlx_gzipped else '') + 'cat "' + input_tlx_filename + '" | python3 {}/yyx_tlx2bed.20181223.py '.format(scripts_root_dirname) + part + ' >' + tmp_bed_filename
			check_file_then_exec_command(tmp_bed_filename, command, should_time=True)
			command = 'bedtools intersect -a ' + tmp_bed_filename + ' -b ' + bed_filename + ' -wao >' + out_bed_filename
			check_file_then_exec_command(out_bed_filename, command, should_time=True)
			check_final_file_then_remove_intermediate_file(out_bed_filename, tmp_bed_filename)

		stop_if_file_not_exist('{}/yyx_sequence_segment_tlx.20181221.py'.format(scripts_root_dirname))
	#	list_assign(output_filenames, 4, '{}.sequence_segmented.{}.tlx'.format(output_prefix, date_postfix))
		command = ('z' if is_tlx_gzipped else '') + 'cat "' + input_tlx_filename + '" | python3 {}/yyx_sequence_segment_tlx.20181221.py >{}'.format(scripts_root_dirname, output_filenames[4])
		check_file_then_exec_command(output_filenames[4], command, should_time=True)

		stop_if_file_not_exist('{}/yyx_annotate_tlx_with_intersectBedResults.20181223.py'.format(scripts_root_dirname))
	#	list_assign(output_filenames, 6, '{}.intersectBed_annotated.{}.tlx'.format(output_prefix, date_postfix))
		command = 'python3 {}/yyx_annotate_tlx_with_intersectBedResults.20181223.py {} {} {} {}'.format(scripts_root_dirname, output_filenames[4], output_filenames[0], output_filenames[1], output_filenames[2]) + " | perl -pe 'BEGIN{$r=-1;} $r++; if($r==0){ s/anno1/junction/g; s/anno2/prey/g; s/anno3/bait/g; }' >" + output_filenames[6]
		check_file_then_exec_command(output_filenames[6], command, should_time=True)
		check_final_file_then_remove_intermediate_file(output_filenames[6], output_filenames[0:5])

	if not exist_file_or_dir(output_filenames[8], 'So skip miD_LCS part...', mode='any'):
		stop_if_file_not_exist('{}/yyx_annotate_tlx_midD_LCS.20190116.py'.format(scripts_root_dirname))
	#	list_assign(output_filenames, 8, '{}.annotate_tlx_midD_LCS.{}.tlx'.format(output_prefix, date_postfix) )
		command = ('z' if is_tlx_gzipped else '') + 'cat "' + input_tlx_filename + '" | python3 {}/yyx_annotate_tlx_midD_LCS.20190116.py {} {} >{}'.format(scripts_root_dirname, D_ref_filename, 3, output_filenames[8])
		check_file_then_exec_command(output_filenames[8], command, should_time=True)

	stop_if_file_not_exist('{}/yyx_uniq_count_and_merge.20190111.py'.format(scripts_root_dirname))
#	list_assign(output_filenames, 10, '{}.HTGTS_annotate_merged.{}'.format(output_prefix, date_postfix) )
#	list_assign(output_filenames, 11, output_filenames[10] + '.tsv' )
#	list_assign(output_filenames, 12, output_filenames[10] + '.log' )
	command = 'python3 {}/yyx_uniq_count_and_merge.20190111.py - 0 {}  {}:1-{} {}:1-{}'.format(scripts_root_dirname, output_filenames[10], output_filenames[6], original_colnum+5, output_filenames[8], original_colnum+5)
	check_file_then_exec_command(output_filenames[11], command, should_time=True)

	stop_if_file_not_exist('{}/yyx_show_or_skip_or_retrieve_columns.20190122.py'.format(scripts_root_dirname))
#	list_assign(output_filenames, 15, '{}.HTGTS_VDJ_annotate.{}.tsv'.format(output_prefix, date_postfix) )
	command = 'cat {} | python3 {}/yyx_show_or_skip_or_retrieve_columns.20190122.py skip '.format(output_filenames[11], scripts_root_dirname) + "'^count_in_file_.*'" + ' >' + output_filenames[15]
	check_file_then_exec_command(output_filenames[15], command, should_time=True)
	check_final_file_then_remove_intermediate_file(output_filenames[15], output_filenames[9:15])


_print('[PYTHON-END] ' + time.ctime(), file=sys.stderr)
check_elapsed_time(start_time)





