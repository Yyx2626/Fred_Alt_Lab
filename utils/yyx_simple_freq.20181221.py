#!/usr/bin/env python3
#### Usage: cat <input> | python3 this.py [should_skip_blank(default:True)]
#### Input: STDIN
#### Output: STDOUT   3 columns
####    level   count   proportion


import sys
import re

should_skip_blank = True
blank_pattern = re.compile('^\s*$')

def str2bool(x):
	try:
		x = int(x)
		return x > 0
	except ValueError:
		return x.lower() in ('true', 't', 'yes', 'y')

if len(sys.argv) > 1:
	should_skip_blank = str2bool(sys.argv[1])


hash = {}
total = 0
for line in sys.stdin:
	line = line.rstrip()
	if should_skip_blank:
		m = blank_pattern.search(line)
		if m:
			continue
	if line not in hash:
		hash[line] = 0
	hash[line] += 1
	total += 1

sorted_by_value = sorted(hash.items(), key=lambda kv: -kv[1])
for kv in sorted_by_value:
	print('{}\t{}\t{}'.format(kv[0], kv[1], kv[1]/total))
