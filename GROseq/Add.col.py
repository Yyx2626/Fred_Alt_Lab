import sys
if len(sys.argv) != 2:
    sys.exit('sorted gmean file')

corlist = {0: '0,51,17', 1: '0,77,26', 2: '0,102,34', 3: '0,128,43', 4: '0,153,51', 
           5: '0,179,60', 6: '0,204,68', 7: '0,230,77', 8: '0,255,85', 9: '26,255,102'}

total_line = 0
for line in open(sys.argv[1]):
    total_line += 1
total_line += 10

num_col = {}
for c in range(0, 10):
    color = corlist[c]
    for g in range(int(total_line/10*c), int(total_line/10*(c+1))):
        num_col[g] = color

n = 0
for line in open(sys.argv[1]):
    print('%s\t%s' % (line.strip(), num_col[n]))
    n += 1

