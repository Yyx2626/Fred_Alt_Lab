import sys
region_pos = {}
region_neg = {}
for line in open('%s.positive.cnt' % sys.argv[1]):
    l = line.strip().split()
    region_pos[l[0]] = float(l[-1])

for line in open('%s.negative.cnt' % sys.argv[1]):
    l = line.strip().split()
    region_neg[l[0]] = float(l[-1])

for region in region_pos:
    r = region.split('_')
    value = (region_pos[region] * region_neg[region]) ** 0.5
    print('%s\t%f\t0\t.\t%s\t%s' % ('\t'.join(r), value, r[1], r[2]))

