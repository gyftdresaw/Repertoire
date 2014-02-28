
# general poking around
import numpy as np
import csv

data_names = ['O1-SP', 'O1-THY', 'Y1-SP', 'Y1-THY', 'O2-SP', 'O2-THY', 'Y2-SP', 'Y2-THY', 'O3-SP', 'O3-THY', 'Y3-SP', 'Y3-THY', 'O4-SP', 'O4-THY', 'Y4-SP', 'Y4-THY', 'O5-SP', 'O5-THY', 'Y5-SP', 'Y5-THY']

data_dir = 'IMGT Condensed/'

# simply complete list of V, J, and CDR3 categories
Vset = set()
Jset = set()
Cset = set()
combset = set()

for data_file in data_names:
    from_file = open(data_dir + data_file + '_condensed.txt','r')
    for line in from_file:
        d = line.strip().split(',')
        Vset.add(d[1])
        Jset.add(d[2])
        Cset.add(d[3])
        combset.add((d[1],d[2],d[3]))
    from_file.close()

# there are only ~530,000 unique V,J and CDR3 combinations
# let's save the complete data in a matrix format based on these combinations
to_dir = 'IMGT Combined/'
combos = sorted(combset)

# make dictionary to hold combination index
combo_index = {}
for i in xrange(len(combos)):
    combo_index[combos[i]] = i

ages = ['Y','O']
tissue = ['SP','THY']
samples = 5 # 5 samples in each of the 4 conditions
total_samples = 20

# put complete data in monster array
# generate ordered list of data file names
data_ordered = []
for i in xrange(len(ages)):
    for j in xrange(len(tissue)):
        for k in xrange(1,samples+1):
            data_ordered.append(ages[i]+str(k)+'-'+tissue[j])

# preallocate standard array to hold count data
counts = [[0 for i in xrange(total_samples)] for j in xrange(len(combos))]

# go through each data file and update counts
for i in xrange(len(data_ordered)):
    from_file = open(data_dir + data_ordered[i] + '_condensed.txt','r')
    for line in from_file:
        d = line.strip().split(',')
        current_combo = (d[1],d[2],d[3])
        current_count = int(d[0])
        
        counts[combo_index[current_combo]][i] += current_count

    from_file.close()

# write count data to file
to_file = open(to_dir + 'IMGT_combined.txt', 'w')
writer = csv.writer(to_file,delimiter='\t')

# write header
writer.writerow(['V_GENE','J_GENE','CDR3'] + data_ordered)
for i in xrange(len(combos)):
    writer.writerow(list(combos[i]) + counts[i])

to_file.close()

