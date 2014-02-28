#will generate a dictionary for each sample, where each V gene is a key, and values are count, then relative count of each J_gene


from bisect import bisect_left

def CDR3_usage(f):
    '''file --> list of 
takes a file open for reading, reads through each line of IMGT sequence, finds the V-gene used, and increments the value of that gene
in the V_genes_dict by one, then returns that dictionary'''
    V_genes_dict = {'count': 0, 'IGKV10-96': 0, 'IGKV10-95': 0, 'IGKV10-94': 0, 'IGHV1S32': 0, 'IGHV1S31': 0, 'IGHV1S37': 0, 'IGKV16-104': 0, 'IGKV1-35': 0, 'IGKV8-27': 0, 'IGKV3-10': 0, 'IGKV3-12': 0, 'IGHV3-3': 0, 'IGKV12-38': 0, 'IGKV19-93': 0, 'IGLV3': 0, 'IGLV2': 0, 'IGLV1': 0, 'IGKV8-28': 0, 'IGKV5-39': 0, 'IGKV8-21': 0, 'IGKV8-24': 0, 'IGKV8-26': 0, 'IGKV5-37': 0, 'IGKV13-84': 0, 'IGKV17-121': 0, 'IGKV17-127': 0, 'IGKV6-17': 0, 'IGKV2-a': 0, 'IGKV18-36': 0, 'IGKV6-14': 0, 'IGKV14-126': 0, 'IGKV1-117': 0, 'IGKV8-30': 0, 'IGKV8-34': 0, 'IGKV9-120': 0, 'IGKV9-123': 0, 'IGHV9-2': 0, 'IGHV9-1': 0, 'IGKV7-33': 0, 'IGKV9-129': 0, 'IGKV4-80': 0, 'IGKV4-81': 0, 'IGKV4-86': 0, 'IGKV1-110': 0, 'IGHV12-2-1': 0, 'IGHV1-9': 0, 'IGKV14-130': 0, 'IGHV1S65': 0, 'IGKV2-137': 0, 'IGHV3S1': 0, 'IGKV4-91': 0, 'IGKV4-90': 0, 'IGHV9-2-1': 0, 'IGKV12-46': 0, 'IGKV12-44': 0, 'IGKV12-41': 0, 'IGKV14-100': 0, 'IGHV1S96': 0, 'IGKV8-16': 0, 'IGKV8-19': 0, 'IGKV8-18': 0, 'IGKV1-99': 0, 'IGHV1S75': 0, 'IGKV4-92': 0, 'IGKV9-124': 0, 'pr IGKV10-96': 0, 'IGKV6-32': 0, 'IGHV15-2': 0, 'IGKV14-111': 0, 'IGKV4-69': 0, 'IGHV13-1': 0, 'IGKV4-61': 0, 'IGKV4-62': 0, 'IGKV4-63': 0, 'IGKV1-88': 0, 'IGKV20-101-2': 0, 'IGKV6-d': 0, 'IGKV6-c': 0, 'IGKV6-b': 0, 'IGKV5-48': 0, 'IGKV6-13': 0, 'IGKV3-7': 0, 'IGKV3-4': 0, 'IGKV3-5': 0, 'IGKV3-2': 0, 'IGKV3-3': 0, 'IGKV6-15': 0, 'IGKV3-1': 0, 'IGKV3-9': 0, 'IGKV4-79': 0, 'IGKV4-78': 0, 'IGKV4-73': 0, 'IGKV4-72': 0, 'IGKV4-71': 0, 'IGKV4-70': 0, 'IGKV4-74': 0, 'IGHV8S9': 0, 'pr IGLV2': 0, 'pr IGLV3': 0, 'IGKV12-98': 0, 'IGKV6-25': 0, 'IGKV6-23': 0, 'IGKV6-20': 0, 'IGKV6-29': 0, 'IGKV11-125': 0, 'IGKV2-112': 0, 'IGKV1-122': 0, 'IGHV1S44': 0, 'IGHV10-1': 0, 'IGKV1-132': 0, 'IGKV1-133': 0, 'IGKV12-89': 0, 'IGKV15-103': 0, 'IGKV5-43': 0, 'IGKV4-51': 0, 'IGKV4-50': 0, 'IGKV4-53': 0, 'IGKV4-52': 0, 'IGKV4-55': 0, 'IGKV2-109': 0, 'IGKV4-57': 0, 'IGKV4-59': 0, 'IGKV4-58': 0, 'IGKV1-131': 0, 'IGKV4-57-1': 0, 'IGKV1-135': 0}
    f.readline()
    List = []
    # count = 0
    for line in f:
        a = line.split('\t')
        if a[2] == 'productive':                              #if it is a productive rearrangement
            reads = int(a[1].split('-')[1])
            vs = a[3].split('*')[0]
            V_gene = vs[vs.find('IG'):] # HGEDIT explicitly find IG to start V gene ID
            J_gene = a[9].split()[1].split('*')[0]       #gets J_gene segment minus additional nomenclature
            # HGEDIT spliting CDR3 to get rid of parentheticals sometime
            CDR3 = a[20].split()[0]   #gets the IgV region, before the allele nomenclature '*' and after the species nomenclature 'Musmus'
            # HGEDIT getting rid of counts in csv
            List.append((reads, V_gene, J_gene, CDR3))
            # count = count + 1
    f.close()     #closes the file 
    return List


def compare_overlap(a, b):
    end = len(b)
    count = 0
    for item in a:
        pos = bisect_left(b, item)
        if pos < end:
            if b[pos] == item:
                count = count + 1
    return count            

import csv

# HGEDIT condense each of these data files to just
# count, V, J, and CDR3
data_names = ['O1-SP', 'O1-THY', 'Y1-SP', 'Y1-THY', 'O2-SP', 'O2-THY', 'Y2-SP', 'Y2-THY', 'O3-SP', 'O3-THY', 'Y3-SP', 'Y3-THY', 'O4-SP', 'O4-THY', 'Y4-SP', 'Y4-THY', 'O5-SP', 'O5-THY', 'Y5-SP', 'Y5-THY']

from_dir = 'compile IMGT Summary/'
to_dir = 'IMGT Condensed/'

for data_file in data_names:
    from_file = open(from_dir + data_file + '.txt', 'r')
    CSVList = CDR3_usage(from_file)

    to_file = open(to_dir + data_file + '_condensed.txt', 'w')
    writer = csv.writer(to_file)

    for item in CSVList:
	writer.writerow(item)

    from_file.close()
    to_file.close()
