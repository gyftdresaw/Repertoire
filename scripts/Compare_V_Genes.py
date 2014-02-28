#will generate a dictionary for each sample, where each V gene is a key, and values are count, then relative count of each J_gene


#location = '/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/IMGT 10K/O2_SP/1_Summary_O2_SP_051113.txt'            #test_doc
location = '/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/compile IMGT Summary/'                 #compressed sequence files

#a list of all Vgenes found in our sequencing data
V_genes_list = ['IGHV1-9', 'IGHV10-1', 'IGHV12-2-1', 'IGHV13-1', 'IGHV15-2', 'IGHV1S31', 'IGHV1S32', 'IGHV1S37', 'IGHV1S44', 'IGHV1S65', 'IGHV1S75', 'IGHV1S96', 'IGHV3-3', 'IGHV3S1', 'IGHV8S9', 'IGHV9-1', 'IGHV9-2', 'IGHV9-2-1', 'IGKV1-110', 'IGKV1-117', 'IGKV1-122', 'IGKV1-131', 'IGKV1-132', 'IGKV1-133', 'IGKV1-135', 'IGKV1-35', 'IGKV1-88', 'IGKV1-99', 'IGKV10-94', 'IGKV10-95', 'IGKV10-96', 'IGKV11-125', 'IGKV12-38', 'IGKV12-41', 'IGKV12-44', 'IGKV12-46', 'IGKV12-89', 'IGKV12-98', 'IGKV13-84', 'IGKV14-100', 'IGKV14-111', 'IGKV14-126', 'IGKV14-130', 'IGKV15-103', 'IGKV16-104', 'IGKV17-121', 'IGKV17-127', 'IGKV18-36', 'IGKV19-93', 'IGKV2-109', 'IGKV2-112', 'IGKV2-137', 'IGKV2-a', 'IGKV20-101-2', 'IGKV3-1', 'IGKV3-10', 'IGKV3-12', 'IGKV3-2', 'IGKV3-3', 'IGKV3-4', 'IGKV3-5', 'IGKV3-7', 'IGKV3-9', 'IGKV4-50', 'IGKV4-51', 'IGKV4-52', 'IGKV4-53', 'IGKV4-55', 'IGKV4-57', 'IGKV4-57-1', 'IGKV4-58', 'IGKV4-59', 'IGKV4-61', 'IGKV4-62', 'IGKV4-63', 'IGKV4-69', 'IGKV4-70', 'IGKV4-71', 'IGKV4-72', 'IGKV4-73', 'IGKV4-74', 'IGKV4-78', 'IGKV4-79', 'IGKV4-80', 'IGKV4-81', 'IGKV4-86', 'IGKV4-90', 'IGKV4-91', 'IGKV4-92', 'IGKV5-37', 'IGKV5-39', 'IGKV5-43', 'IGKV5-48', 'IGKV6-13', 'IGKV6-14', 'IGKV6-15', 'IGKV6-17', 'IGKV6-20', 'IGKV6-23', 'IGKV6-25', 'IGKV6-29', 'IGKV6-32', 'IGKV6-b', 'IGKV6-c', 'IGKV6-d', 'IGKV7-33', 'IGKV8-16', 'IGKV8-18', 'IGKV8-19', 'IGKV8-21', 'IGKV8-24', 'IGKV8-26', 'IGKV8-27', 'IGKV8-28', 'IGKV8-30', 'IGKV8-34', 'IGKV9-120', 'IGKV9-123', 'IGKV9-124', 'IGKV9-129', 'IGLV1', 'IGLV2', 'IGLV3', 'pr IGKV10-96', 'pr IGLV2', 'pr IGLV3']



def V_gene_usage(f):
    '''file --> list
takes a file open for reading, reads through each line of IMGT sequence, finds the V-gene used, and increments the value of that gene
in the V_genes_dict by one, then returns that dictionary'''
    V_genes_dict = {'count': 0, 'IGKV10-96': 0, 'IGKV10-95': 0, 'IGKV10-94': 0, 'IGHV1S32': 0, 'IGHV1S31': 0, 'IGHV1S37': 0, 'IGKV16-104': 0, 'IGKV1-35': 0, 'IGKV8-27': 0, 'IGKV3-10': 0, 'IGKV3-12': 0, 'IGHV3-3': 0, 'IGKV12-38': 0, 'IGKV19-93': 0, 'IGLV3': 0, 'IGLV2': 0, 'IGLV1': 0, 'IGKV8-28': 0, 'IGKV5-39': 0, 'IGKV8-21': 0, 'IGKV8-24': 0, 'IGKV8-26': 0, 'IGKV5-37': 0, 'IGKV13-84': 0, 'IGKV17-121': 0, 'IGKV17-127': 0, 'IGKV6-17': 0, 'IGKV2-a': 0, 'IGKV18-36': 0, 'IGKV6-14': 0, 'IGKV14-126': 0, 'IGKV1-117': 0, 'IGKV8-30': 0, 'IGKV8-34': 0, 'IGKV9-120': 0, 'IGKV9-123': 0, 'IGHV9-2': 0, 'IGHV9-1': 0, 'IGKV7-33': 0, 'IGKV9-129': 0, 'IGKV4-80': 0, 'IGKV4-81': 0, 'IGKV4-86': 0, 'IGKV1-110': 0, 'IGHV12-2-1': 0, 'IGHV1-9': 0, 'IGKV14-130': 0, 'IGHV1S65': 0, 'IGKV2-137': 0, 'IGHV3S1': 0, 'IGKV4-91': 0, 'IGKV4-90': 0, 'IGHV9-2-1': 0, 'IGKV12-46': 0, 'IGKV12-44': 0, 'IGKV12-41': 0, 'IGKV14-100': 0, 'IGHV1S96': 0, 'IGKV8-16': 0, 'IGKV8-19': 0, 'IGKV8-18': 0, 'IGKV1-99': 0, 'IGHV1S75': 0, 'IGKV4-92': 0, 'IGKV9-124': 0, 'pr IGKV10-96': 0, 'IGKV6-32': 0, 'IGHV15-2': 0, 'IGKV14-111': 0, 'IGKV4-69': 0, 'IGHV13-1': 0, 'IGKV4-61': 0, 'IGKV4-62': 0, 'IGKV4-63': 0, 'IGKV1-88': 0, 'IGKV20-101-2': 0, 'IGKV6-d': 0, 'IGKV6-c': 0, 'IGKV6-b': 0, 'IGKV5-48': 0, 'IGKV6-13': 0, 'IGKV3-7': 0, 'IGKV3-4': 0, 'IGKV3-5': 0, 'IGKV3-2': 0, 'IGKV3-3': 0, 'IGKV6-15': 0, 'IGKV3-1': 0, 'IGKV3-9': 0, 'IGKV4-79': 0, 'IGKV4-78': 0, 'IGKV4-73': 0, 'IGKV4-72': 0, 'IGKV4-71': 0, 'IGKV4-70': 0, 'IGKV4-74': 0, 'IGHV8S9': 0, 'pr IGLV2': 0, 'pr IGLV3': 0, 'IGKV12-98': 0, 'IGKV6-25': 0, 'IGKV6-23': 0, 'IGKV6-20': 0, 'IGKV6-29': 0, 'IGKV11-125': 0, 'IGKV2-112': 0, 'IGKV1-122': 0, 'IGHV1S44': 0, 'IGHV10-1': 0, 'IGKV1-132': 0, 'IGKV1-133': 0, 'IGKV12-89': 0, 'IGKV15-103': 0, 'IGKV5-43': 0, 'IGKV4-51': 0, 'IGKV4-50': 0, 'IGKV4-53': 0, 'IGKV4-52': 0, 'IGKV4-55': 0, 'IGKV2-109': 0, 'IGKV4-57': 0, 'IGKV4-59': 0, 'IGKV4-58': 0, 'IGKV1-131': 0, 'IGKV4-57-1': 0, 'IGKV1-135': 0}
    f.readline()
    total = 0
    for line in f:
        a = line.split('\t')
        if a[2] == 'productive':                              #if it is a productive rearrangement
   #        reads = int(a[1].split('-')[1])                   #gets the number of reads from the sequence ID header 
            V_gene = a[3].split('*')[0].strip('Musmus ')      #gets the IgV region, before the allele nomenclature '*' and after the species nomenclature 'Musmus'
            V_genes_dict[V_gene] += 1                         ###if you want to disregard identical sequences, use this line 
   #        V_genes_dict[V_gene] += reads                     ###if you want to account for identical sequences, use this line
            total += 1                                        #keeps a running count of how many V genes were read (Unique)
   #        total += reads                                    #keeps a running tally of how many V genes were read (including identical)
    V_genes_dict['count'] = total
    f.close()                                         #closes the file
    return V_genes_dict



names = ['O1-SP', 'O1-THY', 'Y1-SP', 'Y1-THY', 'O2-SP', 'O2-THY', 'Y2-SP', 'Y2-THY', 'O3-SP', 'O3-THY', 'Y3-SP', 'Y3-THY', 'O4-SP', 'O4-THY', 'Y4-SP', 'Y4-THY', 'O5-SP', 'O5-THY', 'Y5-SP', 'Y5-THY']
names2 = ['SP-O1', 'SP-O2', 'SP-O3', 'SP-O4', 'SP-O5', 'SP-Y1', 'SP-Y2', 'SP-Y3', 'SP-Y4', 'SP-Y5', 'THY-O1', 'THY-O2', 'THY-O3', 'THY-O4', 'THY-O5', 'THY-Y1', 'THY-Y2', 'THY-Y3', 'THY-Y4', 'THY-Y5']

V_usage = {}
for item in names:
    file1 = open(location+item+'.txt', 'r')
    new_dict = V_gene_usage(file1)
    V_usage[str(item[3:] + '-' + item[:2])] = new_dict       #rearranges file names to a more easily sortable name
    print (item + ' is complete')


def write_txt():
    location = ('/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/Prism and Excel/Compare_V_Genes.txt')
    file1 = open(location, 'w')
    file1.



##file1.write('Sample Name' + '\t')                   Strategy for getting files into excel format
##for item in V_genes_list:
##	file1.write(item + '\t')
##
##
##
##for organ in names2:
##	file1.write(organ + '\t')
##	for family in V_genes_list:
##		file1.write(str(V_usage[organ][family]) + '\t')
##	file1.write('\n')
