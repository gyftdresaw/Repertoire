

#location = '/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/IMGT 10K/O2_SP/1_Summary_O2_SP_051113.txt'            #test_doc
location = '/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/compile IMGT Summary/'

names = ['O1-SP', 'O1-THY', 'Y1-SP', 'Y1-THY', 'O2-SP', 'O2-THY', 'Y2-SP', 'Y2-THY', 'O3-SP', 'O3-THY', 'Y3-SP', 'Y3-THY', 'O4-SP', 'O4-THY', 'Y4-SP', 'Y4-THY', 'O5-SP', 'O5-THY', 'Y5-SP', 'Y5-THY']


def V_genes_list(f):
    '''file --> list
takes a file open for reading, and returns a list of all the V genes that are used in that file.'''
    V_genes_list = []
    f.readline()
    for line in f:
        a = line.split('\t')
        if a[2] == 'productive':                      #if it is a productive rearrangement
            V_gene = a[3].split('*')[0].strip('Musmus ')      #gets the IgV region, before the allele nomenclature '*' and after the species nomenclature 'Musmus'
            if V_gene not in V_genes_list:
                V_genes_list.append(V_gene) 
    f.close()                                         #closes the file
    return V_genes_list






V_list = []
new_list = []

for item in names:
    file1 = open(location+item+'.txt', 'r')
    new_list = V_genes_list(file1)
    for gene in new_list:
        if gene not in V_list:
            V_list.append(gene)
    print (item +' complete')
    
    

