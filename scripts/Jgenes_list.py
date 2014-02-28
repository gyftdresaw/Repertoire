

#test_doc = '/Users/jasonperera/Documents/3H9 Project/IMGT10K Summary/1_Summary_O1_SP.txt'
location = '/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/DC IMGT Summary/'

names = ['O1-SP', 'O1-THY', 'Y1-SP', 'Y1-THY', 'O2-SP', 'O2-THY', 'Y2-SP', 'Y2-THY', 'O3-SP', 'O3-THY', 'Y3-SP', 'Y3-THY', 'O4-SP', 'O4-THY', 'Y4-SP', 'Y4-THY', 'O5-SP', 'O5-THY', 'Y5-SP', 'Y5-THY']


def J_genes(f):
    J = {}
    f.readline()
    count = 0
    for line in f:
        a = line.split('\t')
        if a[2] == 'productive':                      #if it is a productive rearrangement
            name = a[9].split()                       #gets the J gene data from xls
            name1 = name[0] + ' ' + name[1]           #scrapes just the inital J gene call, removing any additional guesses
            if name1 not in J:                        #if that J gene not in dict, then add with count = 1
                J[name1] = 1
            if name1 in J:                            #if that J gene is already there, increment by 1
                J[name1] += 1                           
            count += 1                                #keeps an overall count, in case you want to calculate percentage at end
    f.close()                                         #closes the file
    return_list = []
    for key in J:
        return_list.append((key, J[key]))           #returns a list of tup so that it can be more easily sorted
    return return_list

J = []

    
for item in names:
    file1 = open(location+item+'.txt', 'r')
    J.append((item, J_genes(file1)))
    print (item +' complete')
    


