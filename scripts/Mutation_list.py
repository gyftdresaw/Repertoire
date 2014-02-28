

#test_doc = '/Users/jasonperera/Documents/3H9 Project/IMGT10K Summary/1_Summary_O1_SP.txt'
location = '/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/compile IMGT Mutation/'

names = ['O1-SP', 'O1-THY', 'Y1-SP', 'Y1-THY', 'O2-SP', 'O2-THY', 'Y2-SP', 'Y2-THY', 'O3-SP', 'O3-THY', 'Y3-SP', 'Y3-THY', 'O4-SP', 'O4-THY', 'Y4-SP', 'Y4-THY', 'O5-SP', 'O5-THY', 'Y5-SP', 'Y5-THY']


def Mutations(f):
    J = {}
    f.readline()
    count = 0
    total = 0
    nonsilent = 0
    silent = 0
    for line in f:
        a = line.split('\t')
        if a[2] == 'productive':                      #if it is a productive rearrangement
            count = count + 1
            total = total + int(a[7].split()[0])      #takes the # of mutations minus the weird parenthetical term that I can't figure out (N-addition?)
            nonsilent = nonsilent + int(a[8].split()[0])
            silent = silent + int(a[9].split()[0])
    f.close()                                      #closes the file
    return_list = (count, total, nonsilent, silent)
    return return_list

J = []

    
for item in names:
    file1 = open(location+item+'.txt', 'r')
    J.append((item[3:] + '-' +item[:2], Mutations(file1)))
    print (item +' complete')
    


