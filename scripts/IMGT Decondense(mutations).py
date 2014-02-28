''' this program will take a IMGT file, and will uncompress the FASTA compression introduced earlier'''


list1 = ['O1-SP', 'O1-THY', 'Y1-SP', 'Y1-THY', 'O2-SP', 'O2-THY', 'Y2-SP', 'Y2-THY', 'O3-SP', 'O3-THY', 'Y3-SP', 'Y3-THY', 'O4-SP', 'O4-THY', 'Y4-SP', 'Y4-THY', 'O5-SP', 'O5-THY', 'Y5-SP', 'Y5-THY']
#list1 = ['O1-SP']
for item in list1:
    file1 = open('/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/DC IMGT Mutation/' + item +'.txt', 'w')
    file2 = open('/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/compile IMGT Mutation/' + item + '.txt', 'r')
    file2.readline()
    for line in file2:
        if line != '':
            a = (line.split()[1]).split('-')[1]
            for i in range(int(a)):
                file1.write(line)
    file2.close()    
    file1.close()
    print(item + ' is finished')
