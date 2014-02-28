'''compressed FASTA files had to be split because IMGT only takes 150,000 sequences at a time
this program will rejoin all of the files that come back from a given sample using a
user generated list of sample names, and number of files they were split into'''



list1 = [('O1-SP', 8), ('O1-THY', 8), ('Y1-SP', 7), ('Y1-THY', 8), ('O2-SP', 8), ('O2-THY', 3), ('Y2-SP', 8), ('Y2-THY', 7), ('O3-SP', 5), ('O3-THY', 5), ('Y3-SP', 7), ('Y3-THY', 4), ('O4-SP', 12), ('O4-THY', 7), ('Y4-SP', 8), ('Y4-THY', 7), ('O5-SP', 8), ('O5-THY', 7), ('Y5-SP', 9), ('Y5-THY', 7)]
#list1 = [('O1-SP' , 8)]
for item in list1:
    file1 = open('/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/compile IMGT Summary/' + item[0] +'.txt', 'w')
    file2 = open('/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/ALL IMGT/' + item[0] + '_C-' + str(1) + '/1_Summary_'+ item[0] + '_C-' + str(1) + '_091113.txt', 'r')
    for line in file2:
        file1.write(line)
    file2.close()    
    for i in range(2, item[1]+1):
        file2 = open('/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/ALL IMGT/' + item[0] + '_C-' + str(i) + '/1_Summary_'+ item[0] + '_C-' + str(i) + '_091113.txt', 'r')
        file2.readline()                #removes header from write loop
        for line in file2:
            file1.write(line)
        file2.close()
    file1.close()        
    print(item[0] + ' is finished')
