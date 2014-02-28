name1 = '/Users/haochuhuang/Documents/Jason Perera/3H9 BCR HTS Sequencing/FASTA/O1-SP.fasta'

file1 = open(name1, 'r')


x = file1.readlines()
Fasta_tup = []
for i in range(0, len(x)-1, 2):
    a = x[i]
    b = x[i+1]
    Fasta_tup.append((a,b))

Unique_list = []
Unique_tup = []
for item in Fasta_tup:
    if item[1] not in Unique_list:
        Unique_list.append(item[1])
        Unique_tup.append(item)
        
