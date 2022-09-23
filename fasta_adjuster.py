import os
import sys
data = []
def adjuster(fasta):
    #print(fasta)
    for x in fasta.strip().split('>'):
        content = x.split('\n')
        #Since after splitting by '>' first element will be blank so len(element) >1 is required
        if len(content) > 1:
            header = content[0]
            sequence = str(''.join(content[1:]))
            adjusted_fasta = ''.join(header + '\t'+ sequence)
            data.append(adjusted_fasta)
    return data

fasta = open(sys.argv[1])
fasta = fasta.read()
os.remove('adjusted_fasta.txt')
for i in adjuster(fasta):
    print(i)
    with open('adjusted_fasta.txt', 'a') as myfasta:
        myfasta.write(i+'\n')
        myfasta.close()