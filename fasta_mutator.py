import os,sys
# python3 fasta_mutator.py myfasta.fasta snp.bed check_output.fasta
def read_input(myfasta):
    with open(myfasta, "r") as myinput:
        myinput = myinput.readlines()
        sequences = {}
        chr_id = None
        chr_seq = []
        for line in myinput:
            line = line.strip()
            # get the header here
            if line.startswith('>'):
                # print(line)
                if chr_id != None:
                    sequences[chr_id] = "".join(chr_seq)
                chr_id = line.split(' ')[0]
                #print(chr_id)
                chr_seq = []
            # get the sequence here
            else:
                chr_seq.append(line)
                # print(chr_seq)

        '''
        Since all > containing lines finished the last chr element could not be inserted into sequences. 
        However it got stored in chr_id and chr_seq
        Therefore we have to manually add the last chromosome and the respective sequence
        '''
        if chr_id != None:
            sequences[chr_id] = "".join(chr_seq)

    return(sequences)

def reads_snp_file(snpfile):
    snps_dict = {}
    with open(snpfile) as snps:
        snps = snps.readlines()
        for line in snps:
            line = line.strip().split('\t')
            # print(line[0], line[1])
            sites = str(''.join(line[0] + '|' + line[1]))
            snps_dict[sites] = [line[2], line[3]]
    return(snps_dict)

def fasta_mutator(sequences, snps_dict):
    new_chr = {}
    count = 0
    for ids in sequences:
        chrm = ids.replace('>', '')
        '''
        Let's convert string to bytearray. Byte array are easily mutable and useful very long sequences
        '''
        seqs = bytearray(sequences[ids], 'utf-8')
        for snps in snps_dict:
            chrp = snps.split('|')[0]
            pos = int(snps.split('|')[1])
            ref = snps_dict[snps][0]
            alt = snps_dict[snps][1]
            if chrp == chrm:
                if chr(seqs[pos-1]) == ref:
                    seqs[pos-1] = ord(alt)
                    count += 1
                    new_seqs = seqs.decode('utf-8')
                    new_chr[chrp] = new_seqs
                    # print(chrp,pos,ref,alt, chrm, pos, chr(seqs[pos-1]), new_seqs)
    print('=>  total sites changed: ', count)
    return(new_chr)

def get_output(new_chr, myoutput):
    with open(myoutput, 'w') as output:
        for nlines in new_chr:
            # print(nlines)
            output.write('>' + nlines + '\n' + new_chr[nlines] +  '\n')
        output.close()

print('reading fasta ...')
fasta_file  = read_input(sys.argv[1])
print('reading snp file ...')
snp_file  = reads_snp_file(sys.argv[2])
print('creating desired single nt mutations')
new_fasta_file = fasta_mutator(fasta_file, snp_file)
print('fetching mutated file')
get_output(new_fasta_file, sys.argv[3])
