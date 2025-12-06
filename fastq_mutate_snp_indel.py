import os, sys

# python3 fasta_mutator.py myfasta.fasta snps.bed indels.bed output.fasta

def read_input(myfasta):
    with open(myfasta, "r") as myinput:
        myinput = myinput.readlines()
        sequences = {}
        chr_id = None
        chr_seq = []
        for line in myinput:
            line = line.strip()
            if line.startswith('>'):
                if chr_id != None:
                    sequences[chr_id] = "".join(chr_seq)
                chr_id = line.split(' ')[0]
                chr_seq = []
            else:
                chr_seq.append(line)

        if chr_id != None:
            sequences[chr_id] = "".join(chr_seq)

    return sequences


def reads_snp_file(snpfile):
    snps_dict = {}
    with open(snpfile) as snps:
        for line in snps:
            line = line.strip().split('\t')
            key = line[0] + "|" + line[1]
            snps_dict[key] = [line[2], line[3]]   # REF, ALT
    return snps_dict


def read_indel_file(indel_file):
    """
    chr   pos   INS/DEL   size_or_seq
    """
    indels = []
    with open(indel_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrid  = parts[0]
            pos    = int(parts[1])
            itype  = parts[2].upper()

            if itype == "DEL":
                size = int(parts[3])
                indels.append((chrid, pos, "DEL", size))
            elif itype == "INS":
                seq = parts[3]
                indels.append((chrid, pos, "INS", seq))
    return indels


def fasta_mutator(sequences, snps_dict, indels):
    new_chr = {}
    total_snps = 0
    total_indels = 0

    # sort indels descending so coordinate shifts don't break positions
    indels_sorted = sorted(indels, key=lambda x: x[1], reverse=True)

    for ids in sequences:
        chrm = ids.replace('>', '')
        seq = bytearray(sequences[ids], 'utf-8')

        # --------------------
        # Apply SNPs
        # --------------------
        for snp in snps_dict:
            chrp, pos = snp.split('|')
            pos = int(pos)
            ref, alt = snps_dict[snp]

            if chrp == chrm:
                if chr(seq[pos-1]) == ref:
                    seq[pos-1] = ord(alt)
                    total_snps += 1

        # --------------------
        # Apply indels
        # --------------------
        for (chrid, pos, itype, val) in indels_sorted:
            if chrid != chrm:
                continue

            if itype == "DEL":
                size = val
                # Python slice deletion
                del seq[pos-1 : pos-1+size]
                total_indels += 1

            elif itype == "INS":
                insseq = val.encode('utf-8')
                seq[pos-1:pos-1] = insseq
                total_indels += 1

        new_chr[chrm] = seq.decode('utf-8')

    print("=> SNPs applied:   ", total_snps)
    print("=> Indels applied: ", total_indels)
    return new_chr


def get_output(new_chr, myoutput):
    with open(myoutput, 'w') as output:
        for chrid in new_chr:
            output.write('>' + chrid + '\n')
            output.write(new_chr[chrid] + '\n')


if __name__ == "__main__":
    print("Reading FASTA ...")
    fasta_file = read_input(sys.argv[1])

    print("Reading SNP file ...")
    snp_dict = reads_snp_file(sys.argv[2])

    print("Reading INDEL file ...")
    indel_list = read_indel_file(sys.argv[3])

    print("Applying mutations ...")
    new_fasta = fasta_mutator(fasta_file, snp_dict, indel_list)

    print("Writing output FASTA...")
    get_output(new_fasta, sys.argv[4])

    print("Done.")
