#from Bio import SeqIO

#for seq_record in SeqIO.parse("parse_fasta.fasta", "fasta"):
    #print(seq_record.id)
    #print(seq_record.seq)


num_lines = sum(1 for line in open("rosalind_gc(1).fasta", "r"))
fasta_file=open("rosalind_gc(1).fasta", "r")
seq_n, count = "", 0
for s in fasta_file:
    s = s.strip()
    count += 1
    if s[0]==">":
        if seq_n != "":
            print(seq_n, ((seq.count("G") + seq.count("C"))/float(len(seq)))*100, seq, sep="\n")
        seq_n, seq = s, ""
    else:
        seq += s
        if count == num_lines:
            print(seq_n, ((seq.count("G") + seq.count("C"))/float(len(seq)))*100, seq, sep="\n")

#parse_fasta("rosalind_gc.fasta")

