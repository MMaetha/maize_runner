from Bio import SeqIO
for seq in SeqIO.parse("./sample_seq/IRGSP-1.0_genome.fasta", "fasta"):
	print(seq.id)
	print(len(seq))
	print(repr(seq))
#dt2 = data[76527-50:76527+51]