from Bio import SeqIO
from Bio.Seq import Seq
count=0
gap=100
place=76527
for import_seq in SeqIO.parse("./sample_seq/IRGSP-1.0_genome.fasta", "fasta"):
	count = count+1
	if(count==1):
		print(import_seq.seq[place-gap:place+gap+1])
#dt2 = data[76527-50:76527+51]