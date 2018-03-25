from Bio import SeqIO
from Bio.Seq import Seq
for import_seq in SeqIO.parse("./sample_seq/IRGSP-1.0_genome.fasta", "fasta"):
	print(repr(import_seq))
#dt2 = data[76527-50:76527+51]