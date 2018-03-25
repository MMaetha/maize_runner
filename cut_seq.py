from Bio import SeqIO

count=0
f = open('./sample_seq/IRGSP-1.0_genome.fasta','r')
data = f.readlines()
if(data=='>'):
    count = count +1
while(count<2):
    break;
f.close()
dt2 = data[76527-50:76527+51]

for i in range(len(data)):
    if(data[i]=='a'):
        data[i]=='A'
    if(data[i]=='t'):
        data[i]=='T'
    if(data[i]=='c'):
        data[i]=='C'
    if(data[i]=='g'):
        data[i]=='G'
        
print(dt2)