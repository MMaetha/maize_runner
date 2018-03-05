from collections import Counter;

f = open('TE_search07_chr1.csv','r')
data = f.read()
f.close()
num = Counter(data);
print(num);