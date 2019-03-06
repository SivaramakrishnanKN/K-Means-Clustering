import numpy as np
from Bio import SeqIO
from Bio import pairwise2 as pw
from Bio.pairwise2 import format_alignment
import numpy as np
import random
from collections import defaultdict
import edlib
import parasail

lists = []

for record in SeqIO.parse("DNASequences.fasta", "fasta"):
  lists.append(record)

X = "xy"
Y = "abdf"

# Parameters
penalty = -1
reward = 1

# Match or not
def score(n1,n2):
    if(n1==n2):
        return reward
    else:
        return penalty

k = 4

a = random.sample(range(0, len(lists)), k)

clusters = np.full(shape = (len(lists)), fill_value = -1, dtype = int)
for i in range(len(a)):
  clusters[a[i]] = i


t = pw.align.globalxx(X, Y)
print(t[0][2])


distances = np.zeros(shape=(len(lists),len(lists)))
for i in range(len(lists)):
  for j in range(i+1,len(lists)):
    t = pw.align.globalxx(lists[i].seq, lists[j].seq)
    #t = parasail.nw(str(lists[i].seq), str(lists[j].seq), 10, 1, parasail.blosum62)
    distances[i,j] = t[0][2]
    distances[j,i] = t[0][2]
    print(j)
  print(i)
    
f = open("temp.bin",'wb')
np.save(f,distances)
f.close()
#
#t = pw.align.globalxx(lists[1].seq, lists[127].seq)
#print(t[0][2])

#for i in range(len(lists)):
#  mins = 9999999
 # index = -1
  
  #for j in range(len(a)):
   # t = pw.align.globalxx(lists[i].seq, lists[a[j]].seq)
    #print(t[0][2])
    #if t[0][2]<mins:
     # mins = t[0][2]
      #index = j
  #print(index)
  #print('\n')
  #clusters[i] = index

