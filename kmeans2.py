# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 22:54:21 2019

@author: siva9
"""

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
#for i in range(len(lists)):
#  for j in range(i+1,len(lists)):
#    t = parasail.nw(str(lists[i].seq), str(lists[j].seq), 10, 1, parasail.blosum62)
#    distances[i,j] = t.score
#    distances[j,i] = t.score
#    print(j)
#  print(i)
#    
#f = open("temp1.bin",'wb')
#np.save(f,distances)
#f.close()

#t = pw.align.globalxx(lists[1].seq, lists[127].seq)
#print(t[0][2])

f = open("temp1.bin", 'rb')
distances = np.load(f)


for i in range(len(lists)):
  mins = 9999999
  index = -1
  
  for j in range(len(a)):
    t = distances[i, a[j]]
    print(t)
    if t<mins:
      mins = t
      index = j
  print(index)
  print('\n')
  clusters[i] = index

def min_avg_distance_index(k):
    for i in range(len(lists)):
        min = 99999
        index = -1
        dist = 0
        count = 0
        if clusters[i]==k:
            for j in range(len(lists)):
                if clusters[j]==k and i!=j:
                    dist+=distances[i,j]
                    count = count + 1
            if count == 0       :
              avg_dist = 0
            else:
              avg_dist = dist/count
            if avg_dist<min:
                min = avg_dist
                index = i
    return index

new_clusters = np.zeros(shape=(k))
# Average distances in cluster q
for q in range(k):
    new_clusters[q] = min_avg_distance_index(q)