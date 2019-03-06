import numpy as np
from Bio import SeqIO
import random
from Levenshtein import distance

lists = []

for record in SeqIO.parse("DNASequences.fasta", "fasta"):
  lists.append(record)

# Number of CLusters
k = 4

# Randomly Initialize the center of each cluster
a = random.sample(range(0, len(lists)), k)

clusters = np.full(shape = (len(lists)), fill_value = -1, dtype = int)
for i in range(len(a)):
  clusters[a[i]] = i

distances = np.zeros(shape=(len(lists),len(lists)))
d2 = np.zeros(shape=(len(lists),len(lists)))

for i in range(len(lists)):
  for j in range(i+1,len(lists)):
    t = distance(str(lists[i].seq), str(lists[j].seq))
    distances[i][j] = t
    distances[j][i] = t
  print(i)
    
f = open("temp1.bin",'wb')
np.save(f,distances)
f.close()


f = open("temp1.bin", 'rb')
distances = np.load(f)
f.close()

f = open("prox_first.bin", 'rb')
d2 = np.load(f)
f.close()

for i in range(len(lists)):
    
  for j in range(len(a)):
    t = distances[i, a[j]]
    

#for i in range(len(lists)):
#  mins = 9999999
#  index = -1
#  
#  for j in range(len(a)):
#    t = distances[i, a[j]]
#    print(t)
#    if t<mins:
#      mins = t
#      index = j
#  print(index)
#  print('\n')
#  clusters[i] = index
#
#def min_avg_distance_index(k):
#    for i in range(len(lists)):
#        min = 99999
#        index = -1
#        dist = 0
#        count = 0
#        if clusters[i]==k:
#            for j in range(len(lists)):
#                if clusters[j]==k and i!=j:
#                    dist+=distances[i,j]
#                    count = count + 1
#            if count == 0       :
#              avg_dist = 0
#            else:
#              avg_dist = dist/count
#            if avg_dist<min:
#                min = avg_dist
#                index = i
#    return index
#
#new_clusters = np.zeros(shape=(k))
## Average distances in cluster q
#for q in range(k):
#    new_clusters[q] = min_avg_distance_index(q)