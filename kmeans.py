import numpy as np
from Bio import SeqIO
import numpy as np
from collections import defaultdict

lists = defaultdict(str)

for record in SeqIO.parse("DNASequences.fasta", "fasta"):
  lists[record.id] += record

X = "GTCCATACA"
Y = "TCATATCAG"

# Parameters
penalty = -1
reward = 1

# Match or not
def score(n1,n2):
    if(n1==n2):
        return reward
    else:
        return penalty

# Function to find the similarity
def similarity_score(X,Y):
    # Making the matrix
    score_matrix = np.zeros(shape=(len(X)+1,len(Y)+1))
    for i in range(len(X)+1):
        score_matrix[i,0] = i*penalty
    for j in range(len(Y)+1):
        score_matrix[0,j] = j*penalty

    # Filling up the rest of the values
    for i in range(1,len(X)+1):
        for j in range(1,len(Y)+1):
            score_matrix[i,j] = max(score_matrix[i-1,j-1]+score(X[i-1],Y[j-1]),score_matrix[i-1,j]+penalty,score_matrix[i,j-1]+penalty)

    # Getting the similarity value
    i = len(X)
    j = len(Y)
    value = 0
    while i>0 or j>0:
        if i>0 and j>0 and X[i-1]==Y[j-1]:
            value = value + 1
            i = i-1
            j = j-1
        elif i>0 and score_matrix[i,j]==score_matrix[i-1,j]+penalty:
            i = i-1
        else:
            j = j-1

    return value

print(similarity_score(X,Y))
