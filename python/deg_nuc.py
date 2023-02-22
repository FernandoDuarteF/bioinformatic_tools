primer1 = 'GCTANATCGATCGGACTACY'

deg = {'R': ['A', 'G'],
'Y': ['C', 'T'],
'S': ['G', 'C'],
'W': ['A', 'T'],
'K': ['G', 'T'],
'M': ['A', 'C'],
'B': ['C', 'G', 'T'],
'D': ['A', 'G', 'T'],
'H': ['A', 'C', 'T'],
'V': ['A', 'C', 'G'],
'N': ['A', 'C', 'G', 'T'],}


append1 = []
for a in range(len(primer1)):
    if primer1[a] not in deg:
        append1.append(primer1[a])
    else:
        append1.append(deg[primer1[a]])

for seq1 in seqs1:
    print(seq1)