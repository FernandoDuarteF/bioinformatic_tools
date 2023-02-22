comp = {"A": "T", "C": "G", "G": "C", "T" :"A"}

seq = "ATTT"

rev = ""

for i in seq:
        rev += i.replace(i, comp[i])
        rev_com = rev[::-1]

print(rev_com)
