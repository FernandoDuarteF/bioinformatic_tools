file = open("multiqc_general_stats.txt", "r")
file.readline()
for line in file:
    col = float(line.split("\t")[2])
    if col > 59.0 or col < 57.0:
        print(line.strip())
