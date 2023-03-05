from itertools import product
import os, sys
import argparse

try:
	import regex
except ImportError:
	print("regex module required ($ pip install regex)")

#import regex


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-s","--sequence", help="sequence in fasta format")
  parser.add_argument("-f","--forward", help="forward primer", type = str)
  parser.add_argument("-r","--reverse", help="reverse primer", type = str)
  parser.add_argument("-m","--mismatches", default = 0, help="number of mismatches", type = str)
  parser.add_argument("--minlength", default = 100, help="minimum product length", type = int)
  parser.add_argument("--maxlength", default = 1000, help="maximum product length", type = int)
  parser.add_argument("-p","--primerfile", help="file with primers in csv format")
  return parser.parse_args()

#parser.add_argument("-p","--primerfile", help="file with primers in csv format")
#return parser.parse_args()

#regex.search would look for any of the words in [] in a string (e.g. IN [0]: re.search(r'[a-z]+', 'Hello')
																#     OUT [0]: <re.Match object; span=(1, 5), match='ello'> (this is case sensitive))
def replace_ambiguity_codes(sequence):
	codes={"A":"A", "C":"C", "G":"G", "T":"T", "R":"[AG]", "S":"[GC]", "B":"[CGT]", "Y":"[CT]", "W":"[AT]", "D":"[AGT]", "K":"[GT]", "N":"[ACGT]", "H":"[ACT]", "M":"[AC]", "V":"[ACG]", "X": "[ACGT]"}
	regexlist=[]
	for base in sequence:
		if not base in codes:
			print("Unidentified nucleotide code in primer:", base)
			sys.exit()
		else:
			regexlist.append(codes[base])
	
	return ''.join(regexlist)


def revcom(sequence):
  comp = {"A": "T", "C": "G", "G": "C", "T" :"A", "[":"]", "]":"[", "(":")", ")":"("}
  rev = ""
  for i in sequence:
    rev += i.replace(i, comp[i])
    rev_com = rev[::-1]
  return rev_com

if __name__ == "__main__":
  
  args = main()

  mis = args.mismatches
  mismatches = f"{{s<={mis}}}" #Fuzzy matching (mismatch). From regex module f"{{s<={mis}}}" substitution

  print(f"{mis} mismatches")

  #Read sequence file

  print("Reading sequence file...")

  sequence = args.sequence

  fasta_file = open(sequence, "r")
  num_lines = sum(1 for line in fasta_file)
  fasta_file.close() #Look for a shortcut

  fasta_file = open(sequence, "r")
  seq_n, count = "", 0
  seqs = []

  for s in fasta_file:
    s = s.strip()
    count += 1
    if s[0]==">":
      if seq_n != "":
        seqs.append([seq_n, seq])
      seq_n, seq = s, ""
    else:
      seq += s
      if count == num_lines:
        seqs.append([seq_n, seq])
        
  print("Done")

  if args.primerfile is None:

    print("Checking primers...")

    #Forward primer

    primerF = args.forward
    primerF = primerF.replace(" ", "")
    primerF = [*primerF]
    primerF = replace_ambiguity_codes(primerF)
    primerF = f"({primerF})" #Brackets are necessary for fuzzy matching

    print(f"Forward primer: {primerF}")

    patternF = regex.compile(primerF + mismatches)

	  #Forward primer reverse complement

    primerFRv = revcom(primerF)
    patternFRv = regex.compile(primerFRv + mismatches)

	  #Reverse primer

    primerR = args.reverse
    primerR = primerR.replace(" ", "")
    primerR = [*primerR]
    primerR = replace_ambiguity_codes(primerR)
    primerR = f"({primerR})"

    print(f"Reverse primer: {primerR}")

    patternR = regex.compile(primerR + mismatches)

	  #Reverse primer reverse complement

    primerRRv = revcom(primerR)
    patternRRv = regex.compile(primerRRv + mismatches)

    print("Done")

    minlen = args.minlength
    maxlen = args.maxlength

    print(f"{minlen}bp minimum product length\n{maxlen}bp maximum product length")

    print("Starting analysis...")

    for a,b in seqs:
      try:
        matchF = patternF.finditer(b)
        matchF = [[i.end() + 1, i.fuzzy_counts[0]] for i in matchF]
        matchRRV = patternRRv.finditer(b)
        matchRRV = [[i.start() - 1, i.fuzzy_counts[0]] for i in matchRRV]
        matchFR = list(product(matchF, matchRRV))
        for d,e in matchFR:
          start = d[0]
          end = e[0]
          mis_F = d[1]
          mis_R = e[1]
          length = end - start
          length = abs(length)
          if length < maxlen and length > minlen:
            print(a, start, mis_F, end, mis_R, length, sep = "\t")

	  		#For reverse complement, needs testing

        matchFRv = patternFRv.finditer(b)
        matchFRv = [[i.end() + 1, i.fuzzy_counts[0]] for i in matchFRv]
        matchR = patternR.finditer(b)
        matchR = [[i.start() - 1, i.fuzzy_counts[0]] for i in matchR]
        matchFRRv = list(product(matchFRv, matchR))
        for f,g in matchFRRv:
          start = f[0]
          end = g[0]
          mis_F = f[1]
          mis_R = g[1]
          length = end - start
          length = abs(length)
          if length < maxlen and length > minlen:
            print(a, start, mis_F, end, mis_R, length, sep = "\t")
      except:
        pass
      
  print("Done")

  if args.primerfile is not None:

    with open(args.primerfile) as primers:
      results = open("pcr_results.txt","w")

      for lines in primers:

        lines = lines.replace(" ","")
        lines = lines.split(",")

        experiment = lines[0]

        print (f"Experiment {experiment}")

        #Forward primer

        primerF = lines[1]
        primerF = [*primerF]
        primerF = replace_ambiguity_codes(primerF)
        primerF = f"({primerF})" #Brackets are necessary for fuzzy matching

        print(f"Forward primer: {primerF}")

        patternF = regex.compile(primerF + mismatches)

	      #Forward primer reverse complement

        primerFRv = revcom(primerF)
        patternFRv = regex.compile(primerFRv + mismatches)

        #Reverse primer

        primerR = lines[2]
        primerR = [*primerR]
        primerR = replace_ambiguity_codes(primerR)
        primerR = f"({primerR})"

        print(f"Reverse primer: {primerR}")

        patternR = regex.compile(primerR + mismatches)

	      #Reverse primer reverse complement

        primerRRv = revcom(primerR)
        patternRRv = regex.compile(primerRRv + mismatches)

        print("Done")

        minlen = lines [3]
        minlen = int(minlen)
        maxlen = lines [4]
        maxlen = int(maxlen)

        print(f"{minlen}bp minimum product length")
        print(f"{maxlen}bp maximum product length")

        print("Starting analysis...")

        for a,b in seqs:
          try:
            matchF = patternF.finditer(b)
            matchF = [[i.end() + 1, i.fuzzy_counts[0]] for i in matchF]
            matchRRV = patternRRv.finditer(b)
            matchRRV = [[i.start() - 1, i.fuzzy_counts[0]] for i in matchRRV]
            matchFR = list(product(matchF, matchRRV))
            for d,e in matchFR:
              start = d[0]
              end = e[0]
              mis_F = d[1]
              mis_R = e[1]
              length = end - start
              length = abs(length)
              if length < maxlen and length > minlen:
                print(experiment, a, start, mis_F, end, mis_R, length, sep = "\t", file=results)
                #results.write(experiment, a, start, mis_F, end, mis_R, length, sep = "\t")

	  		    #For reverse complement, needs testing

            matchFRv = patternFRv.finditer(b)
            matchFRv = [[i.end() + 1, i.fuzzy_counts[0]] for i in matchFRv]
            matchR = patternR.finditer(b)
            matchR = [[i.start() - 1, i.fuzzy_counts[0]] for i in matchR]
            matchFRRv = list(product(matchFRv, matchR))
            for f,g in matchFRRv:
              start = f[0]
              end = g[0]
              mis_F = f[1]
              mis_R = g[1]
              length = end - start
              length = abs(length)
              if length < maxlen and length > minlen:
                print(experiment, a, start, mis_F, end, mis_R, length, sep = "\t", file=results)
                #results.write(experiment, a, start, mis_F, end, mis_R, length, sep = "\t")
          except:
            pass
      
      results.close()

      print("Done")


        #print(experiment, primerR, primerF, minlen, maxlen)
  
  fasta_file.close()