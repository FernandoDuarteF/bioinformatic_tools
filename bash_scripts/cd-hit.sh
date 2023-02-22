#!/usr/bin/bash
#Cluster sequences into one representative sequence. Remove redundant sequences

catalogue=$1
out=$2
outdir=$PWD/$out

mkdir -p $outdir

cd $outdir

cd-hit -i $catalogue -o `basename $catalogue .aa.fa`_prot.fasta -c 0.95 -aL 0.9 -M 0 -T 0

#Extract ids

gawk -F'[#;=]' '/^>/ {print ">" $6; next}{print}' $catalogue > `basename $catalogue _prot.fasta`_ID_prot.fasta


exit

#Non redundant gene catalogue
#for f1 in ./annotations/*.aa.fa
#do
#  cd-hit -i $catalogue -o ./cd-hit/`basename $f1 .aa.fa`_prot.fasta -c 0.95 -aL 0.9 -M 0 -T 0
#  echo $f1
#done
#De la descripción de las secuencias, se extraen únicamente los IDs para posteriormente
#poder asociarlos a la abundancia
#cd cd-hit
#for f2 in ../cd-hit/*.prot.fa
#do
#  gawk -F'[#;=]' '/^>/ {print ">" $6; next}{print}' $f2 > `basename $f2 .prot.fa`_ID_prot.fasta
#  echo $f
#done
#cd ..
