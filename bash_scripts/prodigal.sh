#!/usr/bin/bash
#                       ----------------
#                       |   PRODIGAL   |
#                       ----------------
#Fast, reliable protein-coding gene prediction for prokaryotic genomes.

contig=$1
out=$2
outdir=$PWD/$out

mkdir -p $outdir

cd $outdir

prodigal -a `basename $contig .contigs.fa`.aa.fa \
         -d `basename $contig .contigs.fa`.nuc.fa \
         -i $contig -f gff -p meta > `basename $contig .contigs.fa`.gff

exit