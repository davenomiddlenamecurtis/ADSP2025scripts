#!/bin/bash
# for single chromosome files e.g. ukb23158_c9_b0_v1.bim
# I think these can be run in parallel
# file names use X but internally use 23

set -x

root=$1

mult=no # set --pick_allele_gene

if [ -z "$mult" ]
then
  mult=yes
# annotations for multiple transcripts
fi

if [ -z "$root" ]
then
  echo Usage: $0 vcfRoot [ target is vcfRoot.vars.vcf.gz, can also first export mult=no for only one transcript and X=no if no X chromosome data ]
  exit
fi

PICK=
MULT=.mult

if [ "$mult" == no ]
then
  PICK="--pick_allele_gene"
  MULT=
fi

# export PERL5LIB=/usr/lib/perl64:$PERL5LIB
OLDPERL5LIB=$PERL5LIB
export PERL5LIB=/usr/lib/perl64:/share/apps/ensembl-vep-97/plugins/
inputFile=$root.vars.vcf
outputFile=$root.annot.vcf
zcat $root.vcf.gz | grep -v '^#' | cut -f 1-5 > $inputFile
if [ ! -e $outputFile.done ]
then
  perl /share/apps/ensembl-vep-97/vep \
    --dir /share/apps/ensembl-vep-97/plugins \
	--input_file $inputFile \
    --synonyms ~/vep/chr_synonyms.txt \
	--cache --dir /cluster/project9/bipolargenomes/vepcache --merged --force_overwrite \
	--sift b --polyphen b --assembly GRCh38 --format vcf \
	--fasta /cluster/project9/bipolargenomes/vepcache/homo_sapiens_merged/97_GRCh38 \
	--canonical --regulatory \
	--plugin AlphaMissense,file=/share/ref/VEP/AlphaMissense/AlphaMissense_hg38.tsv.gz,cols=all \
	--vcf --output_file $outputFile $PICK
  if [ -s ${outputFile}_summary.html ]
  then
    bgzip $outputFile
	tabix -p vcf $outputFile.gz
	# rm $inputFile
	echo done > $outputFile.done
  else
    echo  ${outputFile}_summary.html is zero length
  fi
fi

