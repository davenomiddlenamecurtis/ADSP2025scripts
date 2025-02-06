#!/bin/bash

# run one GVA analysis on one gene for ADSP 2025

model=$1
gene=$2
test=$3

wd=/cluster/project9/bipolargenomes/ADSP2025/genes
tmpdir=/cluster/project9/bipolargenomes/ADSP2025/genes/$gene
resultsdir=/cluster/project9/bipolargenomes/ADSP2025/$model/results

cd $wd
if [ ! -e ADSP2025.exons.$gene.vcf.gz.tbi ]
then
  echo ADSP2025.exons.$gene.vcf.gz has not been created yet, exiting...
  exit
  # bash /home/rejudcu/ADSP2025/ADSP2025scripts/getGeneExons.sh $gene
fi

if [ ! -e ADSP2025.exons.$gene.annot.vcf.gz.tbi ]
then
  bash /home/rejudcu/ADSP2025/ADSP2025scripts/annotate.geneVCF.with.AlphaMissense.sh ADSP2025.exons.$gene
fi

if [ ! -e $resultsdir/$model.$gene.sco ]
then
mkdir $tmpdir
cd $tmpdir

ln -s /cluster/project9/bipolargenomes/ADSP2025/genes/ADSP2025.exons.$gene.annot.vcf.gz myAnnots.gz
ln -s /cluster/project9/bipolargenomes/ADSP2025/genes/ADSP2025.exons.$gene.annot.vcf.gz.tbi myAnnots.gz.tbi

if [ .$3 -eq . ]
  geneVarAssoc --arg-file /home/rejudcu/ADSP2025/ADSP2025scripts/pars/gva.$model.arg --gene $gene --case-file /cluster/project9/bipolargenomes/ADSP2025/genes/ADSP2025.exons.$gene.vcf.gz
else
  geneVarAssoc --arg-file /home/rejudcu/ADSP2025/ADSP2025scripts/pars/gva.$model.arg --gene $gene --case-file /cluster/project9/bipolargenomes/ADSP2025/genes/ADSP2025.exons.$gene.vcf.gz --testfile $test 
fi
cp *.s?o $resultsdir
  if [ -e $resultsdir/$model.$gene.sco -a -s $resultsdir/$model.$gene.sco ]
  then
    cd .. 
	rm -r $tmpdir
  fi
fi
