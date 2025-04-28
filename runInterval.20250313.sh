#!/bin/bash

# run one GVA analysis on one gene for ADSP 2025

set +e

model=$1
chr=$2
start=$3
end=$4
test=$5

interval=$chr.$start.$end

wd=/cluster/project9/bipolargenomes/ADSP2025/intervals
tmpdir=/cluster/project9/bipolargenomes/ADSP2025/intervals/$interval
resultsdir=/cluster/project9/bipolargenomes/ADSP2025/$model/results

cd $wd
if [ ! -e ADSP2025.interval.$interval.vcf.gz.tbi ]
then
  echo ADSP2025.interval.$interval.vcf.gz has not been created yet, exiting...
  exit
fi

# if [ ! -e ADSP2025.interval.$interval.annot.vcf.gz.tbi ]
# then
#   bash /home/rejudcu/ADSP2025/ADSP2025scripts/annotate.geneVCF.with.AlphaMissense.sh ADSP2025.interval.$interval
# fi

if [ ! -e $resultsdir/$model.$interval.sco ]
then
mkdir $tmpdir
cd $tmpdir

# ln -s /cluster/project9/bipolargenomes/ADSP2025/intervals/ADSP2025.interval.$interval.annot.vcf.gz myAnnots.gz
# ln -s /cluster/project9/bipolargenomes/ADSP2025/intervals/ADSP2025.interval.$interval.annot.vcf.gz.tbi myAnnots.gz.tbi

echo $chr:$start-$end > interval.$interval.lst

if [ .$test == . ]
then
  intVarAssoc --arg-file /home/rejudcu/ADSP2025/ADSP2025scripts/pars/iva.$model.arg --interval-list-file interval.$interval.lst --case-file /cluster/project9/bipolargenomes/ADSP2025/intervals/ADSP2025.interval.$interval.vcf.gz --test-name $model.$interval
else
  intVarAssoc --arg-file /home/rejudcu/ADSP2025/ADSP2025scripts/pars/iva.$model.arg --interval-list-file interval.$interval.lst --case-file /cluster/project9/bipolargenomes/ADSP2025/intervals/ADSP2025.interval.$interval.vcf.gz --test-name $model.$interval --testfile $test 
fi
cp *.s?o $resultsdir
  if [ -e $resultsdir/$model.$interval.sco -a -s $resultsdir/$model.$interval.sco ]
  then
    cd .. 
	rm -r $tmpdir
  fi
fi
