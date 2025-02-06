#!/bin/bash
# run this on little so cluster does not get jammed

set -x
wd=/cluster/project9/bipolargenomes/ADSP2025/genes
TopTestsFile=/home/rejudcu/ADSP2024/genes/ADSP.topTests.20250120.txt
maxToRun=10

pushd $wd

running=0
tail -n +2 $TopTestsFile | while read gene rest
do
  exonsVCF=ADSP2025.exons.$gene.vcf.gz
  if [ ! -e $exonsVCF ]
  then
    if [ -e temp.$gene ] # previous run still running
	then
	  continue
	fi
    running=$(( running + 1 ))
	nohup bash ~/ADSP2025/ADSP2025scripts/getGeneExons.20250124.sh $gene  &> download.$gene.log &
	if [ $running -ge $maxToRun ]
	then
	  exit
	fi
  fi
done
