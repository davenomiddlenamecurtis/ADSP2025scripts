#!/bin/bash
set -x
# just download needed vcfs and their indices

gene=$1
wd=/cluster/project9/bipolargenomes/ADSP2025/genes
tmpDir=/cluster/project9/bipolargenomes/ADSP2025/genes/temp.$gene
downloadedVCF=ADSP2025downloaded.$gene.vcf.gz # this should be OK to use for analyses
exonsVCF=ADSP2025.exons.$gene.vcf.gz 
fileList=todownload.$gene.lst

refGeneFile=~/reference38/refseqgenes.hg38.20191018.sortedAndCleaned.onePCDHG.txt
index=~/ADSP2025/vcfs.index.r5.txt
bucketPrefix=s3://wanglab-dss-tier0/distribution/adsp/genotype/ALL/fsa000116/preview
bcftools=/share/apps/genomics/bcftools-1.9/bin/bcftools
argFile=/home/rejudcu/pars/gva.ADSP2025.extractExons.arg

extractGeneCoords=' BEGIN { start=300000000; end=0 } { chr= $3; if ($5<start) start=$5; if ($6>end) end=$6 } END { print chr, start, end }'
getFiles='{ if ($1==chr) { if ($2>start) { print last } ; last=$4; if ($2 > end) { exit } }  } END { print last }'

pushd $wd

if [ ! -e $exonsVCF ]
then

if [ ! -e $downloadedVCF ]
then
  mkdir $tmpDir
  cd $tmpDir
  geneArgs=`grep -w $gene $refGeneFile | awk "$extractGeneCoords"`
  name=$gene
  coords=($geneArgs)
  chr=${coords[0]}
  start=${coords[1]}
  end=${coords[2]}
  awk -v chr=$chr -v start=$start -v end=$end "$getFiles" $index > $fileList
  allDownloaded=yes
  cat $fileList | while read f
  do
    if [ ! -e $f ]
    then
      allDownloaded=no
      /share/apps/aws-cli-tools/v2/current/bin/aws s3 cp $bucketPrefix/$f .
	  /share/apps/aws-cli-tools/v2/current/bin/aws s3 cp $bucketPrefix/$f.csi .
    fi
  done
fi
fi

