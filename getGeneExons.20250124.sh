#!/bin/bash
set -x

set +e
# I am baffled as to why the script silently finishes after successfully downloading a file
# Clutching at straws I will guess that maybe aws returns an error

# make a local vcf containing all the exons of a specified ADSP gene
# need about 8 hours to tabix the downloaded vcf and more time to download it and do other stuff

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
  mkdir $tmpDir
  cd $tmpDir
  geneArgs=`grep -w $gene $refGeneFile | awk "$extractGeneCoords"`
  name=$gene
  coords=($geneArgs)
  chr=${coords[0]}
  start=${coords[1]}
  end=${coords[2]}
  awk -v chr=$chr -v start=$start -v end=$end "$getFiles" $index > $fileList
  cat $fileList | while read f
  do
    if [ ! -e $f ]
    then
      /share/apps/aws-cli-tools/v2/current/bin/aws s3 cp $bucketPrefix/$f .
      /share/apps/aws-cli-tools/v2/current/bin/aws s3 cp $bucketPrefix/$f.csi .
    fi
  done
  cat $fileList | while read f
  do
    for ff in $f $f.csi
	do
	  if [ ! -e $ff ]
	  then
	    echo $ff did not download
		exit
	  fi
	  if [ ! -s $ff ]
	  then
	    echo $ff was zero size, deleting...
		rm $ff
		exit
	  fi
	done
  done

  cat $fileList | while read f
  do
      geneVarAssoc --arg-file $argFile --case-file $f --gene $gene 
      mv gva.$gene.case.1.vcf $f.$gene.exons.vcf
  done
  firstFile=`head -n 1 $fileList`
  wco=`wc -l $fileList` 
  words=($wco)
  nFiles=${words[0]}	
  if [ $nFiles -gt 1 ]
  then
    grep '^#' $firstFile.$gene.exons.vcf > ADSP2025.exons.$gene.vcf
    cat $fileList | while read f
    do
      grep -v '^#' $f.$gene.exons.vcf >> ADSP2025.exons.$gene.vcf
    done
  else
    mv $firstFile.$gene.exons.vcf ADSP2025.exons.$gene.vcf
  fi

  bgzip ADSP2025.exons.$gene.vcf
  mv ADSP2025.exons.$gene.vcf.gz ../$exonsVCF
  cd $wd
  tabix -p vcf $exonsVCF
fi

if [ -e $exonsVCF.tbi -a -s $exonsVCF.tbi ]
then
  rm -r $tmpDir
fi