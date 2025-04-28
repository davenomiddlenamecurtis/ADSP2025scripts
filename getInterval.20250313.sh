#!/bin/bash
set -x

set +e
# I am baffled as to why the script silently finishes after successfully downloading a file
# Clutching at straws I will guess that maybe aws returns an error

# make a local vcf containing a specific interval from the ADSP WGS data

chr=$1
start=$2
end=$3
interval=$chr.$start.$end

wd=/cluster/project9/bipolargenomes/ADSP2025/intervals
tmpDir=/cluster/project9/bipolargenomes/ADSP2025/intervals/temp.$interval
downloadedVCF=ADSP2025downloaded.$interval.vcf.gz # this should be OK to use for analyses
intervalVCF=ADSP2025.interval.$interval.vcf.gz 
fileList=todownload.$interval.lst

index=~/ADSP2025/vcfs.index.r5.txt
bucketPrefix=s3://wanglab-dss-tier0/distribution/adsp/genotype/ALL/fsa000116/preview
bcftools=/share/apps/genomics/bcftools-1.9/bin/bcftools

getFiles='{ if ($1==chr) { if ($2>start) { print last } ; last=$4; if ($2 > end) { exit } }  } END { print last }'

pushd $wd

if [ ! -e $intervalVCF ]
then
  mkdir $tmpDir
  cd $tmpDir
  name=$interval
  awk -v chr=chr$chr -v start=$start -v end=$end "$getFiles" $index > $fileList
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
      tabix -h $f chr$chr:$start-$end > $f.$interval.interval.vcf
  done
  firstFile=`head -n 1 $fileList`
  wco=`wc -l $fileList` 
  words=($wco)
  nFiles=${words[0]}	
  if [ $nFiles -gt 1 ]
  then
    grep '^#' $firstFile.$interval.interval.vcf > ADSP2025.interval.$interval.vcf
    cat $fileList | while read f
    do
      grep -v '^#' $f.$interval.interval.vcf >> ADSP2025.interval.$interval.vcf
    done
  else
    mv $firstFile.$interval.interval.vcf ADSP2025.interval.$interval.vcf
  fi

  bgzip ADSP2025.interval.$interval.vcf
  mv ADSP2025.interval.$interval.vcf.gz ../$intervalVCF
  cd $wd
  tabix -p vcf $intervalVCF
fi

if [ -e $intervalVCF.tbi -a -s $intervalVCF.tbi ]
then
  rm -r $tmpDir
fi