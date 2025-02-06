!#/bin/bash

# make a local vcf containing a specified ADSP gene

gene=$1
wd=/cluster/project9/bipolargenomes/ADSP2025/genes
tmpDir=/cluster/project9/bipolargenomes/ADSP2025/genes/temp.$gene
downloadedVCF=ADSP2025downloaded.$gene.vcf.bgz # this should be OK to use for analyses
exonsVCF=ADSP2025.exons.$gene.vcf.bgz # potentially write code to make this if needed
fileList=todownload.$gene.lst

refGeneFile=~/reference38/refseqgenes.hg38.20191018.sorted.onePCDHG.txt
index=~/ADSP2025/vcfs.index.r5.txt
bucketPrefix=s3://wanglab-dss-tier0/distribution/adsp/genotype/ALL/fsa000116/preview
bcftools=/share/apps/genomics/bcftools-1.9/bin/bcftools

extractGeneCoords=' BEGIN { start=300000000; end=0 } { chr= $3; if ($5<start) start=$5; if ($6>end) end=$6 } END { print chr, start, end }'
getFiles='{ if ($1==chr) { if ($2>start) { print last } ; last=$4; if ($2 > end) { exit } } }'

pushd $wd

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
    fi
  done
  firstFile=`head -n 1 $fileList`
  if [ `wc -l $fileList` -gt 1 ]
  then
    $bcftools view --header-only $first > downloadedVCF.$gene.vcf
    cat $fileList | while read f
    do
      $bcftools view --no-header $f >> downloadedVCF.$gene.vcf
    done
    bgzip downloadedVCF.$gene.vcf
    mv downloadedVCF.$gene.vcf.gz ../$downloadedVCF
  else
    mv $firstFile ../$downloadedVCF
  fi
  cd $wd
  tabix -p vcf $downloadedVCF
  rm -rf $tmpDir
fi
