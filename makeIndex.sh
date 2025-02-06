!#/bin/bash

# make index of vcff regions from csv file with locations

CSVFile=10819-non-sample-S3-locations.r5.csv
IndexFile=vcfs.index.r5.txt

MakeIndex='{ split($1,w,"."); pos=w[11]; split(pos,c,":"); split(c[2],p,"-"); print c[1],p[1],p[2],$1 }'

cat $CSVFile | grep bgz | grep -v csi | awk "$MakeIndex" | sort -k1,1 -k2,2n> $IndexFile

