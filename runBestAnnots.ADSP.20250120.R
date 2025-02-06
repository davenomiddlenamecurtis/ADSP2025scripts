#!/share/apps/R-3.6.1/bin/Rscript

# script to take summary file, select most significant gene/annotation results
# then run scoreassoc on second dataset with those pairs

# submit this script to run all analyses

TopTestsFile="/home/rejudcu/ADSP2024/genes/ADSP.topTests.20250120.txt"
NewModel="ADSP.best.annot.20250120"
CountModel="ADSP.raw.counts.20250131"
NewArgFile=sprintf("/home/rejudcu/ADSP2024/ADSPscripts.2024/pars/gva.%s.arg",NewModel)
AllResultsFile="ADSP2025.AllResults.txt"
BestResultsFile="ADSP2025.BestResults.txt"

wd="/home/rejudcu/ADSP2025"
setwd(wd)

TopTests=data.frame(read.table(TopTestsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
AllDone=TRUE
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sco",NewModel,NewModel,gene))) {
		if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/genes/ADSP2025.exons.%s.vcf.gz.tbi",gene))) {
			print(sprintf("No exon vcf for %s",gene))
			next
		}
		commStr=sprintf("subComm.sh bash /home/rejudcu/ADSP2025/ADSP2025scripts/runOneGene.sh %s %s %s",NewModel,gene,TopTests[r,2])
		print(commStr)
		system(commStr)
		AllDone=FALSE
	}
}

if (!AllDone) {
	quit()
}

for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	lines=readLines(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sao",NewModel,NewModel,gene))
	for (ll in 1:length(lines)) {
		words=strsplit(lines[ll],"\\s+")[[1]]
		if (length(words>0)) {
			if (words[1]=="tMLP") {
				TopTests$WGSMLP[r]=words[3]
				break
			}
		}
	}
}
write.table(TopTests, AllResultsFile, row.names=FALSE, quote=FALSE, col.names = TRUE,sep="\t")

first=TRUE
TopTests=data.frame(read.table(AllResultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
for (r in 1:nrow(TopTests)) {
	if (TopTests$WGSMLP[r]>-log10(0.05/nrow(TopTests))) {
		if (first) {
			Best=TopTests[r,]
			first=FALSE
		} else {
			Best=rbind(Best,TopTests[r,])
		}
	}
}

write.table(Best, BestResultsFile, row.names=FALSE, quote=FALSE, col.names = TRUE,sep="\t")

AllDone=TRUE
Best=data.frame(read.table(BestResultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
for (r in 1:nrow(Best)) {
	gene=Best[r,1]
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sco",CountModel,CountModel,gene))) {
		AllDone=FALSE
		commStr=sprintf("subComm.sh bash /home/rejudcu/ADSP2025/ADSP2025scripts/runOneGene.sh %s %s %s",CountModel,gene,Best[r,2])
		print(commStr)
		system(commStr)
	}
	
}

