#!/share/apps/R-3.6.1/bin/Rscript

# script to take summary file, select most significant gene/annotation results
# then run scoreassoc on second dataset with those pairs

# submit this script to run all analyses

TopTestsFile="/home/rejudcu/ADSP2024/genes/ADSP.topTests.20250120.txt"
NewModel="ADSP.best.annot.20250120"
CountsModel="ADSP.raw.counts.20250131"
NewArgFile=sprintf("/home/rejudcu/ADSP2024/ADSPscripts.2024/pars/gva.%s.arg",NewModel)
AllResultsFile="ADSP2025.AllResults.txt"
BestResultsFile="ADSP2025.BestResults.txt"

TREM2NoSNPsModel="ADSP.TREM2NoSNPs.20250416"
TREM2JustSNPsModel="ADSP.TREM2JustSNPs.20250416"
NoSNPsArgFile=sprintf("/home/rejudcu/ADSP2024/ADSPscripts.2024/pars/gva.%s.arg",TREM2NoSNPsModel)
JustSNPsArgFile=sprintf("/home/rejudcu/ADSP2024/ADSPscripts.2024/pars/gva.%s.arg",TREM2JustSNPsModel)

CategoryNames=c("Five prime UTR","InDel etc","Intronic etc","LOF","Protein altering","Splice region","Synonymous","Three prime UTR")
MainWeights=c("VEPWeight","FivePrime","InDelEtc","IntronicEtc","LOF","ProteinAltering","SpliceRegion","Synonymous","ThreePrime")
ExtraWeightsFile="/home/rejudcu/ADSP2025/ADSP2025scripts/pars/extraWeights.20250122.txt"
ExtraWeights=data.frame(read.table(ExtraWeightsFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))[,1]
DefaultTestFile="/home/rejudcu/ADSP2025/ADSP2025scripts/pars/test.counts.20250110.tst"

ScoreThreshold=0.35

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

# Now run TREM2 without associated SNPs
AllDone=TRUE
for (r in 1:nrow(TopTests)) {
	gene=TopTests[r,1]
	if (gene!="TREM2") {
		next
	}
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sco",TREM2NoSNPsModel,TREM2NoSNPsModel,gene))) {
		commStr=sprintf("subComm.sh bash /home/rejudcu/ADSP2025/ADSP2025scripts/runOneGene.sh %s %s %s",TREM2NoSNPsModel,gene,TopTests[r,2])
		print(commStr)
		system(commStr)
		AllDone=FALSE
	}
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sco",TREM2JustSNPsModel,TREM2JustSNPsModel,gene))) {
		commStr=sprintf("subComm.sh bash /home/rejudcu/ADSP2025/ADSP2025scripts/runOneGene.sh %s %s /home/rejudcu/ADSP2025/ADSP2025scripts/pars/test.counts.for.TREM2.JustSNPs.tst",TREM2JustSNPsModel,gene)
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
	Best[r,2]=gsub("/home/rejudcu/ADSP2024/ADSPscripts.2024/pars/test.annot.","",Best[r,2])
	Best[r,2]=gsub(".tst","",Best[r,2])
	Best[r,2]=gsub(".NoLOF","",Best[r,2]) # may or may not be there
	test=Best[r,2]
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sco",CountsModel,CountsModel,gene))) {
		AllDone=FALSE
		TestFile=sprintf("/home/rejudcu/ADSP2025/ADSP2025scripts/pars/test.counts.for.%s.tst",gene)
		CommStr=sprintf("cp %s %s",DefaultTestFile,TestFile)
		system(CommStr)
		write(sprintf("%s 0.0 0 1",test),TestFile,append=TRUE)
		commStr=sprintf("subComm.sh bash -x /home/rejudcu/ADSP2025/ADSP2025scripts/runOneGene.sh %s %s %s",CountsModel,gene,TestFile)
		print(commStr)
		system(commStr)
	}
}

if (!AllDone) {
	quit()
}

for (r in 1:nrow(Best)) {
	gene=Best[r,1]
	test=Best[r,2]
	SaoFile=sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sao",CountsModel,CountsModel,gene)
		SaoTable=na.omit(data.frame(read.table(SaoFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE)))
	colnames(SaoTable)[16:24]=MainWeights
	colnames(SaoTable)[25]="AM_prediction"
	colnames(SaoTable)[26]="AM_score"
	colnames(SaoTable)[27:(27+length(ExtraWeights)-1)]=ExtraWeights
	GlmTable=data.frame(read.table(SaoFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))
	First=which("L1"==GlmTable[,1])[[1]]
	GlmTable=GlmTable[-First:-1,1:4]
	colnames(GlmTable)=GlmTable[1,]
	Cats=c(MainWeights[-1],test)
	GeneResults=data.frame(matrix(ncol=8,nrow=length(Cats),0))
	colnames(GeneResults)=c("Category","NumVars","TotCountControls","MeanCountControls","TotCountCases","MeanCountCases","OR","SLP")
	for (c in 1:length(Cats)) {
		Cat=Cats[c]
		GeneResults[c,1]=Cat
		if (c==length(Cats)) {
			Tab=SaoTable[as.numeric(SaoTable[,Cat])>ScoreThreshold,]
		} else {
			Tab=SaoTable[as.numeric(SaoTable[,Cat])!=0,]
		}
		for (cc in c(2,4,6,8,10,12)) {
			Tab[,cc]=as.numeric(Tab[,cc])
		}
		GeneResults$NumVars[c]=nrow(Tab)
		if (nrow(Tab)>0) {
			GeneResults$TotCountControls[c]=sum(Tab$contAB)+2*sum(Tab$contBB)
			GeneResults$MeanCountControls[c]=GeneResults$TotCountControls[c]/((sum(Tab$contAA)+sum(Tab$contAB)+sum(Tab$contBB))/nrow(Tab))
			GeneResults$TotCountCases[c]=sum(Tab$caseAB)+2*sum(Tab$caseBB)
			GeneResults$MeanCountCases[c]=GeneResults$TotCountCases[c]/((sum(Tab$caseAA)+sum(Tab$caseAB)+sum(Tab$caseBB))/nrow(Tab))
			gr=which(Cat==GlmTable$beta)[[1]]
			b=as.numeric(GlmTable$value[gr])
			SE=as.numeric(GlmTable$SE[gr])
			GeneResults$OR[c]=sprintf("%.2f (%.2f-%.2f)",exp(b),exp(b-2*SE),exp(b+2*SE))
			GeneResults$SLP[c]=log10(2*pnorm(abs(as.numeric(GlmTable$z[gr])),lower.tail=FALSE))
			if (b>0) {
				GeneResults$SLP[c]=GeneResults$SLP[c]* -1
			}
		}
	}
	print(GeneResults)
	GeneResults[1:length(CategoryNames),1]=CategoryNames
	write.table(GeneResults,sprintf("%s.counts.txt",gene),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}


