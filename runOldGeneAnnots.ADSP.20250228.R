#!/share/apps/R-3.6.1/bin/Rscript

# submit this script to run all analyses

SummaryFile="/home/rejudcu/ADSP2024/genes/ADSP20240526.coding.summ.txt"
OldGenesTestsFile="/home/rejudcu/ADSP2025/ADSP.OldGenesTests.20250216.txt"

NewModel="ADSP.best.annot.20250120"
CountsModel="ADSP.raw.counts.20250131"
NewArgFile=sprintf("/home/rejudcu/ADSP2024/ADSPscripts.2024/pars/gva.%s.arg",NewModel)
AllResultsFile="ADSP2025.AllOldGenesResults.txt"
BestResultsFile="ADSP2025.BestOldGenesResults.txt"

CategoryNames=c("Five prime UTR","InDel etc","Intronic etc","LOF","Protein altering","Splice region","Synonymous","Three prime UTR")
MainWeights=c("VEPWeight","FivePrime","InDelEtc","IntronicEtc","LOF","ProteinAltering","SpliceRegion","Synonymous","ThreePrime")
ExtraWeightsFile="/home/rejudcu/ADSP2025/ADSP2025scripts/pars/extraWeights.20250122.txt"
ExtraWeights=data.frame(read.table(ExtraWeightsFile,header=TRUE,stringsAsFactors=FALSE,fill=TRUE))[,1]
DefaultTestFile="/home/rejudcu/ADSP2025/ADSP2025scripts/pars/test.counts.20250110.tst"

OldGenes=c("GRN", "GBA", "ADAM10", "FRMD8", "DDX1", "DNMT3L", "MORC1", "TGM2","ATP8B4","ABCA1" )
TestsFile="/home/rejudcu/ADSP2024/ADSPscripts.2024/pars/ADSP.annot.allTests.20240529.txt"

ScoreThreshold=0.35

wd="/home/rejudcu/ADSP2025"
setwd(wd)

Summary=data.frame(read.table(SummaryFile,header=TRUE,stringsAsFactors=FALSE))
Summary=Summary[Summary[,1] %in% OldGenes,]
# I need to add code to find best annotation

Summary$MaxMLP=apply(Summary[,11:ncol(Summary)], 1, max)

Tests=data.frame(read.table(TestsFile,header=FALSE,stringsAsFactors=FALSE))

OldGenesTests=data.frame(matrix("",nrow=nrow(Summary),ncol=3),stringsAsFactors=FALSE)
colnames(OldGenesTests)=c("Gene","TestFile","MaxMLP")
OldGenesTests[,1]=Summary[,1]
for (r in 1:nrow(Summary)) {
	for (t in 1:nrow(Tests)) {
		if (Summary[r,t+10]==Summary$MaxMLP[r]) {
			if (abs(Summary[r,6])>1.3) { # LOF
				OldGenesTests$TestFile[r]=Tests[t,1]
			} else {
# 				OldGenesTests$TestFile[r]=gsub(".tst",".NoLOF.tst",Tests[t,1])
				OldGenesTests$TestFile[r]=Tests[t,1] # for OldGenes, always include LOF
			}
			OldGenesTests$MaxMLP[r]=Summary$MaxMLP[r]
			break
		}
	}
}

write.table(OldGenesTests,OldGenesTestsFile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

AllDone=TRUE
for (r in 1:nrow(OldGenesTests)) {
	gene=OldGenesTests[r,1]
	if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sco",NewModel,NewModel,gene))) {
		if (!file.exists(sprintf("/cluster/project9/bipolargenomes/ADSP2025/genes/ADSP2025.exons.%s.vcf.gz.tbi",gene))) {
			print(sprintf("No exon vcf for %s",gene))
			next
		}
		commStr=sprintf("subComm.sh bash /home/rejudcu/ADSP2025/ADSP2025scripts/runOneGene.sh %s %s %s",NewModel,gene,OldGenesTests[r,2])
		print(commStr)
		system(commStr)
		AllDone=FALSE
	}
}

if (!AllDone) {
	quit()
}

for (r in 1:nrow(OldGenesTests)) {
	gene=OldGenesTests[r,1]
	lines=readLines(sprintf("/cluster/project9/bipolargenomes/ADSP2025/%s/results/%s.%s.sao",NewModel,NewModel,gene))
	for (ll in 1:length(lines)) {
		words=strsplit(lines[ll],"\\s+")[[1]]
		if (length(words>0)) {
			if (words[1]=="tMLP") {
				OldGenesTests$WGSMLP[r]=words[3]
				break
			}
		}
	}
}
write.table(OldGenesTests, AllResultsFile, row.names=FALSE, quote=FALSE, col.names = TRUE,sep="\t")

first=TRUE
OldGenesTests=data.frame(read.table(AllResultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
for (r in 1:nrow(OldGenesTests)) {
	if (OldGenesTests$WGSMLP[r]>-log10(0.05/nrow(OldGenesTests))) {
		if (first) {
			Best=OldGenesTests[r,]
			first=FALSE
		} else {
			Best=rbind(Best,OldGenesTests[r,])
		}
	}
}

write.table(Best, BestResultsFile, row.names=FALSE, quote=FALSE, col.names = TRUE,sep="\t")

AllDone=TRUE
Best=data.frame(read.table(AllResultsFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
# For OldGenes genes, I am going to get counts for all genes, not just best ones
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


