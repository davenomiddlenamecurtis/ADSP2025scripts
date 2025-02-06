#!/share/apps/R-3.6.1/bin/Rscript
InCSV="ADSPIntegratedPhenotypes_DS_2024.11.22.csv"
APOEFile="ADSPIntegratedPhenotypes_DS.withAPOE.2024.11.22.txt"
wd="/home/rejudcu/ADSP2025"

setwd(wd)
phenos=data.frame(read.csv(InCSV,header=TRUE, as.is=TRUE, sep=",", stringsAsFactors=FALSE))
phenos$APOE=phenos$APOE_WGS
phenos$APOE[is.na(phenos$APOE)]=phenos$APOE_reported[is.na(phenos$APOE)]
write.table(phenos, APOEFile, sep = "\t", row.names = FALSE, quote = FALSE)

# need to get doses for e3 and e4
APOE=as.character(phenos$APOE)
phenos$DoseAPOE3=lengths(regmatches(APOE, gregexpr("3", APOE)))
phenos$DoseAPOE4=lengths(regmatches(APOE, gregexpr("4", APOE)))
phenos=na.omit(phenos[c("SampleID","Sex","DX_harmonized","APOE","DoseAPOE3","DoseAPOE4")])

colnames(phenos)[1]="IID"
write.table(phenos[c("IID","DX_harmonized")], "ADSP.PrevAD.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(phenos[c("IID","Sex")], "ADSP.Sex.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(phenos[c("IID","DoseAPOE3","DoseAPOE4")], "ADSP.APOE.txt", sep = "\t", row.names = FALSE, quote = FALSE)
