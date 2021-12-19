library(readxl)
library(matlabr)

source(<path to the file "counts_to_TPM.R">)
source(<path to the file "EMT_score_func.R">) 

## Folder with all the raw read counts
setwd("F:/Jolly_Lab_IISc/EMT_Scoring_RNA_seq/RawCounts_Data")

fileList = list.files(pattern = "*.tsv")
num = length(fileList)
gseIDs = sapply(strsplit(fileList, split = "_"), function(x) x[1])

for(dataNum in 1:num){ 
	counts = read.delim(fileList[dataNum], header = TRUE, sep = "\t")
    log2_TPM = countToTpm(counts,gseIDs[dataNum])
	MA_val = rnaToMA(log2_TPM, gseIDs[dataNum])
    EMT76GS_Score = EMT76GS(MA_val, gseIDs[dataNum])
	KS_score=KSScore(MA_val, gseIDs[dataNum])
}