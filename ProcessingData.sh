### ****Survival Analysis of TCGA survival data ( Overall survival) within the context of metabolism****


### Download survival data of TCGA cancer cohort from UCSC Xena browser (https://xenabrowser.net/)
### Match the samples (patient ID) in the ssGSEA scores with that of survival data
### Divide the samples into EMT high and EMT low, Glycolysis high and Glycolysis low, OXPHOS high and OXPHOS low groups based on the median of the respective ssGSEA scores of the samples
### For convenience use GNU datamash - a nice tool for "command-line statistical operations"  to calculate the median easily
### Plot KM plots using "survival" R package as described in survival.r
### Calculate Hazard ratio using Cox regression using the “coxph” function and plot forest plots and log2 Hazard ratio plots using ggplots

### With the combination of EMT, OXPHOS, and Glycolysis, the samples were divided into eight possible combinations: E+G+O+, E+G+O-, E+G-O+, E+G-O-, E-G+O+, E-G+O-, E-G-O+ and E-G-O-

for i in $(ls *_survData.tsv | sed -e 's/_survData.tsv//' | sort -u)

do
	awk '{ if ($4=="high" && $5=="high" && $6=="high") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E+G+O+""\t""1"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E+G+O+_survData.tsv
	awk '{ if ($4=="high" && $5=="high" && $6=="low") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E+G+O-""\t""2"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E+G+O-_survData.tsv
	awk '{ if ($4=="high" && $5=="low" && $6=="high") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E+G-O+""\t""3"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E+G-O+_survData.tsv
	awk '{ if ($4=="high" && $5=="low" && $6=="low") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E+G-O-""\t""4"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E+G-O-_survData.tsv
	awk '{ if ($4=="low" && $5=="high" && $6=="high") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E-G+O+""\t""5"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E-G+O+_survData.tsv
	awk '{ if ($4=="low" && $5=="high" && $6=="low") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E-G+O-""\t""6"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E-G+O-_survData.tsv
	awk '{ if ($4=="low" && $5=="low" && $6=="high") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E-G-O+""\t""7"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E-G-O+_survData.tsv
	awk '{ if ($4=="low" && $5=="low" && $6=="low") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""E-G-O-""\t""8"}' ${i}_survData.tsv > EMT_Glyco_OXPHOS/${i}_E-G-O-_survData.tsv
	
done


### Plot KM plots using "survival" R package as described in survival.r
### Calculate Hazard ratio using Cox regression using the “coxph” function and plot forest plots using ggplots
### Plot Heatmap plots using "heatmap" function and using RdBu color scale from the RcolorBrewer package.
