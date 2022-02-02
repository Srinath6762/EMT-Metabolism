library("survival")
library("ggplot2")
library("survminer")
library("ggforce")
library("ggfortify")
library("ggpubr")
x = read.table("../ACC_survData.tsv", header = T, sep="\t")
y = x[,c("sample", "OS", "OS.time", "EMT")]    ## change EMT to Glycolysis, OXPHOS and FAO as required
fit <- survfit(Surv(OS.time, OS) ~ EMT_, data = y)
png(file = "EMT_/ACC_EMT_KMplot.png",
	width = 800,
	height = 600,
	type = "cairo")
ggsurvplot(fit, data = y, pval = T, xlab ="Time in days", ggtheme = theme_classic())
dev.off()

fit.coxph <- coxph(Surv(OS.time, OS) ~ EMT_, data = y)
png(file = "EMT_/ACC_EMT_HRplot.png",
	width = 800,
	height = 600,
	type = "cairo")
ggforest(fit.coxph, data = y)
dev.off()
