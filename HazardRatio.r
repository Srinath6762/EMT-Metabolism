library("survival")
library("ggplot2")
library("survminer")
library("ggforce")
library("ggfortify")
library("ggpubr")

x = read.table("../EMT_OXPHOS_FAO/ACC_EMT_OXPHOS_FAO_survData.tsv", header = T, sep="\t")
y = x[,c("sample", "OS", "OS.time", "Type")]
y$Type = factor(y$Type, levels = c("E-O-F-", "E+O+F+", "E+O+F-", "E+O-F+", "E+O-F-", "E-O+F+", "E-O+F-", "E-O-F+")) ## to make E-G-O- as reference
fit.coxph <- coxph(Surv(OS.time, OS) ~ Type, data = y)
png(file = "EMT_OXPHOS_FAO/E-O-F-/ACC_EMT_OXPHOS_FAO_HRplot.png",
	width = 800,
	height = 600,
	type = "cairo")
ggforest(fit.coxph, data = y)
dev.off()

rex.cox = summary(coxph(Surv(OS.time, OS) ~ Type, data = y))
res = rex.cox$conf.int
write.table(t(res[,1]), file="/media/csb/NewVolume/susmita/EMT_Metabolism/survival/KMplots/EMT_OXPHOS_FAO/E-O-F-/ACC_E-O-F-_HR.tsv", col.names = T, row.names = F, sep="\t")




