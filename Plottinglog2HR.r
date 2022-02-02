data = read.table("log2HR_Data.tsv", sep="\t", header = T)
ggplot(data=data, aes(CancerType, Hazard_Ratio)) +
geom_pointrange(aes(ymin=lower95, ymax=upper95, colour = ssGSEA), position = position_dodge(0.5), shape=20, size= 0.5) +
theme_classic() +
geom_hline(yintercept=1, lty=2) + coord_flip()
