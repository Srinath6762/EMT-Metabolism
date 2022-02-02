library("RcolorBrewer")
x = read.table("E-O-F-_HR.tsv", header = T, sep="\t", row.names = 1, fill = T)
y = as.matrix(-log2(x))
col= colorRampPalette(brewer.pal(8, "RdBu"))(25)
png(file = "E-O-F-_heatmap.png",
    width = 800,
    height = 600)
heatmap(y, col= colorRampPalette(brewer.pal(8, "RdBu"))(25))
dev.off()
