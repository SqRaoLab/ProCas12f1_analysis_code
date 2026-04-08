library(qs)
qsave(mouse,'crop0.qs',nthreads = 8) 
mouse <- qread('~/project/crop/sgrna_R/crop0.qs',nthreads = 8) 
DimPlot(mouse, reduction="umap", group.by="celltype",label = T)+scale_color_npg()#scale_color_d3(palette ="category20")

umap_coords <- Embeddings(mouse, reduction = "umap")

mouse@meta.data$UMAP1 <- umap_coords[, 1]
mouse@meta.data$UMAP2 <- umap_coords[, 2]
 
gg <- c("STC","GATA1","SPI1","CEBPE","ELK1","ZNF555","SMYD3","POU2AF1","IRF5","SCMH1")
gg <- c("STC")

#gg <- setdiff(unique(mouse$target_gene), c(NA, "NFE2"))
#gg <- c("STC", sort(setdiff(gg, "STC"))) 
#gg <- setdiff(unique(mouse$target_gene), c(NA, "NFE2","HLX","PCGF6",
                                           "STC","GATA2","SPI1","RUNX3","ELK1","ZNF555","SMYD3","POU2AF1","IRF5","SCMH1"))
gg <- c(sort(gg)) 
pDT.top <- mouse@meta.data[mouse@meta.data$target_gene %in% c(gg), ]

pDT.ntc<-mouse@meta.data
library(data.table)  
setDT(pDT.ntc)
pDT.ntc[, diff := "Background"]

pDT.final <- data.table()
setDT(pDT.top)
#for(x in c("Kmt2a","Kmt2d","Smarcd2","Smarcd1","Brd9")){
for(x in unique(pDT.top$target_gene)){
  pDTg <- rbind(pDT.ntc, pDT.top[target_gene == x], fill=TRUE)
  pDTg$plot <- x
  pDT.final <- rbind(pDT.final, pDTg, fill=TRUE)
}

sDT <- pDT.final[,c("UMAP1", "UMAP2", "plot", "diff","target_gene"),with=F]
sDT[,.N, by=c("plot", "diff")]
sDT[is.na(diff), diff := "target_gene"]
#names(sDT) <- c("UMAP1", "UMAP2", "Plot", "Diff")
#sDT[,.N, by=c("Plot", "Diff")]

# Plot
n.bins = 50
ggplot(pDT.final[is.na(diff)], aes(x=UMAP1, y=UMAP2)) + 
  #themeNF() +
  stat_binhex(
    data=pDT.final[diff == "Background"], 
    mapping = aes(fill="x"), 
    bins=n.bins, 
    fill="lightgrey") +
  geom_hex(data=pDT.final[is.na(diff)], bins=n.bins) +
  scale_fill_gradientn(colours=c("#ffff99", "#e31a1c")) +
  facet_wrap(~plot, ncol=3) 

pDT.final$plot <- factor(pDT.final$plot, levels = gg)
ggplot(pDT.final[is.na(diff)], aes(x=UMAP1, y=UMAP2)) + 
  theme_minimal()  +theme(panel.grid=element_blank())+
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(colour = "black",face = "bold", size = 12))+
  stat_binhex(
    data=pDT.final[diff == "Background"], 
    mapping = aes(fill="x"), 
    bins=n.bins, 
    fill="lightgrey") +
  stat_binhex(
    data=pDT.final[is.na(diff)], 
    mapping = aes(fill=(..count..)),#log10(..count..)
    bins=n.bins) +
  labs(fill="cell count")+#\n(log 10)
  scale_fill_gradientn(
    colours=c("#FCFDBF", "#FD9567","#F1605D","#CC0000"),#custom_mako,
    limits = c(0, 6),
    oob = squish
    ) +#viridis::scale_fill_viridis(option = "mako", direction = -1)+
  facet_wrap(~plot, ncol=6) + coord_fixed()
