library(data.table)
annList<-mouse@meta.data
annList <- annList[!is.na(annList$type), ]

setDT(annList)
clDT <- dcast.data.table(annList[,.N, by=c("type", "celltype")], 
                         celltype ~ type, value.var = "N")
clDT[, frac := GUIDES/STC]

(clDT.remove <- clDT[(STC < 5 | is.na(STC)) & (frac > 25 | is.na(frac))])
annList[celltype %in% clDT.remove$celltype, celltype := "remove"]
annList.allClusters <- copy(annList)
annList <- annList[celltype != "remove"]

annList$tissue <- "invivo"
fish.test.sets <- list()
TISSUES <- unique(annList$tissue)
celltype <- unique(annList$celltype)

for(tx in TISSUES){
  fish.test.sets[[paste("basic", tx, sep="_")]] <- annList[tissue == tx]
}

for(tx in TISSUES){
  x <- annList.allClusters[tissue == tx]
  x[, Clusters := paste("cl", celltype)]
  fish.test.sets[[paste("numeric", tx, sep="_")]] <- x
}

sapply(fish.test.sets, nrow)
fish.test.sets <- fish.test.sets[sapply(fish.test.sets, nrow) > 0]

(fish.test.x <- names(fish.test.sets)[1])
res <- data.table()
pDT1 <- fish.test.sets[[fish.test.x]]
#write.tsv(pDT1[,"rn",with=F], out("Guides_Fisher_Mixscape_",fish.test.x,"_Cells.tsv"))

for(gx in unique(pDT1[target_gene != "STC"]$target_gene)){
      cx <- pDT1$celltype[1]
      for(cx in unique(pDT1$celltype)){
        ntc <- unique(pDT1[type == "STC"][,.N, by="gRNA"][N > 3]$gRNA)[1] #N>20
        for(ntc in unique(pDT1[type == "STC"][,.N, by="gRNA"][N > 3]$gRNA)){
          mx <- as.matrix(with(pDT1[target_gene == gx | gRNA == ntc], table(celltype == cx, target_gene == gx)))
          if(all(dim(mx) == c(2,2))){
            fish <- fisher.test(mx)
            res <- rbind(res, data.table(
              celltype=cx, 
              target_gene=gx, 
              ntc=ntc, 
              p=fish$p.value, 
              OR=fish$estimate, 
              total.cells=sum(mx),
              guide.cells=nrow(pDT1[target_gene  == gx])
            ))
          }
        }
      }
    }

if(nrow(res) < 3) next
res[,padj := p.adjust(p, method="BH")]
res[, log2OR := pmax(-5, pmin(5, log2(OR)))]
#write.tsv(res, out("Guides_Fisher_Mixscape_",fish.test.x,".tsv"))
  
  
res <- hierarch.ordering(res, toOrder = "target_gene", orderBy = "celltype", value.var = "log2OR", aggregate = TRUE)
library(ggplot2)
ggplot(res[1:50, ], aes(
  x=celltype,
  y=ntc,
  color=log2OR,
  size=pmin(-log10(padj), 5))) +
  geom_point(shape=16) +
  scale_color_gradient2(name="log2OR", low="blue", high="red") +
  scale_size_continuous(name="padj") +
  facet_grid(target_gene ~ ., space = "free", scales = "free") +
  theme_bw(12) +
  theme(strip.text.y = element_text(angle=0)) #+
  xRot()
ggsave(out("Guides_Fisher_Mixscape_",fish.test.x,".pdf"), w=10, h=length(unique(res$grp)) * 0.50 + 1, limitsize = FALSE)

#dotplot------
fish <- res
fish <- fish[, .(
  log2OR=mean(log2OR), 
  dir=length(unique(sign(log2OR[padj < 0.05]))) <= 1, 
  #dir=length(unique(sign(log2OR)))==1, 
  padj=sum(padj < 0.05),
  N=.N,
  Ncells=sum(unique(guide.cells))), by=c( "celltype", "target_gene")]
fish[dir == FALSE, padj := 0]
fish[dir == FALSE, log2OR := NA]

fish[, gene := target_gene]
fish[padj == 0 | is.na(log2OR), log2OR := 0]
fish[, sig.perc := padj / N]
fish[,log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
#fish <- hierarch.ordering(dt = fish, toOrder = "gene", orderBy = "celltype", value.var = "log2OR")
#fish[, Clusters := gsub("^Gran$", "Gran.", Clusters)]
#fish[, Clusters := cleanCelltypes(Clusters)]
#fish <- hierarch.ordering(dt = fish, toOrder = "Clusters", orderBy = "gene", value.var = "log2OR")

genes_all_zero <- fish %>%
  group_by(gene) %>%
  summarise(all_zero = all(sig.perc == 0)) %>%
  filter(all_zero) %>%
  pull(gene)
fish_filtered <- fish[!fish$gene %in% genes_all_zero, ]
fish_filtered <- fish[!fish$gene %in% "NFE2", ]
library(scales)
ggplot(fish_filtered, aes(x=gene, y=celltype, size=sig.perc, color=log2OR_cap)) +
  theme_linedraw()+
  theme(axis.text.x=element_text(hjust = 0,vjust= 0,angle = 90),text = element_text(size = 16))+
  #coord_flip()+
  theme(panel.grid=element_blank())+
  theme(strip.text = element_text(colour = "black",face = "bold", size = 12))+
  scale_color_gradientn(name="log2(OR)",
    colours=c("#3300cc","#4157A4","white","#FF3333","#CC0000"),
    limits = c(-2, 2),
    oob = squish
  ) +
  #scale_color_gradient2(name="log2(OR)",low="#3300cc",, midpoint = 0, high="#E02925",limits = c(-2, 2),oob= squish) +#
  scale_size_continuous(name="% sign.", range = c(0,8)) +
  scale_x_discrete(position="top") +
  geom_point() +
  geom_point(shape=1, color="lightgrey") +
  xlab("") + ylab("")#+ scale_color_gradientn(colours = c('#2571B2','#A6CDE1','#FCCDB5','#BF3E3F'))
"#4157A4"
"#3300cc","#3399FF","white","#FF3333","#CC0000"

order <- c("Mk","Ery","MEP","Mast","EBM","Mono","Gran","LMPP")
fish$celltype <- factor(fish$celltype,levels = order)

library(forcats) 
priority_genes <- gg
fish_filtered$gene <- factor(fish_filtered$gene,
                             levels = c(priority_genes, 
                                        setdiff(levels(factor(fish_filtered$gene)), priority_genes)))

fish_filtered <- fish[!fish$gene %in% c("NFE2", "HLX","PCGF6")]

gg <-c("HOXA7","SPI1",
       "CEBPE", "LEF1", "TAL1","GATA2",
       "GATA1","KLF1", "RUNX3","RIOK2")
fish_filtered <- fish_filtered %>%
  mutate(group = ifelse(gene %in% gg, "positive control", "Other"))

fish_filtered$group <- factor(fish_filtered$group, levels = c("positive control","Other"))
ordered_genes <- c(gg,sort(setdiff(unique(fish_filtered$gene), gg)))
fish_filtered$gene <- factor(fish_filtered$gene, levels = ordered_genes)

ggplot(fish_filtered, aes(x = gene, y = celltype)) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(hjust = 0, vjust = 0.5, angle = 90),
    text = element_text(size = 15),
    panel.grid = element_blank(),
    strip.text = element_text(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines")
  ) +
  geom_point(
    aes(
      size = sig.perc,
      color = ifelse(sig.perc > 0, log2OR_cap, NA)  
    ),
    shape = 16
  ) +
  scale_color_gradientn(
    name = "log2(OR)",
    colours = c("#3300cc", "#4157A4", "white", "#FF3333", "#CC0000"),
    limits = c(-2, 2),
    oob = squish,
    na.value = "lightgrey"  
  ) +
  scale_size_continuous(name = "% sign.", range = c(0.6, 8)) +
  scale_x_discrete(position = "top") +
  xlab("") + ylab("")+
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  xlab("") + ylab("")

ggplot(fish_filtered, aes(x = gene, y = celltype)) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(hjust = 0, vjust = 0, angle = 90),
    text = element_text(size = 16),
    panel.grid = element_blank(),
    strip.text = element_text(colour = "black", face = "bold", size = 12)
  ) +
 
  geom_point(
    aes(
      size = sig.perc,
      color = ifelse(sig.perc > 0, log2OR_cap, NA)  
    ),
    shape = 16
  ) +
  scale_color_gradientn(
    name = "log2(OR)",
    colours = c("#3300cc", "#4157A4", "white", "#FF3333", "#CC0000"),
    limits = c(-2, 2),
    oob = squish,
    na.value = "lightgrey" 
  ) +
  scale_size_continuous(name = "% sign.", range = c(0.6, 8)) +
  scale_x_discrete(position = "top") +
  xlab("") + ylab("")


library(data.table)
annList<-mouse@meta.data
annList <- annList[!is.na(annList$celltype), ]

setDT(annList)
annList <- annList[
  , Clusters := fifelse(
    grepl("Ery|MEP",  celltype, ignore.case = TRUE),  
    "MEP",
    fifelse(grepl("GMP|Gran|Mono|EBM|Mast", celltype, ignore.case = TRUE), "GMP", NA_character_)
  )
]
annList <- annList[!is.na(Clusters)]

fish.test.sets         <- list()
fish.test.sets[["eryVsMye"]] <- annList 

res <- data.table()
pDT1 <- fish.test.sets[["eryVsMye"]]
#write.tsv(pDT1[,"rn",with=F], out("Guides_Fisher_Mixscape_",fish.test.x,"_Cells.tsv"))

for(gx in unique(pDT1[target_gene != "STC"]$target_gene)){
  #cx <- pDT1$seurat_clusters[1]
  for(cx in unique(pDT1$Clusters)){
    ntc <- unique(pDT1[type == "STC"][,.N, by="gRNA"][N > 3]$gRNA)[1] #N>20
    for(ntc in unique(pDT1[type == "STC"][,.N, by="gRNA"][N > 3]$gRNA)){
      mx <- as.matrix(with(pDT1[target_gene == gx | gRNA == ntc], table(Clusters == cx, target_gene == gx)))
      if(all(dim(mx) == c(2,2))){
        fish <- fisher.test(mx)
        res <- rbind(res, data.table(
          Clusters=cx, 
          target_gene=gx, 
          ntc=ntc, 
          p=fish$p.value, 
          OR=fish$estimate, 
          total.cells=sum(mx),
          guide.cells=nrow(pDT1[target_gene  == gx])
        ))
      }
    }
  }
}

if(nrow(res) < 3) next
res[,padj := p.adjust(p, method="BH")]
res[, log2OR := pmax(-5, pmin(5, log2(OR)))]
#write.tsv(res, out("Guides_Fisher_Mixscape_",fish.test.x,".tsv"))

fish1 <- res
fish1 <- fish1[, .(
  log2OR=mean(log2OR), 
  dir=length(unique(sign(log2OR[padj < 0.05]))) <= 1, 
  #dir=length(unique(sign(log2OR)))==1, 
  padj=sum(padj < 0.05),
  N=.N,
  Ncells=sum(unique(guide.cells))), by=c( "Clusters", "target_gene")]
fish1[dir == FALSE, padj := 0]
fish1[dir == FALSE, log2OR := NA]

fish1[, gene := target_gene]
fish1[padj == 0 | is.na(log2OR), log2OR := 0]
fish1[, sig.perc := padj / N]
fish1[,log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]

fish1 <- fish1[Clusters == "GMP"]
genes_all_zero <- fish %>%
  group_by(gene) %>%
  summarise(all_zero = all(sig.perc == 0)) %>%
  filter(all_zero) %>%
  pull(gene)
fish_filtered1 <- fish1[!fish1$gene %in% c("NFE2","HLX","PCGF6") ]
ggplot(fish_filtered1, aes(x=log2OR, y=gene, fill=sign(log2OR_cap))) + #alpha=sig.perc,
  theme_linedraw()+
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(color = "black",angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    panel.grid = element_blank()
  )+
  theme(axis.text.x=element_text(hjust = 0.5,vjust=0.5,angle = 90))+
  geom_col() +
  scale_fill_gradient2(name="log2(OR)",low="#4157A4", midpoint = 0, high="#E02925") +
  #scale_size_continuous(name="% sign.", range = c(0,5), limits = c(0,1)) +
  #geom_point() +
  scale_y_discrete(position="left") +
  ylab("") + xlab("Mye vs \nery log2(OR)") +
  geom_vline(xintercept = 0)+coord_flip()+
  facet_grid(. ~ group, scales = "free_x", space = "free_x") + #分面
  theme_linedraw() +
  theme(
    axis.text.x = element_text(hjust = 0, vjust = 0.5, angle = 90),
    text = element_text(size = 12),
    panel.grid = element_blank(),
    strip.text = element_text(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines")
  ) 

gg <-c("HOXA7","SPI1",
       "CEBPE", "LEF1", "TAL1","GATA2",
       "GATA1","KLF1", "RUNX3","RIOK2")
fish_filtered1 <- fish_filtered1 %>%
  mutate(group = ifelse(gene %in% gg, "positive control", "Other"))

fish_filtered1$group <- factor(fish_filtered1$group, levels = c("positive control","Other"))
ordered_genes <- c(gg,sort(setdiff(unique(fish_filtered1$gene), gg)))
fish_filtered1$gene <- factor(fish_filtered1$gene, levels = ordered_genes)

library(forcats) 
priority_genes <- gg
fish_filtered1$gene <- factor(fish_filtered1$gene,
                             levels = c(priority_genes, 
                                        setdiff(levels(factor(fish_filtered$gene)), priority_genes)))
