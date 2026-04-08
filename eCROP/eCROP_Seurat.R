library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)

sample_name <- "crop"
input_dir <- "~/project/crop/0610/celescope/crop/outs/filtered/"
assignment_df <- read_tsv(sprintf("~/project/crop/0610/sgrna_R/sgrna_FinalAssignments.tsv"))
data <- Read10X(input_dir)
seurat_obj<- CreateSeuratObject(data, project = "CROP-seq", assay = "RNA")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")

metadata <- seurat_obj@meta.data
metadata$cell_barcode <- rownames(metadata)
assignment_df <- as.data.frame(assignment_df)
rownames(assignment_df) <- assignment_df$cell_barcode
colnames(assignment_df) <- c("cell_barcode", "type", "gRNA", "target_gene")

metadata <- metadata%>%select(-type) %>% select(-gRNA)%>% select(-target_gene)
metadata <- left_join(metadata, assignment_df, by = "cell_barcode")
#assignment_df <- assignment_df[, c("type", "gRNA",)]

# Add back to seurat obj
seurat_obj[["type"]] <- metadata$type
seurat_obj[["gRNA"]] <- metadata$gRNA
seurat_obj[["target_gene"]] <- metadata$target_gene
mouse <-seurat_obj 

# QC plots
umi_hist <- ggplot(metadata, aes(nCount_RNA)) + geom_histogram(binwidth = 1000) + 
  geom_vline(xintercept = 9000, color = "red") + geom_vline(xintercept = 25000, color = "red") + theme_classic() + 
  ggtitle(sprintf("Total UMIs for %s", sample_name)) + xlab("Total UMIs") + ylab("Number of cells")
feature_hist <- ggplot(metadata, aes(nFeature_RNA)) + geom_histogram(binwidth = 100) + 
  geom_vline(xintercept = 500, color = "red") + theme_classic() + 
  ggtitle(sprintf("Total Features for %s", sample_name)) + xlab("Total Features") + ylab("Number of cells")
mt_hist <- ggplot(metadata, aes(percent.mt)) + geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = 10, color = "red") + theme_classic() + 
  ggtitle(sprintf("Mitochondrial Expression for %s", sample_name)) + xlab("% Mt Expression") + ylab("Number of cells")
rb_hist <- ggplot(metadata, aes(percent.rb)) + geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = 40, color = "red") + theme_classic() + 
  ggtitle(sprintf("Ribosomal Expression for %s", sample_name)) + xlab("% Rb Expression") + ylab("Number of cells")

pdf(sprintf("%s_QC.pdf", sample_name), width = 11, height = 8.5)
grid.arrange(umi_hist, feature_hist, mt_hist, rb_hist, ncol = 2)
dev.off()

# Filter data
seurat_obj <- subset(seurat_obj, nCount_RNA > 0 & nCount_RNA < 25000 & nFeature_RNA < 6000 & 
                       percent.mt < 10 & percent.rb < 30)
mouse <- seurat_obj
saveRDS(mouse, sprintf("~/project/crop/sgrna_R/crop_SeuratObject.rds"))

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)

sample_name <- "crop"
input_path <- sprintf("~/project/crop/sgrna_R/crop_SeuratObject.rds")
seurat_obj <- readRDS(input_path)
metadata <- seurat_obj@meta.data
metadata$cell_barcode <- rownames(metadata)
metadata <- metadata %>% mutate(sample = sample_name) %>% select(cell_barcode, sample, type, gRNA, target_gene) 
write_tsv(metadata, sprintf("~/project/crop/sgrna_R/sgrna_PostFilterAssignments.tsv"))


#Seurat------
library(Seurat)
VlnPlot(mouse, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

mouse <- PercentageFeatureSet(mouse, pattern = "^MT-", col.name = "percent.mt")
#mouse <- subset(mouse, subset = nFeature_RNA > 200 & percent.mt < 20)

#options(future.globals.maxSize = 3000 * 1024^2)
#BiocManager::install("glmGamPoi")
#library(glmGamPoi)
#library(future)
#options(future.globals.maxSize = 6 * 1024^3)
mouse <- SCTransform(mouse, 
                     vars.to.regress = NULL) 

mouse <- NormalizeData(mouse)
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
mouse <- ScaleData(mouse, features = rownames(mouse))
mouse <- RunPCA(mouse)
mouse <- RunUMAP(mouse, dims = 1:50)
mouse <- FindNeighbors(mouse, dims = 1:30)
mouse <- FindClusters(mouse, resolution = 0.6)

DimPlot(mouse, reduction="umap",group.by = "seurat_clusters",#ncol=3,split.by = "type","gRNA","seurat_clusters"
        label = T)+ scale_color_d3(palette ="category20")


mouse<- readRDS("~/project/crop/sgrna_R/crop_SeuratObject.rds")
mouse$type[is.na(mouse$type)] <- "Na"
DimPlot(mouse, reduction="umap",group.by = "seurat_clusters",#ncol=3,split.by = "type",#seurat_clusters
        label = T)+ scale_color_d3(palette ="category20")
DimPlot(mouse, reduction="umap",group.by =  "type",
        label = F,pt.size=0.5)+  scale_color_d3(palette ="category20")
split.by = "orig.ident",ncol=4,
ElbowPlot(mouse)
DimPlot(mouse, reduction = "pca")

#cellratio
mouse$type <- factor(mouse$type, level=c("GUIDES","STC","Na"))
cellratio <- prop.table(table(mouse$type,mouse$celltype), margin=2) 
cellratio
cellratio <- as.data.frame(cellratio)
colorcount = length(unique(cellratio$Var1))
ggplot(cellratio)+
  geom_bar(aes(x=Freq*100, y=Var2, fill=Var1),stat = "identity",width = 0.8, size = 0.2, color="black", linetype="blank")+
  theme_linedraw(base_size = 16, base_family = "Arial", base_line_size = 0.5)+
  theme(panel.grid=element_blank())+
  labs(x="Proportion(%)", y="", fill="")+
  scale_fill_d3(palette ="category20")

#barplot
mouse$group <- "invivo"
cellratio <- prop.table(table(mouse$type,mouse$group), margin=2) 
cellratio
cellratio <- as.data.frame(cellratio)
ggplot(cellratio, aes(x = "", y = Freq*100, fill = Var1)) +
  geom_bar(width = 0.8, stat = "identity", color = "white") +
  coord_polar(theta = "y") +theme_void() +
  theme(plot.title = element_text(hjust = 0.5))+ scale_fill_brewer(palette = "Paired")+
  labs(fill="Type")+
  geom_text(aes(label = paste0(round(Freq * 100, 1), "%")), 
            position = position_stack(vjust = 0.5))

meta <- mouse@meta.data
meta$target_gene <- as.character(meta$target_gene)

gene_counts <- meta %>%
  filter(!is.na(target_gene), !target_gene %in% c("STC","NFE2","HLX","PCGF6")) %>%
  group_by(target_gene) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))

ggplot(gene_counts, aes(x = target_gene, y = cell_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black") +  
  labs(x = "", y = "Number of cells", title = "") +
  theme_linedraw(base_size = 16, base_line_size = 0.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()   
  )
ggsave("~/project/crop/0610/merge_R/supp_barplot_v2.pdf", width =12.6, height = 3.28)

#STC corrplot-----
library(data.table)
library(Seurat)
library(ggplot2)
meta <- mouse@meta.data
dt <- as.data.table(meta)
dt_stc <- dt[target_gene == "STC"]
pDT <- dt_stc[, .N, by = .(celltype, gRNA)]
pDT[, total := sum(N), by = .(gRNA)]
pDT <- pDT[total > 100]
pDT[, frac := N / total * 100]
pDT[, guide_i := as.numeric(factor(gRNA))]
pDT[, Clusters := as.factor(celltype)]

ggplot(pDT, aes(x = Clusters, y = frac, fill = factor(guide_i))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Percentage (%)",fill ="STC") +
  scale_fill_brewer(palette = "Paired")+
  theme_linedraw(base_size = 16, base_line_size = 0.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )

pDT_wide <- dcast(pDT, Clusters ~ paste0("STC", guide_i), value.var = "frac")

guide_cols <- grep("^STC", names(pDT_wide), value = TRUE)
guide_nums <- as.numeric(gsub("STC", "", guide_cols))
guide_cols_sorted <- guide_cols[order(guide_nums)]
sub <- pDT_wide[, ..guide_cols_sorted]

library(corrplot)
library(GGally)
ggpairs(sub, title = "")

my_palette <- colorRampPalette(c("#FFFFCC", "#FF6600", "#990000"))

my_upper <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr_val <- cor(x, y, use = "complete.obs")
  corr_rounded <- round(corr_val, 2)
  
  df <- data.frame(x = 0.5, y = 0.5, corr = corr_val)
  
  ggplot(df, aes(x = x, y = y, fill = corr)) +
    geom_tile(color = "black") +
    geom_text(aes(label = corr_rounded), color = "black", size = 4.5) +
    scale_fill_gradientn(
      colors = my_palette(100),
      limits = c(0.95, 1),
      oob = scales::squish 
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "none"
    )
}

my_lower <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(size = 0.7) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.text = element_text(),
      axis.ticks = element_line(),
      axis.title = element_blank()
    )
}

my_diag <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(color = "black", fill = NA) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.text = element_text(),
      axis.ticks = element_line(),
      axis.title = element_blank()
    )
}

ggpairs(
  sub,
  upper = list(continuous = my_upper),
  lower = list(continuous = my_lower),
  diag  = list(continuous = my_diag),
  title = ""
)



library(pheatmap)
cor_mat <- cor(sub, use = "pairwise.complete.obs")
my_colors <- colorRampPalette(c("#ffff99", "#FF3333"))(20)
pheatmap(
  cor_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = my_colors,
  number_color = "black", 
  main = ""
)
pheatmap(cor_mat, cluster_rows=F,cluster_cols=F,
         display_numbers = TRUE, 
         main = "")



library(SingleR)
library(celldex)
library(ggsci)
db <-HumanPrimaryCellAtlasData() 
db1 <- MonacoImmuneData()
human_data <- GetAssayData(mouse, layer ="data")
human_annot <- SingleR(test=human_data, ref=db, labels=db$label.fine, clusters=mouse$seurat_clusters)
mouse$celltype_SingleR <- factor(mouse$seurat_clusters, levels=0:14, labels=human_annot$labels)
DimPlot(mouse, reduction="umap", group.by="celltype_SingleR",
        label = T)+ scale_color_d3(palette ="category20")

VlnPlot(mouse1,features=c("nCount_RNA", "nFeature_RNA", "percent.mt"), split.by = "seurat_clusters")
mouse1 <- subset(mouse,seurat_clusters%in% c(0,1,3,4,5,6,7,8))
mouse <- mouse1

marker <- c("AVP","CCR7","CD79B","GATA1","GZMH","SPI1","DNTT","VWF")
DefaultAssay(mouse)<-"RNA"
FeaturePlot(mouse,features = marker,reduction="umap", ncol=4)

markergenes <- c(#"LGALS1","SPI1", #MD
                 "CSF1R","FCN1","VCAN",#Mono
                 "IRF8","CCR2",#DC-prog "IGKC","SPIB","MPEG1",
                 "HLA-DPA1","HLA-DRA",#DC
                 "HLF","AVP",# HSC,"CRHBP"
                 "CRHBP","HOPX","SPINK2", #MPP
                 "ITGA2B","PF4", #MKP "VWF",
                 "ITGA2B","PLEK","GP9", #Mk
                 "SLC40A1","GATA1","CNRIP1","CTNNBL1",#MEP
                 "SMIM24","SELL","SPINK2",#LMPP
                 "CLC","HDC", #EBM,"RNASE2", "GATA2","MS4A3"
                 "SPINK2","CSF3R",#CMP
                 "ELANE","MPO","CTSG",#"GMP
                 "HBA1","HBB","HBM", #Ery
                 "JCHAIN","HOPX","SPINK2","TCF4","EBF1", #MLP
                 "CD3E","CD3D",#T
                 "FCG3RA", "NCAM1", "NKG7","GNLY", #NK
                 "CD79A","CD19","CD79B", #B
                 "MPO", #Neu1,Mono-prog
                 "NCF1","ITGAM","LYZ" #Neu2
)
markergenes <- c("SPINK2","HOPX","SMIM24","SPNS3",#LMPP
                 "CLC","PRG2","RNASE2","HDC","GATA2", #EBM,"RNASE2", "GATA2","MS4A3"
                 #"HPGDS","KIT","TPSAB1", #MAST
                 "MPO","CTSG","ELANE","AZU1",#Gran
                 "CSF1R","FCN1","VCAN",#Mono
                 "ITGA2B","PLEK","GP9","PF4",#Mk
                 "GATA1","KLF1","CTNNBL1",#MEP
                 "APOC1","AHSP","HBB"#Ery
)
unique_markergenes <- unique(markergenes)
## Dot Plot
DotPlot(mouse, features = unique_markergenes,group.by = "seurat_clusters")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
DotPlot(mouse,features=unique_markergenes,group.by = "celltype")+ theme_linedraw()+
  theme(axis.text.x=element_text(hjust = 1,vjust=0.5,angle = 90),text = element_text(size = 18))+
  scale_size_continuous(range=c(2,8))+
  scale_color_gradientn(colours = c('#2571B2','#A6CDE1','#FCCDB5','#BF3E3F'))+
  coord_flip()+
  guides(size  = guide_legend(title.position = "top",direction= "horizontal"),
         color = guide_colourbar(title.position = "top",direction= "horizontal",
                                 barwidth  = unit(6,  "cm"),barheight = unit(0.4,"cm")))+
  theme(legend.position  = "top",legend.direction = "horizontal", legend.box = "horizontal")
ggsave("~/project/crop/0610/merge_R/marker_v2.pdf", width =6.10, height = 9.17)

order <- c("LMPP","EBM","Gran","Mono","Mk","MEP","Ery")
mouse$celltype <- factor(mouse$celltype,levels = order)

## Feature Plot
FeaturePlot(mouse,features = unique_markergenes,reduction="umap", ncol=5)

marker <- c("AVP","CCR7","CD79B","GATA1","GZMH","SPI1","DNTT","VWF")
marker <- c("CD14","CD33","PTPRC","ITGA2B","CD34")
DefaultAssay(mouse)<-"RNA"
FeaturePlot(mouse,features = marker,reduction="umap", ncol=3)

cluster <- mouse$seurat_clusters
celltype <- rep("cell", length(cluster))
celltype[cluster %in% c(7)] <- "LMPP"
celltype[cluster %in% c(9,12,13)] <- "Mono"
celltype[cluster %in% c(8)] <- "Gran"
celltype[cluster %in% c(5,11)] <- "Mk"
celltype[cluster %in% c(3)] <- "MEP"
celltype[cluster %in% c(2)] <- "Ery"
celltype[cluster %in% c(4,6)] <- "EBM"
celltype[cluster %in% c(0,1,10)] <- "Mast"
mouse$celltype <- celltype

mouse <- subset(mouse, seurat_clusters %in% setdiff(unique(seurat_clusters), c(9, 11, 16)))
mouse <- subset(mouse, seurat_clusters %in% setdiff(unique(seurat_clusters), c(13)))
cluster <- mouse$seurat_clusters
celltype <- rep("cell", length(cluster))
celltype[cluster %in% c(8)] <- "LMPP"
celltype[cluster %in% c(0,9)] <- "Mono"
celltype[cluster %in% c(6)] <- "Gran"
celltype[cluster %in% c(7)] <- "Mk"
celltype[cluster %in% c(1,2,3,11)] <- "MEP"
celltype[cluster %in% c(12)] <- "Ery"
celltype[cluster %in% c(4,5,10)] <- "EBM"
mouse$celltype <- celltype

DimPlot(mouse, reduction="umap", group.by="celltype",
        pt.size = 0.1,label.size = 8,alpha = 1,
        label = T)+scale_color_npg()+#+scale_color_d3(palette ="category20") 
        guides(color = guide_legend(override.aes = list(size = 3, shape = 16, stroke = 0))) 
#split.by = "sample",label = F, ncol=2
order <- c("LMPP","EBM","Gran","Mono","Mk","MEP","Ery")
order <- c("EBM","Ery","LMPP","Gran","MEP","Mono","Mk")
mouse$celltype <- factor(mouse$celltype,levels = order)