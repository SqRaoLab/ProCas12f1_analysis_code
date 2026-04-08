library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(readr)


data <- Read10X("~/project/crop/celescope/count/GFP/outs/filtered")
data <- Read10X("~/project/crop/0610/celescope/crop/outs/filtered")
seurat_obj <- CreateSeuratObject(data, project = "CROP-seq", assay = "RNA")

# Grab list of gRNAs
grna_list <- grep("-gene", rownames(seurat_obj), value = TRUE)
counts <- seurat_obj[["RNA"]]$counts
grna_matrix <- counts[grna_list, ]

# For each cell, identify most abundant guide, discard if there are other molecules
transcriptome_assignments <- apply(grna_matrix, 2, function(x){
  if(any(x > 0)){
    max_val <- max(x)
    max_guide <- names(x)[which(x == max_val)]
    return(max_guide)
  } else{
    return(NA)
  }
})

transcriptome_df <- as.data.frame(unlist(transcriptome_assignments))
transcriptome_df$cell_barcode <- rownames(transcriptome_df)
rownames(transcriptome_df) <- NULL
colnames(transcriptome_df) <- c("gRNA", "cell_barcode")
transcriptome_df <- transcriptome_df[, c("cell_barcode", "gRNA")]
#write_tsv(transcriptome_df, sprintf("~/project/crop/sgrna_R/sgrna_TranscriptAssignments.tsv"))

#sgrna------
library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(Seurat)


grna_count_df <- read_tsv("~/project/crop/star/sgrna_count.filter.txt")
grna_count_df <- read_tsv("~/project/crop/sgrna_R/filtered_crop1.txt")

grna_count_df <- read_tsv("~/project/crop/sgrna_R/crop_filtered_output.txt", col_names = FALSE)
grna_count_df <- read_tsv("~/project/crop/sgrna_R/crop-1_filtered_output.txt", col_names = FALSE)
grna_count_df <- read_tsv("~/project/crop/sgrna_R/crop-2_filtered_output.txt", col_names = FALSE)
grna_count_df <- read_tsv("~/project/crop/sgrna_R/crop-mix_filtered_output.txt", col_names = FALSE)

grna_count_df <- read_tsv("~/project/crop/0610/sgrna_count/cropseq.output.txt")
grna_count_df <- read_tsv("~/project/crop/0610/sgrna_R/filtered_output.txt")

colnames(grna_count_df) <- c("cell", "barcode", "sgrna", "gene", "read_count", "umi_count")

summarized_reads <- grna_count_df %>% group_by(cell, barcode) %>% summarize(total_reads = sum(umi_count))
summarized_reads <- summarized_reads %>% separate(barcode, "_", into = c("type", "gRNA", "target"))
summarized_reads$gRNA[which(summarized_reads$type == "STC")] <- paste0("STC", summarized_reads$target[which(summarized_reads$type == "STC")])
summarized_reads$target[which(summarized_reads$type == "STC")] <- "STC"
summarized_reads <- ungroup(summarized_reads)

selectGuide1 <- function(x, y){
  if (nrow(x) > 1){
    # 1. If all of them are 1, we will not assign it something
    if (all(x["total_reads"] == 1)){
      output <- data.frame(cell = y, gRNA = NA)
    } 
    # 2. Otherwise, if all guides are equal, we will not assign it something
    else if (n_distinct(x["total_reads"]) == 1){
      output <- data.frame(cell = y, gRNA = NA)
    }
    # 3. If one is larger than the rest, and it is more than (sum(others) * 3) we will keep it
    else{
      # Find largest value and retrieve rows with largest value
      max_val <- max(x["total_reads"])
      max_row <- x[which(x["total_reads"] == max_val), ]
      
      # If there are multiple rows that are tied, check if they target the 
      # same gene
      # If they do, run standard check
      if (nrow(max_row) > 1){
        if (n_distinct(max_row["target"]) == 1){
          max_row <- max_row[1, ]
          other_rows <- x[which(!(x["gRNA"] %in% max_row["gRNA"])), ]
          if (max_row["total_reads"] > (sum(other_rows["total_reads"]))*3){
            output <- data.frame(cell = y, gRNA = max_row["gRNA"])
          } else{
            output <- data.frame(cell = y, gRNA = NA)
          }
        } else{
          output <- data.frame(cell = y, gRNA = NA)
        }        
      } else{
        # Find largest value and retrieve rows with largest value
        other_rows <- x[which(!(x["gRNA"] %in% max_row["gRNA"])), ]
        if (max_row["total_reads"] > (sum(other_rows["total_reads"]))*3){
          output <- data.frame(cell = y, gRNA = max_row["gRNA"])
        } else{
          output <- data.frame(cell = y, gRNA = NA)
        }
      }
    }} 
  else{
    output <- data.frame(cell = y, gRNA = x["gRNA"])
  }
  return(output)
}

selectGuide <- function(x, y){
  if (nrow(x) > 1){  
    if (all(x$total_reads == 1)){
      output <- data.frame(cell = y, gRNA = NA, type = "Ambiguous")
    }
    else if (n_distinct(x$total_reads) == 1){
      output <- data.frame(cell = y, gRNA = NA, type = "Ambiguous")
    }
    else {
      max_val <- max(x$total_reads)
      max_row <- x[which(x$total_reads == max_val), ]

      if (nrow(max_row) > 1){
        if (n_distinct(max_row$target) == 1){
          max_row <- max_row[1, ]
          other_rows <- x[!(x$gRNA %in% max_row$gRNA), ]
          if (max_row$total_reads > (sum(other_rows$total_reads)) * 3){
            output <- data.frame(cell = y, gRNA = max_row$gRNA, type = "Assigned")
          } else {
            output <- data.frame(cell = y, gRNA = NA, type = "Ambiguous")
          }
        } else {
          output <- data.frame(cell = y, gRNA = NA, type = "NP")
        }
      } else {
        other_rows <- x[!(x$gRNA %in% max_row$gRNA), ]
        if (max_row$total_reads > (sum(other_rows$total_reads)) * 3){
          output <- data.frame(cell = y, gRNA = max_row$gRNA, type = "Assigned")
        } else {
          output <- data.frame(cell = y, gRNA = NA, type = "Ambiguous")
        }
      }
    }
  } else {
    output <- data.frame(cell = y, gRNA = x$gRNA, type = "Single")
  }
  return(output)
}
# Remove null values
summarized_reads <- summarized_reads %>% group_by(cell)
assignment_list <- group_map(summarized_reads, selectGuide1)
assignments_df <- bind_rows(assignment_list)
assignments_df <- left_join(assignments_df, summarized_reads[, c("cell", "type", "gRNA", "target")], by = c("cell", "gRNA"))
colnames(assignments_df)[colnames(assignments_df) == "cell"] <- "cell_barcode"
write_tsv(assignments_df, sprintf("sgrna_ReadAssignments.tsv"))

#sgrna assign
library(readr)
library(tidyr)
library(dplyr)
library(purrr)

transcriptome_df <- read_tsv(sprintf("~/project/crop/sgrna_R/sgrna_TranscriptAssignments.tsv"))
read_df <- read_tsv(sprintf("~/project/crop/0610/sgrna_count/sgrna_ReadAssignments.tsv"))
read_df <- merged_df
transcriptome_df <- separate(transcriptome_df, gRNA, "-", into = c("type", "gRNA", "gene", "feature"))
transcriptome_df[which(transcriptome_df$gRNA == "Human"), "gRNA"] <- paste0("STC",transcriptome_df$gene[which(transcriptome_df$gRNA == "Human")])
transcriptome_df[which(transcriptome_df$type == "STC"), "gene"] <- "STC"
transcriptome_df <- transcriptome_df[, c("cell_barcode", "type", "gRNA", "gene")]

# Combine assignments
combined_df <- left_join(transcriptome_df, read_df, by = "cell_barcode", suffix = c(".UMI", ".read"))

# Filter conflicts and all NA values
filtered_combined_df <- combined_df %>% filter_at(vars(gRNA.UMI, gRNA.read), any_vars(!is.na(.)))
filtered_combined_df <- filtered_combined_df %>% filter_at(vars(gene, target), any_vars(is.na(gene) | is.na(target) | (gene == target)))

mergeResult <- function(x){
  if (any(is.na(x))){
    if (all(is.na(x[2:4]))){
      x <- x %>% select(cell_barcode = cell_barcode, type = type.read, gRNA = gRNA.read, target = target)
    } else{
      x <- x %>% select(cell_barcode = cell_barcode, type = type.UMI, gRNA = gRNA.UMI, target = gene)
    }
  } else{
    x <- x %>% select(cell_barcode = cell_barcode, type = type.UMI, gRNA = gRNA.UMI, target = gene)
  }
  return(x)
}

filtered_rows <- lapply(1:nrow(filtered_combined_df), function(x) mergeResult(filtered_combined_df[x, ]))
final_assignments <- bind_rows(filtered_rows)
write_tsv(final_assignments, sprintf("~/project/crop/sgrna_R/sgrna_FinalAssignments.tsv"))
write_tsv(final_assignments, sprintf("~/project/crop/0610/sgrna_R/sgrna_FinalAssignments.tsv"))
