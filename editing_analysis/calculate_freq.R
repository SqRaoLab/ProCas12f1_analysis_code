suppressPackageStartupMessages({
  library(data.table)
  library(pbapply)
  library(dplyr)
  library(glue)
  library(tidyr)
})


LABELs = c(
  "OsCas12f1"="enOsCas12f1", 
  "enAsCas12f1"="AsCas12f-HKRA", "enAsCas12f"="AsCas12f-HKRA",
  "RhCas12f1"="enRhCas12f1", "CasMINI"="CasMINI",
  "SpaCas12f1"="SpaCas12f1"
)


read_edit = function(path) {
  genome = as.data.frame(fread(path))

  colnames(genome)[colnames(genome) == "CasProtein"] = "cas"
  
  # rename and format cas protein name
  genome$cas = sapply(genome$cas, function(x) {LABELs[[x]]})
  
  # generate unique id
  # genome$id = pbapply(genome[, c("cas", "spacer", "before", "pam", "after")], 1, function(row){paste(row, collapse = "|")})
  # 直接合并这几列，生成新列，并移除原列
  genome$id <- tidyr::unite(genome, col = "temp", 
                            cas, spacer, before, pam, after, 
                            sep = "|", remove = FALSE)$temp

  genome$match = ifelse(genome$target == genome$edit, "match", "mismatch")
  
  # count number of different edits
  diff_kind = genome %>%
    filter(match == "mismatch") %>%
    dplyr::select(id, edit) %>%
    unique() %>%
    group_by(id) %>%
    add_tally() %>%
    dplyr::select(id, n) %>%
    unique() %>%
    as.data.frame()
  
  colnames(diff_kind)[colnames(diff_kind) == "n"] = "number_of_edit"
  
  # stats inel and match
  genome$value = 1
  genome = reshape2::dcast(genome, id+cas+spacer+before+pam+after~match, fun.aggregate = sum)
  genome$total = genome$match + genome$mismatch
  genome$freq = genome$mismatch / genome$total * 100
  genome = merge(genome, diff_kind, by = "id")
  rownames(genome) = genome$id
  genome
}



calculate_freq = function(root_path) {
  genome = glue("{root_path}/genome")
  plasmid = glue("{root_path}/plasmid")
  
  df = pblapply(Sys.glob(glue("{genome}/*.tsv.gz")), function(path) {
    cat(glue("reading {path}\n\n"))
    genome = read_edit(path)
    
    plasmid = glue("{plasmid}/{basename(path)}")
    cat(glue("reading {plasmid}\n\n"))
    plasmid = read_edit(plasmid)

    genome$bgMatch = plasmid[genome$id, "match"]
    genome$bgMisMatch = plasmid[genome$id, "mismatch"]
    genome$bgFreq = plasmid[genome$id, "freq"]
    genome$bgTotal = plasmid[genome$id, "total"]
    
    genome[is.na(genome)] = 0
    
    # calculate corrected_efficiency
    genome$estimate_bg_indel = genome$total * (genome$bgFreq / 100)
    genome$corrected_efficiency = (genome$mismatch - genome$estimate_bg_indel) / (genome$total - genome$estimate_bg_indel) * 100
    
    # filtering
    df = genome[genome$corrected_efficiency > 0 & !is.na(genome$corrected_efficiency), ]
    df = df[df$mismatch > 10 & df$bgFreq < 8 & df$number_of_edit > 1 & df$bgTotal > 0, ]
    df
  }, cl = 5)
  # cl here for multiprocessing, bigger cl bigger memory usage
  
  df = do.call(rbind, df)
  
  saveRDS(df, glue("{root_path}/freq.rds"))
  df
}


## Batch 1 PoolA
root_path = "res/batch1/decode/A"
df = calculate_freq(root_path)


## Batch 2 PoolA
root_path = "res/batch2/decode/A"
df = calculate_freq(root_path)


## Batch 1 PoolC
root_path = "res/batch1/decode/C"
df = calculate_freq(root_path)


## Batch 2 PoolC
root_path = "res/batch2/decode/C"
df = calculate_freq(root_path)
