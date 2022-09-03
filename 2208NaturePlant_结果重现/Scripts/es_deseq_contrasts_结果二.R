library(DESeq2)
library(biobroom)
library(tidyverse)

otu <- readRDS("./Data/otu_pers.RDS")
map <- readRDS("./Data/drought_map.RDS")

map <- filter(map, Compartment == "ES") 

map <- map %>% 
  mutate(Group = paste(Treatment2,Time, sep = "."))

otu <- otu[, colnames(otu) %in% map$SampleID]
otu <- otu[, match(map$SampleID, colnames(otu))]
otu <- otu[rowSums(otu) > 0, ]

dds<- DESeqDataSetFromMatrix(countData = otu,
                             colData = map,
                             design = ~ Group)

dds <- DESeq(dds)

# contrasts <- vector(mode = "list")

#--创建两两分组
tem = resultsNames(dds) %>% strsplit("_") %>% sapply(`[`, 2)
tem.1 = resultsNames(dds) %>%
  strsplit("p_") %>% sapply(`[`, 2) %>%
  strsplit("_vs_") %>% sapply(`[`, 1)

tem.2 = resultsNames(dds) %>%
  strsplit("_vs_") %>% sapply(`[`, 2) 

tem3 = paste("Group",tem.1[-1],tem.2[-1],sep = "=") %>%
  strsplit("=") 
names(tem3) = tem.1[-1]
contrasts = tem3

#--去掉作者的分析方式--不适用于所有数据
# for(i in 1:13) {
#   
#   contrasts[[paste("WC_TRN", i, sep = ".")]] <- c("Group", paste("WC_TRN", i, sep = "."), paste("WC", i, sep = "."))
#   contrasts[[paste("D1", i, sep = ".")]] <- c("Group", paste("D1", i, sep = "."), paste("WC", i, sep = "."))
#   contrasts[[paste("D2", i, sep = ".")]] <- c("Group", paste("D2", i, sep = "."), paste("WC", i, sep = "."))
#   contrasts[[paste("D3", i, sep = ".")]] <- c("Group", paste("D3", i, sep = "."), paste("WC", i, sep = "."))
#   
# }


results <- vector(mode = "list")
shrinkFC <- vector(mode = "list")


library("apeglm")
for(i in seq_along(contrasts)) {
  results[[names(contrasts)[i]]] <- tidy(results(dds, contrast = contrasts[[i]]) %>% as.data.frame()) %>% mutate(Day = i)
  shrinkFC[[names(contrasts)[i]]] <- tidy(lfcShrink(dds, coef = c(i + 1)) %>% as.data.frame()) %>% mutate(Day = i)
  print(contrasts[[i]])
}



results.df <- plyr::ldply(results, function(x) x)
names(results.df)[1] <- "Contrast"
names(results.df)[2] <- "OTU_ID"

shrinkFC.df <- plyr::ldply(shrinkFC, function(x) x)
names(shrinkFC.df)[1] <- "Contrast"
names(shrinkFC.df)[2] <- "OTU_ID"

saveRDS(results.df, "./Data/es_dab_bal.RDS")
saveRDS(shrinkFC.df, "./Data/es_shlfc.RDS")
saveRDS(dds, "./Data/es_dds.RDS")

