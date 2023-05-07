

## 



# #--安装R包
# BiocManager::install("lefser")
# devtools::install_github('cafferychen777/ggpicrust2')


#--实例

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
map <-
  read_csv(
    "./map_16s.csv"

  )
map$...1 = NULL

group <- "Group"
library(R.utils)
gunzip("./picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", remove = TRUE)

daa_results_list <-
  ggpicrust2(
    file = "./picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv",
    metadata = map,
    group = "Group",
    pathway = "KO",
    daa_method = "LinDA",
    p_values_bar = TRUE,
    p.adjust = "BH",
    ko_to_kegg = TRUE,
    order = "pathway_class",
    select = NULL,
    reference = NULL 
  )
#--这个函数出现错误，我们开始用作者的逐条代码
#If you want to analysis kegg pathway abundance instead of ko within the pathway. You should turn ko_to_kegg to TRUE.
#The kegg pathway typically have the more explainable description.
#metadata should be tibble.
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
metadata <-
 map

kegg_abundance <-
  ko2kegg_abundance(
    "./picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv"
  )

group <- "Group"
#--进行差异分析
daa_results_df <-
  pathway_daa(
    abundance = kegg_abundance,
    metadata = map,
    group = group,
    daa_method = "ALDEx2",
    select = NULL,
    reference = NULL
  )

#if you are using LinDA, limme voom and Maaslin2, please specify a reference just like followings code.
# daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = group, daa_method = "LinDA", select = NULL, reference = "Harvard BRI")

daa_sub_method_results_df <-
  daa_results_df[daa_results_df$method == "ALDEx2_Kruskal-Wallace test", ]

daa_annotated_sub_method_results_df <-
  pathway_annotation(pathway = "KO",
                     daa_results_df = daa_sub_method_results_df,
                     ko_to_kegg = TRUE)

Group <-
  metadata$Group # column which you are interested in metadata

# select parameter format in pathway_error() is c("ko00562", "ko00440", "ko04111", "ko05412", "ko00310", "ko04146", "ko00600", "ko04142", "ko00604", "ko04260", "ko04110", "ko04976", "ko05222", "ko05416", "ko00380", "ko05322", "ko00625", "ko00624", "ko00626", "ko00621")


#  展示组合柱状图和误差图#--------
daa_results_list <-
  pathway_errorbar(
    abundance = kegg_abundance,
    daa_results_df = daa_annotated_sub_method_results_df,
    Group = Group,
    p_values_threshold = 0.04,
    order = "pathway_class",
    select = NULL,
    ko_to_kegg = TRUE,
    p_value_bar = TRUE,
    colors = NULL,
    x_lab = "pathway_name"
  )


# 展示热图#---------
# Create example functional pathway abundance data
abundance_example <- matrix(rnorm(30), nrow = 10, ncol = 3)
rownames(abundance_example) <- paste0("Sample", 1:10)
colnames(abundance_example) <- c("PathwayA", "PathwayB", "PathwayC")

# Create example metadata
# Please change your sample id's column name to sample_name
metadata_example <- data.frame(sample_name = rownames(abundance_example),
                               group = factor(rep(c("Control", "Treatment"), each = 5)))



heatmap_plot <- ggpicrust2::pathway_heatmap(t(abundance_example), metadata_example, "group")

#--展示排序图#---------
# Create example functional pathway abundance data
abundance_example <- data.frame(A = rnorm(10), B = rnorm(10), C = rnorm(10))

# Create example metadata
metadata_example <- tibble::tibble(sample_id = 1:10,
                                   group = factor(rep(c("Control", "Treatment"), each = 5)))

pca_plot <- ggpicrust2::pathway_pca(t(abundance_example), metadata_example, "group")






