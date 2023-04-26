
#--核心根系分泌物

library(tidyverse)
library(phyloseq)
library(ggClusterNet)
#---开始分析#--------
ps0 = readRDS("./ps_result3.rds")

#去除NA值#-----
otu = ps0 %>% vegan_otu()
otu[is.na(otu)] = 0
otu_table(ps0) = otu_table(otu,taxa_are_rows = FALSE)

#--整理map文件#----
map = sample_data(ps0)
head(map)

map$zone = map$Group %>% strsplit( "-") %>% 
  sapply(`[`, 1)
map$ste = map$ID %>% strsplit( "-") %>%
  sapply(`[`,6)
map$ste2 = map$ID %>% strsplit( "-") %>%
  sapply(`[`,7)

map$plant = map$Group %>% strsplit( "-") %>% 
  sapply(`[`, 3)
map$Group = paste(map$zone,map$plant,sep = ".")
map$gro2 = paste(map$ste,map$ste2,sep = ".")
map$gro2 = gsub("a","",map$gro2 )
map$gro2 = gsub("b","",map$gro2 )
map$gro2 = gsub("c","",map$gro2 )
sample_data(ps0) = map


# ps = scale.solu(ps = ps0,col.nm = "Group",solu.ck = c("ExCtrl.NA"))
ps1 = scale.black (ps = ps0,
                   black.col1 = "zone",
                   black.nm1 = "ExCtrl",
                   black.col2 = "gro2",
                   black.nm2 = c("3.R","3.S","4.R","4.S"))

otu = ps1 %>% vegan_otu()
otu[is.na(otu)] = 0
otu_table(ps1) = otu_table(otu,taxa_are_rows = FALSE)



ps.root = ps1
#--昼夜根系分泌物#-------
ps0 = readRDS("./ps_result4.rds")

#去除NA值#-----
otu = ps0 %>% vegan_otu()
otu[is.na(otu)] = 0
otu_table(ps0) = otu_table(otu,taxa_are_rows = FALSE)

#--整理map文件#----
map = sample_data(ps0)
head(map)

map$exu = map$Group %>% strsplit( "-") %>% 
  sapply(`[`, 1)
map$EO = map$ID %>% strsplit( "-") %>%
  sapply(`[`, 2)
map$plant = map$Group %>% strsplit( "-") %>% 
  sapply(`[`, 3)
map$Group = paste(map$exu,map$EO,map$plant,sep = ".")

# map$gro2 = paste(map$EO,map$plant ,sep = ".")

sample_data(ps0) = map

sample_data(ps)
ps = scale.solu(ps = ps0,col.nm = "exu",solu.ck = c("InjBl"))
ps1 = scale.black (ps = ps,
                   black.col1 = "plant",
                   black.nm1 = "TxCtrl",
                   black.col2 = "EO",
                   black.nm2 = c("EON","EOD"))

sample_data(ps1)

ps.exu = ps1

#--合并茎根，根系分泌物的数据#------

ps3 = merge.GC(ps1 = ps.exu,
               ps2 = ps.root)


#--1 代谢物注释HMDB和KEGG数据库#------
id = ps3 %>% ggClusterNet::vegan_otu() %>% t() %>%
  as.data.frame() %>% row.names()

#-HMDB数据库注释
source("E:\\Shared_Folder\\Function_local\\R_function\\micro//ann.HMDB.R")
repath2 = "E:/Shared_Folder/Function_local/R_function/micro/"
tax1 = ann.HMDB (id = id,repath  = repath2 )
colnames(tax1)

tax1 = tax1 %>% distinct(id,.keep_all = TRUE) %>%
  column_to_rownames("id")
head(tax1)
phyloseq::tax_table(ps3) = as.matrix(tax1)

#去除NA值#-----
otu = ps3 %>% vegan_otu()
otu[is.na(otu)] = 0
otu_table(ps3) = otu_table(otu,taxa_are_rows = FALSE)



#--合并后计算核心代谢物#----
repath = "./GCMS_result_and_plot_R5/"
dir.create(repath)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven.Upset.gg.R")

# BiocManager::install("ggupset")

  Venpath = paste(repath,"/Ven_Upset_super/",sep = "")
  dir.create(Venpath)
  library(ggVennDiagram)
  
  map = sample_data(ps3)
  map$Group = map$Group %>% strsplit( "[.]") %>%
    sapply(`[`, 1)
  sample_data(ps3) = map
  
  
  res = Ven.Upset(ps =  ps3 ,#%>% subset_samples.wt("Group","shoot",TRUE),
                  group = "Group",
                  N =0.5,
                  size = 3)
  
  p1 = res[[1]]
  
  p1
  p2 = res[[2]]
  
  filename3 <- paste(Venpath,"Ven_gg.pdf", sep = "")
  ggsave(filename3, p1, width = 8, height = 8)
  filename3 <- paste(Venpath,"Upset_gg.pdf", sep = "")
  ggsave(filename3, p2, width = 8, height = 8)


# 2核心代谢物分析#---
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/CoreSeper.GC.R")
source("E:/Shared_Folder/Function_local/R_function/micro/barMainplot.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven-Upset.R")

Venpath = paste(repath,"/Core_metabolite_super_3/",sep = "")
dir.create(Venpath)
# ps1 = rm.low.area(ps = ps,threshold = 100000)
#---每个部分
result =CoreSeper.GC(ps = ps3,#%>% subset_samples.wt("Group","shoot",TRUE),
                     N= 0.8,
                     path = Venpath,
                     group = "Group",
                     Top = 10,
                     j = "Class"
                     
)
# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]
#每个部分序列的数量占比，并作差异
p8 <- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]


FileName <- paste(Venpath,j,"count_Facet_ven", ".pdf", sep = "")
ggsave(FileName, p7_1, width = 15, height = 12)

FileName <- paste(Venpath,j,"diff_count_box", ".pdf", sep = "")
ggsave(FileName, p8, width = 15, height = 12)

FileName <- paste(Venpath,j,"count_Facet_ven_flow", ".pdf", sep = "")
ggsave(FileName, p7_2, width = 15, height = 12)

FileName <- paste(Venpath,j,"count_Facet_ven", ".jpg", sep = "")
ggsave(FileName, p7_1, width = 15, height = 12)

FileName <- paste(Venpath,j,"diff_count_box", ".jpg", sep = "")
ggsave(FileName, p8, width = 15, height = 12)

FileName <- paste(Venpath,j,"count_Facet_ven_flow", ".jpg", sep = "")
ggsave(FileName, p7_2, width = 15, height = 12)






