
#--代谢组学数据
#-非靶向代谢组学数据类似环境数据
#但是由于代谢物数量较多，所以需要单独分析

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


ps = ps1



map$Group  %>% unique()
#--分作物进行排序分析#----
# ps = ps0 %>% subset_samples.wt("Group","ExCtrl.NA",TRUE)

# 设定排序顺序
axis_order = phyloseq::sample_data(ps)$Group %>% unique()
# axis_order =c("Arabidopsis-00h","Arabidopsis-02h","Arabidopsis-04h",
#   "Arabidopsis-p5h",
#   "Arabidopsis-1d",  "Arabidopsis-4d"  )


#---主题颜色设置#-------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro//total_amplicon.R")
#---扩增子环境布置

res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = brewer.pal(9,"Set1")
colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

#--提取有多少个分组
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum#-设定结果保存路径
repath = "./GCMS_result_and_plot_R3/"
fs::dir_create(repath)


#--这部分代码正在调试--尚不可使用-这部分问题先不解决
#-本来设计的使用爬虫实时调用，但是网络问题很多
# -所以就使用下载数据库来进行，现在正在权衡中

#--1 代谢物注释HMDB和KEGG数据库#------
id = ps %>% ggClusterNet::vegan_otu() %>% t() %>%
  as.data.frame() %>% row.names()

#-HMDB数据库注释
source("E:\\Shared_Folder\\Function_local\\R_function\\micro//ann.HMDB.R")
repath2 = "E:/Shared_Folder/Function_local/R_function/micro/"
tax1 = ann.HMDB (id = id,repath  = repath2 )
colnames(tax1)

tax1 = tax1 %>% distinct(id,.keep_all = TRUE) %>%
  column_to_rownames("id")
head(tax1)
phyloseq::tax_table(ps) = as.matrix(tax1)
ps


#--0基本表格保存#----------
tabpath = paste(repath,"/report_table/",sep = "")
dir.create(tabpath)
#--raw otu tab
otu = as.data.frame(t(ggClusterNet::vegan_otu(ps)))
head(otu)
FileName <- paste(tabpath,"/otutab.csv", sep = "")
write.csv(otu,FileName,sep = "")
# tax table
tax = as.data.frame((ggClusterNet::vegan_tax(ps)))
head(tax)
FileName <- paste(tabpath,"/tax.csv", sep = "")
write.csv(tax,FileName,sep = "")

ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
#--norm otu tab
otu_norm = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
FileName <- paste(tabpath,"/otutab_norm.csv", sep = "")
write.csv(otu_norm,FileName,sep = "")

otutax <- cbind(as.data.frame(t(ggClusterNet::vegan_otu(ps_rela))),as.data.frame((ggClusterNet::vegan_tax(ps_rela))))
FileName <- paste(tabpath,"/otutax_norm.csv", sep = "")
write.csv(otutax,FileName,sep = "")


#2代谢物ggplot升级版本韦恩图和Upset#---
ps1 = rm.low.area(ps = ps,threshold = 100000)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven.Upset.gg.R")

# BiocManager::install("ggupset")

if (gnum < 6) {
  Venpath = paste(repath,"/Ven_Upset_super/",sep = "")
  dir.create(Venpath)
  
  library(ggVennDiagram)
  res = Ven.Upset(ps =  ps1,
                  group = "Group",
                  N =1,
                  size = 3)
  
  p1 = res[[1]]
  p2 = res[[2]]
  
  filename3 <- paste(Venpath,"Ven_gg.pdf", sep = "")
  ggsave(filename3, p1, width = 8, height = 8)
  filename3 <- paste(Venpath,"Upset_gg.pdf", sep = "")
  ggsave(filename3, p2, width = 8, height = 8)
}


# 2核心代谢物分析#---
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/CoreSeper.GC.R")
source("E:/Shared_Folder/Function_local/R_function/micro/barMainplot.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven-Upset.R")

Venpath = paste(repath,"/Core_metabolite_super/",sep = "")
dir.create(Venpath)
ps1 = rm.low.area(ps = ps,threshold = 100000)
#---每个部分
result =CoreSeper.GC(ps = ps1,
                     N= 1,
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


#--3 分类堆叠柱状图#--------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot_GC.R")
barpath = paste(repath,"/Microbial_composition/",sep = "")
dir.create(barpath)

phyloseq::rank_names(ps)
j = "Class"

strbar = c("Super_class", "Sub_class" , "Class" )
# strbar = c("Superclass","Class"  )
# ps  = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)
for (j in strbar) {
  result = barMainplot(ps = ps,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 12)
  p4_1 <- result[[1]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset3) +
    scale_x_discrete(limits = axis_order) +
    mytheme2
  p4_1
  p4_2  <- result[[3]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset3) +
    scale_x_discrete(limits = axis_order) +
    mytheme2
  p4_2
  
  databar <- result[[2]] %>% 
    dplyr::group_by(Group,aa) %>%
    dplyr::summarise(sum(Abundance)) %>% as.data.frame()
  head(databar)
  colnames(databar) = c("Group",j,"Abundance(%)")
  
  
  FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
  ggsave(FileName1, p4_2, width = (5+ gnum), height =8 )
  FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
  ggsave(FileName2, p4_2, width = (5+ gnum), height =8 )
  
  FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
  ggsave(FileName1, p4_1, width = (5+ gnum), height =8 )
  FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
  ggsave(FileName2, p4_1, width = (5+ gnum), height =8 )
  
  FileName <- paste(barpath,"/a2_",j,"_bar_data",".csv", sep = "")
  write.csv(databar,FileName)
}



#--5-分类化合物分组差异#-------
library(EasyStat)
library(ggClusterNet)
barpath = paste(repath,"/Different_Class_EasyStat/",sep = "")
dir.create(barpath)

map = sample_data(ps)
head(map)
map = map[,1:2]
sample_data(ps) = map
for (j in strbar) {
  dat <- ps %>% scale_micro(method = "rela") %>%
    tax_glom_wt(ranks = j) %>%
    vegan_otu() %>% 
    as.data.frame()
  head(dat)
  
  dat$id = row.names(dat)
  
  dat2 = dat %>% 
    dplyr::left_join(as.tibble(sample_data(ps)),by = c("id" = "ID")) %>%
    # dplyr::filter(Group != "qiao") %>%
    dplyr::rename(group = Group) %>%
    select(id,group,everything())
  # dat2 %>%
  #   dim()
  
  dat2$group = as.factor(dat2$group)
  head(dat2)
  dat2$id = gsub("[-]",".",dat2$id)
  dat2$group = gsub("[-]",".",dat2$group)
  result = MuiKwWlx2(data = dat2,num = c(3:dim(dat2)[2]))
  
  FileName <- paste(barpath,"/",j,"_classification_different_label.csv", sep = "")
  write.csv(result,FileName,sep = "")
  FileName <- paste(barpath,"/",j,"_classification_data.csv", sep = "")
  write.csv(dat2,FileName,sep = "")
  
  result1 = EasyStat::FacetMuiPlotresultBox(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",
                                            ncol = 4 )
  p1_1 = result1[[1]] + 
    # scale_x_discrete(limits = axis_order) + 
    mytheme2 +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = colset1)
  p1_1
  
  res = FacetMuiPlotresultBar(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",
                              ncol = 4)
  p1_2 = res[[1]]+
    # scale_x_discrete(limits = axis_order) + 
    guides(color = FALSE) +
    mytheme2 + 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_2
  
  res = FacetMuiPlotReBoxBar(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",ncol = 4)
  p1_3 = res[[1]]+ 
    # scale_x_discrete(limits = axis_order) + 
    mytheme2 + 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_3
  
  h = dim(dat2)[2]%/%4
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_box", ".pdf", sep = "")
  ggsave(FileName, p1_1, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_bar", ".pdf", sep = "")
  ggsave(FileName, p1_2, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_boxbar", ".pdf", sep = "")
  ggsave(FileName, p1_3, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_box", ".jpg", sep = "")
  ggsave(FileName, p1_1, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_bar", ".jpg", sep = "")
  ggsave(FileName, p1_2, width = 12, height =3*h,limitsize = FALSE)
  
  FileName <- paste(barpath,"/",j,"_classification_Facet_boxbar", ".jpg", sep = "")
  ggsave(FileName, p1_3, width = 12, height =3*h,limitsize = FALSE)
  
  
}


#---4-分组热图#-----------
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\GC_ggheatmap_buplot.R")
rank.names(ps)

barpath = paste(repath,"/Different_Class_heatmap_EasyStat/",sep = "")
dir.create(barpath)



for (j in strbar) {
  ps_rela <- ps %>% scale_micro(method = "rela") %>%
    tax_glom_wt(ranks = "Class")
  
  result <- GCheatmap (ps_rela,
                       label =  F,
                       col_cluster = F,
                       row_cluster = F)
  p1 <- result[[1]] 
  p1
  # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
  p2 <- result[[2]]
  p2
  
  h = taxa_names(ps_rela) %>% length()
  w = sample_names(ps_rela)  %>% length()
  
  
  filename = paste(barpath,"/",j,"_classification_","ggheatmap.pdf",sep = "")
  ggsave(filename,p1,width = w/1.7,height = h/3)
  
  filename = paste(barpath,"/",j,"_classification_","ggbubble.pdf",sep = "")
  ggsave(filename,p2,width = w/1.7,height = h/3)
  
  filename = paste(barpath,"/",j,"_classification_","ggheatmap.png",sep = "")
  ggsave(filename,p1,width = w/1.7,height = h/3)
  
  filename = paste(barpath,"/",j,"_classification_","ggbubble.png",sep = "")
  ggsave(filename,p2,width = w/1.7,height = h/3)
  
}



#---6 排序分析PCA等#-------------
betapath = paste(repath,"/beta_orda_total/",sep = "")
dir.create(betapath)


# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")


methodlist = c("PCoA")
for (method in methodlist) {
  result = BetaDiv(ps = ps %>% scale_micro(), group = "Group", dist = "bray",
                   method = method, Micromet = "anosim",
                   pvalue.cutoff = 0.05,pair = F)
  p3_1 = result[[1]] + 
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) +
    mytheme1 
  p3_1
  #带标签图形出图
  p3_2 = result[[3]] +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_2
  
  FileName <- paste(betapath,"/a2_",method,"bray.pdf", sep = "")
  ggsave(FileName, p3_1, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"",method,"bray.jpg", sep = "")
  ggsave(FileName1 , p3_1, width = 12, height = 12)
  
  FileName <- paste(betapath,"/a2_",method,"bray_label.pdf", sep = "")
  ggsave(FileName, p3_2, width = 12, height = 12)
  FileName1 <- paste(betapath,"/a2_",method,"bray_label.jpg", sep = "")
  ggsave(FileName1 , p3_2, width = 12, height = 12)
  
  # 提取出图数据
  plotdata = result[[2]]
  FileName <-  paste(betapath,"/a2_",method,"bray.csv", sep = "")
  write.csv(plotdata,FileName)
  #---------排序-精修图
  plotdata =result[[2]]
  head(plotdata)
  # 求均值
  cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
  cent
  # 合并到样本坐标数据中
  segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
                by = 'Group', sort = FALSE)
  
  # p2$layers[[2]] = NULL
  # library(ggcor)
  library(ggsci)
  p3_3 = p3_1 +geom_segment(data = segs,
                            mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
    geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_3
  
  FileName <- paste(betapath,"/a2_",method,"bray_star.pdf", sep = "")
  ggsave(FileName, p3_3, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"bray_star.jpg", sep = "")
  ggsave(FileName1 , p3_3, width = 8, height = 8)
  
}

map
#提取总体比较
TResult =result[[5]]
head(TResult)

# 提取两两检测结果
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
write.csv(pair,FileName)
FileName <- paste(betapath,"Total_anosim.csv", sep = "")
write.csv(TResult,FileName)

#--换用adonis差异分析

# title1 = MicroTest(ps = ps, Micromet = "adonis", dist = "bray")
# title1
# FileName <- paste(betapath,"Total_adonis.csv", sep = "")
# write.csv(title1,FileName)
# pairResult = pairMicroTest(ps = ps, Micromet = "adonis", dist = "bray")
# FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
# write.csv(pair,FileName)

#---7 层次聚类#--------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/cluster_plot.R")
clupath = paste(repath,"/cluster_plot/",sep = "")
dir.create(clupath)
res = cluster_plot (ps= ps,hcluter_method = "complete",
                    dist = "bray",cuttree = gnum,row_cluster = T,col_cluster =  T)

p0 = res[[1]]
p0



FileName <- paste(clupath,"cluster", ".jpg", sep = "")
ggsave(FileName, p0, width = 6, height =8,limitsize = FALSE)
FileName <- paste(clupath,"cluster", ".pdf", sep = "")
ggsave(FileName, p0, width = 6 , height = 8,limitsize = FALSE)

p1 = res[[2]]
p2 = res[[3]]

FileName <- paste(clupath,"heatmap_cluster", ".jpg", sep = "")
ggsave(FileName, p1, width = 8, height =8,limitsize = FALSE)
FileName <- paste(clupath,"heatap_cluster", ".pdf", sep = "")
ggsave(FileName, p1, width = 8 , height = 8,limitsize = FALSE)

FileName <- paste(clupath,"bubble_cluster", ".jpg", sep = "")
ggsave(FileName, p2, width = 8, height =8,limitsize = FALSE)
FileName <- paste(clupath,"bubble_cluster", ".pdf", sep = "")
ggsave(FileName, p2, width = 8 , height = 8,limitsize = FALSE)

dat = res[4]
FileName <- paste(clupath,"clu_data.csv", sep = "")
write.csv(dat,FileName)


#-----差异代谢物#----------
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")
alppath = paste(repath,"/All_different_metabolites/",sep = "")
dir.create(alppath)

#--非参数检验
result = statSuper(ps = ps,group  = "Group",artGroup = NULL,method = "wilcox")
head(result)
FileName <- paste(alppath,"/data_wlx_all.compounds.csv", sep = "")
write.csv(result,FileName,sep = "")
#--t检验检验--建议四个重复以上
result = statSuper(ps = ps,group  = "Group",artGroup = NULL,method = "ttext")
head(result)
FileName <- paste(alppath,"/data_ttest_all.compounds.csv", sep = "")
write.csv(result,FileName,sep = "")





source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper2.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Mui.cluster-group.volcano.R")

# 聚类差异火山图指定分组#-------
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Plot.CompareWithCK.R",encoding = "utf-8")
id = sample_data(ps)$Group %>% unique()
# group1 = c("Group1","Group2")
# b= data.frame(group1)
repath = "./230424/"
diffpath.1 = paste(repath,"/Mui.cluster.v/",sep = "")
fs::dir_create(diffpath.1)
map
i = 1
for (i in 1:length(id)) {
  aaa = combn(id,2)
  diffpathv = paste(diffpath.1,"/Mui.cluster.v",paste(aaa[,i][1],aaa[,i][2],sep = "_"),sep = "")
  dir.create(diffpathv)
  # sample_data(ps)
  ps2 = ps %>% subset_samples.wt("Group",aaa[,i])
  res = EdgerSuper(ps = ps2,group  = "Group",artGroup = NULL,
                   j = "OTU",
                   path = diffpathv
  )
  
  head(res)
  
  result = Mui.cluster.volcano(res = res,rs.k = 5)
  p = result[[1]]
  # p
  
  p1 = result[[2]]
  # p
  filename = paste(diffpathv,"/","Mui.cluster.volcano.label.pdf",sep = "")
  ggsave(filename,p,width = 8,height = 4)
  filename = paste(diffpathv,"/","Mui.cluster.volcano.pdf",sep = "")
  ggsave(filename,p1,width = 8,height = 4)
}


#----多组差异分析火山图#------
library(ggrepel)
library(
  tidyverse
)
library(phyloseq)
library(ggClusterNet)
otupath = "./GCMS_result_and_plot/"
diffpath.1 = paste(otupath,"/Mui.Group.v/",sep = "")
dir.create(diffpath.1)
res = EdgerSuper (ps = ps,group  = "Group",artGroup =NULL,
                  j = "OTU",
                  path = diffpath.1
)
head(res)
res$id = row.names(res)


tax = ps %>% vegan_tax() %>% as.data.frame() %>% rownames_to_column("id")
head(tax)


res2 = res %>% left_join(tax,by = "id")

result = Mui.Group.volcano(res = res)
p = result[[1]]
p
filename = paste(diffpath.1,"/","Mui.group.volcano.pdf",sep = "")
ggsave(filename,p,width = 12,height = 6,limitsize = FALSE)

FileName <- paste(diffpath.1,"/data_diff_all.compounds.csv", sep = "")
write.csv(res2,FileName,sep = "")



#---8 单变量统计分析-箱线图等可视化#--------

#--提取差异代谢物标签
# head(result)

alppath = paste(repath,"/summary_stat_plot/",sep = "")
dir.create(alppath)

dat = ps %>% 
  ggClusterNet::vegan_otu() %>% 
  as.data.frame()
head(dat)
# colnames(map)
map = sample_data(ps) %>% as.tibble() %>%
  select(ID,Group)
data = cbind(map[,c(1,2)],dat)
head(data)
colnames(data)[2] = "group"
num = c(3:ncol(data))

# num = 18#--为了减少运行压力，修改为18个化合物
num

#--分割数据
n.fac = length(num)/ 25 
n.fac2 = ceiling(n.fac)
A = list()
# j  =2
for (j in 1:n.fac2) {
  
  if (j == 1) {
    A[[j]] = num[1:25]
  } else if(j != n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:y]
  }else if (j == n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:length(num)]
    
  }
  
}



for (i in 1:n.fac2) {
  result = EasyStat::MuiaovMcomper2(data = data,num = A[[i]])
  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = A[[i]],
                                            result = result,
                                            sig_show ="abc",ncol = 5 )
  p1_1 = result1[[1]] + 
    ggplot2::scale_x_discrete(limits = axis_order) + 
    mytheme2 +
    ggplot2::guides(fill = guide_legend(title = NULL)) +
    ggplot2::scale_fill_manual(values = colset1)
  p1_1
  
  res = EasyStat::FacetMuiPlotresultBar(data = data,num = A[[i]],
                                        result = result,sig_show ="abc",ncol = 5)
  p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) +
    mytheme2+ 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_2
  
  res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = A[[i]],
                                       result = result,sig_show ="abc",ncol = 5)
  p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + 
    mytheme2 + 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_3
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_box", ".pdf", sep = "")
  ggsave(FileName, p1_1, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_bar", ".pdf", sep = "")
  ggsave(FileName, p1_2, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_boxbar", ".pdf", sep = "")
  ggsave(FileName, p1_3, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_box", ".jpg", sep = "")
  ggsave(FileName, p1_1, width = 18, height =16,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_bar", ".jpg", sep = "")
  ggsave(FileName, p1_2, width = 18, height =16,limitsize = FALSE)
  
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"Facet_boxbar", ".jpg", sep = "")
  ggsave(FileName, p1_3, width = 18, height =16,limitsize = FALSE)
  # result
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"_data_aov_abc.csv", sep = "")
  write.csv(result,FileName,sep = "")
  FileName <- paste(alppath,paste("part_",i,sep = ""),"_data.csv", sep = "")
  write.csv(data,FileName,sep = "")
  
  res = EasyStat::MuiHeatmapBubplot(
    data = data,
    i =A[[i]],
    col_cluster = F,
    row_cluster = F,
    label = TRUE,
    result = result,
    sample = TRUE,
    scale = TRUE
  )
  p1 = res[[1]]
  p1
  h = sample_names(ps) %>% length()
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatmap", ".jpg", sep = "")
  ggsave(FileName, p1, width = 12, height =h/2,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatap", ".pdf", sep = "")
  ggsave(FileName, p1, width = 12, height =h/2,limitsize = FALSE)
  
  p2 = res[[2]]
  p2
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble", ".jpg", sep = "")
  ggsave(FileName, p2, width = 12, height =h/2,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble", ".pdf", sep = "")
  ggsave(FileName, p2, width = 12 , height = h/2,limitsize = FALSE)
  
  
  res = EasyStat::MuiHeatmapBubplot(
    data = data,
    i =A[[i]],
    result = result,
    col_cluster = F,
    row_cluster = F,
    label = TRUE,
    sample = FALSE,
    scale = TRUE
    
    
  )
  
  p1 = res[[1]]
  p1
  
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatmap_group", ".jpg", sep = "")
  ggsave(FileName, p1, width = gnum*1.5, height =8,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"heatap_group", ".pdf", sep = "")
  ggsave(FileName, p1, width = gnum*1.5, height =8,limitsize = FALSE)
  
  
  p2 = res[[2]]
  p2
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble_group", ".jpg", sep = "")
  ggsave(FileName, p2, width = gnum*1.5, height =6,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"buble_group", ".pdf", sep = "")
  ggsave(FileName, p2, width = gnum*1.5, height =6,limitsize = FALSE)
  
  
  res = EasyStat::value_stackBar(
    data = data,
    i =A[[i]],
    result = result,
    add_abc = TRUE)
  
  
  p1 = res[[1]]
  p1
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"sample_relative_abundacne", ".jpg", sep = "")
  ggsave(FileName, p1, width = 6, height =5,limitsize = FALSE)
  
  FileName <- paste(alppath,paste("part_",i,sep = ""),"sample_relative_abundacne", ".pdf", sep = "")
  ggsave(FileName, p1, width = 6, height = 5,limitsize = FALSE)
  
}

#---9 载荷矩阵挑选重要代谢物#------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/loadingPCA.R")
pcapath = paste(repath,"/loadingPCA/",sep = "")
dir.create(pcapath)
res = loadingPCA(ps = ps,Top = 20)

otu = ps %>% vegan_otu() %>% t() %>%
  as.data.frame()
head(otu)
p = res[[1]]
p
dat = res[[2]]

filemane = paste(pcapath,"/PCALoading.pdf",sep = "")
ggsave(filemane, p, width = 8, height = 6)
filemane = paste(pcapath,"/PCALoading.jpg",sep = "")
ggsave(filemane, p, width = 8, height = 6)
FileName <- paste(pcapath,"/Loadsing_pca.csv", sep = "")
write.csv(dat,FileName,sep = "")


#--10 机器学习#-------
matpath = paste(repath,"/Machine_learing/",sep = "")
dir.create(matpath )
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\MicroMachine_learning.R")
ROC  = F
rfcv = F
# library(randomForest)
# library(caret)
# library(ROCR) ##用于计算ROC
# library(e1071)

if (ROC ) {
  #--三种机器学习方法评测
  result = MicroRoc( ps = ps,group  = "Group")
  #--提取roc曲线
  p <- result[[1]] + 
    mytheme1
  p
  #提取AUC值
  data <- result[[2]]
  
  filename = paste(matpath,"/three_method_AUCvalue.csv",sep = "")
  write.csv(data,filename,quote = F)
  
  data <- result[[3]]
  filename = paste(matpath,"/three_method_AUCdata.csv",sep = "")
  write.csv(data,filename,quote = F)
  
  filename = paste(matpath,"/three_method_AUC_plot.pdf",sep = "")
  ggsave(filename,p,width = 8,height = 8)
  filename = paste(matpath,"/three_method_AUC_plot.jpg",sep = "")
  ggsave(filename,p,width = 8,height = 8)
  
}


mapping = as.data.frame(phyloseq::sample_data(ps))
optimal = 40
#--随机森林全套-如果圈图尚未显示前面几个，就设定max大一点
result = MicroRF_GC(ps = ps,group  = "Group",optimal = 40,
                    rfcv =F,nrfcvnum = 5,
                    min = -1,max = 5)
#火柴图展示前二十个重要的OTU
p <- result[[1]] + 
  mytheme1
p

filename = paste(matpath,"/randonforest_loading.pdf",sep = "")
ggsave(filename,p,width = 8,height = optimal/3)
filename = paste(matpath,"/randonforest_loading.jpg",sep = "")
ggsave(filename,p,width = 8,height = optimal/3)
# 圈图展示
p <- result[[2]]
p
filename = paste(matpath,"/randonforest_loading_circle.pdf",sep = "")
ggsave(filename,p,width = 8,height = 10)
filename = paste(matpath,"/randonforest_loading_circle.jpg",sep = "")
ggsave(filename,p,width = 8,height = 10)

p <- result[[6]]
p
filename = paste(matpath,"/Show_model.pdf",sep = "")
ggsave(filename,p,width = 8,height = 4)
filename = paste(matpath,"/Show_model.jpg",sep = "")
ggsave(filename,p,width = 8,height = 4)


if (rfcv) {
  # 展示交叉验证结果
  p <- result[[3]]
  filename = paste(matpath,"/randonforest_cross_check.pdf",sep = "")
  ggsave(filename,p,width = 8,height = 12)
  data <- result[[4]]
  filename = paste(matpath,"/randomforest_cross_data.csv",sep = "")
  write.csv(data,filename,quote = F)
}

data <- result[[5]]
filename = paste(matpath,"/randomforest_data.csv",sep = "")
write.csv(data,filename,quote = F)




## 附件#-------

#---0 构建phyloseq对象#------

dat = readxl::read_excel("./data_raw.xlsx",sheet = 1) %>% as.data.frame()
head(dat)

dat$ID = paste("M",1:length(row.names(dat)),sep = "")
dat = dat %>% distinct(ID, .keep_all = TRUE)  
row.names(dat) = dat$ID
dat$ID = NULL
dat = as.matrix(dat)
dat[is.na(dat)] = 0



map = readxl::read_excel("./data_raw.xlsx",sheet = 2) %>% as.data.frame()
head(map)
colnames(map)[1] = "ID"
row.names(map) = map$ID

#--注释文件
tax = readxl::read_excel("./data_raw.xlsx",sheet = 3) %>% as.data.frame()
head(tax)
tax$KEGG.ID
tax$ID = paste("M",1:length(row.names(dat)),sep = "")
# row.names(tax) = tax$Name
row.names(tax) = tax$ID
tax$ID = NULL

tax = tax %>% filter(!is.na(KEGG.ID) ) %>%
  filter(KEGG.ID!= "NA" )
tax = as.matrix(tax)



ps = phyloseq::phyloseq(
  phyloseq::otu_table(as.matrix(dat),taxa_are_rows = TRUE),
  phyloseq::sample_data(map),
  tax_table(as.matrix(tax))
  
)

head(tax)
ps = cg.Rownm.ps(ps,id = "Name")


saveRDS(ps,"./ps_GC.rds")



