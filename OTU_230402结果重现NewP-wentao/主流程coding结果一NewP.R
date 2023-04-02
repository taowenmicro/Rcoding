
library(fs)

# # 建立结构保存一级目录#--------
result_path <- paste("./","/result_and_plot/",sep = "")
dir_create(result_path)

#导入R包#-------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(EasyStat)
# library(EasyMicrobiome)

# 设置主题#-----
library(ggthemes)
mytheme1 = theme_bw() + theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  legend.position="right",
  
  legend.title = element_blank(),
  legend.background=element_blank(),
  legend.key=element_blank(),
  # legend.text= element_text(size=7),
  # text=element_text(),
  # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
  plot.title = element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
  axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
  axis.text = element_text(size = 20,face = "bold"),
  axis.text.x = element_text(colour = "black",size = 14),
  axis.text.y = element_text(colour = "black",size = 14),
  legend.text = element_text(size = 15,face = "bold")
)

mytheme2 = theme_bw() + theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  legend.position="right",
  
  legend.title = element_blank(),
  legend.background=element_blank(),
  legend.key=element_blank(),
  # legend.text= element_text(size=7),
  # text=element_text(),
  # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
  plot.title = element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
  axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
  axis.text = element_text(size = 20,face = "bold"),
  axis.text.x = element_text(colour = "black",size = 14,angle = 90),
  axis.text.y = element_text(colour = "black",size = 14),
  legend.text = element_text(size = 15,face = "bold")
)

#设定颜色#------------
library(RColorBrewer)#调色板调用包
#调用所有这个包中的调色板
RColorBrewer::display.brewer.all()
#提取特定个数的调色板颜色，会出图显示
RColorBrewer::display.brewer.pal(9,"Set1")
colset1 <- c(brewer.pal(12,"Set1"))
colset2 <- brewer.pal(12,"Paired")
colset3 <- c(brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"))
colset4 = colset3



# 微生物组数据输入#------
ps0 = readRDS("./ps.rds")

res1path <- paste(result_path,"/Base_diversity_16s",sep = "")
dir.create(res1path)

#--检查map文件，是否需要修改#----------
map  =sample_data(ps0)
head(map)
map$Group = gsub("N","",map$Group)
map$Group = gsub("B0_0","Control",map$Group)
head(map)
sample_data(ps0) = map
map
#--样本筛选-根据分组和ID等map文件中的信息#-------
ps_sub <- subset_samples(ps0,!zone %in% c("BN"));ps_sub
map  =sample_data(ps_sub)
head(map)
map$Group
ps0 <- ps_sub 
# 序列筛选-根际七大门类信息或者OTU的ID信息#---------
ps0 <- ps0 %>%
  subset_taxa(
    # Kingdom == "Fungi"
    # Kingdom == "Bacteria"
    # Genus  == "Genus1"
    # Species %in%c("species1") 
    # row.names(tax_table(ps0))%in%c("OTU1")
  )
ps0

#--过滤OTU#----------
ps0 = filter_taxa(ps0, function(x) sum(x ) > 0 , TRUE);ps0

#--最终确定的phyloseq对象定义为ps#--------
ps = ps0


#--提取有多少个分组#-----------
map = sample_data(ps0)
gnum <- unique(map$Group) %>% length()
gnum
# 设定排序顺序
map$Group %>%
  unique()
axis_order = c("Control","RC1","RC5", "RC8", "RF1","RF5","RF8")
# axis_order = c("Control","BC1","BC5", "BC8", "BF1","BF5","BF8")



# 堆叠柱状图展示TOp前是个门,j 展示的水平为门水平#-----
Top = 10
rank.names(ps0)
j = "Phylum"
# j = "Genus"

# -类韦恩图的二分网络设置过滤阈值#-------
# bionum = 200
ps_biost = filter_OTU_ps(ps = ps,Top = 500)




# 差异分析 edger设定分组#------
# group1 = c("Gro1","Gro2")
# b= data.frame(group1)
b = NULL# 如果每两个组之间都做差异，那就指定

# 设置CK，用于双向柱状图绘制#------
CK = "Control"

# 热图和气泡图对丰度最高的前多少个OTU#--------
heatnum　=　20



# 用python做lefse#-----
lefsenum = 0
tax = vegan_tax(ps)
head(tax)
ps_lefse <- ps %>%
  subset_taxa(
    # Kingdom == "Fungi"
    Kingdom == "k__Bacteria"
    # Genus  == "Genus1"
    # Species %in%c("species1") 
    # row.names(tax_table(ps0))%in%c("OTU1")
  )

ps_lefse = filter_OTU_ps(ps = ps_lefse,Top = 400)

# #文件预处理
# format_input.py LEFSE_to_run_G_level.txt pri_lefse.in -c 1 -u 2 -o 1000000
# # 注意这里 –c用来指定分组信息-u 1指定样品信息
# #文件分析,这里-l设置LDA阈值，默认为2，我们使用4 会更加严格
# ~/src/nsegata-lefse/run_lefse.py pri_lefse.in pri_lefse_2.res  -l 2
# #柱状图绘制
# plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
# #树状图绘制
# plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf
# #做每个差异的柱状图
# mkdir biomarkers_raw_images
# plot_features.py pri_lefse.in pri_lefse_2.res biomarkers_raw_images/


#--R 语言做lefse法分析-过滤#----------
# ps_sub = filter_taxa(ps0, function(x) sum(x ) > 1000 , TRUE);ps_sub
ps_Rlefse = filter_OTU_ps(ps = ps,Top = 400)


#---机器学习部分#------
# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
ROC = FALSE
# 是否做交叉检验
rfcv = FALSE
# 选择多少个重要变量
optimal = 40




# 网络
# 过滤多少丰度的做网络
N = 200
zipi = FALSE






#---附件#-----------

# #--修改OTUID
# otu = t(vegan_otu(ps0))
# head(otu)
# row.names(otu) = paste("OTU_",1:length(row.names(otu)),sep = "")
# 
# tax = vegan_tax(ps0)
# head(tax)
# row.names(tax) = paste("OTU_",1:length(row.names(otu)),sep = "")
# 
# ps0 = phyloseq::phyloseq(otu_table(otu,taxa_are_rows=  T),
#                          tax_table(tax),
#                          sample_data(ps0)
#                          )

# #--分类等级水平合并
# rank_names(ps)
# psP  <- tax_glom_wt(ps = ps,ranks = "Phylum")
# psC  <- tax_glom_wt(ps = ps,ranks = "Class")
# psO  <- tax_glom_wt(ps = ps,ranks = "Order")
# psF  <- tax_glom_wt(ps = ps,ranks = "Family")
# psG  <- tax_glom_wt(ps = ps,ranks = "Genus")




ps 


#---统计基本微生物数据#------------

#step_1：
where <- "where"
method <- "using 16S rRNA gene amplicon sequencing."
rep = 8
step_1 <- paste("We analyzed the composition of microbial communities in the",where,method)


#step_2 统计测序样本总量每个样品中的序列数量
a = sum(sample_sums(ps))
b = "high-quality sequences were obtained"
# 统计样本数量
b1 = paste("across the" ,length(sample_sums(ps)),"samples",sep = " ")
b1
#统计重复数量
repead <- paste("For this analysis, we collected",rep,"repeats for each samples.",sep = " ")
repead 
#合并句子
each_count <- paste(repead,b1,a,b," and an average read count per sample of ",round(mean(sample_sums(ps)),0),"(standard deviation (SD)",round(sd(sample_sums(ps)),2),").",sep = " ")
each_count 

# step3 
aa = vegan_otu(ps)
otu_table = as.data.frame(t(aa))
otu_table = as.matrix(otu_table)
otu_table [otu_table > 0] <-1

otu_table = t(otu_table)

OTU_sum <- colSums(otu_table)


d = length(OTU_sum[OTU_sum > 0])


c = paste("All sequences were clustered into", d, "operational taxonomic units (OTUs).",sep = " ")

sample_tax <- paste(c,"the numbers of OTU, generally ranged between ",
                    min(OTU_sum)," and ",max(OTU_sum)," per sample with an average of ",
                    round(mean(OTU_sum),0),"(SD ",round(sd(OTU_sum)),")",sep = "")
sample_tax 

###step 4 统计门水平的总体丰度信息
library("tidyverse")
Taxonomies <- ps %>%
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)} )%>% 
  psmelt() %>%
  #filter(Abundance > 0.05) %>%
  arrange(Phylum)
iris_groups<- group_by(Taxonomies, Phylum)
ps0_sum <- dplyr::summarise(iris_groups, mean(Abundance), sd(Abundance))
ps0_sum[is.na(ps0_sum)] <- 0
colnames(ps0_sum) = c("ID","mean","sd")
ps0_sum <- dplyr::arrange(ps0_sum,desc(mean))
ps0_sum$mean <- ps0_sum$mean *100
ps0_sum <- as.data.frame(ps0_sum)
a = paste(ps0_sum[1,1],"(",round(ps0_sum[1,2],2),"%"," with sd ",round(ps0_sum[1,3],3),")",sep = " ")
b = paste(ps0_sum[2,1],"(",round(ps0_sum[2,2],2),"%"," with sd ",round(ps0_sum[2,3],3),")",sep = " ")
c = paste(ps0_sum[3,1],"(",round(ps0_sum[3,2],2),"%"," with sd ",round(ps0_sum[3,3],3),")",sep = " ")
d = paste(ps0_sum[4,1],"(",round(ps0_sum[4,2],2),"%"," with sd ",round(ps0_sum[4,3],3),")",sep = " ")
e = paste(ps0_sum[5,1],"(",round(ps0_sum[5,2],2),"%"," with sd ",round(ps0_sum[5,3],3),")",sep = " ")
f = paste(ps0_sum[6,1],"(",round(ps0_sum[6,2],2),"%"," with sd ",round(ps0_sum[6,3],3),")",sep = " ")
g =  paste(ps0_sum[7,1],"(",round(ps0_sum[7,2],2),"%"," with sd ",round(ps0_sum[7,3],3),")",sep = " ")
h = paste(ps0_sum[8,1],"(",round(ps0_sum[8,2],2),"%"," with sd ",round(ps0_sum[8,3],3),")",sep = " ")
i = paste(ps0_sum[9,1],"(",round(ps0_sum[9,2],2),"%"," with sd ",round(ps0_sum[9,3],3),")",sep = " ")
j =  paste(ps0_sum[10,1],"(",round(ps0_sum[10,2],2),"%"," with sd ",round(ps0_sum[10,3],3),")",sep = " ")

tax_sum = paste("The majority of OTU belonged to the phyla",a,b,c,d,e,f,g,h,i,"and",j,".",sep = " ")
tax_sum
##all_first 
line = paste(step_1,each_count ,sample_tax ,tax_sum,sep = "")
line

#--开始分析#--------
#-分析共四个部分，这是第二个部分，代号2，第一部分是扩增子原始序列处理
#---result1 base_diversity analyses-------------
otupath = paste(res1path,"/OTU_230402/",sep = "");otupath
dir.create(otupath)


#--1--基础多样性分析#——-----

#--alpha多样性#---------
alppath = paste(otupath,"/alpha/",sep = "")
dir.create(alppath)

#---多种指标alpha多样性分析加出图-标记显著性
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")
index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )

#--多种组合alpha分析和差异分析出图
alp = alpha(ps = ps,inde="Shannon",group = "Group",Plot = TRUE )
index= alp
head(index)
colnames(index)[1] = "ID"
data = as.tibble(map) %>%inner_join(index,by = "ID")
head(data)
colnames(data)[2] = "Group"
data$time = gsub("RC","",data$Group)
data$time = gsub("RF","",data$time)

tem = data %>% group_by(treat,time) %>%
  summarise(mean(Shannon)) 
head(tem)

tem2 = tibble(treat = "C",time = "Control",`mean(Shannon)` = 7.40)
tem2.2 = tibble(treat = "F",time = "Control",`mean(Shannon)` = 7.40)

tem3 = rbind(tem[-1,],tem2,tem2.2)
# tem[1,2] = "0"
# tem[1,1] = "C"

mi = brewer.pal(11,"BrBG")[c(2,10)]
head(data)
library(ggnewscale)
p = ggplot(data) + geom_point(aes(x = time,y = Shannon,color = treat),pch = 21,size = 4) +
  ggplot2::scale_x_discrete(limits = c("Control","1","5","8")) +
  new_scale_fill() +
  geom_point(data = tem3,aes(x = time,y = `mean(Shannon)`,group = treat,fill = treat),
             pch = 21,size = 6) +
  geom_line(data = tem3,aes(x = time,y = `mean(Shannon)`,group = treat,color = treat),size = 3) +
  ggplot2::scale_fill_manual(values = mi) +
  ggplot2::scale_color_manual(values = c("Black",mi)) + 
  theme_classic()
p
FileName <- paste(alppath,"line_shannon", ".pdf", sep = "")
ggsave(FileName, p, width = 5, height =4,limitsize = FALSE)



#----基于时间序列的beta排#------

betapath = paste(otupath,"/beta/",sep = "")
dir.create(betapath)


# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

method = "NMDS"

result = BetaDiv(ps = ps, group = "Group", dist = "bray",
                 method = method, Micromet = "anosim", pvalue.cutoff = 0.05,
                 pair = F)
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
map = sample_data(ps) %>% as.tibble()
head(map)
plotdata = plotdata %>% inner_join(map)
head(plotdata)

plotdata$time = gsub("RC","",plotdata$Group)
plotdata$time = gsub("RF","",plotdata$time)


# 求均值
cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
cent
cent$treat = c("Control","C","C","C","F","F","F")
cent$time = gsub("RC","",cent$Group)
cent$time = gsub("RF","",cent$time)
head(cent)
cent$Group = NULL
tem2 = tibble(x = -1.15906116,y = 0.56810010,time = "Control",treat = "C")
tem2.2 = tibble(x = -1.15906116,y = 0.56810010,time = "Control",treat = "F")

tem3 = rbind(cent[-1,],tem2,tem2.2)


library(ggnewscale)
head(plotdata)

p = ggplot(plotdata) + geom_point(aes(x = time,y = x,color = treat),pch = 21,size = 3) +
  ggplot2::scale_x_discrete(limits = c("Control","1","5","8")) +
  # ggplot2::scale_color_manual(values = c("Black",mi)) + 
  new_scale_fill() +
  geom_point(data = tem3,aes(x = time,y = x,group = treat,fill = treat),
             pch = 21,size = 4) +
  geom_line(data = tem3,aes(x = time,y = x,group = treat,color= treat),size = 2) +
  ggplot2::scale_fill_manual(values = mi) +
  ggplot2::scale_color_manual(values = c("Black",mi)) +
  theme_classic()
p



FileName <- paste(betapath,"/a2_",method,"bray_star.pdf", sep = "")
ggsave(FileName, p, width = 3, height = 4)

p = ggplot(plotdata) + geom_point(aes(x = time,y = y,color = treat),pch = 21,size = 3) +
  ggplot2::scale_x_discrete(limits = c("Control","1","5","8")) +
  # ggplot2::scale_color_manual(values = c("Black",mi)) + 
  new_scale_fill() +
  geom_point(data = tem3,aes(x = time,y = y,group = treat,fill = treat),
             pch = 21,size = 4) +
  geom_line(data = tem3,aes(x = time,y = y,group = treat,color= treat),size = 2) +
  ggplot2::scale_fill_manual(values = mi) +
  ggplot2::scale_color_manual(values = c("Black",mi)) +
  theme_classic()
p



FileName <- paste(betapath,"/a2_",method,"bray_stary.pdf", sep = "")
ggsave(FileName, p, width = 3, height = 4)






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

#----差异分析结合聚类分析-流动图形#------


#-----微生物组数据DEsep2——矫正效应#-----------
library(DESeq2)
library(biobroom)
library(tidyverse)
library(ggClusterNet)

otu <- ps %>% vegan_otu() %>% t() %>%
  as.data.frame()
head(otu)
map <- ps %>% sample_data() 
map$ID = row.names(map)
map = map %>% as.tibble() %>%
  select(ID,Group,everything())
map$Group = as.factor(map$Group)

dds<- DESeqDataSetFromMatrix(countData = otu,
                             colData = map,
                             design = ~ Group)

dds <- DESeq(dds)
dds
# contrasts <- vector(mode = "list")

#--创建两两分组
# tem = resultsNames(dds) %>% strsplit("_") %>% sapply(`[`, 2)
# tem.1 = resultsNames(dds) %>%
#   strsplit("p_") %>% sapply(`[`, 2) %>%
#   strsplit("_vs_") %>% sapply(`[`, 1)
# 
# tem.2 = resultsNames(dds) %>%
#   strsplit("_vs_") %>% sapply(`[`, 2) 
# 
# 
# tem3 = paste("Group",tem.1[-1],tem.2[-1],sep = "=") %>%
#   strsplit("=") 
# names(tem3) = paste(tem.1[-1],tem.2[-1],sep = "_")
# contrasts = tem3


Desep_group <- ps %>% 
  phyloseq::sample_data() %>%
  .$Group %>%
  as.factor() %>%
  levels() %>%
  as.character()
aaa = combn(Desep_group,2)

tem.1 = aaa[1,]
tem.2 = aaa[2,]
tem3 = paste("Group",tem.1,tem.2,sep = "=") %>%
  strsplit("=") 
names(tem3) = paste(tem.1,tem.2,sep = "_")
contrasts = tem3


results <- vector(mode = "list")
shrinkFC <- vector(mode = "list")


library("apeglm")

for(i in seq_along(contrasts)) {
  tem = (DESeq2::results(dds, contrast = contrasts[[i]])%>% as.data.frame()) %>% mutate(compareG = i)
  tem$OTU_ID = row.names(tem)
  results[[names(contrasts)[i]]] <- tem
  tem = (lfcShrink(dds, contrast = contrasts[[i]],type = "normal") %>% as.data.frame() ) %>% mutate(compareG = i)
  tem$OTU_ID = row.names(tem)
  shrinkFC[[names(contrasts)[i]]] <- tem
  print(contrasts[[i]])
}

results.df <- plyr::ldply(results, function(x) x)
names(results.df)[1] <- "Contrast"


shrinkFC.df <- plyr::ldply(shrinkFC, function(x) x)
names(shrinkFC.df)[1] <- "Contrast"

rs.fc <- shrinkFC.df
write.csv(rs.sig,"图1D全部差异微生物表格.csv")


#--基于desep2的差异微生物数据取子集
rs.sig <- rs.fc %>% 
  select(-compareG) %>% 
  # separate(Contrast, c("Treatment", "Time"), sep = "\\.") %>% 
  # filter(Treatment != "WC_TRN") %>% 
  filter(!is.na(padj)) %>% 
  ungroup() %>% 
  mutate(p.adjusted2 = p.adjust(pvalue, method = "fdr")) %>%
  ungroup() %>% 
  filter(p.adjusted2 < 0.05)
rs.sig$OTU_ID %>% unique() %>% length()

rs.fc.mtx <- rs.fc %>% 
  # separate(Contrast, c("Treatment2", "Time"), sep = "\\.") %>% 
  # mutate(Treatment = Treatment2) %>% 
  # filter(Treatment2 != "WC_TRN") %>% 
  filter(OTU_ID %in% rs.sig$OTU_ID) %>% 
  # mutate(Contrast = paste(Treatment2, Time, sep = ".")) %>% 
  select(Contrast, OTU_ID, lfcSE) %>% 
  spread(key = Contrast, value = lfcSE)
rs.fc.mtx <- as.data.frame(rs.fc.mtx) 
rownames(rs.fc.mtx) <- rs.fc.mtx$OTU_ID 
rs.fc.mtx <- rs.fc.mtx[,-1] 

rs.k <- 10
#-就算微生物矩阵
rs.dist <- dist(as.matrix(rs.fc.mtx)) 
rs.clust <- hclust(rs.dist, method = "ward.D") 
rs.ord.names <- rs.clust$labels[rs.clust$order] 
rs.ord <- data.frame(OTU_ID = rs.ord.names, order = 1:length(rs.ord.names))

rs.cut <- cutree(rs.clust[c(1,2,4)], k = rs.k)
rs.ord$Cluster <- as.factor(rs.cut[rs.ord$OTU_ID])

rs.ord <- rs.ord %>% 
  mutate(Cluster = paste("C", Cluster, sep = "")) %>% 
  group_by(Cluster) %>% 
  mutate(nOTU = n()) %>% 
  ungroup() %>% 
  mutate(Cluster2 = paste(Cluster, " (", nOTU, " OTUs)", sep = ""))


source("E:/Shared_Folder/Study_project/21年干旱natureplant/RiceDroughtRecovery-master/General/rmb_functions.R")
### with relative abundances

#--相对丰度转化
otu <- rel_ab(otu)
# 转化为长数据-根际数据
rs.otu.tidy <- tidy_otu(otu) %>% 
  mutate(Count = Count/100) %>% 
  filter(!is.na(Count)) 
head(rs.otu.tidy )
head(map)

colnames(map)[1] = "SampleID"
rs.zs.tidy <- rs.otu.tidy %>%
  inner_join(map, by = "SampleID") %>%
  # inner_join(rs.ord,by = "OTU_ID") %>%
  filter(OTU_ID %in% rs.sig$OTU_ID) %>%
  group_by(OTU_ID) %>%
  mutate(zscore = (Count - mean(Count))/sd(Count)) %>%
  group_by(SampleID,Group, OTU_ID) %>%
  summarise(MeanZS = mean(zscore)) %>%
  ungroup()



# trt.lines <- data.frame(Treatment = c("KO", "KO", "OE"),
#                         Treatment2 = c("OE", "WT", "WT"),
#                         Contrast = c("KO vs OE", "KO vs WT", "OE vs WT"),
#                         IniTreatment = c(4.5, 4.5, 4.5),
#                         EndTreatment = c(5.5, 6.5, 7.5))

rs.hm.p <- rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  mutate(MeanZS2 = ifelse(MeanZS > 2,2,MeanZS)) %>%
  # mutate(Treatment = fct_relevel(Treatment, "WC")) %>%
  ggplot(aes(SampleID, reorder(OTU_ID, order), fill = MeanZS2)) +
  geom_tile() +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment*10), linetype = 3, color = "white") +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment*10), linetype = 3, color = "white") +
  scale_fill_viridis_c(name = "Mean Rel. Abund.\n(z-score)") +
  ylab("Differentially Abundant OTU") +
  xlab("Group") +
  # scale_fill_gradientn(name = "log2FC",
  #                      colors = RColorBrewer::brewer.pal(11, "BrBG")) +
  facet_grid(Cluster ~ Group, scales = "free", space = "free") +
  theme_void() +
  theme(text = element_text(size = 15),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title.align = 1)

rs.hm.p
ggsave("4.2.pdf",width = 8,height = 20)

head(rs.ord)

tax = ps %>% vegan_tax() %>%
  as.data.frame()
head(rs.ord)
tax$Genus %>% unique()
tem = tax[tax$Genus == "g__Bacillus",]
id = row.names(tem)
tem = tax[tax$Genus == "g__Sphingomonas",]
id2 = row.names(tem)
tem = tax[tax$Genus == "g__Kaistobacter",]
id3 = row.names(tem)
id3

rs.ord %>% filter(OTU_ID  %in% id3)
tem = rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  filter(Cluster == "C2" )
tem$OTU_ID %>% unique()



rs.hm.p <- rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  filter(Cluster == "C2" ) %>%
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>%
  # mutate(Treatment = fct_relevel(Treatment, "WC")) %>%
  ggplot(aes(SampleID, reorder(OTU_ID, order), fill = MeanZS)) +
  geom_tile() +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment*10), linetype = 3, color = "white") +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment*10), linetype = 3, color = "white") +
  scale_fill_viridis_c(name = "Mean Rel. Abund.\n(z-score)") +
  ylab("Differentially Abundant OTU") +
  xlab("Group") +
  # scale_fill_gradientn(name = "log2FC",
  #                      colors = RColorBrewer::brewer.pal(11, "BrBG")) +
  facet_grid(Cluster ~ Group, scales = "free", space = "free") +
  theme_classic() +
  theme(text = element_text(size = 15),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title.align = 1)

# rs.hm.p
ggsave("5.pdf",width = 8,height = 6)

tem = rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  group_by(Cluster,Group) %>%
  summarise(mean(MeanZS)) %>%
  filter(Cluster == "C2")
head(tem)

write.csv(tem,"图1D中C2模块差异微生物表格.csv")

p = rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  group_by(Cluster,Group) %>%
  summarise(mean(MeanZS)) %>%
  ggplot(aes(x = Group,y = `mean(MeanZS)`, group = Cluster)) +
  geom_point(pch = 21,alpha = .3) +geom_line(alpha = .5) +
  geom_line(data = tem,aes(x = Group,y = `mean(MeanZS)`, group = Cluster),
            color = "#01665E",size = 3) +
  theme_classic()

ggsave("cluster_change.pdf",width = 6,height = 5)




p = rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  group_by(Cluster,Group) %>%
  summarise(mean(MeanZS)) %>%
  ggplot(aes(x = Group,y = `mean(MeanZS)`, group = Cluster)) +
  geom_point(aes(fill = Cluster),pch = 21) +geom_line(aes(color = Cluster)) +
  # scale_fill_manual(values = RColorBrewer::brewer.pal(11, "BrBG")) +
  # scale_color_manual(values = RColorBrewer::brewer.pal(11, "BrBG")) +
  theme_classic()

ggsave("cluster.pdf",width = 6,height = 5)




rs.hm.p <- rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  filter(Cluster == "C2" ) %>%
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>%
  # mutate(Treatment = fct_relevel(Treatment, "WC")) %>%
  ggplot(aes(SampleID, reorder(OTU_ID, order), fill = MeanZS)) +
  geom_tile() +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment*10), linetype = 3, color = "white") +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment*10), linetype = 3, color = "white") +
  scale_fill_viridis_c(name = "Mean Rel. Abund.\n(z-score)") +
  ylab("Differentially Abundant OTU") +
  xlab("Group") +
  # scale_fill_gradientn(name = "log2FC",
  #                      colors = RColorBrewer::brewer.pal(11, "BrBG")) +
  facet_grid(Cluster ~ Group, scales = "free", space = "free") +
  theme_classic() +
  theme(text = element_text(size = 15),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title.align = 1)

# rs.hm.p
ggsave("C2_cluster.pdf",width = 8,height = 20)

tax = ps %>% vegan_tax() %>%
  as.data.frame()
tax$OTU_ID = row.names(tax)

otu = ps %>% tax_glom_wt(6) %>%
  vegan_otu() %>% #t() %>%
  as.data.frame()
head(otu)

tem = colMeans(otu) %>% as.data.frame()
colnames(tem) = "mean"
tem$ID = row.names(tem)

tem2 = tem %>% arrange(desc(mean))
head(tem2)

tem3 = rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  inner_join(tax,by = "OTU_ID") %>%
  filter(Cluster == "C2") %>%
  filter(!Genus %in% c("","1.00","0.67","3","g__")) %>%
  group_by(Genus,Group) %>%
  summarise(mean(MeanZS)) %>%
  inner_join(tem,by = c("Genus" = "ID")) %>%
  arrange(Group,desc(mean)) 
head(tem3)

tem3$Genus = factor(tem3$Genus,levels = unique(tem3$Genus )[49:1])
tem3 %>%
  ggplot(aes(x= Group, y = Genus, fill = `mean(MeanZS)`)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Rel. Abund.\n(z-score)") +
  ylab("Differentially Abundant OTU") +
  xlab("Group") +
  # scale_fill_gradientn(name = "log2FC",
  #                      colors = RColorBrewer::brewer.pal(11, "BrBG")) +
  # facet_grid(. ~ Group, scales = "free", space = "free") +
  theme_classic() +
  theme(text = element_text(size = 15),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title.align = 1)

ggsave("C2_Genus.pdf",width = 8,height = 10)

