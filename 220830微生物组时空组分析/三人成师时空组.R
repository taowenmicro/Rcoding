
# 时空数据处理

#这部分数据主要包括时间序列样本和空间区域样本

# 这是高级的多分组模式

# 1-alpha多样性:时空序列和多分组微生物组#-----

ps = readRDS("./data_all/ps.rds")

ps_meta = readRDS("./data_all/ps_meta.rds")

#--对不同空间的基于处理的差异微生物进行聚类的结果
all.clust.summary <- readRDS("./data_all//drought_clusters_updated.RDS")
#---表型数据
phenotypes <- readRDS("./data_all//phenotypes.RDS")
#筛菌菌种序列比对结果
seqs <- readRDS("./data_all//actino_seqs.RDS")
# #--叶片定殖后微生物组数据
ps_verify = readRDS("./data_all/ps_verify.rds")

#干旱影响的植物生育期的变化数据
stages <- read.table("./data_all//growth_stages.txt", header = T, sep = "\t")

# source("./General/rmb_functions.R")
# source("./General/parameters.R")
library(vegan)
library(broom)
library(tidyverse)
library(phyloseq)
library(ggClusterNet)

otu = ps %>% vegan_otu() %>% t() %>%
  as.data.frame()
map = sample_data(ps) %>% as.tibble()
head(map)

# 定义三个分组，处理组，事件组，空间组
group = c("Treatment","Compartment","Age")

#--定义颜色
phy.pal <- c("gray35",
             RColorBrewer::brewer.pal(8, "Set2")[1:4],
             "indianred3",
             RColorBrewer::brewer.pal(8, "Set2")[5:8],
             RColorBrewer::brewer.pal(11, "RdYlBu")[7:10])

trt.pal <- RColorBrewer::brewer.pal(11, "BrBG")[c(10,3,2,1)]
#-排序颜色
cmp.pal <- RColorBrewer::brewer.pal(12, "Set3")[c(5,4)]
#-仅仅计算shannon多样性
alpha <- data.frame(ID = colnames(otu),
                    Shannon = diversity(t(otu))) %>% 
  inner_join(map, by = "ID") 

#-时间序列展示：计算多样性均值
#-构造时间和生育期表格
time.age <- map %>% group_by(Time, Age) %>% count() %>% select(-n)

alpha.rs.means <- alpha %>% 
  group_by(Treatment, Age,Compartment) %>% 
  mutate(Mean = mean(Shannon)) %>% 
  ungroup()

rs.plot <- alpha %>% 
  ggplot(aes(x = Age,y =  Shannon)) +
  # geom_vline(xintercept = ws.line, linetype = 3, color = "gray20", size = 0.5) +
  geom_point(aes(color = Treatment), alpha = 1, shape = 1, size = 2) +
  geom_line(data = alpha.rs.means, aes(y = Mean, color = Treatment), size = 1) +
  xlab("Plant age (days)") +
  ylab("Rhizosphere\nShannon index") +
  scale_fill_manual(values = trt.pal) +
  scale_color_manual(values = trt.pal) +
  # scale_x_continuous(breaks = seq(0,140,20)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 15)) + facet_wrap(.~ Compartment)

rs.plot


# 时间和处理进行模型-差异分析
alpha$TimeFctr # 时间梯度做了一个因子变量，为了进行各种建模操作

#  (1|Row) 注意写法的含义
alpha.rs.lmm <- lmerTest::lmer(Shannon ~ Treatment * TimeFctr * Compartment + (1|Row), data = alpha)
summary(alpha.rs.lmm)

#对建模结果进行检测，处理和时间都对微生物产生了极显著的影响，交互作用不显著
# 分析结果转化为tidy格式，也就是一个tibble。
alpha.rs.aov <- anova(alpha.rs.lmm) %>% 
  tidy()

# 适用于多种方差分析模型，包括重复测量和嵌套设计，其中初始建模将使用‘aov’、‘lm’、‘ez’或‘lme4’(混合模型)输入，进行两两比较，或者时事后比较
alpha.rs.contrasts <- emmeans::emmeans(alpha.rs.lmm, specs = trt.vs.ctrl ~ Treatment|TimeFctr|Compartment, adjust = "none") %>% 
  .$contrasts %>% 
  as.data.frame() %>% 
  separate(contrast, c("Treatment", "Control")) %>% 
  rename("Time" = "TimeFctr") %>% 
  mutate(Time = as.integer(Time))

# 得到两两检测差异表格
alpha.rs.contrasts

# 可视化检测结果，时间 * 空间正好是一个矩阵，所以，显著性可以通过标记展示在图片面板上，然后使用颜色区分是否显著



rs.sig <- alpha.rs.contrasts %>% 
  ungroup() %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% 
  inner_join(time.age, by = "Time") %>% 
  ggplot(aes(x = Age, y = Treatment)) +
  geom_point(aes(color = Treatment, alpha = p.value < 0.05), shape = "*", size = 7) +
  # scale_y_discrete(limits = c("DS3", "DS2", "DS1")) +
  # scale_x_continuous(breaks = seq(0,140,20)) +
  scale_color_manual(values = trt.pal[-1]) +
  scale_alpha_manual(values = c(0,1)) +
  xlab("Plant age (days)") +
  theme_classic() + facet_wrap(.~ Compartment) +
  theme(text = element_text(size = 13),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") 

#-横坐标是时间，纵坐标有两个，一个是alpha多样性指标，另一个就是处理对了
rs.final <- cowplot::plot_grid(rs.sig, rs.plot, ncol = 1, rel_heights = c(4,17), align = "v")
rs.final

head(alpha.rs.contrasts)

alpha.rs.contrasts <- alpha.rs.contrasts %>% 
  ungroup() %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% 
  mutate(Contrast = paste(Treatment, Control, sep = "-")) %>% 
  select(Contrast, Time:p.adj)



# 2-beta多样性:时空序列和多分组微生物组#-----

# source("./General/rmb_functions.R")
# source("./General/parameters.R")
library(vegan)
library(broom)
library(contrast)
library(cowplot)
library(ggConvexHull)
# devtools::install_github("cmartin/ggConvexHull")
library(tidyverse)


dist = "bray"
wuf = ps %>% phyloseq::distance(method=dist)

map = sample_data(ps) %>% as.tibble()
head(map)
pmanova <- adonis(wuf ~ Compartment * as.factor(Time) * Treatment,
                  strata = as.factor(map$Row),  data = map)
pmanova


pcoa_axes <- function(dist, map) {
  require(ape)
  pcoa <- pcoa(dist)
  as.data.frame(pcoa$vectors) %>%
    mutate(ID = rownames(.)) %>%
    inner_join(map, by = "ID")
}

pcoa_eigval <- function (dist, map) {
  require(ape)
  # map <- filter(map, SampleID %in% rownames(dist))
  # dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- pcoa(as.dist(dist))
  eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)
  data.frame( PC = 1:length(eigval), Eigval = eigval, CumEigval = cumsum(eigval))
}

#--基于距离矩阵的pcoa排序分析
pcoa.axes <- pcoa_axes(dist = as.matrix(wuf), map)
#-提取解释度
pcoa.eig <- pcoa_eigval(as.matrix(wuf), map)

head(pcoa.axes)

tem = pcoa.axes[,c("Axis.1","Axis.2",group)]
head(tem)

all.pcoa.p <- pcoa.axes %>% 
  # mutate(Compartment = fct_recode(Compartment,
  #                                 "Rhizosphere" = "RS",
  #                                 "Endosphere" = "ES")) %>% 
  # mutate(Compartment = fct_relevel(Compartment, "Rhizosphere")) %>% 
  ggplot(aes(x= Axis.1,y =  Axis.2)) +
  geom_point(aes(color = Compartment), size = 2, shape = 16) +
  xlab(paste("PCo1 (", pcoa.eig$Eigval[1], "%)", sep = "")) +
  ylab(paste("PCo2 (", pcoa.eig$Eigval[2], "%)", sep = "")) +
  scale_color_manual(values = cmp.pal) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "top")

all.pcoa.p


dist.tidy <- wuf %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  mutate(SampleID.x = row.names(.)) %>% 
  gather(key = "SampleID.y", value = "Distance", -SampleID.x) %>% 
  filter(!is.na(Distance)) %>% 
  filter(Distance > 0) %>% 
  inner_join(select(map, ID, Compartment, Treatment, Time), 
             by = c("SampleID.x" = "ID")) %>% 
  inner_join(select(map, ID, Compartment, Treatment, Time), 
             by = c("SampleID.y" = "ID"))

dist.within <- dist.tidy %>% 
  filter(Compartment.x == Compartment.y) %>% 
  filter(Treatment.x == Treatment.y) %>% 
  filter(Time.x == Time.y)

head(dist.within)
#--组内和组间距离展示
cmp.box <- dist.within %>% 
  mutate(Compartment = Compartment.x) %>%
  # mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>%
  ggplot(aes(Compartment, Distance, color = Compartment)) +
  geom_boxplot(size = 1) +
  geom_rug() +
  ylab("Within-Compartment\nDistance") +
  xlab("Compartment") +
  scale_color_manual(values = cmp.pal) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "top")

cmp.box

#-空间展示
plot_grid(all.pcoa.p,
          cmp.box,
          nrow = 1,
          align = "v",
          labels = c("A", "B"),
          label_size = 20)




#-根际的样本单独展示排序图，按照时间上色
rs.pc.time <- pcoa.axes %>% 
  ggplot(aes(Axis.1, Axis.2, color = Age)) +
  geom_point(size = 2) +
  scale_color_viridis_c(name = "Plant Age\n(days)") +
  xlab(paste("PCo1 (", pcoa.eig$Eigval[1], "%)", sep = "")) +
  ylab(paste("PCo2 (", pcoa.eig$Eigval[2], "%)", sep = "")) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "top",
        legend.title.align = 1)

# 按照处理上色 看不出来规律
rs.pc.trt <- pcoa.axes %>% 
  # mutate(Treatment = fct_recode(Treatment,
  #                               "DS1" = "D1",
  #                               "DS2" = "D2",
  #                               "DS3" = "D3")) %>% 
  # mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  ggplot(aes(Axis.1, Axis.2, color = Treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = trt.pal,
                     guide = guide_legend(title.position = "top",
                                          title.hjust = 0.5)) +
  xlab(paste("PCo1 (", pcoa.eig$Eigval[1], "%)", sep = "")) +
  ylab(paste("PCo2 (", pcoa.eig$Eigval[2], "%)", sep = "")) +
  theme_classic() 



rs.expanded <- pcoa.axes %>% 
  ggplot(aes(Axis.1, Axis.2, color = Treatment)) +
  geom_point() +
  geom_convexhull(aes(fill = Treatment), alpha = 0.7) +
  xlab(paste("PCo1 (", pcoa.eig$Eigval[1], "%)", sep = "")) +
  ylab(paste("PCo2 (", pcoa.eig$Eigval[2], "%)", sep = "")) +
  scale_color_manual(values = trt.pal) +
  scale_fill_manual(values = trt.pal) +
  facet_wrap(Compartment ~ Age, nrow = 2,scales="free_y") +
  theme_bw() +
  theme(legend.position = c(0.9,0.15),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

supp.legend <- plot_grid(get_legend(all.pcoa.p), get_legend(rs.pc.time))

supp.top <- plot_grid(all.pcoa.p + theme(legend.position = "none"),
                      cmp.box + theme(legend.position = "none"), 
                      rs.pc.time + theme(legend.position = "none"),
                      rs.pc.trt +theme(legend.position = "none"),
                      nrow = 1,
                      align = "vh",
                      axis = "l",
                      rel_widths = c(3,3,3,3),
                      labels = c("A","B", "C", "D"),
                      label_size = 20)

supp.bottom <- plot_grid(rs.expanded,
                         nrow = 2,
                         align = "vh",
                         axis = "l",
                         labels = c("E", "F"),
                         label_size = 20)

plot_grid(supp.legend, supp.top, supp.bottom,
          nrow = 3,
          rel_heights = c(1,2,5))


#--微生物组差异分析#--------

#-----微生物组数据DEsep2——矫正效应#-----------
library(DESeq2)
library(biobroom)
library(tidyverse)
library(ggClusterNet)

rel_ab <- function(otu, total = 100) {
  t(t(otu)/colSums(otu)) * 100
}

tidy_otu <- function(otu) {
  as.data.frame(otu) %>%
    mutate(OTU_ID = row.names(otu)) %>%
    tidyr::gather(key = "ID", value = "Count", -OTU_ID)
}



#--这里仅仅选择一个区域的四个时间点做差异分析
# 关注同一个时间点的不同处理做差异分析
otu <- ps %>% 
  subset_samples(Compartment == "RS") %>%
  subset_samples(Time %in% c(1,2,5,7)) %>%
  vegan_otu() %>% 
  t() %>%
  as.data.frame()

head(otu)


map <- ps %>% 
  subset_samples(Compartment == "RS") %>%
  subset_samples(Time %in% c(1,2,5,7)) %>%
  sample_data() 
# map$ID = row.names(map)
head(map)
map = map %>% as.tibble() %>%
  mutate(Group = paste(Treatment2,Time, sep = ".")) %>%
  select(ID,Group,everything())
map$Group = as.factor(map$Group)
map$Group %>% table()
dds<- DESeqDataSetFromMatrix(countData = otu,
                             colData = map,
                             design = ~ Group)

dds <- DESeq(dds)
dds

#--1-时间空间梯度如果太多需要人工指定分组
contrasts =  NULL
for(i in c(1,3,5,7)) {
  
  contrasts[[paste("WC_TRN", i, sep = ".")]] <- c("Group", paste("WC_TRN", i, sep = "."), paste("WC", i, sep = "."))
  contrasts[[paste("D1", i, sep = ".")]] <- c("Group", paste("D1", i, sep = "."), paste("WC", i, sep = "."))
  contrasts[[paste("D2", i, sep = ".")]] <- c("Group", paste("D2", i, sep = "."), paste("WC", i, sep = "."))
  contrasts[[paste("D3", i, sep = ".")]] <- c("Group", paste("D3", i, sep = "."), paste("WC", i, sep = "."))
  
}

# #--2-自动生成两两分组
# Desep_group <- ps %>% 
#   phyloseq::sample_data() %>%
#   .$Group %>%
#   as.factor() %>%
#   levels() %>%
#   as.character()
# aaa = combn(Desep_group,2)
# 
# tem.1 = aaa[1,]
# tem.2 = aaa[2,]
# tem3 = paste("Group",tem.1,tem.2,sep = "=") %>%
#   strsplit("=") 
# names(tem3) = paste(tem.1,tem.2,sep = "_")
# contrasts = tem3


results <- vector(mode = "list")
shrinkFC <- vector(mode = "list")
library("apeglm")
for(i in seq_along(contrasts)[1:2]) {
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

# saveRDS(results.df, "./Data/es_dab_bal.RDS")
# saveRDS(shrinkFC.df, "./shlfc.RDS")
# saveRDS(dds, "./Data/es_dds.RDS")
# rs.fc <- readRDS("./Data/rs_shlfc.RDS")

rs.fc <- shrinkFC.df

#--基于desep2的差异微生物数据取子集
head(rs.fc)
rs.sig <- rs.fc %>% 
  select(-Contrast) %>% 
  # separate(Contrast, c("Treatment", "Time"), sep = "\\.") %>% 
  # filter(Treatment != "WC_TRN") %>% 
  filter(!is.na(padj)) %>% 
  ungroup() %>% 
  mutate(p.adjusted2 = p.adjust(pvalue, method = "fdr")) %>%
  ungroup() %>% 
  filter(p.adjusted2 < 0.05)
# 查看差异微生物的数量
rs.sig$OTU_ID %>% unique() %>% length()



#----可视化差异微生物在全部微生物中分布#------
otu = ps %>% 
  vegan_otu() %>% t() %>%
  as.data.frame()
head(otu)

otu <- rel_ab(otu)
otu.tidy <- tidy_otu(otu)
head(otu.tidy )


map = sample_data(ps) %>% as.tibble()

tax = ps %>% vegan_tax() %>%
  as.data.frame()
tax$OTU_ID = row.names(tax)
head(oc.ab )

oc.ab <- otu.tidy %>% 
  mutate(Count = Count/100) %>% 
  inner_join(map, by = "ID") %>% 
  group_by(Compartment, OTU_ID) %>% 
  summarise(MeanRelAb = mean(Count),
            Occupancy = sum(Count>0)/dim(otu)[2] * 100) %>% 
  # mutate(Compartment = fct_recode(Compartment,
  #                                 "Rhizosphere" = "RS",
  #                                 "Endosphere" = "ES")) %>% 
  # mutate(Compartment = fct_relevel(Compartment, "Rhizosphere")) %>% 
  inner_join(tax, by = "OTU_ID")
int.ids <- rs.sig$OTU_ID %>% unique()
oc.ab.p <- oc.ab  %>% 
  filter(!OTU_ID %in% int.ids) %>% 
  filter(Occupancy > 0) %>% 
  ggplot(aes(MeanRelAb, Occupancy)) +
  geom_point(shape = 1, color = "gray25", alpha = 0.5) +
  geom_point(data = filter(oc.ab, OTU_ID %in% int.ids),
             shape = 16, aes(color = PhyClass), stroke = 1) +
  # geom_point(data = filter(oc.ab, OTU_ID %in% "1037355"), shape = 21, fill = phy.pal[3], color = "black", size = 3, stroke = 1) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_color_manual(name = "",
  #                    values = phy.pal[3],
  #                    labels = "Semi-persistently enriched\nActinobacteria") +
  xlab("Mean relative\nabundance (log10)") +
  ylab("Occupancy (% samples)") +
  facet_grid(Compartment ~ ., scales = "free") +
  theme_classic() +
  theme(text = element_text(size = 15), 
        legend.position = "bottom",
        legend.title = element_blank())

oc.ab.p

#--基于全局关注特定差异微生物#---------
library(ggalluvial)
library(ggClusterNet)
otu.t = ps %>% 
  scale_micro() %>%
  filter_OTU_ps(100) %>%
  vegan_otu() %>% t() %>%
  as.data.frame()
dim(otu.t)
otu.tidy2 <- tidy_otu(otu.t)
head(otu.tidy2 )

tmp <- otu.tidy2 %>% 
  mutate(Count = Count/100) %>% 
  inner_join(map, by = "ID") %>% 
  group_by(Compartment, Age, Treatment, OTU_ID) %>% 
  summarise(MeanAb = mean(Count)) %>% 
  ungroup() %>% 
  filter(MeanAb > 0)

head(tmp)
dim(tmp)

#--选择特定的微生物做展示
prel.ranked.p <- tmp %>% 
  ggplot(aes(Age, MeanAb, alluvium = OTU_ID)) +
  geom_alluvium(aes(fill = OTU_ID == "1037355",
                    color = OTU_ID == "1037355", alpha = OTU_ID == "1037355"), decreasing = F, size = 0.25) +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment), linetype = 3) +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment), linetype = 3) +
  facet_grid(Compartment ~ Treatment, scales = "free", space = "free") +
  scale_fill_manual(name = expression(paste("OTU: 1037355 (", italic("Streptomyces "), "sp.)", sep = "")), values = c("white", phy.pal[3])) +
  scale_color_manual(name = expression(paste("OTU: 1037355 (", italic("Streptomyces "), "sp.)", sep = "")), values = c("gray50", "black")) +
  scale_alpha_manual(name = expression(paste("OTU: 1037355 (", italic("Streptomyces "), "sp.)", sep = "")), values = c(0.5, 1)) +
  xlab("Plant age\n(days)") +
  ylab("Ranked relative abundance") +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
# ggsave("./2.3.pdf",prel.ranked.p)


#--基于差异分析结果进行聚类分析#----------
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

rs.k <- 3
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

#--相对丰度转化
otu <- rel_ab(otu)
# 转化为长数据-根际数据


rs.otu.tidy <- tidy_otu(otu) %>% 
  mutate(Count = Count/100) %>% 
  filter(!is.na(Count)) 
head(rs.otu.tidy )
head(map)

# colnames(map)[1] = "ID"
rs.zs.tidy <- rs.otu.tidy %>%
  inner_join(map, by = "ID") %>%
  # inner_join(rs.ord,by = "OTU_ID") %>%
  filter(OTU_ID %in% rs.sig$OTU_ID) %>%
  group_by(OTU_ID) %>%
  mutate(zscore = (Count - mean(Count))/sd(Count)) %>%
  group_by(ID,Compartment,Age, OTU_ID) %>%
  summarise(MeanZS = mean(zscore)) %>%
  ungroup()


rs.hm.p <- rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  mutate(MeanZS2 = ifelse(MeanZS > 2,2,MeanZS)) %>%
  # mutate(Treatment = fct_relevel(Treatment, "WC")) %>%
  ggplot(aes(ID, reorder(OTU_ID, order), fill = MeanZS2)) +
  geom_tile() +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment*10), linetype = 3, color = "white") +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment*10), linetype = 3, color = "white") +
  scale_fill_viridis_c(name = "Mean Rel. Abund.\n(z-score)") +
  ylab("Differentially Abundant OTU") +
  xlab("Group") +
  # scale_fill_gradientn(name = "log2FC",
  #                      colors = RColorBrewer::brewer.pal(11, "BrBG")) +
  facet_grid(Cluster ~ Age, scales = "free", space = "free") +
  theme_void() +
  theme(text = element_text(size = 8),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title.align = 1)

rs.hm.p



tem = rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  filter(Cluster == "C1" )
tem$OTU_ID %>% unique()


#-挑选指定的模块展示微生物组差异#--------

rs.hm.p <- rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  filter(Cluster == "C1" ) %>%
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>%
  # mutate(Treatment = fct_relevel(Treatment, "WC")) %>%
  ggplot(aes(ID, reorder(OTU_ID, order), fill = MeanZS)) +
  geom_tile() +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment*10), linetype = 3, color = "white") +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment*10), linetype = 3, color = "white") +
  scale_fill_viridis_c(name = "Mean Rel. Abund.\n(z-score)") +
  ylab("Differentially Abundant OTU") +
  xlab("Group") +
  # scale_fill_gradientn(name = "log2FC",
  #                      colors = RColorBrewer::brewer.pal(11, "BrBG")) +
  facet_grid(Cluster ~ Age, scales = "free", space = "free") +
  theme_classic() +
  theme(text = element_text(size = 10),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.title.align = 1)

rs.hm.p

# 模块平均z值展示，查看模块整体微生物在全部样本中的变化
p = rs.zs.tidy %>%
  inner_join(rs.ord, by = "OTU_ID") %>%
  group_by(Cluster,Age) %>%
  summarise(mean(MeanZS)) %>%
  ggplot(aes(x = Age,y = `mean(MeanZS)`, group = Cluster)) +
  geom_point(pch = 21,alpha = .3) +geom_line(alpha = .5) +
  theme_classic()
p

#--不同分类等级展示差异微生物#------
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
  filter(Cluster == "C1") %>%
  filter(!Genus %in% c("","1.00","0.67","3","g__")) %>%
  group_by(Genus,Age) %>%
  summarise(mean(MeanZS)) %>%
  inner_join(tem,by = c("Genus" = "ID")) %>%
  arrange(Age,desc(mean)) 
head(tem3)

tem3 %>%
  ggplot(aes(x= Age, y = Genus, fill = `mean(MeanZS)`)) +
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

#-组合图表

plot_grid(oc.ab.p, prel.ranked.p, 
          align = "h",
          axis = "tb",
          rel_widths = c(1,3),
          labels = c("A", "B"),
          label_size = 20)


#-4-基于微生物大数据整合特征微生物丰度调查#---------


tax = ps_meta %>% vegan_tax() %>%
  as.data.frame()
head(tax)
tax$OTU_ID = row.names(tax)

otu = ps_meta %>% 
  vegan_otu() %>% t() %>%
  as.data.frame()
head(otu)

tidy_otu <- function(otu) {
  as.data.frame(otu) %>%
    mutate(OTU_ID = row.names(otu)) %>%
    tidyr::gather(key = "ID", value = "Count", -OTU_ID)
}
otu <- otu[!row.names(otu) %in% organelle.id,]
tidy.otu <- otu %>% 
  tidy_otu()

meta.map = sample_data(ps_meta) %>% as.tibble()
head(meta.map )

head(tidy.otu)

#---绘制特定菌在不同研究中的丰度箱线图#-------
drought.p <- tidy.otu %>% 
  inner_join(meta.map , by = "ID") %>% 
  group_by(ID) %>% 
  mutate(RelAb = Count/sum(Count)) %>% 
  ungroup() %>% 
  filter(OTU_ID == "1037355") %>% 
  filter(Compartment %in% c("Rhizosphere", "Endosphere")) %>% 
  mutate(Compartment = fct_relevel(Compartment, "Rhizosphere", "Endosphere")) %>% 
  filter(Compartment != "Bulk Soil") %>% 
  mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  ggplot() +
  geom_boxplot(aes(Treatment, RelAb, color = Treatment), outlier.shape = NA) +
  geom_jitter(aes(Treatment, RelAb, color = Treatment)) +
  # geom_segment(data = wald, aes(x = 1, xend = 2, y = 0.25, yend = 0.25)) +
  # geom_text(data = wald, aes(x = 1.5, y = 0.3, label = "***", hjust = 0.5)) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = trt.pal[c(1,3)]) +
  ylab("Relative abundance (log10)") +
  facet_grid(. ~ Compartment + Soil) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "none",
        panel.border = element_blank())

drought.p


n.samples <- meta.map %>% 
  group_by(Experiment, Compartment, Soil) %>% 
  summarise(nSamples = n())

oc.ab <- tidy.otu %>% 
  group_by(ID) %>% 
  mutate(RelAb = Count/sum(Count)) %>% 
  inner_join(meta.map) %>% 
  filter(Compartment %in% c("Rhizosphere", "Endosphere")) %>% 
  group_by(Experiment, Compartment, Soil, OTU_ID) %>% 
  summarise(MeanRelAb = mean(RelAb),
            Occupancy = sum(Count > 0)) %>% 
  ungroup() %>% 
  inner_join(n.samples, by = c("Experiment", "Compartment", "Soil")) %>% 
  mutate(OccupancyPct = Occupancy/nSamples)

oc.ab.p <- oc.ab %>% 
  filter(Occupancy > 0) %>% 
  ggplot(aes(MeanRelAb, OccupancyPct)) +
  geom_point(shape = 1, color = "gray25", alpha = 0.5) +
  geom_point(data = filter(oc.ab, OTU_ID == "1037355"), shape = 21, fill = phy.pal[3], color = "black", size = 3, stroke = 1) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab("Mean relative\nabundance (log10)") +
  ylab("Occupancy (% samples)") +
  facet_grid(Compartment ~ Soil) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        legend.position = "bottom",
        legend.title = element_blank())

oc.ab.p


oc.ab %>% filter(OTU_ID == "1037355")


cowplot::plot_grid(oc.ab.p, drought.p, nrow = 2,
                   labels = c("A", "B"), label_size = 20)






# 5-筛菌微生物挑选和比对#-----

### 分离株和OTU进行比对

#这里的逻辑我帮大家梳理一下。这里首先用的97的聚类得到的OTU，我们的目标OTU是OTU 1037355，
# 但是所有被分类为这个OTU的序列并不是完全一致的，所有这里选择了最能代表这个OTU的序列作为最可能的序列和分菌的序列比对。
#这里候选了10中相似的链霉菌，有物种之和目标菌株有一个碱基的差异，根际这个分菌的来源，最终选择了177号菌株。

iso.pal <- wesanderson::wes_palette(name = "Zissou1", n = 5, type = "discrete")[c(1,3,5)]


seq.p <- seqs %>%
  select(-(OTU.major.sequence:SLBN.111)) %>%
  gather("Isolate","match",OTU.major.seq:`SLBN-111`) %>%
  mutate(Isolate = fct_recode(Isolate,
                              "OTU: 1037355\nprevalent seq" = "OTU.major.seq")) %>% 
  mutate(Isolate = fct_relevel(Isolate,
                               "OTU: 1037355\nprevalent seq",
                               "SLBN-197",
                               "SLBN-175",
                               "SLBN-134",
                               "SLBN-193",
                               "SLBN-191",
                               "SLBN-186",
                               "SLBN-177",
                               "SLBN-162",
                               "SLBN-161",
                               "SLBN-111")) %>%
  ggplot(.,aes(x=Position,y=1,fill=match))+
  geom_tile()+
  facet_grid(Isolate~.)+
  scale_fill_manual(values=c("black","grey80"))+
  theme_minimal()+
  theme(text = element_text(size = 15),
        strip.text.y = element_text(angle = 0,hjust=0, size = 10),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")



# 绘制用于表型实验的菌株的丰度占用信息。

oc.ab <- otu.tidy %>% 
  #filter(Count > 0) %>% 
  mutate(Count = Count/100) %>% 
  inner_join(map, by = "ID") %>% 
  group_by(Compartment, OTU_ID) %>% 
  summarise(MeanRelAb = mean(Count),
            Occupancy = sum(Count>0)/208 * 100) %>% 
  mutate(Compartment = fct_recode(Compartment,
                                  "Rhizosphere" = "RS",
                                  "Endosphere" = "ES")) %>% 
  mutate(Compartment = fct_relevel(Compartment, "Rhizosphere"))

oc.ab.p <- oc.ab  %>% 
  ggplot(aes(MeanRelAb, Occupancy)) +
  geom_point(shape = 1, color = "gray25", alpha = 0.5) +
  geom_point(data = filter(oc.ab, OTU_ID %in% c("1037355", "1108350")), aes(fill = OTU_ID), shape = 21, color = "white", size = 3) +
  scale_x_log10() +
  scale_fill_manual(name = "OTU",
                    values = iso.pal[3:2]) +
  xlab("Mean relative abundance\n(log10)") +
  ylab("Occupancy\n(% samples)") +
  facet_wrap(~ Compartment) +
  theme_classic() +
  theme(text = element_text(size = 15))

plot_grid(seq.p, oc.ab.p,
          ncol = 1,
          labels = c("A", "B"),
          label_size = 20)


# 6-基于筛菌验证实验分析流程#-----


# source("./General/rmb_functions.R")
library(emmeans)
library(cowplot)
library(tidyverse)







tax = ps_verify %>% vegan_tax() %>%
  as.data.frame()
head(tax)
tax$OTU_ID = row.names(tax)

iso.otu = ps_verify %>% 
  vegan_otu() %>% t() %>%
  as.data.frame()
head(iso.otu)

tidy_otu <- function(otu) {
  as.data.frame(otu) %>%
    mutate(OTU_ID = row.names(otu)) %>%
    tidyr::gather(key = "ID", value = "Count", -OTU_ID)
}


tidy.otu <- iso.otu %>% 
  tidy_otu()

iso.map = sample_data(ps_meta) %>% as.tibble()
head(iso.map)



iso.pal <- wesanderson::wes_palette(name = "Zissou1", n = 5, type = "discrete")[c(1,3,5)]
trt.colors <- RColorBrewer::brewer.pal(11, "BrBG")[c(10,3)]



### 验证实验设计路线展示
diagram.df <- data.frame(Isolate = c(rep("Mock control",6),rep("Streptomyces sp. SLBN-177",6),rep("Microbacterium sp. SLBN-111",6)),
                         Treatment = rep(c(rep("WC",3),rep("DS",3)),3),
                         PrevAge= rep(c(0,10,24),6),
                         Age=rep(c(10,24,31),6),
                         WaterStatus = rep(c("W","W","W","W","S","W"),3),
                         Tube = rep(c("Closed", "Open", "Open"),3)) %>% 
  mutate(Isolate = fct_relevel(Isolate, "Mock control"),
         Treatment = fct_relevel(Treatment, "WC"))
xpt.iso.p <- diagram.df %>% 
  filter(WaterStatus != "S") %>% 
  ggplot(aes(PrevAge, Isolate)) +
  geom_rect(data = diagram.df, aes(xmin = PrevAge, xmax = Age, ymin = 0.5, ymax = 3.5, fill = Tube)) +
  geom_segment(aes(xend = Age, yend = Isolate, color = Isolate), size = 3, linetype = 1) +
  geom_segment(data = filter(diagram.df, WaterStatus == "S"), aes(xend = Age, yend = Isolate, color = Isolate), linetype = 2) +
  scale_color_manual(values = iso.pal, 
                     guide = F) +
  scale_fill_manual(name = "Growth\nsystem",
                    values = c("gray75", "gray25")) +
  scale_y_discrete(limits = c("Streptomyces sp. SLBN-177", "Microbacterium sp. SLBN-111", "Mock control"),
                   labels = c(expression(paste(italic("Streptomyces "), "sp. SLBN-177", sep = "")),
                              expression(paste(italic("Microbacterium "), "sp. SLBN-111", sep = "")),
                              "Mock control")) +
  scale_x_continuous(breaks = seq(0,31,10),
                     position = "bottom") +
  facet_grid(Treatment ~ .) +
  theme_classic() +
  xlab("Plant age (days)") +
  theme(text = element_text(size = 12),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right")
xpt.iso.p



#验证实验表型排序分析#--------
phenotype.mtx <- phenotypes %>%
  mutate(key = paste(isolate,trt,tray,row,col,sep=".")) %>%
  select(key,shoot_len,mass_root_length,shoot_wt,root_wt,n_leaves,n_roots) %>%
  column_to_rownames(.,var = "key") %>%
  na.omit() %>%
  as.matrix()
pca <- prcomp(phenotype.mtx,center = TRUE, scale. = TRUE)

pca.out <- as.data.frame(pca$x) %>%
  add_rownames("key")
pca.rotation <- as.data.frame(pca$rotation) %>%
  add_rownames("trait") %>%
  mutate(trait = case_when(trait=="shoot_len" ~ "Shoot length",
                           trait=="mass_root_length" ~ "Root length",
                           trait=="shoot_wt" ~ "Shoot weight",
                           trait=="root_wt" ~ "Root weight",
                           trait=="n_leaves" ~ "# Leaves",
                           trait=="n_roots" ~ "# Roots")) %>%
  mutate(PC1.1=PC1-.025)
map.pca <- phenotypes %>%
  mutate(key = paste(isolate,trt,tray,row,col,sep=".")) %>%
  select(isolate,trt,tray,row,col,key) %>%
  inner_join(pca.out,by="key")
iso.pca.p <- map.pca %>% 
  mutate(isolate = fct_recode(isolate,
                              "Mock Ctrl" = "ctr",
                              "SLBN-177" = "177",
                              "SLBN-111" = "111")) %>% 
  mutate(isolate = fct_relevel(isolate, "Mock Ctrl", "SLBN-111")) %>% 
  ggplot()+
  geom_point(aes(x=PC1,y=PC2,color=isolate,shape=trt),size=2)+
  geom_label(data=pca.rotation, aes(x=PC1.1*3, y=PC2*3, label=trait), size = 4, hjust=1, alpha = 0.5)+
  geom_segment(data=pca.rotation,aes(x=0, y=0, xend=PC1*3, yend=PC2*3),arrow=arrow(length=unit(.25,"picas")))+
  scale_color_manual(name = "Microbial\ntreatment", values = iso.pal)+
  scale_shape_manual(name = "Watering\ntreatment", values = c(16,1))+
  guides(color = guide_legend(nrow = 3),
         shape = guide_legend(nrow = 2)) +
  theme_bw()+
  labs(x="PC1 (58.68%)",y="PC2 (15.14%)")+
  theme_classic() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
iso.pca.p

### 基于混合线性模型差异统计分析


mrl.model <- lmerTest::lmer(mass_root_length ~ (1|tray) + isolate * trt, phenotypes)
sl.model <- lmerTest::lmer(shoot_len ~ (1|tray) + isolate * trt, phenotypes)
sw.model <- lmerTest::lmer(shoot_wt ~ (1|tray) + isolate * trt, phenotypes)
rw.model <- lmerTest::lmer(root_wt ~ (1|tray) + isolate * trt, phenotypes)
nr.model <- lmerTest::lmer(n_roots ~ (1|tray) + isolate * trt, phenotypes)
anova.res <- rbind(anova(mrl.model) %>% as.data.frame() %>% add_rownames("term") %>% mutate(phenotype="Root Length"),
                   anova(sl.model) %>% as.data.frame() %>% add_rownames("term") %>% mutate(phenotype="Shoot Length"),
                   anova(rw.model) %>% as.data.frame() %>% add_rownames("term") %>% mutate(phenotype="Root Weight"), 
                   anova(sw.model) %>% as.data.frame() %>% add_rownames("term") %>% mutate(phenotype="Shoot Weight"),
                   anova(nr.model) %>% as.data.frame() %>% add_rownames("term") %>% mutate(phenotype="# Roots")) %>%
  select(phenotype,term:`Pr(>F)`) %>%
  mutate(p.adjusted=p.adjust(`Pr(>F)`))
anova.res
# write.table(anova.res, "./Tables/pheno_anova.tsv", sep = "\t", quote = F, row.names = F)


#验证实验挑选出来主要的差异表型展示#-------


mrl.contrasts <- phenotypes %>%
  mutate(group=paste(isolate,trt,sep=".")) %>%
  rstatix::games_howell_test(mass_root_length ~ group)


p.values <- as.vector(mrl.contrasts$p.adj)
names(p.values) <- paste(mrl.contrasts$group1,mrl.contrasts$group2,sep="-")
multcompView::multcompLetters(p.values)

group.letters <- data.frame(isolate=c("Mock ctrl","Mock ctrl","SLBN-177","SLBN-177","SLBN-111","SLBN-111"),
                            trt=rep(c("WC","DS"),3),
                            position=rep(17,6),
                            group.lets=c("ab","a","c","bc","a","a"))
root.length.p <- phenotypes %>%
  mutate(isolate = fct_recode(isolate,
                              "Mock ctrl" = "ctr",
                              "SLBN-177" = "177",
                              "SLBN-111" = "111")) %>%
  mutate(isolate = fct_relevel(isolate, "Mock ctrl", "SLBN-111")) %>% 
  ggplot(aes(x=trt,y=mass_root_length, color = trt)) +
  geom_boxplot(size = 1) +
  geom_jitter(aes(shape = trt), width = .25, size = 2, alpha=1) +
  geom_text(data = group.letters, aes(y = position, label = group.lets), size = 5,color = "black") +
  facet_grid(. ~ isolate) +
  scale_color_manual(name = "Watering\ntreatment", values = trt.colors) +
  scale_shape_manual(name = "Watering\ntreatment", values = c(16,1))+
  guides(color = guide_legend(nrow = 2)) +
  xlab("") +
  ylab("Root length (cm)") +
  theme_classic() +
  theme(text = element_text(size = 12),
        legend.position = "bottom")
root.length.p


bottom <- plot_grid(iso.pca.p, 
                    root.length.p, 
                    ncol = 2,
                    align = "h",
                    axis = "b",
                    rel_widths = c(3,3),
                    labels = c("B","C"),
                    label_size = 20)
bottom
plot_grid(xpt.iso.p, bottom, 
          nrow = 2,
          rel_heights = c(1,2),
          labels = "A",
          label_size = 20)


# 除此之外其他几种指标差异#-----------


pheno.box <- phenotypes %>%
  select(tray,row,col,isolate,trt,shoot_len,shoot_wt,root_wt,n_roots) %>%
  gather("measurement","value",shoot_len:n_roots) %>%
  mutate(measurement = case_when(measurement == "shoot_len" ~ "Shoot Length (cm)",
                                 measurement == "shoot_wt" ~"Shoot Weight (g)",
                                 measurement == "root_wt" ~"Root Weight (g)",
                                 measurement == "n_roots" ~ "# Roots")) %>%
  mutate(isolate = fct_recode(isolate,
                              "Mock Ctrl" = "ctr",
                              "SLBN-177" = "177",
                              "SLBN-111" = "111")) %>% 
  mutate(isolate = fct_relevel(isolate, "Mock Ctrl", "SLBN-111")) %>% 
  ggplot(.,aes(x=isolate,y=value,color=trt)) +
  geom_boxplot(size = 1) +
  facet_wrap(~measurement,nrow=1,scales="free")+
  theme_minimal()+
  scale_color_manual(name = "Watering\nTreatment", values = trt.colors)+
  xlab("Microbial Treatment") +
  ylab("Value") +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
pheno.box




#验证微生物群落中目标微生物丰度#-----
iso.map = sample_data(ps_verify) %>% as.tibble()


colonization.p <- iso.otu %>% 
  rel_ab() %>% 
  tidy_otu() %>% 
  left_join(tax, by = "OTU_ID") %>% 
  mutate(Type = case_when(Assignment != "Microbial" ~ "Organelle",
                          OTU_ID == "1108350" ~ "OTU 1108350\n(SLBN-111)",
                          OTU_ID == "1037355" ~ "OTU 1037355\n(SLBN-177)",
                          TRUE ~ "Other OTU")) %>% 
  mutate(Type = fct_relevel(Type, "OTU 1037355\n(SLBN-177)", after = Inf)) %>% 
  group_by(ID, Type) %>% 
  summarise(TotAb = sum(Count)) %>% 
  ungroup() %>% 
  inner_join(iso.map, by = "ID") %>% 
  mutate(isolate = fct_recode(isolate,
                              "Mock Ctrl" = "ctr",
                              "SLBN-177" = "177",
                              "SLBN-111" = "111")) %>%
  mutate(isolate = fct_relevel(isolate, "Mock Ctrl", "SLBN-111"),
         trt = fct_relevel(trt, "WC")) %>% 
  ggplot(aes(ID, TotAb, fill = Type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("gray25", "gray75", iso.pal[c(2,3)])) +
  xlab("Sample") +
  ylab("Relative Abundance") +
  facet_grid(. ~ isolate + trt, scales = "free") +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")
colonization.p



plot_grid(pheno.box, colonization.p, 
          nrow = 2,
          align = "v",
          axis = "lr",
          labels = c("A", "B"),
          label_size = 20)


# 7 时间空间机器学习#-----

# source("./General/parameters.R")
# source("./General/rmb_functions.R")
library(randomForest)
library(cowplot)
library(tidyverse)


# ps = readRDS("./ps.rds")
log_norm <- function(otu) {
  log(otu + 1)
}


tax = ps %>% vegan_tax() %>%
  as.data.frame()
head(tax)
tax$OTU_ID = row.names(tax)

otu = ps %>% 
  vegan_otu() %>% t() %>%
  as.data.frame()
head(otu)

map = sample_data(ps) %>% as.tibble()
head(map )


## G相对丰度转化和log转化
otu <- rel_ab(otu)
otu <- log_norm(otu)
# tax <- readRDS("./Data/tax.RDS")



#机器学习模型训练集
# 选择WC_TRN作为模型训练数据，这是正常水分管理的样本。这里我们要理解一个逻辑，作者使用了群落发育成熟度的概念，这个概念开始于18年那篇文章，说的是干旱让微生物群落发育缓慢，不够成熟。这里正常水分管理的样本用于构建机器学习模型，作为成熟的微生物群落。随后使用这个群落对干旱的样本进行预测，如果准确度低，则证明后续微生物群落发育发生了变化，也就是这里的发育缓慢，不够成熟。


train.id <- map %>% 
  filter(Treatment2 == "WC_TRN") %>% 
  .$ID

# 区分训练和预测样本
map <- mutate(map, 
              Set = ifelse(ID %in% train.id, "Train", "Test"))

# 合并otu表格和分组文件
rf.data <-  as.data.frame(t(otu)) %>%  
  mutate(ID = rownames(.)) %>% 
  inner_join(select(map, ID, Compartment, Set, Treatment, Age), by = "ID")

# The predictors need to be formatted as a data frame or matrix
get_predictors <- function(x){
  select(x, -(ID:Age))
} 

# The response variable needs to be formatted as a vector
get_response <- function(x) {
  select(x, Age) %>% 
    .$Age
}  

# 整理提取训练数据
# train.set[[1,2]][[1]] %>% colnames() %>% tail()

train.set <- rf.data %>% 
  filter(Set == "Train") %>% 
  group_by(Compartment) %>% 
  nest() %>% 
  mutate(otu = map(data, get_predictors), # 去除ID和Age编号
         age = map(data, get_response)) # 选择Age列作为因变量

# 整理全部数据
whole.set <- rf.data %>% 
  group_by(Set, Compartment, Treatment) %>% 
  nest() %>% 
  mutate(otu = map(data, get_predictors),
         age = map(data, get_response))


### 交叉检验-基于时间序列进行回归

# 这里将错误率的的最低点标记一下，很容易就可以看出来。这位分析工作者细节满满。


get_cv <- function(predictors, response) {
  cv <- rfcv(predictors, response, cv.fold = 10, log = T)
  data.frame(nOTU = cv$n.var,
             Error = cv$error.cv) %>% 
    mutate(minError = Error == min(Error))
}

get_nOTU <- function(cv) {
  filter(cv, minError) %>% 
    .$nOTU
}

#-基于tidy进行交叉检验 基于时间
# train.set [[1,4]][[1]]
train.set <- train.set %>% 
  mutate(cv = map2(otu, age, get_cv),
         nOTU = map(cv, get_nOTU)) 

# 交叉检验结果可视化#--------
cv.plot <- train.set %>% 
  unnest(cv, nOTU) %>% 
  mutate(Compartment = fct_relevel(Compartment, "RS")) %>% 
  mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>% 
  ggplot(aes(nOTU, Error, color = Compartment)) +
  geom_line(size =1) +
  geom_point(aes(alpha = minError), size = 10) +
  geom_text(aes(label = nOTU, alpha = minError), size = 6, color = "black") +
  scale_color_manual(values = rev(cmp.pal)) +
  scale_alpha_manual(values = c(0,1)) +
  guides(alpha = F) +
  scale_x_log10() +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = c(0.5, 0.8),
        legend.title = element_blank()) 

cv.plot


#随机森林和重要变量提取

# 这里按照根际和根内分类建模。随机森林回归用年龄做因变量。

# 随机森林训练#----------

get_imp <- function(predictors, response){
  rf <- randomForest(predictors, response, importance = T)
  imp <- data.frame(OTU_ID = row.names(rf$importance),
                    rf$importance)
  names(imp)[2] <- "PercIncMSE"
  as.tibble(imp)
} 
# 提取重要变量
get_imp_otus <- function(data, imp, nOTU){
  id <- imp %>% 
    top_n(nOTU, PercIncMSE) %>% 
    .$OTU_ID
  data[,colnames(data) %in% id]
}
# 提取重要变量的分类信息
get_imp_tax <- function(imp, nOTU){
  imp %>% 
    top_n(nOTU, PercIncMSE) %>% 
    inner_join(tax, by = "OTU_ID")
}

#基于训练数据进行模型学习构建，分别对根际和根内做两个模型
train.set <- train.set %>% 
  mutate(imp = map2(otu, age, get_imp),
         imp.otu = pmap(list(data, imp, nOTU), get_imp_otus),
         imp.tax = map2(imp, nOTU, get_imp_tax))

rf.ranks <- train.set %>% 
  unnest(imp.tax) %>% 
  mutate(cmpOTU = paste(Compartment, OTU_ID),
         Rank = rank(PercIncMSE)) 

# 提取重要变量的id 可视化#------
rf.ranks %>% 
  ggplot(aes(reorder(cmpOTU, Rank), PercIncMSE, fill = PhyClass2)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  scale_y_log10() +
  coord_flip() +
  scale_fill_manual(values = phy.pal) +
  facet_wrap(~Compartment, scales = "free") +
  theme_light() +
  theme(text = element_text(size = 20),
        axis.text.y = element_blank())



otu.rf <- train.set %>% 
  unnest(imp.tax) %>% 
  mutate(cmpOTU = paste(Compartment, OTU_ID),
         Rank = rank(PercIncMSE)) %>% 
  mutate(Compartment = ifelse(Compartment == "RS", "Rhizosphere", "Endosphere"))

rs.id <- filter(otu.rf, Compartment == "Rhizosphere")$OTU_ID
es.id <- filter(otu.rf, Compartment == "Endosphere")$OTU_ID

tidy_otu <- function(otu) {
  as.data.frame(otu) %>%
    mutate(OTU_ID = row.names(otu)) %>%
    tidyr::gather(key = "ID", value = "Count", -OTU_ID)
}


otu.tidy <- tidy_otu(otu)

# 提取聚类结果的顺序
get_order <- function(df, var) {
  mtrx <- df %>% 
    mutate(Group = paste(Compartment, Treatment, Time, sep = ".")) %>% 
    select(Group, var, OTU_ID) %>% 
    spread(key = Group, value = var)
  mtrx <- as.data.frame(mtrx)
  rownames(mtrx) <- mtrx$OTU_ID
  mtrx <- mtrx[,-1]
  dist.tmp <- dist(as.matrix(mtrx))
  dist.clust <- hclust(dist.tmp, method = "average")
  ord.names <- dist.clust$labels[dist.clust$order]
  clust.ord <- data.frame(OTU_ID = ord.names, order = 1:length(ord.names))
  df <- inner_join(df, clust.ord, by = "OTU_ID")
  df$order
}


#定义聚类函数
get_cluster <- function(df, var, nclust) {
  mtrx <- df %>% 
    mutate(Group = paste(Compartment, Treatment, Time, sep = ".")) %>% 
    select(Group, var, OTU_ID) %>% 
    spread(key = Group, value = var)
  mtrx <- as.data.frame(mtrx)
  rownames(mtrx) <- mtrx$OTU_ID
  mtrx <- mtrx[,-1]
  dist.tmp <- dist(as.matrix(mtrx))
  dist.clust <- hclust(dist.tmp, method = "average")
  ord.names <- dist.clust$labels[dist.clust$order]
  clust.ord <- data.frame(OTU_ID = ord.names, order = 1:length(ord.names))
  clust.cut <- cutree(dist.clust[c(1,2,4)], k = nclust)
  clust.ord$Cluster <- as.factor(clust.cut[clust.ord$OTU_ID])
  df <- inner_join(df, clust.ord, by = "OTU_ID")
  df$Cluster
}

### 根内聚类等
head(otu.tidy)
es.master.tidy <- otu.tidy %>% 
  inner_join(select(map, ID, Time, Compartment, Treatment), by = "ID") %>% 
  filter(Compartment == "ES") %>% 
  group_by(Compartment, Treatment, OTU_ID) %>% 
  mutate(zscore = (Count - mean(Count))/sd(Count)) %>% 
  group_by(Compartment, Treatment, Time, OTU_ID) %>% 
  summarise(MeanZS = mean(zscore)) %>% 
  ungroup()

es.order <- es.master.tidy %>% 
  filter(OTU_ID %in% es.id) %>% 
  filter(Treatment == "WC") %>% 
  mutate(ZS_order = get_order(., "MeanZS"),
         ZS_cluster = get_cluster(., "MeanZS", 3)) 

#定义聚类的微生物的意义
# %>% 
#   mutate(Trend = case_when(
#     ZS_cluster == 1 ~ "Complex",
#     ZS_cluster == 2 ~ "Early Colonizer",
#     ZS_cluster == 3 ~ "Late Colonizer"
#   ))

es.plot.df <- es.master.tidy %>% 
  inner_join(select(es.order, OTU_ID, ZS_order, ZS_cluster), by = "OTU_ID") %>% 
  mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  inner_join(tax, by = "OTU_ID") 

# %>% 
#   mutate(EndTreatment = case_when(
#     Treatment == "WC" ~ 4.5,
#     Treatment == "DS10" ~ 5.5,
#     Treatment == "DS20" ~ 6.5,
#     Treatment == "DS30" ~ 7.5,
#   ))



##区分空间聚类等
rs.master.tidy <- otu.tidy %>% 
  inner_join(select(map, ID, Time, Compartment, Treatment), by = "ID") %>% 
  filter(Compartment == "RS") %>% 
  group_by(Compartment, Treatment, OTU_ID) %>% 
  mutate(zscore = (Count - mean(Count))/sd(Count)) %>% 
  group_by(Compartment, Treatment, Time, OTU_ID) %>% 
  summarise(MeanZS = mean(zscore)) %>% 
  ungroup()

rs.order <- rs.master.tidy %>% 
  filter(OTU_ID %in% rs.id) %>% 
  filter(Treatment == "WC") %>% 
  mutate(ZS_order = get_order(., "MeanZS"),
         ZS_cluster = get_cluster(., "MeanZS", 3)) 

# %>% 
#   mutate(Trend = case_when(
#     ZS_cluster == 3 ~ "Complex",
#     ZS_cluster == 1 ~ "Early Colonizer",
#     ZS_cluster == 2 ~ "Late Colonizer"
#   ))

rs.plot.df <- rs.master.tidy %>% 
  inner_join(select(rs.order, OTU_ID, ZS_order, ZS_cluster), by = "OTU_ID") %>% 
  mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  inner_join(tax, by = "OTU_ID") 

# %>% 
#   mutate(EndTreatment = case_when(
#     Treatment == "WC" ~ 4.5,
#     Treatment == "DS10" ~ 5.5,
#     Treatment == "DS20" ~ 6.5,
#     Treatment == "DS30" ~ 7.5,
#   ))

### 提取聚类OTU，分类可视化

trt.lines <- data.frame(Treatment = c("DS1", "DS2", "DS3"),
                        Treatment2 = c("DS1", "DS2", "DS3"),
                        Contrast = c("WC vs DS1", "WC vs DS2", "WC vs DS3"),
                        IniTreatment = c(4.5, 4.5, 4.5),
                        EndTreatment = c(5.5, 6.5, 7.5),
                        IniTreatment2 = c(41,41,41),
                        EndTreatment2 = c(52,62,74))

all.plot.df <- rbind(rs.plot.df, es.plot.df) %>% 
  mutate(OTU_ID = paste(Compartment, OTU_ID))

all.ord <- rbind(group_by(rs.order, Compartment, OTU_ID, ZS_cluster) %>% dplyr::count(),
                 group_by(es.order, Compartment, OTU_ID, ZS_cluster) %>% dplyr::count())

# all.ord$Trend <- as.factor(all.ord$Trend)
#--可视化多模型随机森林结果#--------
trend.p <- all.plot.df %>% 
  # mutate(Treatment = fct_recode(Treatment, WC = "WC", DS1 = "D1", DS2 = "D2", DS3 = "D3")) %>% 
  # mutate(Compartment = fct_relevel(Compartment, "RS")) %>% 
  # mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>% 
  # mutate(Trend = fct_relevel(Trend, "Complex", after = Inf)) %>% 
  ggplot(aes(Time * 10, reorder(OTU_ID, ZS_order), fill = MeanZS)) +
  geom_tile() +
  geom_point(aes(x = -1, color = ZS_cluster)) +
  # geom_vline(data = trt.lines, aes(xintercept = (IniTreatment * 10)), linetype = 3) +
  # geom_vline(data = trt.lines, aes(xintercept = (EndTreatment* 10)), linetype = 3) +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Greys"),
                       name = "Relative Abundance\n(z-score)") +
  scale_color_manual(values = c("dodgerblue", "gold", "grey50")) +
  ylab("Age Discriminant OTU") +
  xlab("Plant Age (days)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = "bottom",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(Compartment ~ Treatment, scales = "free", space = "free") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

trend.p


# 补充表格生成，用于发表文章附件
rf.otus <- rbind(rs.plot.df, es.plot.df) %>% 
  group_by(Compartment, OTU_ID,  Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count() %>% 
  select(-n) %>% 
  ungroup()

#评估时间序列特征 OTU 和干旱响应 OTU 之间的共有微生物数量#--------

overlap.plot <- all.plot.df  %>% 
  separate(OTU_ID, c("tmp", "OTU_ID")) %>% 
  select(OTU_ID, Compartment,ZS_cluster) %>% 
  mutate(RFTrend = ZS_cluster) %>%
  # select(-Trend) %>% 
  group_by(OTU_ID, Compartment, RFTrend) %>% 
  dplyr::count() %>% 
  dplyr::select(-n) %>% 
  left_join(all.clust.summary, by = c("Compartment", "OTU_ID")) %>% 
  ungroup() %>% 
  # mutate(Compartment = fct_relevel(Compartment, "RS")) %>% 
  # mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>% 
  mutate(Trend = ifelse(is.na(Cluster), "None", Cluster)) %>% 
  mutate(Trend = fct_relevel(Trend, "None", after = Inf)) %>% 
  # mutate(RFTrend = fct_relevel(RFTrend, "Early Colonizer", after = Inf)) %>% 
  ggplot(aes(RFTrend, fill = Trend)) +
  geom_bar() +
  facet_grid(Compartment ~.) +
  ylab("nOTU") +
  xlab("") +
  scale_fill_manual(name = "Drought Module",
                    values = c(RColorBrewer::brewer.pal(5,"Dark2")[-2], "gray69")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "right") +
  coord_flip()

overlap.plot



### 随机森林特征几大类的微生物注释信息可视化#------


rf.tax <- all.plot.df %>% 
  # mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>% 
  # mutate(Compartment = fct_relevel(Compartment, "Rhizosphere")) %>% 
  # mutate(Trend = fct_relevel(Trend, "Complex", after = Inf)) %>% 
  # mutate(Trend = fct_recode(Trend, Early = "Early colonizer", Late = "Late colonizer")) %>% 
  group_by(Compartment,ZS_cluster,OTU_ID) %>% 
  dplyr::count() %>% 
  separate(OTU_ID, c("tmp", "OTU_ID")) %>% 
  select(-n, -tmp) %>% 
  inner_join(tax, by = "OTU_ID") %>% 
  group_by(Compartment, ZS_cluster, PhyClass2, Phylum, Class, Order) %>% 
  dplyr::count() %>%
  ggplot(aes(ZS_cluster, paste(Phylum, Class, Order, sep = " / "), fill = PhyClass2, label = n)) +
  geom_tile(color = "white", size = 1) +
  geom_text(color = "black", size = 3) +
  scale_fill_manual(name = "Taxon",
                    values = phy.pal,
                    limits = levels(tax$PhyClass2)) +
  ylab("Phylum / Class / Order") +
  facet_grid(.~Compartment) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "none",
        panel.grid.major = element_line(colour = "gray90"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

rf.tax

#拼图

supp.left <- plot_grid(cv.plot,
                       overlap.plot,
                       nrow = 1,
                       rel_widths = c(2,3),
                       labels = c("A", "B"),
                       label_size = 20)

supp.left

plot_grid(supp.left,
          trend.p,
          ncol = 1,
          rel_heights = c(1,3),
          labels = c(NA, "C"),
          label_size = 20)


#使用这些OTU进行随机森林训练
get_rf <- function(imp.predictors, response) {
  randomForest(imp.predictors, response, importance = F, keep.forest = T)
}

train.set <- train.set %>% 
  mutate(rf = map2(imp.otu, age, get_rf))


# train.set %>% 
#   filter(Compartment == "RS") %>% 
#   .$rf
# 
# train.set %>% 
#   filter(Compartment == "ES") %>% 
#   .$rf



# 这里我们首先来理解一下作者的逻辑：首先使用这些微生物进行随机森林按照时间序列的回归。其次使用另外一部分样本预测模型，判断模型拟合程度
# 早期定殖微生物表现出逐渐下降的初始高丰度，晚期表现出逐渐增加的初始低丰度，复杂殖民者包括 OTU不符合这两种趋势中的任何一种。



get_predict <- function(model, data) {
  predict(model, data)
} 

# 预测时间
test.set <- whole.set %>% 
  inner_join(select(train.set, Compartment, rf), by = "Compartment") %>% 
  mutate(prediction = map2(rf, otu, get_predict)) 

#-分区域提取真实时间和预测时间进行相关性分析#--------
rs.pred <- test.set %>% 
  unnest(prediction, age) %>% 
  filter(Treatment == "WC") %>% 
  filter(Compartment == "RS") %>% 
  .$prediction

rs.age <- test.set %>%
  unnest(prediction, age) %>% 
  filter(Treatment == "WC") %>% 
  filter(Compartment == "RS") %>% .$age

cor(rs.pred, rs.age)


es.pred <- test.set %>% 
  unnest(prediction, age) %>% 
  filter(Treatment == "WC") %>% 
  filter(Compartment == "ES") %>% 
  .$prediction


es.age <- test.set %>%
  unnest(prediction, age) %>%
  filter(Treatment == "WC") %>% 
  filter(Compartment == "ES") %>% .$age
cor(es.pred, es.age)

#--预测和真实时间相关--以及
test.set %>% 
  unnest(prediction, age) %>% 
  ggplot(aes(age, prediction)) +
  geom_point(aes(shape = Set, color = Treatment), alpha = 0.7) + 
  geom_abline(linetype = 2) + 
  geom_smooth(aes(color = Treatment), linetype = 2, stat = "smooth", se = F) +
  scale_color_manual(values = trt.pal) +
  facet_grid(Treatment + Set ~ Compartment) +
  theme_light() +
  theme(text = element_text(size = 20)) 



# metadata $ID = row.names(metadata )
# metadata  = metadata  %>% dplyr::select(ID,Group,everything())
# sample_data(ps) = metadata

### 正常水分管理的预测曲线

predictions <- test.set %>% 
  unnest(prediction, age, data) %>% 
  ungroup() %>% 
  select(Compartment, Set, Treatment, age, prediction, ID)

predictions.wc <- filter(predictions, Treatment == "WC" & Set == "Test") %>% 
  ungroup()

predictions.wc.plotting <- rbind(select(predictions.wc, -Treatment) %>% mutate(Treatment = "DS1"),
                                 select(predictions.wc, -Treatment) %>% mutate(Treatment = "DS2"),
                                 select(predictions.wc, -Treatment) %>% mutate(Treatment = "DS3")) %>% 
  mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) 

# 在水分充足的植物的预测年龄和宿主的实际年龄之间拟合曲线。 该拟合将用作正常微生物组发育的基线，以计算干旱胁迫样品的相对微生物组成熟度。
# R语言的loess函数：当我们想研究不同sample的某个变量A之间的差异时，往往会因为其它一些变量B对该变量的固有影响，而影响不同sample变量A的比较，
# 这个时候需要对sample变量A进行标准化之后才能进行比较。标准化的方法是对sample 的 A变量和B变量进行loess回归，拟合变量A关于变量B的函数 f(b)
# ，f(b)则表示在B的影响下A的理论取值，A-f(B)（A对f(b）残差）就可以去掉B变量对A变量的影响,此时残差值就可以作为标准化的A值在不同sample之间进行比较


loess.pred <- rbind(
  data.frame(loess = predict(loess(prediction~age,filter(predictions.wc, Compartment == "RS")), unique(predictions.wc$age)),
             age = unique(predictions.wc$age),
             Compartment = "RS"),
  data.frame(loess = predict(loess(prediction~age,filter(predictions.wc, Compartment == "ES")), unique(predictions.wc$age)),
             age = unique(predictions.wc$age),
             Compartment = "ES"))

predictions.wc %>% ggplot(aes(age, prediction, color = Compartment)) +
  geom_smooth(linetype = 2,  size = 1, stat = "smooth", se = F) +
  geom_point(data = loess.pred, aes(age, loess))


# 预测的时间和处理和真实时间建模，显著性#-----
# 没有差异便表明
library(broom)
time.age <- map %>% 
  group_by(Time, Age) %>% 
  dplyr::count() %>% 
  select(-n) 

predictions <- inner_join(predictions, time.age, by = c("age" = "Age")) %>% 
  inner_join(select(map, ID, Row), by = "ID")
predictions$TimeFctr <- as.factor(predictions$Time)

run_lmm <- function(df) {
  lmerTest::lmer(prediction ~ Treatment * TimeFctr + (1|Row), data = df)
}

run_anova <- function(fit) {
  anova(fit) %>% tidy()
}

get_contrasts <- function(fit) {
  emmeans::emmeans(fit, specs = trt.vs.ctrl ~ Treatment|TimeFctr, adjust = "none") %>% 
    .$contrasts %>% 
    as.data.frame() %>% 
    separate(contrast, c("Treatment", "Control")) %>% 
    dplyr::rename("Time" = "TimeFctr") %>% 
    mutate(Time = as.integer(Time))
}


# df = predictions.nested[[1,2]][[1]]
# fit = lmerTest::lmer(prediction ~ Treatment * TimeFctr + (1|Row), data = df)


predictions.nested <- predictions %>% 
  group_by(Compartment) %>% 
  nest() %>% 
  mutate(lm = map(data, run_lmm),
         anova = map(lm, run_anova),
         contrasts = map(lm, get_contrasts))

rf.sig <- predictions.nested %>% 
  unnest(contrasts) %>% 
  group_by(Compartment) %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% 
  mutate(Time = as.integer(Time)) %>% 
  inner_join(time.age, by = "Time") %>% 
  ungroup()


a <- predictions %>% 
  filter(Treatment != "WC") %>% 
  mutate(Compartment = fct_relevel(Compartment, "RS", "ES")) %>%
  mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>%
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>%
  ggplot(aes(age, prediction)) +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment2), linetype = 3) +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment2), linetype = 3) +
  geom_point(aes(color = Treatment), size = 2, shape = 1, alpha = 1) + 
  geom_smooth(data = predictions.wc.plotting, linetype = 2,  size = 1, stat = "smooth", se = F, color = trt.pal[1]) +
  xlab("Plant age (days)") +
  ylab("Microbiome age\n(days)") +
  scale_color_manual(values = trt.pal[-1]) +
  facet_grid(Compartment ~ Treatment) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "none") 

a

### 微生物群落成熟度

# 相对成熟度和局部加权回归值差异。我们要注意，作者这里定义的群落的成熟度含义：预测时间和
# 局部加权回归值差异。

b <- predictions %>%
  filter(Treatment != "WC") %>% 
  inner_join(loess.pred, by = c("Compartment", "age")) %>% 
  mutate(RelMat = (prediction - loess)) %>% # 去除时间的影响
  group_by(Compartment, Treatment, age) %>% 
  summarise(MeanMat = mean(RelMat),
            SDMat = sd(RelMat)) %>%
  mutate(SEMat = SDMat/2) %>% 
  inner_join(rf.sig, by = c("Compartment","Treatment", "age" = "Age")) %>% 
  ungroup() %>% 
  mutate(FDR = ifelse(p.adj < 0.05, "< 0.05", "≥ 0.05")) %>% 
  mutate(Compartment = fct_relevel(Compartment, "RS", "ES")) %>% 
  mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>% 
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>% 
  mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  ggplot(aes(age, MeanMat, color = Treatment, shape = FDR)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = MeanMat - SEMat, ymax = MeanMat + SEMat)) +
  geom_vline(data = trt.lines, aes(xintercept = IniTreatment2), linetype = 3) +
  geom_vline(data = trt.lines, aes(xintercept = EndTreatment2), linetype = 3) +
  ylab("Relative microbiome\nmaturity (days)") +
  xlab("Plant age (days)") +
  facet_grid(Compartment ~ Treatment) +
  scale_color_manual(values = trt.pal[-1]) +
  scale_shape_manual(values = rev(ws.shp)) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "bottom") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5),
         shape = guide_legend(title.position = "top", title.hjust = 0.5))
b

# 特征微生物的时间序列和不同处理之间的变化

c <- otu.tidy %>% 
  mutate(Count = Count/100) %>% 
  inner_join(map, by = "ID") %>% 
  filter(Treatment2 != "WC_TRN") %>% 
  inner_join(all.ord, by = c("Compartment", "OTU_ID")) %>% 
  group_by(Compartment, Treatment, Age, Time, ZS_cluster, ID) %>% 
  summarise(Total = sum(Count)) %>% 
  group_by(Compartment, Treatment, Age, Time, ZS_cluster) %>% 
  mutate(MeanAb = mean(Total)) %>% 
  ungroup() %>% 
  # mutate(Trend = fct_relevel(Trend, "Complex", after = Inf)) %>% 
  mutate(Compartment = fct_relevel(Compartment, "RS", "ES")) %>% 
  mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>%
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>% 
  mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  ggplot(aes(Age, Total, fill = ZS_cluster, color = ZS_cluster)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(se = F) +
  # geom_vline(data = trt.lines, aes(xintercept = IniTreatment2), linetype = 3) +
  # geom_vline(data = trt.lines, aes(xintercept = EndTreatment2), linetype = 3) +
  ylab("Relative abundance") +
  xlab("Plant age (days)") +
  scale_fill_manual(values = c("dodgerblue", "gold", "grey50")) +
  scale_color_manual(values = c("dodgerblue", "gold", "grey50")) +
  facet_grid(Treatment ~ Compartment) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "bottom") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5),
         fill = guide_legend(title.position = "top", title.hjust = 0.5))


c

# 组合图形
left <- plot_grid(a,b + theme(legend.position = "none"), get_legend(b),
                  nrow = 3,
                  rel_heights = c(4,4,1),
                  labels = c("A", "B"),
                  label_size = 20)

right <- plot_grid(c + theme(legend.position = "none"), get_legend(c),
                   nrow = 2,
                   rel_heights = c(8,1),
                   labels = "C",
                   label_size = 20)

##1011:849
plot_grid(left, right,
          rel_widths = c(3,2))


### 微生物成熟度与植物生育期#-------

sfig.a <- predictions %>% 
  filter(Treatment == "WC") %>% 
  filter(Compartment == "ES") %>% 
  mutate(Treatment = str_replace(Treatment, "D", "DS")) %>% 
  mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  ggplot(aes(age, prediction)) +
  geom_point(color = trt.pal[1], size = 2, shape = 1, alpha = 1) + 
  geom_smooth(data = filter(predictions.wc.plotting, Treatment == "DS3" & Compartment == "Endosphere"), linetype = 2,  size = 1, stat = "smooth", se = F, color = trt.pal[1]) +
  xlab("Plant Age (days)") +
  ylab("Microbiome Age\n(days)") +
  scale_color_manual(values = trt.pal[-1]) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "none") 

sfig.a



sfig.b <- predictions %>% 
  filter(Treatment == "D3") %>% 
  filter(Compartment == "ES") %>% 
  inner_join(loess.pred, by = c("Compartment", "age")) %>% 
  mutate(RelMat = (prediction - loess)) %>% 
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>% 
  mutate(Treatment = fct_relevel(Treatment, "WC")) %>% 
  filter(ID == "DRGHT.TC.392") %>% 
  ggplot(aes(age, prediction)) +
  geom_point(size = 3, alpha = 1) + 
  geom_smooth(data = filter(predictions.wc.plotting, Treatment == "DS3" & Compartment == "Endosphere"), linetype = 2,  size = 1, stat = "smooth", se = F, color = trt.pal[1]) +
  geom_hline(aes(yintercept = prediction), linetype = 3) +
  geom_hline(aes(yintercept = loess), linetype = 3) +
  xlab("Plant Age (days)") +
  ylab("Microbiome Age\n(days)") +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "right") 

sfig.b


cowplot::plot_grid(sfig.a, sfig.b, sfig.c, nrow = 1, align = "hv", axis = "tblr", labels = c("A", "B", "C"), label_size = 20)


### 进行统计检验

rs.anova <- predictions.nested %>% 
  unnest(anova) %>% 
  ungroup() %>% 
  filter(Compartment == "RS") %>% 
  select(term:p.value)

es.anova <- predictions.nested %>% 
  unnest(anova) %>% 
  ungroup() %>% 
  filter(Compartment == "ES") %>% 
  select(term:p.value)

rs.contrasts <- rf.sig %>% 
  filter(Compartment == "RS") %>% 
  ungroup() %>% 
  mutate(Contrast = paste(Treatment, Control, sep = "-")) %>% 
  select(Contrast, Time:p.adj)

es.contrasts <- rf.sig %>% 
  filter(Compartment == "ES") %>% 
  ungroup() %>% 
  mutate(Contrast = paste(Treatment, Control, sep = "-")) %>% 
  select(Contrast, Time:p.adj)


#绘制每个植物的发育阶段，并将它们与预测的微生物组年龄进行比较#-----------

stages.tidy <- stages %>% 
  gather(key = "Time", value = "Stage", -Tub) %>% 
  mutate(Time = as.numeric(str_remove(Time, "T"))) %>% 
  inner_join(map, by = c("Tub", "Time")) %>% 
  filter(Treatment2 != "WC_TRN") %>% 
  mutate(Stage = fct_recode(Stage,
                            "Anthesis" = "anthesis",
                            "Grain filling" = "grain filling",
                            "Maturity" = "maturity",
                            "Panicle emergence" = "panicle emergence",
                            "Pre-panicle emergence" = "pre-panicle emergence")) %>% 
  mutate(Stage = fct_relevel(Stage,
                             "Pre-panicle emergence",
                             "Panicle emergence",
                             "Anthesis",
                             "Grain filling",
                             "Maturity")) 

stages.top <- stages.tidy %>% 
  ggplot(aes(Time, as.factor(Tub), fill = Stage)) +
  geom_tile(color = "white", size = 1) +
  #geom_vline(xintercept = c(4.5, 5.5, 6.5, 7.5), linetype = 3) +
  scale_fill_viridis_d(direction = -1) +
  xlab("Collection time point") +
  ylab("Plant replicate") +
  facet_grid(Treatment ~ ., scales = "free", space = "free") +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


stages.bottom <- predictions %>% 
  inner_join(select(stages.tidy, ID, Stage), by = c("ID")) %>% 
  filter(Treatment != "WC") %>% 
  mutate(Compartment = fct_relevel(Compartment, "RS", "ES")) %>% 
  mutate(Compartment = fct_recode(Compartment, Rhizosphere = "RS", Endosphere = "ES")) %>% 
  # mutate(Treatment = str_replace(Treatment, "D", "DS")) %>% 
  ggplot(aes(age, prediction)) +
  geom_vline(data = trt.lines, aes(xintercept = IniTreatment2), linetype = 3) +
  geom_vline(data = trt.lines, aes(xintercept = EndTreatment2), linetype = 3) +
  geom_point(aes(fill = Stage), size = 2, shape = 21, alpha = 0.8) + 
  geom_smooth(data = predictions.wc.plotting, linetype = 2,  size = 1, stat = "smooth", se = F, color = trt.pal[1]) +
  xlab("Plant age (days)") +
  ylab("Microbiome age\n(days)") +
  scale_color_manual(values = trt.pal[-1]) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(Compartment ~ Treatment) +
  theme_classic() +
  theme(text = element_text(size = 15),
        legend.position = "none") 

cowplot::plot_grid(stages.top, stages.bottom, nrow = 2, labels = c("A", "B"), label_size = 20)

##---备忘

# iso.otu <- readRDS("./Data/pheno_otu.RDS")
# head(iso.otu)
# 
# 
# iso.map <- readRDS("./Data/pheno_map.RDS")
# head(iso.map)
# row.names(iso.map) = iso.map$SampleID
# colnames(iso.map)[1] = "ID"
# 
# tax <- readRDS("./Data/tax.RDS") %>% 
#   classify_organelle()
# head(tax)
# row.names(tax) = tax$OTU_ID
# tax$OTU_ID = NULL
# 
# ps_verify = phyloseq(
#   otu_table(as.matrix(iso.otu),taxa_are_rows = T),
#   tax_table(as.matrix(tax)),
#   sample_data(iso.map)
#   
# )
# 
# saveRDS(ps_verify,"ps_verify.rds")


