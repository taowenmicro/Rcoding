
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(igraph)
library(RColorBrewer)

data(ps)
ps



#---全部模块的微生物网络#——------
pst =  ps %>%
  scale_micro("rela") %>%
  phyloseq::subset_samples(Group %in% c("KO","WT","OE")) %>%
  filter_OTU_ps(500)

result = cor_Big_micro(ps = pst,
                       N = 0,
                       r.threshold= 0.6,
                       p.threshold=0.05,
                       method = "spearman")

cor = result[[1]]
head(cor)

# igraph = make_igraph(cor)
#--计算模块信息，部分OTU没有模块，注意去除
# netClu  = modulGroup( cor = cor,cut = NULL,method = "cluster_fast_greedy" )
# head(netClu)
# result2 = model_maptree_group(cor = cor,
#                               nodeGroup = netClu,
# )

result2 = model_maptree2(cor = cor, method = "cluster_fast_greedy")
# result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]
netClu = result2[[2]]

# ---node节点注释
nodes = nodeadd(plotcord =node,
                otu_table = pst %>% 
                  vegan_otu() %>%
                  t() %>% 
                  as.data.frame(),
                tax_table = pst %>% vegan_tax() %>%
                  as.data.frame())
head(nodes)

nodes2 = nodes %>% inner_join(netClu,by = c("elements" = "ID"))
nodes2$group = paste("Model_",nodes2$group,sep = "")

#-----计算边
edge = edgeBuild(cor = cor,node = node)

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet


ggsave("cs.pdf",pnet,width = 8,height = 5)

#---选择模块的微生物网络#——------
select.mod = 3
select.mod = c("model_1","model_2","model_3")

mod1 = result2[[2]]
head(mod1)
tem = mod1$group %>% table() %>% 
  as.data.frame() %>% 
  dplyr::arrange(desc(Freq))
colnames(tem) = c("Model","OTU.num")
head(tem)


if (length(select.mod) == 1 & is.numeric(select.mod)) {
  select.mod.name = tem$Model[1:select.mod]
  mod1 = mod1 %>% filter(!group == "mother_no",
                         group %in%c(select.mod.name)
                         
  ) %>% select(ID,group,degree) 
  
} else if (is.character(select.mod)) {
  select.mod.name = select.mod
  mod1 = mod1 %>% filter(!group == "mother_no",
                         group %in%c(select.mod.name)
                         
  ) %>% select(ID,group,degree) 
  
}

head(mod1)
head(node)
node = result2[[1]] %>% filter(elements %in% mod1$ID)
# ---node节点注释
nodes = nodeadd(plotcord =node,
                otu_table = pst %>% 
                  vegan_otu() %>%
                  t() %>% 
                  as.data.frame(),
                tax_table = pst %>% vegan_tax() %>%
                  as.data.frame())
head(nodes)

nodes2 = nodes %>% inner_join(mod1,by = c("elements" = "ID"))

#-----计算边
edge = edgeBuild(cor = cor[mod1$ID,mod1$ID],node = node)

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = mean),pch = 21, data = nodes2) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

ggsave("cs.pdf",pnet,width = 6,height = 5)

#---选择模块的微生物组成#——------
tem = mod1$group %>% table() %>% 
  as.data.frame() %>% 
  dplyr::arrange(desc(Freq))
colnames(tem) = c("Model","OTU.num")
head(tem)

otu = NULL
map = NULL
for (i in 1:length(tem$Model)) {
  id.s = tem$Model %>% as.character()
  id.t =  mod1 %>% filter(group %in% id.s[i]) %>%.$ID
  ps.tem = subset_taxa(pst, row.names(tax_table(pst)) %in%id.t )
  
  otu = ps.tem %>% vegan_otu() %>% t() %>%
    as.data.frame()
  head(otu)
  colnames(otu) = paste(id.s[i],colnames(otu),sep = "_")
  map = data.frame(row.names = colnames(otu),ID = colnames(otu),Group = id.s[i])
  
  if (i == 1) {
    otu.f = otu
    map.f = map
  } else{
    
    otu$ID = row.names(otu)
    otu.f$ID = row.names(otu.f)
    tem.2 = otu.f %>% full_join(otu) %>%
      select(ID,everything()) %>%
      as.data.frame()
    row.names(tem.2) =  tem.2$ID
    tem.2$ID = NULL
    tem.2[is.na(tem.2)] = 0
    otu.f = tem.2
    map.f = rbind(map.f,map)
  }
}

pst.2 = phyloseq(
  otu_table(as.matrix(otu.f),taxa_are_rows = TRUE),
  sample_data(map.f) ,
  tax_table(pst)
)
# library(RColorBrewer)
# colset1 <- brewer.pal(10,"Paired")
# pst.2 = pst.2  %>%
#   subset_taxa(
#     !Genus %in% "Unassigned"
#   )
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")

j = "Family"
result = barMainplot(ps = pst.2,
                     tran = F,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 10)
  p4_1 <- result[[1]] +
    scale_fill_hue() + theme_classic()
  p4_1
  
  p4_2  <- result[[3]]  +
    scale_fill_hue() + theme_classic()
  p4_2

  result = barMainplot(ps = pst.2,
                        tran = T,
                        j = j,
                        # axis_ord = axis_order,
                        label = FALSE,
                        sd = FALSE,
                        Top = 10)
  p3_1 <- result[[1]] +
    scale_fill_hue() + theme_classic()
  p3_1
  
  p3_2  <- result[[3]]  +
    scale_fill_hue() + theme_classic()
  p3_2

#-选择模块的多样性#-----
  otu = NULL
  map = NULL
  for (i in 1:length(tem$Model)) {
    id.s = tem$Model %>% as.character()
    id.t =  mod1 %>% filter(group %in% id.s[i]) %>%.$ID
    ps.tem = subset_taxa(ps %>% scale_micro(method = "sampling"), row.names(tax_table(pst)) %in%id.t )
    
    otu = ps.tem %>% vegan_otu() %>% t() %>%
      as.data.frame()
    head(otu)
    colnames(otu) = paste(id.s[i],colnames(otu),sep = "_")
    map = data.frame(row.names = colnames(otu),ID = colnames(otu),Group = id.s[i])
    
    if (i == 1) {
      otu.f = otu
      map.f = map
    } else{
      
      otu$ID = row.names(otu)
      otu.f$ID = row.names(otu.f)
      tem.2 = otu.f %>% full_join(otu) %>%
        select(ID,everything()) %>%
        as.data.frame()
      row.names(tem.2) =  tem.2$ID
      tem.2$ID = NULL
      tem.2[is.na(tem.2)] = 0
      otu.f = tem.2
      map.f = rbind(map.f,map)
    }
  }
  
  pst.3 = phyloseq(
    otu_table(as.matrix(otu.f),taxa_are_rows = TRUE),
    sample_data(map.f) ,
    tax_table(pst)
  )
  
  #---多种指标alpha多样性分析加出图-标记显著性
  source("E:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")
  index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
  
  #--多种组合alpha分析和差异分析出图
  alp = alpha(ps = pst.3,inde="Shannon",group = "Group",Plot = TRUE,
              sampling = FALSE
              )
  index= alp
  head(index)
  
  #--提取三个代表指标作图
  sel = c(match("Shannon",colnames(index)),match("Richness",colnames(index)),match("Pielou_evenness",colnames(index)))
  data = cbind(data.frame(ID = 1:length(index$Group),group = index$Group),index[sel])
  head(data)
  
  result = EasyStat::MuiaovMcomper2(data = data,num = c(3:5))
  
  # FileName <- paste(alppath,"/alpha_diversity_different_label.csv", sep = "")
  # write.csv(result,FileName,sep = "")
  # FileName <- paste(alppath,"/alpha_diversity_index.csv", sep = "")
  # write.csv(index,FileName,sep = "")
  
  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:5),
                                            result = result,
                                            sig_show ="abc",ncol = 3 )
  p1_1 = result1[[1]] + 
    ggplot2::guides(fill = guide_legend(title = NULL)) 

  p1_1
  
  #如何升级展示-提取数据用小提琴图展示
  p1_1 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) + 
    geom_violin(alpha=1, aes(fill=group)) +
    geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
    labs(x="", y="")+
    facet_wrap(.~name,scales="free_y",ncol  = 3) +
    # theme_classic()+
    geom_text(aes(x=group , y=y ,label=stat)) +
    
    guides(color=guide_legend(title = NULL),
           shape=guide_legend(title = NULL),
           fill = guide_legend(title = NULL)
    ) 
  p1_1
  
#--基于模块的OTU，计算在不同分组中的总丰度zscore 并统计检验#-------
  head(mod1)

id.s = mod1$group %>% unique()

for (i in 1:length(id.s)) {
  id.t =  mod1 %>% filter(group %in% id.s[i]) %>%.$ID
  ps.t = ps %>% 
    scale_micro() %>%
    subset_taxa(row.names(tax_table(ps)) %in%id.t )
  
  otu = ps.t %>% vegan_otu() %>% t()
  
  
  
  colSD = function(x){
    apply(x,2, sd)
  }
  
  dat = (otu - colMeans(otu))/colSD(otu) 
  head(dat)
  otu_table(ps.t) = otu_table(as.matrix(dat),taxa_are_rows = T)
  
  #--计算总丰度
  
  otu = ps.t %>%  vegan_otu() %>% t()
  
  colSums(otu)
  
  dat = data.frame(id = names(colSums(otu)),abundance.zscore = colSums(otu))
  colnames(dat)[2] = id.s[i]
  
  if (i ==1) {
    tem = dat
  } else{
    dat$id = NULL
    tem = cbind(tem,dat)
  }
}
head(tem)
map =sample_data(ps.t)
map$id = row.names(map)
map = map[,c("id","Group")]
data = map %>%
  as.tibble() %>%
  inner_join(tem,by = "id") %>%
  dplyr::rename(group = Group)

num = 3+ length(id.s) -1
result = EasyStat::MuiaovMcomper2(data = data,num = c(3:num))

# FileName <- paste(alppath,"/alpha_diversity_different_label.csv", sep = "")
# write.csv(result,FileName,sep = "")
# FileName <- paste(alppath,"/alpha_diversity_index.csv", sep = "")
# write.csv(index,FileName,sep = "")

result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:num),
                                          result = result,
                                          sig_show ="abc",ncol = length(id.s) )
p1_1 = result1[[1]] + 
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  theme_classic()

p1_1

#如何升级展示-提取数据用小提琴图展示
p1_1 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) + 
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="Module Abundance (z-score)")+
  facet_wrap(.~name,scales="free_y",ncol  = length(id.s)) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat)) +
  
  guides(color=guide_legend(title = NULL),
         shape=guide_legend(title = NULL),
         fill = guide_legend(title = NULL)
  ) +
  theme_classic()
p1_1


#----多功能性和环境因子等和模块相关#-----
# 第一种是模块特征向量
# 模块丰度
data(env1)
head(env1)

select.env = "env1"

env1$id = row.names(env1)
env1 = env1 %>% dplyr::select(id,everything()) %>% select(id,select.env )
head(data)

tab = data %>% left_join(env1,by = "id")

head(tab)
library(reshape2)
mtcars2 = melt(tab, id.vars=c(select.env,"group","id"))
lab = mean(mtcars2[,select.env])
p1_1 = ggplot2::ggplot(mtcars2,aes(x= value,!!sym(select.env), colour=variable)) +
  ggplot2::geom_point() +
  ggpubr::stat_cor(label.y=lab*1.1)+
  ggpubr::stat_regline_equation(label.y=lab*1.1,vjust = 2) +
  facet_wrap(~variable, scales="free_x") +
  geom_smooth(aes(value,!!sym(select.env), colour=variable), method=lm, se=T)+
  theme_classic()

p1_1

#-多功能性和网络属性相关#---------

igraph = make_igraph(cor)

dat = igraph::V(igraph)
names(dat) %>% length()
#--弄清楚每个样本包含的OTU数量
# pst =  ps %>%
#   scale_micro("rela") %>%
#   phyloseq::subset_samples(Group %in% c("KO","WT","OE")) %>%
#   filter_OTU_ps(500) 
  

otu = ps %>% 
  phyloseq::subset_samples(Group %in% c("KO","WT","OE")) %>%
  # filter_OTU_ps(500) %>%
  subset_taxa(row.names(tax_table(ps)) %in% names(dat)) %>%
  vegan_otu() %>% 
  t() 
dim(otu)

otu[otu > 1] = 1
dim(otu)
A = list()
i = 1
for (i in 1:length(colnames(otu))) {
  tem = otu[,colnames(otu)[i]][otu[,colnames(otu)[i]] > 0 ] %>% names()
  A[[colnames(otu)[i]]] = tem
  #-计算性质
  tem.2 = A[[colnames(otu)[i]]]
  tem.g = igraph::induced_subgraph(igraph,tem.2)
  dat = net_properties.2(tem.g,n.hub = FALSE)
  head(dat,n = 16)
  
  dat[16,1] = 0
  dat = as.data.frame(dat)
  dat$value = as.numeric(dat$value)
  colnames(dat) = colnames(otu)[i]
  if (i == 1) {
    dat.f = dat
  } else {
    dat.f = cbind(dat.f,dat)
  }
}

dat.f = dat.f[1:15,] %>% t() %>% as.data.frame()



data(env1)
head(env1)

select.env = "env1"

env1$id = row.names(env1)
env1 = env1 %>% dplyr::select(id,everything()) %>% select(id,select.env )
head(dat.f)
dat.f$id = row.names(dat.f)
dat.f = dat.f %>% dplyr:: select(id,everything())
tab = dat.f %>% left_join(env1,by = "id")
head(tab)

mtcars2 = melt(tab, id.vars=c(select.env,"id"))

lab = mean(mtcars2[,select.env])
p0_1 = ggplot2::ggplot(mtcars2,aes(x= value,!!sym(select.env), colour=variable)) +
  ggplot2::geom_point() +
  ggpubr::stat_cor(label.y=lab*1.1)+
  ggpubr::stat_regline_equation(label.y=lab*1.1,vjust = 2) +
  facet_wrap(~variable, scales="free_x") +
  geom_smooth(aes(value,!!sym(select.env), colour=variable), method=lm, se=T)+
  theme_classic()

p0_1


