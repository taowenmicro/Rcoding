
#--增加网络稳定性内容------

# 输入使用列表，共三列，OTUid，模块信息id，网络分组id
library(ggClusterNet)
library(phyloseq)
library(tidyverse)
library(igraph)
data(ps)
ps 
map = sample_data(ps)
head(map)
id <- map$Group %>% unique()
i  = 1
#------网络出图#---------

netpath = paste("./network2/",sep = "")
dir.create(netpath)

library(igraph)
library(sna)
library(phyloseq)
library(ggClusterNet)
gnum = 3

result =network.2(ps = ps, 
                                 N = 500,
                                 big = TRUE,
                                 maxnode = 5,
                                 select_layout = TRUE,
                                 layout_net = "model_maptree2",
                                 r.threshold=0.8,
                                 p.threshold=0.05,
                                 label = FALSE,
                                 path = netpath,
                                 zipi = F,
                                 ncol = gnum,
                                 nrow = 1,
                                 # method = "sparcc",
                                 fill = "Phylum"
)

# 全部样本的网络比对
p4_1 = result[[1]] 
# 全部样本网络参数比对
data = result[[2]]
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
gnum = 45
ggsave(plotname1, p4_1,width = 6*gnum,height = 6,limitsize = FALSE)
# plotname1 = paste(netpath,"/network_all.jpg",sep = "")
# ggsave(plotname1, p4_1,width = 16*gnum,height = 16)
tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)
# 全部样本的网络比对
p4_2 = result[[3]] + 
  scale_fill_brewer(palette = "Paired")
plotname1 = paste(netpath,"/network_all_cover.pdf",sep = "")
ggsave(plotname1, p4_2,width = 10*gnum,height = 10,limitsize = FALSE)

#-判断相同模块#----
for (i in 1:length(id)) {
  pst =  ps %>% 
    scale_micro() %>%
    phyloseq::subset_samples(Group %in% c(id[i])) %>%
    filter_OTU_ps(500)
  
  result = cor_Big_micro(ps = pst,
                         N = 0,
                         r.threshold= 0.8,
                         p.threshold=0.05,
                         method = "spearman")
  
  cor = result[[1]]
  head(cor)
  
  # igraph = make_igraph(cor)
  #--计算模块信息，部分OTU没有模块，注意去除
  result2 = model_maptree2(cor = cor,
                           method = "cluster_fast_greedy"
  )
  
  mod1 = result2[[2]]
  head(mod1)
  
  mod1 = mod1 %>% filter(!group == "mother_no") %>% select(ID,group) 
  mod1$group = paste(id[i],mod1$group,sep = "")
  mod1$Group = id[i]
  head(mod1)
  
  if (i == 1) {
    dat = mod1
  } else {
    dat = rbind(dat,mod1)
  }
}

node_table2  = dat
head(node_table2)

library(tidyfst)
XX$Group %>% unique()
XX = node_table2
modc = XX %>% filter(Group %in% c("G0" ,"G15","G28","G36","G5","L15","L2","L25","L40"))
head(modc)
#  
#--模块比对
# source("E:\\Shared_Folder\\Function_local\\R_function\\my_R_packages\\ggClusterNet_Decument\\网络稳定性预封装.R")

dat = model_compare(
    node_table2 = XX,
    n = 3,
    padj = FALSE
    )



dim(dat)

dir.create("./NO1_zhong/")


write.table(dat, "./NO1_zhong/preserved_module_pairs.txt", sep="\t")

#--随即取出任意比例节点-网络鲁棒性#---------
library(ggClusterNet)
library(phyloseq)

##read otu table
otutab<- ps %>% 
  vegan_otu() %>% 
  as.data.frame()
dim(otutab)

id <- sample_data(ps)$Group %>% unique()
i  = 1
#计算每个物种的平均丰度，使用测序深度标准化
sp.ra<-colMeans(otutab)/mean(rowSums(otutab))   #relative abundance of each species
library(igraph)

for (i in 1:length(id)) {
  pst =  ps %>% 
    scale_micro() %>%
    phyloseq::subset_samples(Group %in% c(id[i])) %>%
    filter_OTU_ps(500)
  
  result = cor_Big_micro(ps = pst,
                         N = 0,
                         r.threshold= 0.8,
                         p.threshold=0.05,
                         method = "spearman")
  
  cor = result[[1]]
  head(cor)
  
  #存在某些情况计算不出来相关系数，定义相关为0
  cor[is.na(cor)]<-0
  #-去除自相关的点
  diag(cor)<-0  
  #-查看网络边的数量
  sum(abs(cor)>0)/2
  #网络中节点的数量
  sum(colSums(abs(cor))>0)  # node number: number of species with at least one linkage with others.
  #去除没有任何相关的节点.
  network.raw<-cor[colSums(abs(cor))>0,colSums(abs(cor))>0]
  #对应的删除otu表格otu
  sp.ra2<-sp.ra[colSums(abs(cor))>0]
  sum(row.names(network.raw)==names(sp.ra2))  #check if matched
  
  ## 鲁棒性评估robustness simulation 
  #input network matrix, percentage of randomly removed species, and ra of all species
  #return the proportion of species remained
  
  
  #输入相关矩阵 OTU表格
  Weighted.simu<-rmsimu(netRaw = network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, 
                        abundance.weighted=T,nperm=100)
  head(Weighted.simu)
  Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2,
                          abundance.weighted=F,nperm=100)
  head(Weighted.simu)
  tem = pst %>% sample_data() %>% .$Group %>% unique() %>% as.character()
  
  dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),
                   rbind(Weighted.simu,Unweighted.simu),
                   weighted=rep(c("weighted","unweighted"),each=20),
                   Group=tem)
  
  head(dat1)
  
  library(ggplot2)
  
  p = ggplot(dat1[dat1$weighted=="weighted",], 
             aes(x=Proportion.removed, y=remain.mean, group=Group, color=Group)) + 
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
    scale_color_manual(values=c("blue","red"))+
    xlab("Proportion of species removed")+
    ylab("Proportion of species remained")+
    theme_light()
  
  
  p1 = ggplot(dat1[dat1$weighted=="unweighted",], aes(x=Proportion.removed,
                                                      y=remain.mean , group=Group, color=Group)) + 
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
    scale_color_manual(values=c("blue","red"))+
    xlab("Proportion of species removed")+
    ylab("Proportion of species remained")+
    theme_light()
  
  library(patchwork)
  p3 = p|p1
  dir.create("./NO1_zhong/Robustness_Random_removal/")
  write.csv(dat1,
            paste("./NO1_zhong/Robustness_Random_removal/",id[i],"_random_removal_network.csv",sep = ""))
  
  ggsave(paste("./NO1_zhong/Robustness_Random_removal/",id[i],"_random_removal_network.pdf",sep = ""),
         p3,width = 8,height = 4
  )
  
}


#---去除关键节点-网络鲁棒性#------
library(ggClusterNet)
library(phyloseq)
library(igraph)

i = 1

for (i in 1:length(id)) {
  pst =  ps %>% 
    scale_micro() %>%
    phyloseq::subset_samples(Group %in% c(id[i])) %>%
    filter_OTU_ps(500)
  
  result = cor_Big_micro(ps = pst,
                         N = 0,
                         r.threshold= 0.8,
                         p.threshold=0.05,
                         method = "spearman")
  
  cor = result[[1]]
  head(cor)
  
  
  #存在某些情况计算不出来相关系数，定义相关为0
  cor[is.na(cor)]<-0
  #-去除自相关的点
  diag(cor)<-0  
  #-查看网络边的数量
  sum(abs(cor)>0)/2
  #网络中节点的数量
  sum(colSums(abs(cor))>0)  # node number: number of species with at least one linkage with others.
  #去除没有任何相关的节点.
  network.raw<-cor[colSums(abs(cor))>0,colSums(abs(cor))>0]
  
  ##read otu table
  
  otutab<- pst %>% 
    scale_micro() %>%
    subset_taxa(
      row.names(tax_table(pst)) %in% row.names(network.raw)
    ) %>%
    vegan_otu() %>% 
    t() %>%
    as.data.frame()
  
  
  
  #对应的删除otu表格otu
  sp.ra2<- rowSums(otutab)
  sp.ra2
  sum(row.names(network.raw) %in% names(sp.ra2))  #check if matched
  ## robustness simulation 
  #input network matrix, number of removed keystone species, keystonespecies list,  and ra of all species
  #return the proportion of species remained
  
  #get the keystone species list
  igraph = make_igraph(cor)
  
  degree = TRUE
  zipi = FALSE
  if (degree) {
    ret3 = node_properties(igraph) %>% 
      as.data.frame() %>%
      filter(!is.na(igraph.degree) ) %>% 
      arrange(desc(igraph.degree)) 
    head(ret3)
    tem = round(length(ret3$igraph.degree) * 0.05,0)
    module.hub = row.names(ret3)[1:tem]
  }
  
  
  
  if (zipi ) {
    res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
    p <- res[[1]]
    model = res[[2]] %>% filter(roles == "Module hubs")
    head(model)
    model$roles %>% unique()
    module.hub <- as.character(row.names(model))  
  }
  
  
  
  Weighted.simu<-rmsimu(netRaw=network.raw,
                        rm.p.list=1:length(module.hub),
                        keystonelist=module.hub,
                        sp.ra=sp.ra2,
                        abundance.weighted=T,nperm=100)
  Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),
                          keystonelist=module.hub,
                          sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  
  dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),
                   rbind(Weighted.simu,Unweighted.simu),
                   weighted=rep(c("weighted","unweighted"),
                                each=length(module.hub)),
                   Group= id[i])
  
  head(dat1)
  currentdat = dat1
  
  dir.create("./NO1_zhong/Robustness_Targeted_removal/")
  
  write.csv(dat1,
            paste("./NO1_zhong/Robustness_Targeted_removal/",id[i],"_random_removal_network.csv",sep = ""))
  
  
  
  ##plot
  library(ggplot2)
  
  p = ggplot(currentdat[currentdat$weighted=="weighted",], 
             aes(x=Number.hub.removed, y=remain.mean, group=Group, color=Group)) + 
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd, 
                        ymax=remain.mean+remain.sd),size=0.2)+
    scale_color_manual(values=c("red","blue"))+
    xlab("Number of module hubs removed")+
    ylab("Proportion of species remained")+
    theme_light()
  
  
  p2 = ggplot(currentdat[currentdat$weighted=="unweighted",], 
              aes(x=Number.hub.removed, y=remain.mean, group=Group, color=Group)) + 
    geom_line()+
    geom_pointrange(aes(ymin=remain.mean-remain.sd, 
                        ymax=remain.mean+remain.sd),size=0.2)+
    scale_color_manual(values=c("blue","red"))+
    xlab("Number of module hubs removed")+
    ylab("Proportion of species remained")+
    theme_light()
  p3 = p | p2
  
  ggsave(paste("./NO1_zhong/Robustness_Targeted_removal/",id[i],
               "_Targeted_removal_network.pdf",sep = ""),
         p3,width = 8,height = 4
  )
}

#---网络易损性#------
# library(igraph)
# source("./info.centrality.R")
i = 1
A = c()

for (i in 1:length(id)) {
  pst =  ps %>% 
    scale_micro() %>%
    phyloseq::subset_samples(Group %in% c(id[i])) %>%
    filter_OTU_ps(500)
  
  result = cor_Big_micro(ps = pst,
                         N = 0,
                         r.threshold= 0.8,
                         p.threshold=0.05,
                         method = "spearman")
  
  cor = result[[1]]
  
  vulnerability = function(cor = cor){
    
    cor[abs(cor)>0]<-1 # adjacency matrix
    g = graph_from_adjacency_matrix(as.matrix(cor), 
                                    mode="undirected", 
                                    weighted = NULL, diag = FALSE,
                                    add.colnames = NULL) # note: this graph contains isolated nodes.
    # remove isolated nodes
    iso_node_id = which(degree(g)==0)
    g2 = delete.vertices(g, iso_node_id) # graph without isolated nodes
    
    #check node number and links
    length(V(g2));length(E(g2))   
    
    # calculate vulnerability of each node
    node.vul<-info.centrality.vertex(g2)
    return(max(node.vul))
  }
  
  A[i] = vulnerability(cor = cor)
  
  
}

dat = data.frame(ID = id,Vulnerability = A)

dir.create("./NO1_zhong/Vulnerability/")

write.csv(dat,
          paste("./NO1_zhong/Vulnerability//","_Vulnerability_network.csv",sep = ""))

#--计算负相关的比例#----
B = c()
for (i in 1:length(id)) {
  pst =  ps %>% 
    scale_micro() %>%
    phyloseq::subset_samples(Group %in% c(id[i])) %>%
    filter_OTU_ps(500)
  
  result = cor_Big_micro(ps = pst,
                         N = 0,
                         r.threshold= 0.8,
                         p.threshold=0.05,
                         method = "spearman")
  
  cor = result[[1]]
  igraph = make_igraph(cor)
  ret2 = net_properties.2(igraph) %>% as.data.frame()
  head(ret2)
  a = ret2[1,1] %>% as.numeric()
  n = ret2[3,1] %>% as.numeric()
  B[i] = n/a *100
}


dat = data.frame(ID = id,Vulnerability = B)

dir.create("./NO1_zhong/negative_correlation_ratio/")

write.csv(dat,
          paste("./NO1_zhong/negative_correlation_ratio/","_negative_ratio_network.csv",sep = ""))


# 群落稳定性-必须是pair的样本才能计算#——-------

otutab = ps %>% 
  scale_micro() %>%
  phyloseq::subset_samples(Group %in% c(id[i])) %>%
  vegan_otu() %>% t() %>%
  as.data.frame()
head(otutab)
#分组文件
treatused = ps %>% sample_data()

library(tidyverse)

comm=otutab %>% t()
treat=treatused

#去除NA值
sum(is.na(comm)) # check NA
comm[is.na(comm)]=0# if have, should change to zero
head(treat)
treat$pair = paste( "A",c(rep(1:6,3)),sep = "")
plot.lev=unique(treat$pair)
#-提取时间序列
year.lev = sort(unique(treat$Group))
#-构造序列
zeta.lev =2:length(year.lev)

i = 1
# 构造从2到6的全部这组合，这里使用断棍模型构造全部组合
year.windows=lapply(1:length(zeta.lev),
                    function(i)
                    {zetai=zeta.lev[i]
                    lapply(1:(length(year.lev)-zetai+1),function(j){year.lev[j:(j+zetai-1)]})
                    })
names(year.windows)=zeta.lev
year.windows

# 基于不同分组样本的群落稳定性功能函数物种最小丰度和乘以样本数量，得到的结果除以多组全部微生物丰度的和
comstab<-function(subcom){((nrow(subcom)*sum(apply(subcom,2,min)))/sum(subcom))^0.5}
# subcom = comijk
# 下面写了一个循环计算的全部的两两比对，三个比对，等全部比对
# 但是文章中也仅仅用了两组
stabl=lapply(1:length(year.windows),
             function(i)
             {
               stabi=t(sapply(1:length(plot.lev),
                              function(j)
                              {
                                plotj=plot.lev[j]
                                sapply(1:length(year.windows[[i]]),
                                       function(k)
                                       {
                                         yearwdk=year.windows[[i]][[k]] %>% as.character()
                                         sampijk=rownames(treat)[which((treat$pair==plotj) & (treat$Group %in% yearwdk))]
                                         outijk=NA
                                         if(length(sampijk) < length(yearwdk))
                                         {
                                           warning("plot ",plotj," has missing year in year window ",paste(yearwdk,collapse = ","))
                                         }else if(length(sampijk) > length(yearwdk)){
                                           warning("plot ",plotj," has duplicate samples in at least one year of window ",paste(yearwdk,collapse = ","))
                                         }else{
                                           comijk=comm[which(rownames(comm) %in% sampijk),,drop=FALSE]
                                           outijk=comstab(comijk)
                                         }
                                         outijk
                                       })
                              }))
               if(nrow(stabi)!=length(plot.lev) & nrow(stabi)==1){stabi=t(stabi)}
               rownames(stabi) = plot.lev
               colnames(stabi)=sapply(year.windows[[i]],function(v){paste0("Zeta",zeta.lev[i],paste0(v,collapse = ""))})
               stabi
             }) 
library(ggClusterNet)


stabm=Reduce(cbind,stabl)
head(stabm)

