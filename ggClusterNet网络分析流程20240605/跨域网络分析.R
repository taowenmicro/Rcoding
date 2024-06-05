
library(tidyverse)
library(phyloseq)
library(ggClusterNet)


#  差异基因psr3#----------
psr = base::readRDS("./data/ps.trans.rds")
map = sample_data(psr)
head(map)
map$Group  %>% unique()

tax = as.data.frame(vegan_tax(psr))
head(tax)
tax$KO_id = tax$KO_id  %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(psr)  = tax_table(as.matrix(tax))
psr2 <- psr %>% subset_taxa.wt("KO_id" ,"" ,T) %>% 
  tax_glom_wt("KO_id")
psr2



diffpath.2 = "./tem"
dir.create(diffpath.2)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
res = EdgerSuper(ps = psr2,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)
head(res)
id = res %>% filter(`KO-WTlevel` != "nosig") %>% arrange(desc(`KO-WTlogFC`)) %>% .$species
psr3 = psr2 %>% subset_taxa.wt("OTU",id)


# 转录组-Mkegg psr5#-----
getOption("clusterProfiler.download.method")
R.utils::setOption( "clusterProfiler.download.method",'auto' )
#--基于Mkegg进行分析
Mkegg = clusterProfiler::download_KEGG('ko',keggType = "MKEGG")
PATH2ID <- Mkegg$KEGGPATHID2EXTID
PATH2NAME <- Mkegg$KEGGPATHID2NAME
head(PATH2NAME)
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')
head(PATH_ID_NAME )


ko = psr2 %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)

tem2 = ko %>% left_join(PATH_ID_NAME,by = c("KO" = "KO")) %>%
  tidyfst::filter_dt(KEGGID != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("KEGGID","DESCRPTION")) %>%
  as.data.frame()
head(tem2)
colnames(tem2)

otu = tem2[,sample_names(psr2)]
tax = data.frame(MDESCRPTION = tem2$DESCRPTION,MDESCRPTION2 = tem2$DESCRPTION,row.names = tem2$KEGGID)


row.names(otu) = tem2$KEGGID
head(otu)
psr4 = phyloseq(
  otu_table(as.matrix(otu),taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(psr)
)


otupath = "./tem/"
diffpath.2 = paste(otupath,"/EDgeR.Mkegg/",sep = "")
dir.create(diffpath.2)

res = EdgerSuper(ps = psr4,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)

head(res)


id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`KO-WTlevel` != "nosig") %>%
  arrange(desc(`KO-WTlogFC`))
id = id$ID
psr5 = psr4 %>% subset_taxa.wt("OTU",id)

# 转录组-reaction psr7#-----
dat = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/ko-reaction.txt",
                 header = F
)

colnames(dat) = c("ID","reaction")
head(dat)

ko = psr2 %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)
dat$ID = gsub("ko:","",dat$ID)
dat$reaction = gsub("rn:","",dat$reaction)


tem2 = ko %>% left_join(dat,by = c("KO" = "ID")) %>%
  tidyfst::filter_dt(reaction != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("reaction")) %>%
  as.data.frame()

colnames(tem2)

otu = tem2[,sample_names(psr2)]
row.names(otu) = tem2$reaction
head(otu)


dat2 = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/reaction.txt",
                  header = F
)
head(dat2)
colnames(dat2) = c("reaction","DESCRPTION")

tax = data.frame(row.names = dat2$reaction,DESCRPTION = dat2$DESCRPTION,DESCRPTION2 = dat2$DESCRPTION)
head(tax)

psr6 = phyloseq(
  otu_table(as.matrix(otu),taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(psr)
)

psr6
otupath = "./tem/"
diffpath.2 = paste(otupath,"/EDgeR.reaction/",sep = "")
dir.create(diffpath.2)

res = EdgerSuper(ps = psr6,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)

head(res)


id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`KO-WTlevel` != "nosig") %>%
  arrange(desc(`KO-WTlogFC`))
id = id$ID
psr7 = psr6 %>% subset_taxa.wt("OTU",id)





# 差异代谢物 psG3#----------
psG = readRDS("./data/ps.metabolite.rds")
#---开始分析#--------
ps = psG
map = ps %>%sample_data()
head(map)


tax = ps %>% tax_table() %>% as.data.frame()
head(tax)
tax$KEGG.Compound.ID = tax$KEGG.Compound.ID %>% strsplit( "[;]") %>% 
  sapply(`[`, 1)
tax_table(ps)  = tax_table(as.matrix(tax))

psG2 = ps %>% subset_taxa.wt("KEGG.Compound.ID",c("-"),T) %>%
  tax_glom_wt("KEGG.Compound.ID")

#--t检验检验--建议四个重复以上
source("E:\\Shared_Folder\\Function_local\\R_function\\GC-MS\\wlxSuper_GCMS.R")
result2 = statSuper(ps = psG2,group  = "Group",artGroup = NULL,method = "ttext")
head(result2)

id2  = result2 %>% filter(`KO_WT_fdr`<0.05) %>% arrange(desc(`KO_WT_log2_FC`)) %>%
  .$KEGG.Compound.ID
psG3 = psG2 %>% subset_taxa.wt("OTU",id2)




#  网络图中的标签我们要进行很好的设定
ps.all = merge.ps(ps1 = psr3,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "rna",
                  dat2.lab = "compounds")
ps.all 


# 差异宏基因组基因 psko2 #-----
psko = readRDS("./data/ps.meta.rds")


#  去除空缺值
otu = psko %>% vegan_otu() %>% t()
otu[is.na(otu)] =0
otu_table(psko) = otu_table(otu,taxa_are_rows = TRUE)
psko = psko %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)


map = psko %>%sample_data()
head(map)
tax = psko %>% tax_table() %>% as.data.frame()
head(tax)

otupath = "./tem/"
diffpath = paste(otupath,"/Different.gene.ko/",sep = "")
dir.create(diffpath)
md = c("edgr","t","desep2","wlx")
# 准备脚本
source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/EdgerSuper_Meta.R")
# source("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/Plot.CompareWithCK.R",encoding = "utf-8")
supath = paste0(diffpath,"/EdgeR/") 
dir.create(supath)
res = EdgerSuper(ps = psko,group  = "Group",artGroup =NULL,
                 path = supath
)
head(res)
id = res %>% 
  rownames_to_column("ID") %>% 
  filter(`KO-WTlevel` != "nosig") %>% arrange(desc(`KO-WTlogFC`)) %>% .$ID
id
psko2 = psko %>% subset_taxa.wt("OTU",id)

# 宏基因组差异Mkegg  psko4#------
getOption("clusterProfiler.download.method")
R.utils::setOption( "clusterProfiler.download.method",'auto' )
#--基于Mkegg进行分析
Mkegg = clusterProfiler::download_KEGG('ko',keggType = "MKEGG")
PATH2ID <- Mkegg$KEGGPATHID2EXTID
PATH2NAME <- Mkegg$KEGGPATHID2NAME
head(PATH2NAME)
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')
head(PATH_ID_NAME )


ko = psko %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)

tem2 = ko %>% left_join(PATH_ID_NAME,by = c("KO" = "KO")) %>%
  tidyfst::filter_dt(KEGGID != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("KEGGID","DESCRPTION")) %>%
  as.data.frame()
head(tem2)
colnames(tem2)

otu = tem2[,sample_names(psko)]
tax = data.frame(MDESCRPTION = tem2$DESCRPTION,MDESCRPTION2 = tem2$DESCRPTION,row.names = tem2$KEGGID)


row.names(otu) = tem2$KEGGID
head(otu)
psko3 = phyloseq(
  otu_table(as.matrix(otu),taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(ps)
)


otupath = "./tem/"
diffpath.2 = paste(otupath,"/EDgeR.Mkegg/",sep = "")
dir.create(diffpath.2)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
res = EdgerSuper(ps = psko3,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)

head(res)


id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`KO-WTlevel` != "nosig") %>%
  arrange(desc(`KO-WTlogFC`))
id = id$ID
psko4 = psko3 %>% subset_taxa.wt("OTU",id)

# 宏基因组差异reaction psko6#------


dat = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/ko-reaction.txt",
                 header = F
)

colnames(dat) = c("ID","reaction")
head(dat)



# kegg <- clusterProfiler::download_KEGG('ko')
# PATH2ID <- kegg $KEGGPATHID2EXTID
# PATH2NAME <- kegg$KEGGPATHID2NAME
# head(PATH2NAME)
# PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
# colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')
# head(PATH_ID_NAME )

ko = psko %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)
dat$ID = gsub("ko:","",dat$ID)
dat$reaction = gsub("rn:","",dat$reaction)


tem2 = ko %>% left_join(dat,by = c("KO" = "ID")) %>%
  tidyfst::filter_dt(reaction != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("reaction")) %>%
  as.data.frame()

colnames(tem2)

otu = tem2[,sample_names(psr2)]
row.names(otu) = tem2$reaction
head(otu)


dat2 = read.delim("E:/Shared_Folder/Function_local/R_function/Metagenome_Function/meta.mini.db/db.KEGG/reaction.txt",
                  header = F
)
head(dat2)
colnames(dat2) = c("reaction","DESCRPTION")

tax = data.frame(row.names = dat2$reaction,DESCRPTION = dat2$DESCRPTION,DESCRPTION2 = dat2$DESCRPTION)
head(tax)

psko5 = phyloseq(
  otu_table(as.matrix(otu),taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(ps)
)

psko5
otupath = "./tem/"
diffpath.2 = paste(otupath,"/EDgeR.reaction/",sep = "")
dir.create(diffpath.2)

res = EdgerSuper(ps = psko5,group  = "Group",artGroup = NULL,
                 j = "gene",
                 path = diffpath.2
)

head(res)


id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`KO-WTlevel` != "nosig") %>%
  arrange(desc(`KO-WTlogFC`))
id = id$ID
psko6 = psko5 %>% subset_taxa.wt("OTU",id)



# 根际细菌#--------
# 16s数据
# otu = ps.all %>% vegan_otu()

ps01 = readRDS("./data/ps.micro.rds") 
otu1 = ps01 %>% vegan_otu()
map = ps01 %>%sample_data()
head(map)
ps02 = ps01 %>% tax_glom_wt("Genus")
# map = ps02 %>% sample_data()
otu = ps01 %>% vegan_otu()
otu[1:5,1:5]

source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")

res = EdgerSuper(ps = ps01,
                 group  = "Group",
                 artGroup = NULL,
                 j = "Genus",
                 path = diffpath.2
)
head(res)

id = res %>% 
  rownames_to_column("ID" ) %>%
  filter(`KO-WTlevel` != "nosig") %>%
  arrange(desc(`KO-WTlogFC`))
id = id$ID
ps03 = ps02 %>% subset_taxa.wt("OTU",id)


# match(row.names(otu1),row.names(otu))


# 合并全部ps对象#------

ps.all = merge.ps(ps1 = psr5,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "rna",
                  dat2.lab = "compounds")
ps.all
ps.all2 = merge.ps(ps1 = ps.all,
                   ps2 = psko4,
                   N1 = 0,
                   N2 = 0,
                   scale = TRUE,
                   onlygroup = TRUE,#不进行列合并，只用于区分不同域
                   dat1.lab = "",
                   dat2.lab = "meta")

tax = ps.all2 %>% vegan_tax() %>% as.data.frame()
tax$filed %>% unique()


tax = ps.all2 %>% vegan_tax() %>% as.data.frame()
tax$filed %>% unique()

ps.all3 = merge.ps(ps1 = ps.all2,
                   ps2 = ps03,
                   N1 = 0,
                   N2 = 0,
                   scale = TRUE,
                   onlygroup = TRUE,#不进行列合并，只用于区分不同域
                   dat1.lab = "",
                   dat2.lab = "micro")

tax = ps.all3 %>% vegan_tax() %>% as.data.frame()
tax$filed %>% unique()



#  转录组-代谢组-宏基因组基因--微生物群落网络#-------
# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
library(sna)
library(igraph)

res = corBionetwork.st(
  ps.st= ps.all3,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]



repath = "./Muilti.network/"
path = paste0(repath,"./bio.rna.compounds.meta.micro.network/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 60,height = 30,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 90,height = 40,limitsize = FALSE)


map = ps.all3 %>% sample_data()
map$Group = "one"
sample_data(ps.all3) = map


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all3,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

# path = paste0(repath,"./bio.rna.compounds.network.all.togather.reaction/")
# dir.create(path)
path
ggsave(paste(path,"bionetwork.one.pdf",sep = ""),  p,width = 40,height = 40,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.labelone.pdf",sep = ""),  p0,width = 50,height = 50,limitsize = FALSE)






#  转录组-代谢组-宏基因组基因-#-------
# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all2,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.meta/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 60,height = 30,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 90,height = 40,limitsize = FALSE)


map = ps.all2 %>% sample_data()
map$Group = "one"
sample_data(ps.all2) = map


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all2,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

# path = paste0(repath,"./bio.rna.compounds.network.all.togather.reaction/")
# dir.create(path)
path
ggsave(paste(path,"bionetwork.one.pdf",sep = ""),  p,width = 40,height = 40,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.labelone.pdf",sep = ""),  p0,width = 50,height = 50,limitsize = FALSE)






#  转录组-代谢组-微生物-#-------
# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可

ps.all = merge.ps(ps1 = psr5,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "rna",
                  dat2.lab = "compounds")
ps.all
ps.all2 = merge.ps(ps1 = ps.all,
                   ps2 = ps03,
                   N1 = 0,
                   N2 = 0,
                   scale = TRUE,
                   onlygroup = TRUE,#不进行列合并，只用于区分不同域
                   dat1.lab = "",
                   dat2.lab = "micro")

tax = ps.all2 %>% vegan_tax() %>% as.data.frame()
tax$filed %>% unique()


res = corBionetwork.st(
  ps.st= ps.all2,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds.micro/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork.pdf",sep = ""),  p,width = 60,height = 30,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 90,height = 40,limitsize = FALSE)


map = ps.all2 %>% sample_data()
map$Group = "one"
sample_data(ps.all2) = map


# 写死的filed为分组，也就是我们合并好需需要做跨一网络的对象，指定filed作为分组即可
res = corBionetwork.st(
  ps.st= ps.all2,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14)

p = res[[1]]
dat = res[[2]]

# path = paste0(repath,"./bio.rna.compounds.network.all.togather.reaction/")
# dir.create(path)
path
ggsave(paste(path,"bionetwork.one.pdf",sep = ""),  p,width = 40,height = 40,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.5) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.labelone.pdf",sep = ""),  p0,width = 50,height = 50,limitsize = FALSE)





col1 = c("red","blue")
names(col1) = c("+","-")


# 20240326更新 1 转录组和分泌物的相关网络#----------
ps.all = merge.ps(ps1 = psr5,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "rna",
                  dat2.lab = "compounds")
ps.all 


library(igraph)
library(sna)

tax = ps.all  %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.all ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.rna.compounds_2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 18,height = 10,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  # scale_colour_brewer(palette = "Set1") +
  scale_color_manual(values  = col1)+
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 80,height = 30,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对

#  定义处理组
tre = "WT"

A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))


# 2 分泌物,微生物的相关网络#----------
ps.all = merge.ps(ps1 = ps03,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "micro.16s",
                  dat2.lab = "compounds")
ps.all 


library(igraph)
library(sna)

tax = ps.all  %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.all ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.compounds.micro.16s_2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 18,height = 13,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  # scale_colour_brewer(palette = "Set1") +
  scale_color_manual(values  = col1)+
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 80,height = 33,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对


A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))



# 3 分泌物,微生物,功能的相关网络#----------
ps.all = merge.ps(ps1 = ps03,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "micro.16s",
                  dat2.lab = "compounds")
ps.all 

ps.all2 = merge.ps(ps1 = ps.all ,
                   ps2 = psko4,
                   N1 = 0,
                   N2 = 0,
                   scale = TRUE,
                   onlygroup = TRUE,#不进行列合并，只用于区分不同域
                   dat1.lab = "",
                   dat2.lab = "meta")
ps.all2



library(igraph)
library(sna)

tax = ps.all2  %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.all2 ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.compounds.micro.meta._2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 33,height = 12,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 30,height = 33,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对


A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))



# 3 转录组,微生物的相关网络#----------
ps.all = merge.ps(ps1 = ps03,
                  ps2 = psr5,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "micro.16s",
                  dat2.lab = "trans")
ps.all 




library(igraph)
library(sna)

tax = ps.all  %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.all ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.trans.micro._2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 12,height = 12,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  # scale_colour_brewer(palette = "Set1") +
  scale_color_manual(values  = col1)+
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 33,height = 33,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对


A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))



# 4 功能,微生物的相关网络#----------
ps.all = merge.ps(ps1 = ps03,
                  ps2 = psko4,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "micro.16s",
                  dat2.lab = "meta")
ps.all 




library(igraph)
library(sna)

tax = ps.all  %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.all ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.meta.micro._2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 20,height = 12,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  # scale_colour_brewer(palette = "Set1") +
  scale_color_manual(values  = col1)+
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 33,height = 33,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对


A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))





# 5 分泌物,微生物,功能的相关网络#----------
ps.all = merge.ps(ps1 = ps03,
                  ps2 = psG3,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "micro.16s",
                  dat2.lab = "compounds")
ps.all 

ps.all2 = merge.ps(ps1 = ps.all ,
                   ps2 = psko4,
                   N1 = 0,
                   N2 = 0,
                   scale = TRUE,
                   onlygroup = TRUE,#不进行列合并，只用于区分不同域
                   dat1.lab = "",
                   dat2.lab = "meta")
ps.all2


ps.all3 = merge.ps(ps1 = ps.all2,
                   ps2 = psr5,
                   N1 = 0,
                   N2 = 0,
                   scale = TRUE,
                   onlygroup = TRUE,#不进行列合并，只用于区分不同域
                   dat1.lab = "",
                   dat2.lab = "trans")


library(igraph)
library(sna)

tax = ps.all3 %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.all3 ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.compounds.micro.meta.trans._2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 40,height = 18,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  # scale_colour_brewer(palette = "Set1") +
  scale_color_manual(values  = col1)+
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 3 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 60,height = 20,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对


A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))

#  
# 6 分泌物,功能的相关网络#----------


ps.t = merge.ps(ps1 = psG3 ,
                ps2 = psko4,
                N1 = 0,
                N2 = 0,
                scale = TRUE,
                onlygroup = TRUE,#不进行列合并，只用于区分不同域
                dat1.lab = "rootexudate",
                dat2.lab = "meta")
ps.t



library(igraph)
library(sna)

tax = ps.t %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.t ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.compounds.meta._2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 18,height = 12,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 33,height = 33,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对


A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))

# 7 功能 和转录的相关网络#----------


ps.t = merge.ps(ps1 = psr5 ,
                ps2 = psko4,
                N1 = 0,
                N2 = 0,
                scale = TRUE,
                onlygroup = TRUE,#不进行列合并，只用于区分不同域
                dat1.lab = "trans",
                dat2.lab = "meta")
ps.t



library(igraph)
library(sna)

tax = ps.t %>% vegan_tax() %>% as.data.frame()
head(tax)
nod.gro = data.frame(ID = tax$id,group = tax$filed)
nod.gro$group %>% unique()
res = corBionetwork.st(
  ps.st= ps.t ,# phyloseq对象
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  # layout_net = "PolygonRrClusterG",
  layout_net ="model_filled_circle",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  N= 0,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,maxsize = 14,
  group.node = nod.gro ,
  model.node = FALSE
  
)

#  可视化展示#----
p = res[[1]]
p

dat = res[[2]]

path = paste0(repath,"./bio.trans.meta._2/")
fs::dir_create(path)
ggsave(paste(path,"bionetwork2.pdf",sep = ""),  p,width = 18,height = 12,limitsize = FALSE)

node = dat[[2]]
edge = dat[[3]]
head(edge)
head(node)
p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                              data = edge, size = 0.3,alpha = 0.1) +
  geom_point(aes(x = X1, y = X2,size = igraph.degree,fill = filed),pch = 21, data =  node) +
  ggrepel::geom_text_repel(aes(X1, X2,label= id),size=4, data = node) +
  scale_colour_brewer(palette = "Set1") +
  scale_size(range = c(14, 4)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_")) +
  facet_wrap(.~ label,scales="free_y",ncol = 2 ) +
  theme_void()
# p0
ggsave(paste(path,"bionetwork.label.pdf",sep = ""),  p0,width = 33,height = 33,limitsize = FALSE)


#  网络边的挖掘#------
head(edge)
#  定义两组的网络，对照组有的连线和处理组有的连线进行比对


A = paste0(edge$OTU_2,edge$OTU_1)
tem = table (A) %>% as.data.frame() %>% arrange(desc(Freq)) %>%filter(Freq == 2) %>%
  .$A %>% as.character()
gro = edge$group %>%strsplit( "[.]") %>% 
  sapply(`[`, 3)
gro[A %in% tem] 
gro[A %in% tem] = "both"
edge$link.gro = gro
write_csv(edge,paste(path,"bionetwork.edge.csv",sep = ""))
write_csv(node,paste(path,"bionetwork.node.csv",sep = ""))



