# 
# 
# #  这里要探讨什么问题呢？
# # 重复数量对网络构建的影响
# 
# 同一个土壤的样本基本上都没有什么变化，发现重复在10个左右就可以基本上达到很好的状态了
# 似乎没有必要持续测下去，没必要，增加重复数量网络整体结构和参数都差不多
# 
# 
# # 测序深度对网络构建的影响
# 
# 深度越高，检测出来的相关性会越多，深度越低，可能会造成的假阳性，这可能是0值过多导致的，大量0值直接导致相关性接近1，无限逼近1
# 这里虽然我们已经去除了大量的全零情况，但实际上有许多0也会导致大量假阳性.
# 
# 讨论：三个重复的数据为什么会计算出来大量的相关性为1或者-1的关系呢？
# 相关性为1的大概率为这两个OTU都有大量的0 值，并且是对应的0值，啥意思呢，就是一个样本中这两个OTU都是0；
#  负相关那就正好相反，一个OTU是0，另一个是1,这种完美的匹配。无论是哪种，基本上都是因为大量的小值和0值导致的
#  表现到图上那就是集中的一个网络模块都是点，连线几乎看不见，因为相关性太高了。
#  这样的结果也就是不能使用的，可以按照低丰度和相关性绝对值为1的标准去除这些关系。保证真实相关性。
#  这种情况会随着测序深度增加和样本量的增多逐渐弱化,这里测试的结果基本上10个左右，就不存在这种情况了，
# 但是这是在测序深度保证3w以上的情况下的。
# 
# 
# 不同深度抽样增加重复性是否可以带来较好的网络构建过程
# 这是可以的，不同深度的抽样会并且和原来数据合并起来，会增加相关性，但这些增加的基本上是假阳性，不推荐使用
# 
# 
# 我们对一个土壤测定了63个重复，通过相同的网络分析发现了重复越多，网络能鉴定出来的相关性会越来越多
#  这是为什么呢？这次我们测定的是一个土壤，并且使用搅拌机混匀的，也就意味着基本上微生物混合的很均匀
#  这个时候基本上就算是抽样得到的结果了，也就是说我们的测试结果表明实际测样的生物学重复与我们的的模拟重复
# 实际上是一样的。也就是说我们将样本混匀，只需要测两个重复就够了，因为得到的微生物群落基本上是不变的
# 那么我们平时去测定那么多重复是为什么呢？实际上随着技术发展，深度逐渐增加，需要的样本数量越来越少了。
#  这里就有一个问题，那么如此少的样本还能做网络分析吗？很明显不能，因为我们找到的没什么意思了，基本上都是成比例变化
# 的，而不是真正具有相关关系的关系。
# 
#  所以如果我们真的要寻找土壤中互作关系，需要在不同条件下得到更多具有变化的样本才好说这些可能是真正的关系
# 大家也都知道，网络基本上在生态中的，随着生态随机因素的变化，网络实际上可以求出一个很好的值，但是我们却很少能够说
#  这一网络的某种特征与某个因素，环境，或者其他因素有关系。
#  土壤微生物在土壤中的互作关系用共线性网络并不能解决，我们能尝试得到的无非是外界环境变化下，一同变化的微生物。如此
# 我们也无法确定这些微生物是自己互作关系发生变化，还是受到环境独立的改变？
# 所以我们应该如何寻找微生物之间的互作关系呢？这是一个大挑战！


library(ggClusterNet)
library(tidyverse)
library(phyloseq)
library(patchwork)

# 统计最小的深度
samplesize = min(phyloseq::sample_sums(ps))

# 提取分组信息查看具体分组标签
map = sample_data(ps)
map$Group %>%  table() %>% as.data.frame()


# 构造一个临时的网络分析函数，用于快速构建网络
# ps.f %>% vegan_otu() %>% t() %>%
#   as.data.frame() %>% rownames_to_column("id") %>%
#   filter(id %in% c("ASV_4","ASV_43"))

network.tem = function(ps.f,r = 0.8){
  result = cor_Big_micro2(ps = ps.f,
                          N = 500,
                          r.threshold=0.8,
                          p.threshold=0.05,
                          method = "pearson",
                          scale = FALSE
  )

  #--提取相关矩阵
  cor = result[[1]]
  dim(cor)


  result2 <- model_igraph2(cor = cor,
                           method = "cluster_fast_greedy",
                           seed = 12
  )
  node = result2[[1]]
  dim(node)


  dat = result2[[2]]
  head(dat)
  tem = data.frame(mod = dat$model,col = dat$color) %>%
    dplyr::distinct( mod, .keep_all = TRUE)
  col = tem$col
  names(col) = tem$mod

  #---node节点注释#-----------
  otu_table = as.data.frame(t(vegan_otu(ps)))
  tax_table = as.data.frame(vegan_tax(ps))
  nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
  head(nodes)
  #-----计算边#--------
  edge = edgeBuild(cor = cor,node = node)
  colnames(edge)[8] = "cor"
  head(edge)


  tem2 = dat %>%
    dplyr::select(OTU,model,color) %>%
    dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
    dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
  head(tem2)

  tem3 = dat %>%
    dplyr::select(OTU,model,color) %>%
    dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
    dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
  head(tem3)

  tem4 = tem2 %>%inner_join(tem3)
  head(tem4)

  edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                          manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
  )


  col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>%
    select(color,manual)
  col0 = col_edge$manual
  names(col0) = col_edge$color

  library(ggnewscale)

  p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                                data = edge2, size = 1) +
    scale_colour_manual(values = col0)

  # ggsave("./cs1.pdf",p1,width = 16,height = 14)
  p2 = p1 +
    new_scale_color() +
    geom_point(aes(X1, X2,color =model), data = dat,size = 4) +
    scale_colour_manual(values = col) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  p2
}

# 测序深度对网络构建的影响#-------
# 对微生物组数据进行相对丰度标准化之后提取KO组进行网络构建
ps.f = ps %>%
  scale_micro() %>%
  subset_samples.wt("Group","KO")
p0 = network.tem(ps.f)

ggsave("cs1.pdf",p0,width = 8,height = 7)

#--我们降低测序深度，检测对网络的影响
plots = list()

step = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

for (n in 1:length(step)) {

  ps.f = ps %>%
    # scale_micro() %>%
    subset_samples.wt("Group","KO")
  ps.f2 = ps.f %>%
    phyloseq::rarefy_even_depth(sample.size = min(sample_sums(ps.f))*step[n])

  p1 = network.tem(ps.f2)
  plots[[n]] = p1
}

# 通过查看深度基本上不影响相关，但是可能深度特别低会影响
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2

ggsave("cs2.pdf",p2,width = 8*5,height = 7*2,limitsize = FALSE)

#--极低深度测试
plots = list()
data(ps)
step = c(0.001,0.003,0.005,0.008,0.01,0.1,0.3,0.5,0.7,1)

n = 1

for (n in 1:length(step)) {

  ps.f = ps %>%
    # scale_micro() %>%
    subset_samples.wt("Group","KO")
  ps.f2 = ps.f %>%
    phyloseq::rarefy_even_depth(sample.size = min(sample_sums(ps.f))*step[n])%>%
    filter_taxa(function(x) sum(x ) > 0 , TRUE)

  p1 = network.tem(ps.f2) + labs(title = step[n]) +
    theme(plot.title=element_text(hjust=0.5))
  plots[[n]] = p1
}

# 这里就看到深度剩下三千一下，网络假阳性就逐渐出现了
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2

ggsave("cs3.pdf",p2,width = 8*5,height = 7*2,limitsize = FALSE)


# 调试到多少阈值，假阳性开始突增

plots = list()

step = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.1,1)

for (n in 1:length(step)) {

  ps.f = ps %>%
    # scale_micro() %>%
    subset_samples.wt("Group","KO")
  ps.f2 = ps.f %>%
    phyloseq::rarefy_even_depth(sample.size = min(sample_sums(ps.f))*step[n])

  p1 = network.tem(ps.f2) + labs(title = step[n]) +
    theme(plot.title=element_text(hjust=0.5))
  plots[[n]] = p1
}

# 这里看到基本上两千条序列一下做网络假阳性就会逐渐增加
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2

ggsave("cs4.pdf",p2,width = 8*5,height = 7*2,limitsize = FALSE)




#重复数量对网络构建的影响#------
#-基于网络重复性过少的问题现在

ps.f = ps %>%
  scale_micro() %>%
  subset_samples.wt("Group","KO")
p0 = network.tem(ps.f)

#--随机抽样三个样本的网络
ps.f0 = random.add.ps(ps = ps,add = 3)
ps.f = ps.f0 %>%
  scale_micro() %>%
  subset_samples.wt("Group","KO") %>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)
p1 = network.tem(ps =ps.f)
p0|p1


#--四个重复
ps.f0 = random.add.ps(ps = ps,add = 4)
ps.f = ps.f0 %>%
  scale_micro() %>%
  subset_samples.wt("Group","KO") %>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)
p2 = network.tem(ps.f)
p0|p1|p2

#--五个重复
ps.f0 = random.add.ps(ps = ps,add = 5)
ps.f = ps.f0 %>%
  scale_micro() %>%
  subset_samples.wt("Group","KO")%>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)
p3 = network.tem(ps.f)
p0|p1|p2|p3

#--六个重复
ps.f0 = random.add.ps(ps = ps,add = 6,zoom = 0.3,addid = "")
ps.f2 = ps.f0 %>%
  # scale_micro() %>%
  subset_samples.wt("Group","KO")%>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)
p4 = network.tem(ps.f2)
p0|p1|p2|p3|p4



#-七个重复
ps.f0 = random.add.ps(ps = ps,add = 1)
#--合并ps对象
psf = merge_amp.dat(ps.a = ps,
                    ps.m = ps.f0,
                    rank = "OTU")

ps.f = psf %>% subset_samples.wt("Group","KO") %>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)

o1 = network.tem(ps.f)

p0|o1

#-八个重复
ps.f0 = random.add.ps(ps = ps,add = 2)
#--合并ps对象
psf = merge_amp.dat(ps.a = ps,
                    ps.m = ps.f0,
                    rank = "OTU")

ps.f = psf %>% subset_samples.wt("Group","KO") %>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)

o2= network.tem(ps.f)

p0|o1|o2

#-十个重复
ps.f0 = random.add.ps(ps = ps,add = 4)
#--合并ps对象
psf = merge_amp.dat(ps.a = ps,
                    ps.m = ps.f0,
                    rank = "OTU")

ps.f = psf %>% subset_samples.wt("Group","KO")%>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)

o3 = network.tem(ps.f)

p0|o1|o2|o3

#-十二个重复
ps.f0 = random.add.ps(ps = ps,add = 6)
#--合并ps对象
psf = merge_amp.dat(ps.a = ps,
                    ps.m = ps.f0,
                    rank = "OTU")

ps.f = psf %>% subset_samples.wt("Group","KO")

o4 = network.tem(ps.f)

p0/p0|p1/o1|p2/o2|p3/o3|p4/o4


# 上面增加的重复为模拟的，如果我实测数据，从3个重复到12个重复会怎么样？
# 下面我们测了

ps1 = readRDS("E:/Shared_Folder/Function_local/R_function/my_R_packages/ggClusterNet_Decument/测试网络构建的重复数量/ps.rds")

map = sample_data(ps1)
map$Group = map$SampleType
map = map[paste0(60:126,"A")]
map = data.frame(row.names = map$ID[!is.na(map$ID)],
                 ID = map$ID[!is.na(map$ID)], Group = c("one.soil"))
map$Group[2:30] = "two.soil"
sample_data(ps1) = map
ps1 = ps1 %>%  filter_taxa(function(x) sum(x ) > 0 , TRUE)
saveRDS(ps1,"ps.63.rds")

ps1 = readRDS("ps.63.rds")

#--我们做一下排序看看这批数据分布
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")
result = BetaDiv(ps = ps1 %>%  filter_taxa(function(x) sum(x ) > 0 , TRUE) ,
                 group = "Group", dist = "bray",
                   method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                   pair = FALSE)
p3_1 = result[[3]]
p3_1

#--开始用这批数据，抽取不同的重复数量构建网络
map = sample_data(ps1)
map$Group = "one.soil"
sample_data(ps1) = map
ps1 = ps1 %>%  filter_taxa(function(x) sum(x ) > 0 , TRUE)


ps.f = random.add.ps(ps = ps1,add = 63)
u0 = network.tem(ps.f)


plots = list()

step = c(3,6,9,12,15,18,21,24,28,30,36,42,48,50,58,63)

n = 4

for (n in 1:length(step)) {

  ps.f = random.add.ps(ps = ps1,add = step[n])%>%
    filter_taxa(function(x) sum(x ) > 0 , TRUE)
  u0 = network.tem(ps.f)+ labs(title = step[n]) +
    theme(plot.title=element_text(hjust=0.5))
  plots[[n]] = u0
}

# 这里我们看到3个重复还是会增加许多假阳性，但是6个以上似乎就好多了
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 8,nrow = 2)
p2
ggsave("cs5.pdf",p2,width = 70,height = 15,limitsize = FALSE)



# 这里按照1-12个数量再做一次抽样
plots = list()

step = c(3:12)

for (n in 1:length(step)) {

  ps.f = random.add.ps(ps = ps1,add = step[n])%>%
    filter_taxa(function(x) sum(x ) > 0 , TRUE)
  u0 = network.tem(ps.f)+ labs(title = step[n]) +
    theme(plot.title=element_text(hjust=0.5))
  plots[[n]] = u0
}

# 这里我们看到3个重复还是会增加许多假阳性，但是6个以上似乎就好多了
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2
ggsave("cs6.pdf",p2,width = 40,height = 12,limitsize = FALSE)


# 重复数量多了会增加很多弱相关，如果我们提高阈值

plots = list()

step = c(3:12)

for (n in 1:length(step)) {

  ps.f = random.add.ps(ps = ps1,add = step[n])%>%
    filter_taxa(function(x) sum(x ) > 0 , TRUE)
  u0 = network.tem(ps.f,r = 0.8)+ labs(title = step[n]) +
    theme(plot.title=element_text(hjust=0.5))
  plots[[n]] = u0
}

# 我们发现随着相关阈值改变，网络布局基本没有什么变化，全部是大相关
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2
ggsave("cs7.pdf",p2,width = 40,height = 12,limitsize = FALSE)


# 我们这这批数据在测试一次深度对网络的影响

plots = list()

step = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

for (n in 1:length(step)) {

  ps.f = ps1 #%>%
  # scale_micro() %>%
  # subset_samples.wt("Group","KO")
  ps.f2 = ps.f %>%
    phyloseq::rarefy_even_depth(sample.size = min(sample_sums(ps.f))*step[n])

  p1 = network.tem(ps.f2)
  plots[[n]] = p1
}

# 新数据来看通过查看深度基本上不影响相关，但是可能深度特别低会影响
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2
ggsave("cs8.pdf",p2,width = 40,height = 12,limitsize = FALSE)



# 调整深度细点

plots = list()

step = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.1,1)

for (n in 1:length(step)) {

  ps.f = ps
  ps.f2 = ps.f %>%
    phyloseq::rarefy_even_depth(sample.size = min(sample_sums(ps.f))*step[n])

  p1 = network.tem(ps.f2) + labs(title = step[n]) +
    theme(plot.title=element_text(hjust=0.5))
  plots[[n]] = p1
}

# 深度降低，相关性减弱，也就是说有一部分微生物关系检测不到了
#  其假阳性的区间似乎更低一下，如此多的样本的话
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2

ggsave("cs9.pdf",p2,width = 40,height = 12,limitsize = FALSE)


#--再次运行极低深度测试
plots = list()

step = c(0.001,0.003,0.005,0.008,0.01,0.1,0.3,0.5,0.7,1)

n = 3
for (n in 1:length(step)) {

  ps.f = ps
  ps.f2 = ps.f %>%
    phyloseq::rarefy_even_depth(sample.size = min(sample_sums(ps.f))*step[n]) %>%
    filter_taxa(function(x) sum(x ) > 0 , TRUE)

  p1 = network.tem(ps.f2) + labs(title = step[n]) +
    theme(plot.title=element_text(hjust=0.5))
  plots[[n]] = p1
}

# 这里就看到深度降低到180个左右就可以降低很多的假阳性
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 5,nrow = 2)
p2
ggsave("cs10.pdf",p2,width = 40,height = 12,limitsize = FALSE)

