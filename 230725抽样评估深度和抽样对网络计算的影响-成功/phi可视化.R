
library(ggClusterNet)
library(tidyverse)
library(phyloseq)
# cor_Big_micro2 增加了标准化方法和p值矫正方法
result = cor_Big_micro2(ps = ps,
                        N = 1000,
                        r.threshold=0.85,
                        p.threshold=0.05,
                        method = "pearson",
                        scale = FALSE
)

#--提取相关矩阵
cor = result[[1]]
dim(cor)

# model_igraph2
library(igraph)
result2 <- model_Gephi.2(cor = cor
                        )
node = result2[[1]]
dim(node)


dat = result2[[2]]
head(dat)


#---node节点注释#-----------
otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)

#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)


library(ggnewscale)

p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
                              data = edge, size = 1) +
  new_scale_color() +
  geom_point(aes(X1, X2,color = Phylum), data = nodes,size = 4) +
  theme_void()
p1

#--下面是探索通模块聚到一起的过程#---------

head(edge)
head(nodes)
clutab = node[,1:2]

fit1 <- kmeans(clutab, 8)


dat2 = fit1$cluster %>% as.data.frame() %>% rownames_to_column("id")
colnames(dat2)[2] = "cluster"
node3 = nodes %>% rownames_to_column("id") %>%
  left_join(dat2,by = "id")

node3$cluster = as.factor(node3$cluster )
p1 <- ggplot() +
  # geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
  #                             data = edge, size = 1) +
  geom_point(aes(X1, X2,color = cluster), data = node3,size = 4) +
  theme_void()
p1

dim(node)
#--计算距离矩阵，选择距离最近的一批点#-------


gtab = data.frame(cluster = paste0("module_",c(1,2,3,4,5,6,7)),num = c(100,150,50,40,60,80,20))
layout = c("canberra","euclidean","maximum","manhattan","binary","minkowski")

plots = list()
i = 2
for (i in 1:length(layout)) {

  tem = dist(clutab,layout[i]) %>% as.matrix() %>%
    tidyfst::mat_df()
  head(tem)


  # 我们得到统计结果现在有一个列表

  tem2 = tem %>% filter(row == as.character(tem[1,1])
  ) %>%
    arrange(value) %>% select(-row)
  dim(tem2)
  tem2$group = rep(gtab$cluster,gtab$num)

  head(tem2)

  node3 = nodes %>% rownames_to_column("id") %>%
    left_join(tem2,by = c("id"="col"))
  head(node3)

  p1 <- ggplot() +
    # geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
    #                             data = edge, size = 1) +
    geom_point(aes(X1, X2,color = group), data = node3,size = 4) +
    theme_void()
  p1
  plots[[layout[i]]] = p1

}

# 通过查看深度基本上不影响相关，但是可能深度特别低会影响
p2  = ggpubr::ggarrange(plotlist = plots,
                        common.legend = FALSE, legend="right",ncol = 3,nrow = 2)
p2


#---第一个布局算法#-------

head(nodes)
gtab = nodes$Phylum %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(gtab) = c("cluster","num")
# gtab = data.frame(cluster = paste0("module_",c(1,2,3,4,5,6,7)),num = c(100,150,50,40,60,80,20)) %>%arrange(desc(num))
clutab = node[,1:2]
# 使用k-means算法进行聚类
kmeans_result <- kmeans(clutab, centers = (dim(gtab)[1] +1), nstart = 20, iter.max = 100)
cluster_centers <- kmeans_result$centers [1,]
aa = apply(kmeans_result$centers, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
colnames(aa) = "dis"
aa = aa %>% rownames_to_column("id") %>%
  arrange(dis)
center = kmeans_result$centers[aa$id,]

tem = c(0,0)
tem1 = dist(rbind(tem,center)) %>% as.matrix() %>% tidyfst::mat_df() %>%filter(row == "tem") %>%
  arrange(value)

center = center[-as.numeric(tem1[2,2]),]

for (i in 1:dim(gtab)[1]) {
  # 选择指定数量的数据点作为聚类结果
  cluster_centers <- center [i,]
  # 指定要聚为一类的数量
  num_points_in_cluster <- gtab[i,2]
  distances <- apply(clutab, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
  colnames(distances)[1] = "dis"
  id = distances %>% arrange(dis) %>% slice(1:num_points_in_cluster)%>%row.names()

  tem = data.frame(id = id,group = gtab[i,1])

  clutab= clutab %>% rownames_to_column("ID") %>%
    filter(!ID %in% id) %>% column_to_rownames("ID")

  if (i ==1) {
    tem2 = tem
  } else{
    tem2 = rbind(tem2,tem)
  }
}


node3 = nodes %>% rownames_to_column("id") %>%
  left_join(tem2,by = "id")
head(node3)

p1 <- ggplot() +
  # geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
  #                             data = edge, size = 1) +
  geom_point(aes(X1, X2,color = group), data = node3,size = 4) +
  theme_void()
p1

#--第二个布局算法#------
dim(nodes)
gtab = nodes$Phylum %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(gtab) = c("cluster","num")

gtab = data.frame(cluster = paste0("module_",c(1,2,3,4,5,6,7,8,9,10,11,12,13)),num = c(100,150,50,40,60,80,20,60,140,80,120,50,50)) %>%arrange(desc(num))
clutab = node[,1:2]
# # 使用k-means算法进行聚类
# kmeans_result <- kmeans(clutab, centers = (dim(gtab)[1] +1), nstart = 20, iter.max = 100)
# aa = apply(kmeans_result$centers, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
# colnames(aa) = "dis"
# aa = aa %>% rownames_to_column("id") %>%
#   arrange(dis)
# center = kmeans_result$centers[aa$id,]
#
# tem = c(0,0)
# tem1 = dist(rbind(tem,center)) %>% as.matrix() %>% tidyfst::mat_df() %>%filter(row == "tem") %>%
#   arrange(value)
#
# center = center[-as.numeric(tem1[2,2]),]

i = 11
for (i in 1:dim(gtab)[1]) {

  # # 使用k-means算法进行聚类
  # if (dim(clutab)[1]< 4 ) {
  #
  # }

  if (dim(clutab)[1] > (dim(gtab)[1] +2-i)) {
    kmeans_result <- kmeans(clutab, centers = (dim(gtab)[1] +2-i), nstart = 20, iter.max = 100)
    aa = apply(kmeans_result$centers, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
    colnames(aa) = "dis"
    aa = aa %>% rownames_to_column("id") %>%
      arrange(dis)
    center = kmeans_result$centers[aa$id,]

    tem = c(0,0)
    tem1 = dist(rbind(tem,center)) %>% as.matrix() %>% tidyfst::mat_df() %>%filter(row == "tem") %>%
      arrange(value)

    center = center[-as.numeric(tem1[2,2]),]
  } else if(dim(clutab)[1] == (dim(gtab)[1] +2-i)) {
    kmeans_result <- kmeans(clutab, centers = (dim(gtab)[1] +1-i), nstart = 20, iter.max = 100)
    aa = apply(kmeans_result$centers, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
    colnames(aa) = "dis"
    aa = aa %>% rownames_to_column("id") %>%
      arrange(dis)
    center = kmeans_result$centers[aa$id,]

    tem = c(0,0)
    tem1 = dist(rbind(tem,center)) %>% as.matrix() %>% tidyfst::mat_df() %>%filter(row == "tem") %>%
      arrange(value)

    center = center
  }

  cluster_centers <- center[1,]
  # 指定要聚为一类的数量
  num_points_in_cluster <- gtab[i,2]
  distances <- apply(clutab, 1, function(x) sqrt(sum((x - cluster_centers)^2))) %>% as.data.frame()
  colnames(distances)[1] = "dis"
  id = distances %>% arrange(dis) %>% slice(1:num_points_in_cluster)%>%row.names()

  tem = data.frame(id = id,group = gtab[i,1])

  clutab= clutab %>% rownames_to_column("ID") %>%
    filter(!ID %in% id) %>% column_to_rownames("ID")

  if (i ==1) {
    tem2 = tem
  } else{
    tem2 = rbind(tem2,tem)
  }
}


node3 = nodes %>% rownames_to_column("id") %>%
  left_join(tem2,by = "id")
dim(node3)
node3$size = runif(1000,0.1,3)
test_color5 <- RColorBrewer::brewer.pal(n = 12,name = "Set3")
p1 <- ggplot() +
  # geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2),
  #                             data = edge, size = 1) +
  geom_point(aes(X1, X2,color = group,size = size), data = node3) +
  scale_color_manual(values = c(test_color5,"grey80","grey20")[13:1]) +
  theme_void()
p1

