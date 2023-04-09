
library(ggClusterNet)
library(tidyverse)
library(phyloseq)
data(ps)

# 准备脚本
source("./EdgerSuper.R")
source("./EdgerSuper2.R")
source("./Mui.cluster-group.volcano.R")
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Plot.CompareWithCK.R",encoding = "utf-8")
group1 = c("OE","WT")
b= data.frame(group1)
diffpath.1 = "./"

res = EdgerSuper(ps = ps,group  = "Group",artGroup =b,
                   j = "OTU",
                   path = diffpath.1
)

head(res)


result = Mui.cluster.volcano(res = res)
p = result[[1]]
p

p = result[[2]]
p

ggsave("cs.pdf",p,width = 8,height = 6)

#----多组火山图#------

diffpath.1 = "./result_cs"
dir.create(diffpath.1)
res = EdgerSuper2 (ps = ps,group  = "Group",artGroup =NULL,
                   j = "OTU",
                   path = diffpath.1
)

head(res)


result = Mui.Group.volcano (res = res)
p = result[[2]]
p
ggsave("cs4.pdf",p5,width = 12,height = 6,limitsize = FALSE)



