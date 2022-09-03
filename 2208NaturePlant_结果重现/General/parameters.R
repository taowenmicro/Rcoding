# Palettes for each factor 

trt.pal <- RColorBrewer::brewer.pal(11, "BrBG")[c(10,3,2,1)]
org.pal <- c("gray20", "gray40", "salmon")
cmp.pal <- RColorBrewer::brewer.pal(12, "Set3")[c(5,4)]
lib.pal <- RColorBrewer::brewer.pal(3, "Set2")
time.pal <- viridis::viridis(14) 
phy.pal <- c("gray35",
             RColorBrewer::brewer.pal(8, "Set2")[1:4],
             "indianred3",
             RColorBrewer::brewer.pal(8, "Set2")[5:8],
             RColorBrewer::brewer.pal(11, "RdYlBu")[7:10])
resp.pal <- RColorBrewer::brewer.pal(11,"BrBG")[c(9,3)]

phy.pal.df <- data.frame(PhyClass2 = c("other", 
                           "Acidobacteria", 
                           "Actinobacteria", 
                           "Bacteroidetes", 
                           "Chloroflexi", 
                           "Fibrobacteres", 
                           "Firmicutes", 
                           "Gemmatimonadetes", 
                           "Planctomycetes", 
                           "Verrucomicrobia", 
                           "Alphaproteobacteria", 
                           "Betaproteobacteria",
                           "Deltaproteobacteria", 
                           "Gammaproteobacteria"),
             Color = c("gray50",
                       RColorBrewer::brewer.pal(8, "Set2")[1:4],
                       "indianred3",
                       RColorBrewer::brewer.pal(8, "Set2")[5:8],
                       RColorBrewer::brewer.pal(11, "RdYlBu")[7:10]))

ws.shp <- c(1,16)

ws.line <- c(41,52,62,74)