

p.lvl =0.05
lda.lvl = 2.0



# otu = read.csv("./data_filtered.csv",row.names = 1)
otu = read.csv("./data_filtered.csv",row.names = 1)
head(otu)  

map = read.delim("./map.txt",row.names = 1)
head(map)


map = map[colnames(otu),]
otu = t(otu)
claslbl= map$SampleType
set.seed(56290);
  #KW rank sum test
  dat3t = otu
  head(otu)
  rawpvalues <- apply(dat3t, 2, function(x) kruskal.test(x, claslbl)$p.value);
  #--得到计算后得到的p值
  ord.inx <- order(rawpvalues);
  rawpvalues <- rawpvalues[ord.inx];
  clapvalues <- p.adjust(rawpvalues, method ="fdr");
  # p.adjust
  dat3t <- dat3t[,ord.inx];
  dim(dat3t)


  
  wil_datadf <- as.data.frame(dat3t);
  head(wil_datadf)
  #if no subclass within classes then no wilcoxon rank sum test  
  #Linear Discriminant analysis (LDA)
  library( MASS)
  ldares <- lda(claslbl~ .,data = wil_datadf);
  ldares
  ldamean <- as.data.frame(t(ldares$means));
  ldamean 
  class_no <<- length(unique(claslbl));
  ldamean$max <- apply(ldamean[,1:class_no],1,max);
  ldamean$min <- apply(ldamean[,1:class_no],1,min);
  
  #---计算LDA
  ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);
  
  
  
  signif(log10(1+(ldamean$max-ldamean$min)/2),digits=3);
  
  head(ldamean)
  
  a = rep("A",length(ldamean$max))
  i = 1
  for (i in 1:length(ldamean$max)) {
    name =colnames(ldamean[,1:class_no])
    a[i] = name[ldamean[,1:class_no][i,] %in% ldamean$max[i]]
  }
  ldamean$class = a
  head(ldamean)
  
  ldamean$Pvalues <- signif(rawpvalues,digits=5);
  ldamean$FDR <- signif(clapvalues,digits=5);
  resTable <- ldamean;
  head(ldamean)
  # it seems lda add ` around names containing dash "-", need to strip this off
  rawNms <- rownames(resTable);
  rownames(resTable) <- gsub("[`]", '', rawNms)
  
  
