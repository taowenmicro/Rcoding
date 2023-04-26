

#------构建phyloseq对象
# 从文件读取
library(tidyverse)
library(phyloseq)
metadata = readxl::read_xlsx("./NC.xlsx",sheet =4)%>% as.data.frame() 
head(metadata)

colnames(metadata)[1:2] = c("ID","Group")
# metadata$ID = gsub("\\.","",metadata$ID)
row.names(metadata) = metadata$ID

dim(metadata)

otutab = readxl::read_xlsx("./NC.xlsx",sheet =3) %>% as.data.frame() %>% column_to_rownames("ID")
head(otutab)
# row.names(otutab) = otutab$ID

match(colnames(otutab),metadata$ID)


taxonomy = readxl::read_xlsx("./NC.xlsx",sheet =2)%>% as.data.frame() %>% column_to_rownames("ID")

head(taxonomy )
colnames(otutab)
# metadata$SampleID = as.character(metadata$SampleID )


ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))
              
)
ps

saveRDS(ps,"./ps_result1.rds")


#---result2#---------

#------构建phyloseq对象
# 从文件读取
library(tidyverse)
library(phyloseq)
metadata = readxl::read_xlsx("./NC.xlsx",sheet = 8)%>% as.data.frame() 
head(metadata)

colnames(metadata)[1:2] = c("ID","Group")
# metadata$ID = gsub("\\.","",metadata$ID)
row.names(metadata) = metadata$ID

dim(metadata)

otutab = readxl::read_xlsx("./NC.xlsx",sheet =7) %>% as.data.frame() %>% column_to_rownames("ID")
head(otutab)
# row.names(otutab) = otutab$ID

match(colnames(otutab),metadata$ID)


taxonomy = readxl::read_xlsx("./NC.xlsx",sheet =6)%>% as.data.frame() %>% column_to_rownames("ID")

head(taxonomy )
colnames(otutab)
# metadata$SampleID = as.character(metadata$SampleID )


ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))
              
)
ps

saveRDS(ps,"./ps_result2.rds")



#---result3#---------

#------构建phyloseq对象
# 从文件读取
library(tidyverse)
library(phyloseq)
metadata = readxl::read_xlsx("./NC.xlsx",sheet =12)%>% as.data.frame() 
head(metadata)

colnames(metadata)[1:2] = c("ID","Group")
# metadata$ID = gsub("\\.","",metadata$ID)
row.names(metadata) = metadata$ID

dim(metadata)

otutab = readxl::read_xlsx("./NC.xlsx",sheet =11) %>% as.data.frame() %>% column_to_rownames("ID")
head(otutab)
# row.names(otutab) = otutab$ID

match(colnames(otutab),metadata$ID)


taxonomy = readxl::read_xlsx("./NC.xlsx",sheet = 10)%>% as.data.frame() %>% column_to_rownames("ID")

head(taxonomy )
colnames(otutab)
# metadata$SampleID = as.character(metadata$SampleID )


ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))
              
)
ps

saveRDS(ps,"./ps_result3.rds")




#---result4#---------

#------构建phyloseq对象
# 从文件读取
library(tidyverse)
library(phyloseq)
metadata = readxl::read_xlsx("./NC.xlsx",sheet =16)%>% as.data.frame() 
head(metadata)

colnames(metadata)[1:2] = c("ID","Group")
# metadata$ID = gsub("\\.","",metadata$ID)
row.names(metadata) = metadata$ID

dim(metadata)

otutab = readxl::read_xlsx("./NC.xlsx",sheet =15) %>% as.data.frame() %>% column_to_rownames("ID")
head(otutab)
# row.names(otutab) = otutab$ID

match(colnames(otutab),metadata$ID)


taxonomy = readxl::read_xlsx("./NC.xlsx",sheet = 14)%>% as.data.frame() %>% column_to_rownames("ID")

head(taxonomy )
colnames(otutab)
# metadata$SampleID = as.character(metadata$SampleID )


ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))
              
)
ps

saveRDS(ps,"./ps_result4.rds")









