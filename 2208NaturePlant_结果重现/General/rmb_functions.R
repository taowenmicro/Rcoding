library(tidyverse)

read_tax <- function(x) {
  read.table(x, sep = "\t", header = TRUE, quote = "", colClasses = "character")
}

read_otu <- function(x) {
  read.table(x, sep = "\t", header = TRUE, row.names = 1)
}

read_map <- function(x) {
  read.table(x, sep = "\t", header = TRUE)
} 

rel_ab <- function(otu, total = 100) {
  t(t(otu)/colSums(otu)) * 100
}

log_norm <- function(otu) {
  log(otu + 1)
}

classify_organelle <- function(tax) {
  tax %>% 
    mutate(Assignment = case_when(
      Kingdom == "unassigned" ~ "Unassigned",
      Family == "mitochondria" ~ "Mitochondria",
      Class == "Chloroplast" ~ "Chloroplast",
      TRUE ~ "Microbial"
    )) %>% 
    mutate(Assignment = fct_relevel(Assignment, "Unassigned", "Mitochondria", "Chloroplast", "Microbial"))
}

expand_proteo <- function(tax) {
  mutate(tax, PhyClass = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum)))
}

tidy_otu <- function(otu) {
  as.data.frame(otu) %>%
    mutate(OTU_ID = row.names(otu)) %>%
    tidyr::gather(key = "SampleID", value = "Count", -OTU_ID)
}

filter_prevalent <- function(otu, prevalence = 0.05) {
  otu[rowSums(otu > 0) > ncol(otu) * prevalence,]
} 


get_top_taxa <- function(otu, tax, rank = "Phylum", n = 10) {
  data.frame(OTU_ID = row.names(otu), Total = rowSums(otu)) %>%
    inner_join(tax) %>%
    group_by_(rank) %>% # group_by_ is used instead of group_by to allow for standard evaluation semantics
    summarise(RankTotal = sum(Total)) %>%
    arrange(-RankTotal) %>%
    head(n = n) 
}

beta_div_dist <- function(otu, method = "bray") {
  require(phyloseq)
  physeq <- phyloseq(otu_table(otu, taxa_are_rows = TRUE))
  as.matrix(phyloseq::distance(physeq, method))
} 

pcoa_axes <- function(dist, map) {
  require(ape)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- pcoa(as.dist(dist))
  as.data.frame(pcoa$vectors) %>%
    mutate(SampleID = rownames(.)) %>%
    inner_join(map, by = "SampleID")
}

pcoa_eigval <- function (dist, map) {
  require(ape)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- pcoa(as.dist(dist))
  eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)
  data.frame( PC = 1:length(eigval), Eigval = eigval, CumEigval = cumsum(eigval))
}


cap_axes <- function(cap, map) {
  require(vegan)
  as.data.frame(scores(cap, choices = c(1,2,3,4))$sites) %>% # by default, scores gives you only the first two CAPs, by stating the choices you can explore more
    mutate(SampleID = row.names(.)) %>%
    inner_join(map, by = "SampleID")
}

cap_eigval <- function(cap, map){
  require(vegan)
  cap.eigval <- round(cap$CCA$eig/sum(cap$CCA$eig) * 100, digits = 2)
  data.frame(PC = 1:length(cap.eigval), Eigval = cap.eigval)
}

read_rep_set_tax <- function(tax){
  tax.raw <- read.table(tax, 
                        sep = "\t", 
                        col.names = c("OTU_ID", "Taxonomy", "Score", "Hits"))
  tax.raw %>% 
    separate(col = Taxonomy, 
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = "\\; ", 
             remove = FALSE)
}
  
classify_organelle2 <- function(tax) {
  mutate(tax, Assignment = ifelse(grepl("unassigned", Taxonomy, ignore.case = TRUE), "Unassigned",
                                  ifelse(grepl("mitochondria", Taxonomy, ignore.case = TRUE), "Mitochondria",
                                         ifelse(grepl("chloroplast", Taxonomy, ignore.case = TRUE), "Chloroplast", "Microbial"))))
}

expand_proteo2 <- function(tax) {
  mutate(tax, PhyClass = ifelse(Phylum == "p__Proteobacteria", as.character(Class), as.character(Phylum)))
}

collapse_other <- function(tax, top.tax) {
  tax <- mutate(tax, PhyClass2 = as.factor(ifelse(as.character(PhyClass) %in% top.tax$PhyClass, as.character(PhyClass), "other")))
  proteos <- c("Alphaproteobacteria", "Betaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria", "Gammaproteobacteria", "TA18", "Zetaproteobacteria")
  proteos <- proteos[proteos %in% levels(tax$PhyClass2)]
  nonproteos <- levels(tax$PhyClass2)[levels(tax$PhyClass2) != "other" & !levels(tax$PhyClass2) %in% proteos]
  tax$PhyClass2 <- factor(tax$PhyClass2,
                          levels = c("other", nonproteos, proteos))
  tax
}

subset_otu <- function(otu, map, condition) {
  filt_map <- filter(map, condition) 
  filt_otu <- otu[,colnames(otu) %in% map$SampleID]
  filt_otu
}
