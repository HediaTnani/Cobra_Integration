# author : Hédia Tnani
# Step1. Mean analysis microarray
require(GEOquery)
require(oligo)
require(mogene10sttranscriptcluster.db)
require(RCurl)
require(foreign)
require(tidyverse)
require(affyPLM)
library(genefilter)
library(affycoretools)
library(mogene10sttranscriptcluster.db)
library(annotate)
setwd("~/")
celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
targets <- getURL(url)                
my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",") 
names <- my.targets@data$ShortName[1:3]
rawData <- read.celfiles(celFiles, phenoData = my.targets)
eset_rma <- oligo::rma(rawData)
eset <- getMainProbes(eset_rma)
eset <- annotateEset(eset, mogene10sttranscriptcluster.db)
annotation(eset) <- "mogene10sttranscriptcluster.db"
library(genefilter)
eset <- featureFilter(eset, require.entrez=TRUE, require.GOBP=FALSE, require.GOCC=FALSE, require.GOMF=FALSE, require.CytoBand=FALSE, remove.dupEntrez=TRUE, feature.exclude="^AFFX")
#filtered = nsFilter(eset, require.entrez=FALSE, remove.dupEntrez=FALSE)
eset_filtered <-filtered$eset 
ID <- featureNames(eset_filtered)
exp <- data.frame(exprs(eset_filtered))
expb <- exp[,1:3]
expb <- as.data.frame(expb)
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mogene10sttranscriptcluster.db")
Name <- as.character(lookUp(ID, "mogene10sttranscriptcluster.db", "GENENAME"))
Entrez <- as.character(lookUp(ID, "mogene10sttranscriptcluster.db", "ENTREZID"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name,Entrez = Entrez, stringsAsFactors=F)
tmp1 <- select(mogene10sttranscriptcluster.db, ID, 
              c("SYMBOL","GENENAME","ENTREZID"))
write.csv(expb, "ExpressionValues.csv")
GeneExpr <- read.csv("~/ExpressionValues.csv")
colnames(GeneExpr) <- c("affyID", names)
mean_expr=data.frame(IDs=GeneExpr[,1],Value=apply(X=GeneExpr[,2:4],MARGIN = 1,FUN=mean))
mean_expr$IDs <- as.character(mean_expr$IDs)
tmp1$ENTREZID <- as.character(tmp1$ENTREZID)
mean_expr_entrz <- mean_expr %>% inner_join(., tmp1, by=c("IDs"= "PROBEID")) %>% dplyr::select(!c(SYMBOL,GENENAME))
mean_expr_entrzf <- mean_expr_entrz %>% dplyr::select(!IDs) %>% relocate(ENTREZID, Value)
write.csv(mean_expr_entrzf, "MeanExpressionValues_entrez.csv")
gene_model <- read.csv("genes.csv", h=F)
gene_model$V1 %in% tmp1$ENTREZID
mean_expr_entrzf$ENTREZID %in% gene_model$V1
##########################################################################################
# Step2 Gene annotation
library(biomaRt)
ids = as.character(mean_expr[,1])
httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart.current <- useMart(biomart="ensembl",dataset="mmusculus_gene_ensembl",host="http://uswest.ensembl.org")
att=listAttributes(mart.current)
current.results <- getBM(attributes=c("affy_mogene_1_0_st_v1", "entrezgene_id","ensembl_gene_id"),filter="affy_mogene_1_0_st_v1",values=ids, mart=mart.current)
current.results1 <- getBM(attributes=c("affy_mogene_1_0_st_v1", "entrezgene_id"),filter="affy_mogene_1_0_st_v1",values=ids, mart=mart.current)
ensentrzaffy=inner_join(mean_expr,current.results,by = c("IDs"="affy_mogene_1_0_st_v1"))
ensentrzaffy %>% as.tibble() 
ids_map=ensentrzaffy[!duplicated(ensentrzaffy$IDs), ]
ensentrzaffy %>% as.tibble() %>% distinct(IDs, .keep_all = TRUE) %>% drop_na(entrezgene_id) %>% dplyr::select(entrezgene_id, Value) %>% write_csv("NormalizedAbsoluteExpression.csv")
