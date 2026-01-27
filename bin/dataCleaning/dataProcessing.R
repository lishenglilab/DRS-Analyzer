####bin/dataClean/dataProcessing.R
library(getopt)
library(tidyverse)

spec <- matrix(
  c("quant","q",1,"character","input path of diann's quant result",
    "outpath","o",2,"character","outpath, default will be the input path",
    "help","h",0,"logical","This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$quant) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if (is.null(opt$outpath)) {
  opt$outpath <- strsplit(opt$quant,"/")[[1]] %>% 
    .[-(length(.))] %>% paste(.,sep = "",collapse = "/")
}

library(diann)
library(plyr)
quantPath <- opt$quant
outpath <- opt$outpath

df <- diann_load(quantPath)
precursors <- diann_matrix(df, pg.q = 0.01)
peptides <- diann_matrix(df, id.header="Stripped.Sequence", pg.q = 0.01)
peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                                group.header="Stripped.Sequence", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Precursor.Normalised")
unique.genes <- diann_matrix(df, 
                             id.header="Genes", 
                             quantity.header="Genes.MaxLFQ.Unique", 
                             proteotypic.only = T, 
                             pg.q = 0.01)
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                               group.header="Protein.Group", 
                               id.header = "Precursor.Id", 
                               quantity.header = "Precursor.Normalised")
protein.groups <- cbind(pepID=rownames(protein.groups),protein.groups) %>% as.data.frame()
protein.groups.sub <- protein.groups[grep(";",protein.groups$pepID),]
protein.groups.sub <- ddply(protein.groups.sub,"pepID",function(x){
  pepID <- strsplit(x$pepID,";") %>% 
    unlist()
  x <- x[rep(1,length(pepID)),]
  x$pepID <- pepID
  return(x)
})
protein.groups <- protein.groups[-grep(";",protein.groups$pepID),]
protein.groups <- rbind(protein.groups,protein.groups.sub)
protein.groups <- as.matrix(protein.groups)
rownames(protein.groups) <- protein.groups[,1]
protein.groups <- protein.groups[,-1]
protein.groups <- limma::avereps(protein.groups)

write.table(precursors,paste0(outpath,"/precursors.tsv"),col.names = T,row.names = T,quote = F,sep = "\t")
write.table(peptides,paste0(outpath,"/peptides.tsv"),col.names = T,row.names = T,quote = F,sep = "\t")
write.table(peptides.maxlfq,paste0(outpath,"/peptides.maxlfq.tsv"),col.names = T,row.names = T,quote = F,sep = "\t")
write.table(unique.genes,paste0(outpath,"/unique.genes.tsv"),col.names = T,row.names = T,quote = F,sep = "\t")
write.table(protein.groups,paste0(outpath,"/protein.groups.tsv"),col.names = T,row.names = T,quote = F,sep = "\t")
