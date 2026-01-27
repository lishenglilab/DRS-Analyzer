####bin/expAnalysis/expAnalysis.R
library(getopt)
library(tidyverse)

spec <- matrix(
  c("pepQuant","p",1,"character","input path of protein quant file",
    "transcriptQuant","t",1,"character","input path of transcript quant file",
    "cor","c",1,"character","correlation method,default will be pearson",
    "manifest","m",2,"character","manifest of peptide and transcript",
    "outpath","o",2,"character","outpath, default will be the path of quant file",
    "help","h",0,"logical","This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$pepQuant)  || is.null(opt$transcriptQuant) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if (is.null(opt$outpath)) {
  opt$outpath <- strsplit(opt$quant,"/")[[1]] %>% 
    .[-(length(.))] %>% paste(.,sep = "",collapse = "/")
}

if (is.null(opt$cor)) {
  corMethod <- "pearson"
} else {
  corMethod <- opt$cor
}

library(plyr)
library(ComplexHeatmap)
library(ggsci)
library(Rtsne)
library(Hmisc)

pepQuantPath <- opt$pepQuant
transcriptQuantPath <- opt$transcriptQuant
outpath <- opt$outpath
if (!is.null(opt$manifest) & nchar(opt$manifest) > 1) {
  manifest <- read.table(opt$manifest,header = F,quote = "",sep = "\t")
}

pepQuant <- read.table(pepQuantPath,header = T,row.names = 1,sep = "\t",check.names = F)

#heatmap
Pep_freq <- rowSums(!is.na(pepQuant)) %>% as.data.frame()

names(Pep_freq) <- "Freq"
Pep_freq$newId <- rownames(Pep_freq)

Pep_unique <- pepQuant[rownames(pepQuant) %in% Pep_freq$newId[Pep_freq$Freq == 1],] %>%
  mutate(newId = rownames(.)) %>% 
  reshape2::melt(.,id = "newId") %>% 
  na.omit() %>% 
  group_by(variable) %>% 
  top_n(3,value) %>% 
  .$newId
Pep_half <- sample(Pep_freq$newId[Pep_freq$Freq == ncol(pepQuant)/2 ],5)
pep_tot <- sample(Pep_freq$newId[Pep_freq$Freq == ncol(pepQuant)],5)

ht_data <- pepQuant[c(pep_tot,Pep_half,Pep_unique),]
ht_data <- as.matrix(ht_data)
ht_data[is.na(ht_data)] <- 0

col_fun <- circlize::colorRamp2(c(0, max(ht_data)), c("white", "red"))

pdf(paste0(outpath,"/PepExp.pdf"),height = 16,width = 10)
ComplexHeatmap::Heatmap(ht_data,name = "mat",col = col_fun,
                        column_names_rot = 45,
                        cluster_rows = F,
                        cluster_columns = F,
                        row_names_side = "left",
                        row_names_max_width = unit(20, "cm"))
dev.off()

write.table(ht_data,paste0(outpath,"/htData.xls"),
            col.names = T,row.names = T,quote = F,sep = "\t")

#dimension reduction
pepExp <- t(pepQuant)
pepExp[is.na(pepExp)] <- quantile(pepExp[!(is.na(pepExp))],0.01)

pepExp_scale <- apply(pepExp, 2, function(x) (x - mean(x)) / sd(x) ^ as.logical(sd(x)))
pepExp_scale <- unique(pepExp_scale)

cellLine <- rownames(pepExp)

pca_out <- prcomp(pepExp_scale,center = T,scale. = F)

if (nrow(pepExp_scale) < 80) {
  perplexity <- nrow(pepExp_scale)/4
}else{
  perplexity <- 20
}

tsne_out = Rtsne(
  pepExp_scale,
  dims = 2,
  pca = T,
  max_iter = 1000,
  theta = 0.5,
  perplexity = perplexity,
  verbose = T
)

umap_out <- umap::umap(pepExp_scale)

saveRDS(pca_out,paste0(outpath,"/pca_out.RDS"))
saveRDS(tsne_out,paste0(outpath,"/tsne_out.RDS"))
saveRDS(umap_out,paste0(outpath,"/umap_out.RDS"))

pca_result <- as.data.frame(pca_out$x)
pca_result$sample <- rownames(pepExp_scale)
write.table(pca_result[,c("sample","PC1","PC2")],
            paste0(outpath,"/pca.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
tsne_result$sample <- rownames(pepExp_scale)
write.table(tsne_result[,c("sample","tSNE1","tSNE2")],
            paste0(outpath,"/tsne.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

umap_result <- data.frame(umap_out$layout,tissue = rownames(pepExp_scale))
names(umap_result)[1:2] <- c("umap_1","umap_2")
umap_result$sample <- rownames(pepExp_scale)
write.table(umap_result[,c("sample","umap_1","umap_2")],
            paste0(outpath,"/umap.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

pdf(paste0(outpath,"/pep_dimension_reduction.pdf"),height = 6,width = 10)
ggplot(pca_result,aes(PC1,PC2,color=sample)) +
  geom_point()+
  scale_color_manual(values = c(c(pal_d3("category20")(20),pal_d3("category20b")(8)))) + 
  theme_classic()
ggplot(tsne_result,aes(tSNE1,tSNE2,color=sample)) +
  geom_point()+
  scale_color_manual(values = c(c(pal_d3("category20")(20),pal_d3("category20b")(8)))) + 
  theme_classic()
ggplot(umap_result,aes(umap_1,umap_2,color=sample)) +
  geom_point()+
  scale_color_manual(values = c(c(pal_d3("category20")(20),pal_d3("category20b")(8)))) + 
  theme_classic()
dev.off()

#correlation of protein isoforms and transcript
if (!"manifest" %in% ls()) {
  manifest <- data.frame(V1 = rownames(pepQuant),
                         V2 = substr(rownames(pepQuant),
                                     regexpr("[|]",rownames(pepQuant)) + 1,
                                     regexpr(":",rownames(pepQuant)) - 1))
}

transcriptQuant <- read.table(transcriptQuantPath,header = T,row.names = 1,sep = "\t",check.names = F)

correlation <- apply(manifest, 1, function(x){
  Exp <- rbind(pepQuant[x[[1]],,drop=F],
               transcriptQuant[x[[2]],,drop = F])
  Exp <- t(Exp)
  Exp[is.na(Exp)] <- 0
  cors <- rcorr(Exp,type = corMethod)
  res <- data.frame(pep = x[[1]],
                    transcript=x[[2]],
                    cor=cors[["r"]][1,2],
                    pval=cors[["P"]][1,2])
  return(res)
})

correlation <- bind_rows(correlation)

correlation <- correlation[order(correlation$pval),]
correlation$padj <- p.adjust(correlation$pval,method = "BH")

write.table(correlation,paste0(outpath,"/correlation.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

correlation$p_val <- cut(correlation$pval,breaks = c(-Inf,0.01,0.05,Inf),labels = c("<0.01","0.01-0.05",">0.05"))
correlation$p_val <- factor(correlation$p_val,levels = c(">0.05","0.01-0.05","<0.01"))

pdf(paste0(outpath,"/correlationDot.pdf"),height = 10,width = 10)
ggplot(data = correlation[1:50,]) + 
  geom_point(mapping = aes(x = pep,y = transcript,color = cor,size = p_val)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_classic() + theme(axis.title = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
dev.off()
