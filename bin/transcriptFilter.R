####bin/transcriptFilter.R
library(getopt)

spec <- matrix(
  c("flairQuant","t",1,"character","Input path of flair TPM quant file.",
    "gtf","a",1,"character","Input path of transcript genome annotation file exported by flair software.",
    "tpm","q",2,"numeric","The minimum TPM value expressed in at least one sample, default will be 0.1.",
    "outpath","o",2,"character","Outpath, default will be the path of genome annotation gtf file",
    "help","h",0,"logical","This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$flairQuant)  || is.null(opt$gtf) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if (is.null(opt$outpath)) {
  opt$outpath <- strsplit(opt$quant,"/")[[1]] %>% 
    .[-(length(.))] %>% paste(.,sep = "",collapse = "/")
}

if (is.null(opt$tpm)) {
  opt$tpm <- 0.1
}

library(tidyverse)

tpm <- read.table(opt$flairQuant,header = T,row.names = 1,quote = "",check.names = F,sep = "\t")
rownames(tpm) <- substr(rownames(tpm),1,regexpr("_",rownames(tpm)) - 1)
tpm <- tpm[rowSums(tpm > opt$tpm) >= 1,]

gtf <- rtracklayer::import(opt$gtf) %>% as.data.frame()
gtf <- gtf[grep("^chr",gtf$seqnames),] %>% filter(transcript_id %in% rownames(tpm))

tpm <- tpm[rownames(tpm) %in% gtf$transcript_id,]

rtracklayer::export(object = gtf,con = paste0(opt$outpath,"/isoforms_filter.gtf"),format = "gtf")
write.table(tpm,paste0(opt$outpath,"/flair_quantify_filter.tpm.tsv"),
            col.names = T,row.names = T,quote = F,sep = "\t")