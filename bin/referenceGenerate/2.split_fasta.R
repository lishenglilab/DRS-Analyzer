####bin/referenceGenerate/2.split_fasta.R
library(getopt)
library(tidyverse)

spec <- matrix(
  c("input","i",2,"character","input fasta file",
    "split_num","s",2,"numeric","input number you wish to split the fasta file",
    "outpath","o",1,"character","outpath, default will be the input path",
    "name","n",1,"character","output file prefix, default will be the input name",
    "help","h",0,"logical","This is Help!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$input) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if (is.null(opt$outpath)) {
  opt$outpath <- strsplit(opt$input,"/")[[1]] %>% 
    .[-(length(.))] %>% paste(.,sep = "",collapse = "/")
}

if (is.null(opt$name)) {
  opt$name <- strsplit(opt$input,"/")[[1]] %>% 
    .[length(.)] %>% substr(.,1,regexpr(".fa",.) - 1)
}

library(Biostrings)

seqdata <- readBStringSet(opt$input,format = "fasta",nrec = -1L,
                          skip=0L, seek.first.rec=TRUE, use.names=TRUE) %>% 
  as.data.frame()

seqdata_list <- list()
j <- 1
for (i in seq(1,nrow(seqdata),ceiling(nrow(seqdata) / opt$split_num))) {
  seqdata_list[[j]] <- seqdata[i:(i + ceiling(nrow(seqdata) / opt$split_num) - 1),,drop = F]
  j <- j + 1
}

seqdata_list[[opt$split_num]] <- na.omit(seqdata_list[[opt$split_num]])

for (i in 1:(opt$split_num)) {
  fa_path <- paste0(opt$outpath,"/",opt$name,"_",i,".fa")
  writeLines(paste0(">",rownames(seqdata_list[[i]]),"\n",seqdata_list[[i]]$x),fa_path)
}
