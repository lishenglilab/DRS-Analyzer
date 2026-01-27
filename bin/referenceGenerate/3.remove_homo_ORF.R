####bin/referenceGenerate/3.remove_homo_ORF.R
library(getopt)

spec <- matrix(
  c("fasta","f",2,"character","input fasta file",
    "blast","b",2,"character","input blast result",
    "outpath","o",1,"character","outpath of filterd fasta file, default will be the input fasta path",
    "name","n",1,"character","name of fasta file, default is 'ORF_unmatch_uniprot.fa'",
    "identity","d",1,"numeric","identity cut off, default is 100",
    "help","h",0,"logical","This is Help!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$fasta) || is.null(opt$blast) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

library(tidyverse)
library(Biostrings)

if (is.null(opt$outpath)) {
  opt$outpath <- dirname(opt$fasta)
}

if (is.null(opt$name)) {
  opt$name <- "ORF_unmatch_uniprot.fa"
}

if (is.null(opt$identity)) {
  opt$identity <- 100
}

identity <- opt$identity

seqdata <- readBStringSet(opt$fasta,format = "fasta",nrec = -1L,
                          skip=0L, seek.first.rec=TRUE, use.names=TRUE) %>% 
  as.data.frame()
seqdata$`query acc.ver` <- rownames(seqdata)
seqdata$ORF_length <- nchar(seqdata$x)

blast <- data.table::fread(opt$blast, header = F, data.table = F, check.names = F)
colnames(blast) <- c("query acc.ver", "subject acc.ver", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
blast$`% identity` <- as.numeric(blast$`% identity`)
blast <- left_join(blast,seqdata[,c("query acc.ver","ORF_length")],by = "query acc.ver")
blast$query_coverage <- (blast$`q. end` - blast$`q. start` + 1)*100/blast$ORF_length

homo <- blast[blast$evalue < 1e-5 & blast$`% identity` >= identity & blast$ORF_length >= identity,]

seqdata <- seqdata[!(rownames(seqdata) %in% homo$`query acc.ver`),,drop = F]

writeLines(paste0(">",rownames(seqdata),"\n",seqdata$x),
           paste0(opt$outpath,"/",opt$name))
