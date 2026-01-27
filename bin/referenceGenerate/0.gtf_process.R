####bin/referenceGenerate/0.gtf_process.R
library(getopt)
library(tidyverse)

spec <- matrix(
  c("gtf","g",1,"character","input gtf file's path",
    "outpath","o",2,"character","outpath, default will be the input path",
    "name","n",1,"character","output file prefix, default will be the gtf file's name",
    "help","h",0,"logical","This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$gtf) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if (is.null(opt$outpath)) {
  opt$outpath <- strsplit(opt$gtf,"/")[[1]] %>% 
    .[-(length(.))] %>% paste(.,sep = "",collapse = "/")
}

if (is.null(opt$name)) {
  opt$name <- strsplit(opt$gtf,"/")[[1]] %>% 
    .[length(.)] %>% substr(.,1,regexpr(".gtf",.) - 1)
}

library(rtracklayer)
library(Biostrings)

gtf <- opt$gtf
outpath <- opt$outpath
outname <- opt$name

gencode <- import(gtf) %>% as.data.frame()
protein_coding <- c("IG_C_gene","IG_D_gene","IG_V_gene","TR_C_gene","TR_J_gene",
                    "TR_V_gene","IG_J_gene","protein_coding","TR_D_gene")
gencode.sub <- filter(gencode,transcript_type %in% protein_coding)

export(gencode.sub,paste0(outpath,"/",outname,"_proteinCoding.gtf"),format = "gtf")
