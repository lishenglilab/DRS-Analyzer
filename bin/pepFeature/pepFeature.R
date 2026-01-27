####bin/pepFeature/pepFeature.R
library(getopt)
library(tidyverse)

spec <- matrix(
  c("quant","q",1,"character","input path of quant file",
    "pepRef","p",1,"character","input path of peptide reference fasta file",
    "transcriptRef","t",1,"character","input path of transcript reference fasta file",
    "outpath","o",2,"character","outpath, default will be the path of quant file",
    "help","h",0,"logical","This is Help!"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$quant)  || is.null(opt$pepRef)   || is.null(opt$transcriptRef) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if (is.null(opt$outpath)) {
  opt$outpath <- strsplit(opt$quant,"/")[[1]] %>% 
    .[-(length(.))] %>% paste(.,sep = "",collapse = "/")
}

quantPath <- opt$quant
pepRef <- opt$pepRef
transcriptRef <- opt$transcriptRef
outpath <- opt$outpath

library(Biostrings)
library(stringr)
library(ggsci)

peptide <- readBStringSet(pepRef,format = "fasta",nrec = -1L,skip=0L, seek.first.rec=TRUE, 
                          use.names=TRUE) %>% as.data.frame()
names(peptide) <- "pepSeq"
peptide$pepID <- rownames(peptide)
peptide$transcriptID <- substr(peptide$pepID,regexpr("[|]",peptide$pepID) + 1,regexpr(":",peptide$pepID) - 1)
cutpoint <- gregexpr(":",rownames(peptide))
peptide$orfStart <- substr(rownames(peptide),
                           unlist(lapply(cutpoint,function(x)x[1])) + 1,
                           unlist(lapply(cutpoint,function(x)x[2])) - 1) %>% as.numeric()
peptide$orfStart <- peptide$orfStart + 1
peptide$orfEnd <- substring(rownames(peptide),
                            unlist(lapply(cutpoint,function(x)x[2])) + 1) %>% as.numeric()

transcript <- readBStringSet(transcriptRef,format = "fasta",nrec = -1L,skip=0L, seek.first.rec=TRUE, 
                             use.names=TRUE) %>% as.data.frame()
names(transcript) <- "rnaSeq"
transcript$transcriptID <- rownames(transcript)

peptide <- left_join(peptide,transcript,by = "transcriptID")
peptide$rnaSeq <- substr(peptide$rnaSeq,peptide$orfStart,peptide$orfEnd)
peptide$iniCodon <- substr(peptide$rnaSeq,1,3)
peptide$orfLen <- nchar(peptide$rnaSeq)
peptide$GC <- str_count(peptide$rnaSeq,"G|C")/nchar(peptide$rnaSeq)

quant <- read.table(quantPath,header = T,row.names = 1,quote = "",check.names = F,sep = "\t")
recognized_pep_info <- peptide[peptide$pepID %in% rownames(quant),]

write.table(peptide,paste0(outpath,"/reference_info.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")
write.table(recognized_pep_info,paste0(outpath,"/recognized_pep_info.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

#initiation codon stat
ini_codon <- as.data.frame(table(peptide$iniCodon))
names(ini_codon) <- c("ini_codon","Freq")
ini_codon$ratio <- ini_codon$Freq / sum(ini_codon$Freq)

write.table(ini_codon,paste0(outpath,"/ini_codon_ratio.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

myLabel <- as.vector(ini_codon$ini_codon)
myLabel <- paste0(myLabel, "(",ini_codon$Freq," ",round(ini_codon$ratio * 100, 2), "%)") 

pdf(paste0(outpath,"/ini_codon_ratio.pdf"),height = 4,width = 7)
ggplot(data = ini_codon, mapping = aes(x = 'Content', y = Freq, fill = ini_codon)) + 
  geom_bar(stat = 'identity', position = 'stack',width = 1) + 
  coord_polar(theta = 'y')+theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        line = element_blank())+scale_fill_d3(labels = myLabel)
dev.off()

#ORF length stat
ORF_length <- cut(peptide$orfLen,
                  breaks = c(-Inf,50,100,150,200,250,Inf),
                  labels = c("0-50","50-100","100-150","150-200","200-250",">250")) %>% 
  table() %>% as.data.frame()
names(ORF_length) <- c("ORF_length","ORF_counts")

write.table(ORF_length,paste0(outpath,"/ORF_length_distribution.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

pdf(paste0(outpath,"/ORF_length_distribution.pdf"),height = 4,width = 7)
ggplot(data = ORF_length,mapping = aes(x=ORF_length,y = (ORF_counts)/1000,fill=ORF_length))+
  geom_col() + scale_fill_nejm()+theme_classic()+
  geom_text(aes(label=ORF_counts), 
            position = position_dodge2(width = 0.9, preserve = 'single'), 
            vjust = -0.2, hjust = 0.5)+xlab("ORF length (nt)")
dev.off()

#GC content
pdf(paste0(outpath,"/GC_ratio.pdf"),height = 6,width = 12)
ggplot(peptide,aes(x=GC))+
  geom_density(fill="#69b3a2", 
               color="#e9ecef", 
               alpha=0.8)+
  theme_classic()
dev.off()

#amino acid stat
amino_stat <- str_split(peptide$pepSeq,"") %>% 
  unlist() %>% 
  table() %>% 
  as.data.frame()

names(amino_stat)[1] <- "AA"
amino_stat$ratio <- (amino_stat$Freq / sum(amino_stat$Freq))*100

write.table(amino_stat,paste0(outpath,"/amino_acid_stat.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

pdf(paste0(outpath,"/amino_acid_ratio.pdf"),height = 4,width = 7)
ggplot()+
  geom_bar(mapping = aes(x = AA, y = ratio),
           data = amino_stat,
           stat = "identity")+
  theme_classic() + 
  ylab("amino acid ratio (%)") + 
  xlab(NULL)
dev.off()

#recognized peptide initiation codon stat
ini_codon <- as.data.frame(table(recognized_pep_info$iniCodon))
names(ini_codon) <- c("ini_codon","Freq")
ini_codon$ratio <- ini_codon$Freq / sum(ini_codon$Freq)

write.table(ini_codon,paste0(outpath,"/recognized_ini_codon_ratio.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

myLabel <- as.vector(ini_codon$ini_codon)
myLabel <- paste0(myLabel, "(",ini_codon$Freq," ",round(ini_codon$ratio * 100, 2), "%)") 

pdf(paste0(outpath,"/recognized_ini_codon_ratio.pdf"),height = 4,width = 7)
ggplot(data = ini_codon, mapping = aes(x = 'Content', y = Freq, fill = ini_codon)) + 
  geom_bar(stat = 'identity', position = 'stack',width = 1) + 
  coord_polar(theta = 'y')+theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        line = element_blank())+scale_fill_d3(labels = myLabel)
dev.off()

#recognized ORF length stat
ORF_length <- cut(recognized_pep_info$orfLen,
                  breaks = c(-Inf,50,100,150,200,250,Inf),
                  labels = c("0-50","50-100","100-150","150-200","200-250",">250")) %>% 
  table() %>% as.data.frame()
names(ORF_length) <- c("ORF_length","ORF_counts")

write.table(ORF_length,paste0(outpath,"/recognized_ORF_length_distribution.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

pdf(paste0(outpath,"/recognized_ORF_length_distribution.pdf"),height = 4,width = 7)
ggplot(data = ORF_length,mapping = aes(x=ORF_length,y = (ORF_counts)/1000,fill=ORF_length))+
  geom_col() + scale_fill_nejm()+theme_classic()+
  geom_text(aes(label=ORF_counts), 
            position = position_dodge2(width = 0.9, preserve = 'single'), 
            vjust = -0.2, hjust = 0.5)+xlab("ORF length (nt)")
dev.off()

#recognized peptide GC content
pdf(paste0(outpath,"/recognized_GC_ratio.pdf"),height = 6,width = 12)
ggplot(recognized_pep_info,aes(x=GC))+
  geom_density(fill="#69b3a2", 
               color="#e9ecef", 
               alpha=0.8)+
  theme_classic()
dev.off()

#recognized peptide amino acid stat
amino_stat <- str_split(recognized_pep_info$pepSeq,"") %>% 
  unlist() %>% 
  table() %>% 
  as.data.frame()

names(amino_stat)[1] <- "AA"
amino_stat$ratio <- (amino_stat$Freq / sum(amino_stat$Freq))*100

write.table(amino_stat,paste0(outpath,"/recognized_amino_acid_stat.xls"),
            col.names = T,row.names = F,quote = F,sep = "\t")

pdf(paste0(outpath,"/recognized_amino_acid_ratio.pdf"),height = 4,width = 7)
ggplot()+
  geom_bar(mapping = aes(x = AA, y = ratio),
           data = amino_stat,
           stat = "identity")+
  theme_classic() + 
  ylab("amino acid ratio (%)") + 
  xlab(NULL)
dev.off()
