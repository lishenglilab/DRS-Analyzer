#!/usr/bin/env Rscript

####bin/combine_qc_stats.R
library(tidyverse)
library(data.table)

# Read length statistics from pass files
samples <- list.files(".", pattern = "_Length_pass.txt$") %>% 
  gsub("_Length_pass.txt", "", .)
names(samples) <- samples

len <- lapply(samples, function(x){
  fread(paste0(x, "_Length_pass.txt"))
})

# Calculate median lengths
median_len <- lapply(len, function(x) median(x$length))
for (i in names(median_len)) {
  median_len[[i]] <- data.frame(sample_id = i, Median_pass_read_length = median_len[[i]])
}
median_len <- bind_rows(median_len)

# Calculate mean lengths
mean_len <- lapply(len, function(x) mean(x$length))
for (i in names(mean_len)) {
  mean_len[[i]] <- data.frame(sample_id = i, Mean_pass_read_length = mean_len[[i]])
}
mean_len <- bind_rows(mean_len)

# Calculate read numbers
read_num <- lapply(len, nrow)
for (i in names(read_num)) {
  read_num[[i]] <- data.frame(sample_id = i, Total_reads = read_num[[i]])
}
read_num <- bind_rows(read_num)

# Read failed reads
# Read failed reads
if (file.exists("read_number_fail.txt")) {
  # Check if file has content before reading
  if (file.info("read_number_fail.txt")$size > 0) {
    read_num_fail <- read.table("read_number_fail.txt")
    names(read_num_fail) <- c("sample_id", "read_num")
    read_num_fail$read_num <- read_num_fail$read_num / 4
    
    read_num <- left_join(read_num, read_num_fail, by = "sample_id")
    read_num$Total_reads <- read_num$Total_reads + read_num$read_num
    read_num$pass_percent <- round((1 - read_num$read_num/read_num$Total_reads) * 100, 2)
  } else {
    # File exists but is empty - all reads passed
    read_num$pass_percent <- 100
  }
} else {
  read_num$pass_percent <- 100
}

# Combine basic data
dataInfo <- read_num[,c("sample_id", "Total_reads", "pass_percent")]
dataInfo <- left_join(dataInfo, median_len, by = "sample_id")
dataInfo <- left_join(dataInfo, mean_len, by = "sample_id")

# Read file sizes
size_files <- list.files(".", pattern = "_size.txt$")
fileSize <- data.frame()
for (file in size_files) {
  sample <- gsub("_size.txt", "", file)
  size_data <- read.table(file)
  fileSize <- rbind(fileSize, data.frame(sample_id = sample, Size_GB = size_data$V2))
}

dataInfo <- merge(fileSize, dataInfo, by = "sample_id")

# Read mapping statistics
mapping_files <- list.files(".", pattern = "_mapping.txt$")
mapping_read <- data.frame()
for (file in mapping_files) {
  mapping_data <- read.table(file)
  mapping_read <- rbind(mapping_read, data.frame(sample_id = mapping_data$V1, Mapped_read = mapping_data$V2))
}

dataInfo <- merge(dataInfo, mapping_read, by = "sample_id")
dataInfo$Mapped_ratio <- round((dataInfo$Mapped_read / (dataInfo$Total_reads * dataInfo$pass_percent / 100)) * 100, 2)

# Write output
write.table(dataInfo, "QC.xls", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")