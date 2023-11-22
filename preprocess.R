rm(list = ls())
gc()
library(data.table)
# Set seed and working directory
set.seed(1234)
# setwd('/home/chenggg/BIOSTAT625/project')
# Read data
GPL13534 <- fread('data/GPL13534-11288.txt')
# Efficient string splitting and data manipulation
annotated_probe <-
  strsplit(GPL13534$UCSC_RefGene_Name, ';', fixed = TRUE)
names(annotated_probe) <- substr(GPL13534$ID, 1, 10)
# Unlist and create a data frame
annotated_probe <-
  data.frame(id = rep(names(annotated_probe), sapply(annotated_probe, length)),
             gene = unlist(annotated_probe, use.names = FALSE))
# Remove duplicates
annotated_probe <- unique(annotated_probe)
# Optional: Convert to matrix if necessary
annotated_probe <- as.matrix(annotated_probe)
rm(GPL13534)
gc()
data <- read.delim('data/GSE72556_series_matrix.txt', sep = '!')
names(data) <- c("cg", "SampleCharacteristics")

# 定义一个函数来处理和分割数据
process_data <- function(data, pattern, remove_pattern, is_numeric = TRUE) {
  # 提取符合模式的数据
  relevant_data <- data[grepl(pattern, data$SampleCharacteristics), "SampleCharacteristics"]
  # 移除不需要的文本
  relevant_data <- gsub(remove_pattern, "", relevant_data)
  # 额外的清洗步骤，去除制表符等
  relevant_data <- gsub("\\t", "", relevant_data) # 去除制表符
  relevant_data <- gsub("Sample_characteristics_ch1", "", relevant_data) # 去除特定字符串
  # 分割数据
  split_data <- unlist(strsplit(relevant_data, " "))
  # 根据 is_numeric 参数决定是否转换为数值
  if (is_numeric) {
    split_data <- as.numeric(split_data)
    # 移除NA值
    split_data <- split_data[!is.na(split_data)]
  }
  return(split_data)
}

# 应用函数
child_bmi <- process_data(data, "child bmi:", "child bmi:")
adult_bmi <- process_data(data, "adult bmi:", "adult bmi:")
child_waist <- process_data(data, "child waist:", "child waist:")
adult_waist <- process_data(data, "adult waist:", "adult waist:")
gender <- process_data(data, "child gender:", "child gender:", is_numeric = FALSE)[-1]
age <- process_data(data, "child age:", "child age:")

# 处理accession number
accession <- process_data(data, "GSM", "", is_numeric = FALSE)[1:96]
accession[1] <- 'GSM1865141'
accession <- subset(accession, accession != "GSM1865188")
accession <- subset(accession, accession != "GSM1865221")
accession <- subset(accession, accession != "GSM1865231")

# 创建数据框
adult_child <- data.frame(accession, adult_bmi, child_bmi, adult_waist, child_waist, gender, age)

library(data.table)
# 读取数据
y1 <- fread('data/GSE72556_series_matrix.txt', header = FALSE, fill = TRUE)
y2 <- y1
# 移除前65行的注释
y1 <- tail(y2, -65)

# 转置数据
y1 <- transpose(y1)

# 设置第一列为列名,第一行为行名,然后删除
colnames(y1) <- as.character(y1[1, ])
y1 <- y1[-1, ]
