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
data <- read.delim('data/GSE72556_series_matrix.txt')
child_bmi <- data[709:757, ]
child_bmi <-
  c(
    child_bmi$X.Series_title,
    child_bmi$Maternal.Weight.Phenotype.as.a.Predictor.of.Methylation.of.Obesity.Related.CpG.Sites.in.Saliva.Samples.from.Latino.Children
  )
child_bmi <- child_bmi[grepl("child bmi:.", child_bmi)]
child_bmi <- gsub("child bmi:", "", child_bmi)
child_bmi = as.numeric(child_bmi)
child_bmi <- child_bmi[complete.cases(child_bmi)]

adult_bmi <- data[517:565, ]
adult_bmi <-
  c(
    adult_bmi$X.Series_title,
    adult_bmi$Maternal.Weight.Phenotype.as.a.Predictor.of.Methylation.of.Obesity.Related.CpG.Sites.in.Saliva.Samples.from.Latino.Children
  )
adult_bmi <- adult_bmi[grepl("adult bmi:.", adult_bmi)]
adult_bmi <- gsub("adult bmi:", "", adult_bmi)
adult_bmi = as.numeric(adult_bmi)
adult_bmi <- adult_bmi[complete.cases(adult_bmi)]

child_waist <- data[758:806, ]
child_waist <-
  c(
    child_waist$X.Series_title,
    child_waist$Maternal.Weight.Phenotype.as.a.Predictor.of.Methylation.of.Obesity.Related.CpG.Sites.in.Saliva.Samples.from.Latino.Children
  )
child_waist <- child_waist[grepl("child waist:.", child_waist)]
child_waist <- gsub("child waist:", "", child_waist)
child_waist <- as.numeric(child_waist)
child_waist <- child_waist[complete.cases(child_waist)]

adult_waist <- data[566:614, ]
adult_waist <-
  c(
    adult_waist$X.Series_title,
    adult_waist$Maternal.Weight.Phenotype.as.a.Predictor.of.Methylation.of.Obesity.Related.CpG.Sites.in.Saliva.Samples.from.Latino.Children
  )
adult_waist <- adult_waist[grepl("adult waist:.", adult_waist)]
adult_waist <- gsub("adult waist:", "", adult_waist)
adult_waist <- as.numeric(adult_waist)
adult_waist <- adult_waist[complete.cases(adult_waist)]

gender <- data[615:661, ]
gender <-
  c(
    gender$X.Series_title,
    gender$Maternal.Weight.Phenotype.as.a.Predictor.of.Methylation.of.Obesity.Related.CpG.Sites.in.Saliva.Samples.from.Latino.Children
  )
gender <- gender[grepl("child gender:.", gender)]
gender <- gsub("child gender: ", "", gender)

age <- data[622:708, ]
age <-
  c(
    age$X.Series_title,
    age$Maternal.Weight.Phenotype.as.a.Predictor.of.Methylation.of.Obesity.Related.CpG.Sites.in.Saliva.Samples.from.Latino.Children
  )
age <- age[grepl("child age:.", age)]
age <- gsub("child age:", "", age)
age <- as.numeric(age)

accession <- data[76:124, ]
accession <-
  c(
    accession$X.Series_title,
    accession$Maternal.Weight.Phenotype.as.a.Predictor.of.Methylation.of.Obesity.Related.CpG.Sites.in.Saliva.Samples.from.Latino.Children
  )
accession <- accession[grepl("GSM.", accession)]
accession <- subset(accession, accession != "GSM1865188")
accession <- subset(accession, accession != "GSM1865221")
accession <- subset(accession, accession != "GSM1865231")
adult_child <-
  data.frame(accession,
             adult_bmi,
             child_bmi,
             adult_waist,
             child_waist,
             gender,
             age)

data <- read.delim('data/GSE72556_series_matrix.txt', comment.char = '!')
y1 <- as.matrix(t(data))
rm(data)
gc()
colnames(y1) <- y1[1, ]
y1 <- y1[-1, ]
y1 <- subset(y1, rownames(y1) != "GSM1865188")
y1 <- subset(y1, rownames(y1) != "GSM1865221")
y1 <- subset(y1, rownames(y1) != "GSM1865231")
adult_child <-
  subset(adult_child, adult_child$accession != "GSM1865172")
y1 <- subset(y1, rownames(y1) != "GSM1865172")
rownames_y1 <- rownames(y1)
y1 <- apply(y1, 2, function(x)
  as.numeric(as.character(x)))
rownames(y1) <- rownames_y1
# library(pheatmap)
# a <- t(y1)
# distance_matrix <- dist(y1)
# clustering <- hclust(distance_matrix)
# pheatmap(a, cluster_rows = clustering)

#not overweight or fat child
nchild <-
  rbind(adult_child[which(
    adult_child$child_bmi < 17.69 &
      adult_child$child_bmi > 0 &
      adult_child$age == 3 & adult_child$gender == "M"
  ), ],
  #36.5
  adult_child[which(
    adult_child$child_bmi < 17.47 &
      adult_child$child_bmi > 0 &
      adult_child$age == 4 & adult_child$gender == "M"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 17.45 &
      adult_child$child_bmi > 0 &
      adult_child$age == 5 & adult_child$gender == "M"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 17.20 &
      adult_child$child_bmi > 0 &
      adult_child$age == 5 & adult_child$gender == "F"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 17.19 &
      adult_child$child_bmi > 0 &
      adult_child$age == 4 & adult_child$gender == "F"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 17.40 &
      adult_child$child_bmi > 0 &
      adult_child$age == 3 & adult_child$gender == "F"
  ), ])


#overweight child
fchild <-
  rbind(adult_child[which(
    adult_child$child_bmi < 19.39 &
      adult_child$child_bmi > 17.69 &
      adult_child$age == 3 & adult_child$gender == "M"
  ), ],
  #36.5
  adult_child[which(
    adult_child$child_bmi < 19.26 &
      adult_child$child_bmi > 17.47 &
      adult_child$age == 4 & adult_child$gender == "M"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 19.47 &
      adult_child$child_bmi > 17.45 &
      adult_child$age == 5 & adult_child$gender == "M"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 19.34 &
      adult_child$child_bmi > 17.20 &
      adult_child$age == 5 & adult_child$gender == "F"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 19.12 &
      adult_child$child_bmi > 17.19 &
      adult_child$age == 4 & adult_child$gender == "F"
  ), ],
  adult_child[which(
    adult_child$child_bmi < 19.23 &
      adult_child$child_bmi > 17.40 &
      adult_child$age == 3 & adult_child$gender == "F"
  ), ])

y1 <- y1[, which(colSums(is.na(y1)) <= 92 * 0.2)]
library(impute)
y1 <- impute.knn(y1)
y1 <- unlist(y1[["data"]])
y1 <- y1[, intersect(colnames(y1), annotated_probe[, 1])]

normal_child <- y1[nchild$accession, ]
obesity_child <- y1[fchild$accession, ]

all <-
  rbind(cbind("normal_child", normal_child),
        cbind("obesity_child", obesity_child))
all <-
  cbind(data.frame(all[, 1]), data.frame(apply(all[, -1], 2, as.numeric)))
colnames(all)[1] <- 'obesity_status'
dim(all)
all$obesity_status <-
  ifelse(all$obesity_status == "obesity_child", 1, 0)
all <-
  cbind(data.frame(all[, 1]), data.frame(apply(all[, -1], 2, as.numeric)))
colnames(all)[1] <- 'obesity_status'
# all[, -1] <- scale(all[, -1])

library(glmnet)
X <- as.matrix(all[, -1])
X <- scale(X)
y <- as.factor(all$obesity_status)
set.seed(123)
# 使用交叉验证进行弹性网回归
# alpha值在0（岭回归）和1（LASSO回归）之???
# 这里我们尝试一个中间值，例如0.5
cv_fit <- cv.glmnet(X, y, alpha = 0.3, family = "binomial")

# 查看最佳的lambda???
best_lambda <- cv_fit$lambda.min

# 拟合最终模???
final_model <-
  glmnet(X,
         y,
         alpha = 0.5,
         lambda = best_lambda,
         family = "binomial")

# 提取模型系数
coefficients <- coef(final_model, s = best_lambda)

# 查看选入模型的变???
# 系数不为零的变量被认为是选入模型???
# 提取非零系数及其对应的变量名
non_zero_coefficients <-
  coefficients[coefficients[, 1] != 0, , drop = FALSE]
variable_names <- rownames(non_zero_coefficients)

# 创建数据???
selected_genes_df <- data.frame(Gene = variable_names,
                                Coefficient = non_zero_coefficients[, 1])
dim(selected_genes_df)
# print(selected_genes_df)

# First, remove the intercept item from selected_genes_df
selected_genes_df <-
  selected_genes_df[selected_genes_df$Gene != "(Intercept)",]

# merge
merged_data <-
  merge(selected_genes_df,
        annotated_probe,
        by.x = "Gene",
        by.y = "id")

print(merged_data)

# enrichment analysis
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
genes <- merged_data$gene
# 将基因名转换为Entrez ID
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

ego <- enrichGO(
  gene         = entrez_ids,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  # "BP"(biological process), "CC"(cell components)???"MF"(molecular function)
  ont          = "MF",
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable     = TRUE
)
barplot(ego)


kk <- enrichKEGG(gene = entrez_ids,
                 organism = 'hsa',  # 人类的KEGG代码
                 pAdjustMethod = "BH",  # 调整p值的方法
                 pvalueCutoff = 0.1)  # 显著性阈???
barplot(kk)
