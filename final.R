rm(list = ls())
gc()
library(data.table)
# Set seed and working directory
set.seed(1234)
# setwd('/home/chenggg/BIOSTAT625/project')
start.time <- Sys.time()
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
end.time <- Sys.time()
cat('Time for fread annotated data: ', end.time - start.time, '\n')

rm(GPL13534)
gc()

data <- read.delim('data/GSE72556_series_matrix.txt', sep = '!')
names(data) <- c("cg", "SampleCharacteristics")

# Define a function to process and split the data
process_data <-
  function(data,
           pattern,
           remove_pattern,
           is_numeric = TRUE) {
    # Extract data that matches the pattern
    relevant_data <-
      data[grepl(pattern, data$SampleCharacteristics), "SampleCharacteristics"]
    # Remove unwanted text
    relevant_data <- gsub(remove_pattern, "", relevant_data)
    # Additional cleaning steps, TAB removal, etc
    relevant_data <- gsub("\\t", "", relevant_data) # Remove the TAB
    relevant_data <-
      gsub("Sample_characteristics_ch1", "", relevant_data) # Remove specific string
    # Split data
    split_data <- unlist(strsplit(relevant_data, " "))
    # Determine whether to convert to a numeric value based on the is numeric parameter
    if (is_numeric) {
      split_data <- as.numeric(split_data)
      # remove NA
      split_data <- split_data[!is.na(split_data)]
    }
    return(split_data)
  }

# apply the function
child_bmi <- process_data(data, "child bmi:", "child bmi:")
adult_bmi <- process_data(data, "adult bmi:", "adult bmi:")
child_waist <- process_data(data, "child waist:", "child waist:")
adult_waist <- process_data(data, "adult waist:", "adult waist:")
gender <-
  process_data(data, "child gender:", "child gender:", is_numeric = FALSE)[-1]
age <- process_data(data, "child age:", "child age:")

# Processing accession number
accession <- process_data(data, "GSM", "", is_numeric = FALSE)[1:96]
accession[1] <- 'GSM1865141'
accession <- subset(accession, accession != "GSM1865188")
accession <- subset(accession, accession != "GSM1865221")
accession <- subset(accession, accession != "GSM1865231")

# Create data frame
adult_child <-
  data.frame(accession,
             adult_bmi,
             child_bmi,
             adult_waist,
             child_waist,
             gender,
             age)

data <-
  read.delim('data/GSE72556_series_matrix.txt', comment.char = '!')
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
temptForCompareInpute = y1
tempt = y1
library(impute)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(glmnet)
library(missForest)
#### KNN imputation
y1 <- impute.knn(y1)
##################################################### compare for knn and weighted-knn

#sum(is.na(temptForCompareInpute))
k_values <- 1:15
resultsForCompare <- numeric(length(k_values))
for (k in k_values) {
  tempt <- temptForCompareInpute
  tempt <- impute.knn(tempt, k = k)
  tempt <- unlist(tempt[["data"]])
  tempt <- tempt[, intersect(colnames(tempt), annotated_probe[, 1])]
  normal_child <- tempt[nchild$accession, ]
  obesity_child <- tempt[fchild$accession, ]
  normal_child_labeled <- cbind(0, normal_child)
  obesity_child_labeled <- cbind(1, obesity_child)
  all <- rbind(normal_child_labeled, obesity_child_labeled)
  all <-
    cbind(data.frame(all[, 1]), data.frame(apply(all[, -1], 2, as.numeric)))
  colnames(all)[1] <- 'obesity_status'
  X <- as.matrix(all[, -1])
  X <- scale(X)
  y <- as.factor(all$obesity_status)
  cv_fit <-
    cv.glmnet(X, y, alpha = 0.8, family = "binomial", set.seed(1234))
  best_lambda <- cv_fit$lambda.min
  final_model <-
    glmnet(
      X,
      y,
      alpha = 0.8,
      lambda = best_lambda,
      family = "binomial",
      set.seed(1234)
    )
  coefficients <- coef(final_model, s = best_lambda)
  non_zero_coefficients <-
    coefficients[coefficients[, 1] != 0, , drop = FALSE]
  variable_names <- rownames(non_zero_coefficients)
  selected_genes_df <- data.frame(Gene = variable_names,
                                  Coefficient = non_zero_coefficients[, 1])
  resultsForCompare[k] <- nrow(selected_genes_df)
  print(nrow(selected_genes_df))
}

library(ggplot2)

k_values <- 1:15
# best k value,generate from the previous code
resultsForCompare <-
  c(39, 42, 42, 43, 43, 41, 41, 41, 43, 43, 43, 41, 41, 41, 41)

data <- data.frame(K = k_values, Genes = resultsForCompare)
#plot the result for the best k-values
ggplot(data, aes(x = K, y = resultsForCompare)) +
  geom_line(color = 'blue', size = 1) +
  geom_point(color = 'black', size = 3) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = 'aliceblue', color = NA),
    panel.grid.major = element_line(color = 'gray', size = 0.5),
    panel.grid.minor = element_line(color = 'lightgray', size = 0.25)
  ) +
  labs(
    title = 'Number of Genes Selected by Elastic Net for Different K Values',
    subtitle = 'Visual representation of gene count variability',
    x = 'K Value',
    y = 'Number of Genes Selected'
  )

######
################################################
yy = y1[1:10, 1:20]
y1 <- unlist(y1[["data"]])
y1 <- y1[, intersect(colnames(y1), annotated_probe[, 1])]

normal_child <- y1[nchild$accession, ]
obesity_child <- y1[fchild$accession, ]

# 'normal_child' labeled 0
normal_child_labeled <- cbind(0, normal_child)

# 'obesity_child' labeled 1
obesity_child_labeled <- cbind(1, obesity_child)


all <- rbind(normal_child_labeled, obesity_child_labeled)

all <-
  cbind(data.frame(all[, 1]), data.frame(apply(all[, -1], 2, as.numeric)))

# change column name to 'obesity_status'
colnames(all)[1] <- 'obesity_status'
# all <-
#   rbind(cbind("normal_child", normal_child),
#         cbind("obesity_child", obesity_child))
# all <-
#   cbind(data.frame(all[, 1]), data.frame(apply(all[, -1], 2, as.numeric)))
# colnames(all)[1] <- 'obesity_status'
#dim(all)
# all$obesity_status <-
#   ifelse(all$obesity_status == "obesity_child", 1, 0)
# all <-
#   cbind(data.frame(all[, 1]), data.frame(apply(all[, -1], 2, as.numeric)))
# colnames(all)[1] <- 'obesity_status'
# all[, -1] <- scale(all[, -1])

X <- as.matrix(all[, -1])
X <- scale(X)
y <- as.factor(all$obesity_status)
set.seed(123)

# The alpha value is between 0 (ridge regression) and 1 (LASSO regression)

results <- list()


##########  'alpha' = alpha_val, 'enriched_genes' = ego,'enriched_kegg' = kk, 'dim_of_lasso' = dim(selected_genes_df))

# start to find the best alpha_value for the ego and kk method
for (alpha_val in seq(0.1, 1, by = 0.1)) {
  #alpha_val = 0.5
  cv_fit <-
    cv.glmnet(X, y, alpha = alpha_val, family = "binomial", set.seed(1234))
  
  # Check the best lambda
  best_lambda <- cv_fit$lambda.min
  
  # Fits the final mode
  final_model <-
    glmnet(
      X,
      y,
      alpha = alpha_val,
      lambda = best_lambda,
      family = "binomial",
      set.seed(1234)
    )
  
  #plot(cv_fit)
  # summary(final_model)[beta]
  #abline(v = best_lambda, col = "red")
  
  ######
  # library(glmnet)
  #
  # # final_model <- glmnet(X, y, alpha = alpha_val, lambda = best_lambda, family = "binomial")
  #
  # plot(final_model, xvar = "lambda", label = TRUE)
  
  #plot(final_model, xvar = "lambda", label = TRUE)
  
  # Extraction model coefficient
  coefficients <- coef(final_model, s = best_lambda)
  
  # View the changes selected for the model
  # Variables whose coefficients are not zero are considered to be selected into the model
  # Extract non-zero coefficients and their corresponding variable names
  non_zero_coefficients <-
    coefficients[coefficients[, 1] != 0, , drop = FALSE]
  variable_names <- rownames(non_zero_coefficients)
  
  # Create data frame
  selected_genes_df <- data.frame(Gene = variable_names,
                                  Coefficient = non_zero_coefficients[, 1])
  
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
  
  #print(merged_data)
  
  # enrichment analysis
  # BiocManager::install("clusterProfiler")
  # BiocManager::install("org.Hs.eg.db")
  library(knitr)
  install.packages("kableExtra")
  library(kableExtra)
  
  kable(head(merged_data, 10), "html", caption = "Merged Data Table") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
  
  
  #BiocManager::install("org.Hs.eg.db")
  genes <- merged_data$gene
  # Convert the gene name to an Entrez ID
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = genes,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  # enriched_genes method
  ego <- enrichGO(
    gene         = entrez_ids,
    OrgDb        = org.Hs.eg.db,
    keyType      = "ENTREZID",
    # BP(biological process), CC(cell components)MF(molecular function)
    ont          = "MF",
    pvalueCutoff = 0.1, # pvaluecutoff set by reference
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable     = TRUE
  )
  # enriched_kegg method
  kk <- enrichKEGG(
    gene = entrez_ids,
    organism = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1 # pvaluecutoff set by reference
  )
  # barplot(ego)
  
  # save the data into the results
  results[[as.character(alpha_val)]] <-
    list(
      'alpha' = alpha_val,
      'enriched_genes' = ego,
      'enriched_kegg' = kk,
      'dim_of_lasso' = dim(selected_genes_df)
    )
}

max_genes_count <- 0
best_alpha <- NA

#
for (alpha_val in names(results)) {
  enriched_result <- summary(results[[alpha_val]]$enriched_genes)
  enriched_terms_count <- nrow(enriched_result)
  
  if (enriched_terms_count > max_genes_count) {
    max_genes_count <- enriched_terms_count
    best_alpha <- alpha_val
  }
}

print(paste(
  "Best alpha:",
  best_alpha,
  "with",
  max_genes_count,
  "enriched terms"
))

results

library(ggplot2)

alpha_values <- seq(0.1, 1, by = 0.1)
term_counts <- numeric(length(alpha_values))

for (i in seq_along(alpha_values)) {
  alpha_val <- as.character(alpha_values[i])
  if (alpha_val %in% names(results)) {
    term_counts[i] <- nrow(summary(results[[alpha_val]]$enriched_genes))
  }
}

# plot the figure x is Alpha Value,y is Number of Enriched Terms.
plot_data <- data.frame(alpha = alpha_values, terms = term_counts)

ggplot(plot_data, aes(x = alpha, y = terms)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = 'aliceblue', color = NA),
    panel.grid.major = element_line(color = 'gray', size = 0.5),
    panel.grid.minor = element_line(color = 'lightgray', size = 0.25)
  ) +
  labs(title = "Number of Enriched Terms vs Alpha Value",
       x = "Alpha Value",
       y = "Number of Enriched Terms")


## plot the barplot for the ego and kk 

barplot(best_ego <- results[[best_alpha]]$enriched_genes)

## kk part
term_counts_kegg <- numeric(length(alpha_values)) # *

for (i in seq_along(alpha_values)) {
  alpha_val <- as.character(alpha_values[i])
  if (alpha_val %in% names(results)) {
    term_counts_kegg[i] <-
      nrow(summary(results[[alpha_val]]$enriched_kegg)) # *
  }
}


plot_data_kegg <-
  data.frame(alpha = alpha_values, terms = term_counts_kegg) # *

ggplot(plot_data_kegg, aes(x = alpha, y = terms)) + # *
  geom_line() +
  geom_point() +
  theme_minimal() +
  
  theme(
    panel.background = element_rect(fill = 'aliceblue', color = NA),
    panel.grid.major = element_line(color = 'gray', size = 0.5),
    panel.grid.minor = element_line(color = 'lightgray', size = 0.25)
  ) +
  labs(title = "Number of Enriched KEGG Terms vs Alpha Value", # *
       x = "Alpha Value",
       y = "Number of Enriched KEGG Terms") # *
#
# barplot
barplot(best_kk <-
          results[[which.max(term_counts_kegg)]]$enriched_kegg)

# kk <- enrichKEGG(gene = entrez_ids,
#                  organism = 'hsa',  # Human KEGG code
#                  pAdjustMethod = "BH",  # Method of adjusting p value
#                  pvalueCutoff = 0.1)  # Significance threshold
# barplot(kk)
