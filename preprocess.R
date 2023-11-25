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

start.time <- Sys.time()

data <- read.delim('data/GSE72556_series_matrix.txt', sep = '!')
names(data) <- c("cg", "SampleCharacteristics")

# Define a function to process and split the data
process_data <- function(data, pattern, remove_pattern, is_numeric = TRUE) {
  # Extract data that matches the pattern
  relevant_data <- data[grepl(pattern, data$SampleCharacteristics), "SampleCharacteristics"]
  # Remove unwanted text
  relevant_data <- gsub(remove_pattern, "", relevant_data)
  # Additional cleaning steps, TAB removal, etc
  relevant_data <- gsub("\\t", "", relevant_data) # Remove the TAB
  relevant_data <- gsub("Sample_characteristics_ch1", "", relevant_data) # Remove specific string
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
gender <- process_data(data, "child gender:", "child gender:", is_numeric = FALSE)[-1]
age <- process_data(data, "child age:", "child age:")

# Processing accession number
accession <- process_data(data, "GSM", "", is_numeric = FALSE)[1:96]
accession[1] <- 'GSM1865141'
accession <- subset(accession, accession != "GSM1865188")
accession <- subset(accession, accession != "GSM1865221")
accession <- subset(accession, accession != "GSM1865231")

# Create data frame
adult_child <- data.frame(accession, adult_bmi, child_bmi, adult_waist, child_waist, gender, age)

end.time <- Sys.time()
cat('Time for optimization: ', end.time - start.time, '\n')

# library(data.table)
# # read data
# y1 <- fread('data/GSE72556_series_matrix.txt', header = FALSE, fill = TRUE)
# y2 <- y1
# # Remove comments from the first 65 lines
# y1 <- tail(y2, -65)
# 
# # Transposed data
# y1 <- transpose(y1)
# 
# # Set the first column name, the first row name, and then delete
# colnames(y1) <- as.character(y1[1, ])
# y1 <- y1[-1, ]
