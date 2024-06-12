#USAGE: Rscript DEseq_softcode.R [1]rawcounts_filename.txt [2]group-samples.csv [3]list-contrastsCSV.csv [4]out-directory(include '/' at the beginning and end of the file path)

# set library path

libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

#install.packages("vctrs") 
#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("dplyr")
#install.packages("EnhancedVolcano")
#install.packages("DESeq2")
#install.packages("purrr")
#install.packages("geneplotter")
#library(tidyverse)
#BiocManager::install("geneplotter")
#BiocManager::install("DESeq2")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(EnhancedVolcano)
library(DESeq2)
library(purrr)

## Read command line arguments ------------------------
args <- commandArgs(TRUE)
counts_file <- args[1]
grp_samples_file <- args[2]
list_contrasts_file <- args[3]
outs_directory <- as.character(args[4])
# Read counts data -----------------------
# to be changed to counts_data object

counts_data <- read.table(counts_file, header=TRUE, sep="\t", row.names = 1) # non-hardcoded script should contain {file_name} as the input


# make a copy for reinserting Gene Symbol column later
counts_copy <- counts_data
# remove the gene symbol as a column 
counts_data <- counts_data %>% select(-Gene_Symbol)


# turn character values into numeric values
counts_data <- counts_data %>% mutate(across(everything(),as.numeric))

# ColData -------------------------------
# to be changed to ColData object
ColData <- read.csv(grp_samples_file, header=TRUE, sep=",") 
# ColData should be read in as args[2]

# now DEseq processes can begin ---------------------------

counts_data_colfilter <- counts_data[colnames(counts_data) %in% ColData$sample] 
# code is selecting columns in the counts_data object that are present in the ColData$sample vector. The result of this code will be a new object, counts_data_filtered, which contains only the columns that were present in the ColData$sample vector.
ColData_rowfilter <- ColData[ColData$sample %in% colnames(counts_data),]
# code is selecting the rows in the ColData object that have the same column names as the columns in the counts_data object. The result of this code will be a new object, ColData_test, which contains only the rows that have the same column names as the columns in the counts_data object.

# this is like a venn diagram middle section between these two dataframes that relate to each other
# not needed - take out the line
ColData_rowfilter <- ColData_rowfilter %>% distinct(ColData_rowfilter$sample, .keep_all = TRUE)
#  code is removing any duplicate rows from the ColData_test object based on the values in the V2 column. The result of this code will be a new object, ColData_rowfilter, which does not contain any duplicate rows based on the values in the V2 column.
# order the counts data and the ColData to be the same arrangement
counts_data_colfilter = counts_data_colfilter[,order(colnames(counts_data_colfilter))]

ColData_rowfilter = ColData_rowfilter[order(ColData_rowfilter$sample),]

all(colnames(counts_data_colfilter) %in% ColData_rowfilter$sample)
all(colnames(counts_data_colfilter) == ColData_rowfilter$sample)

# 2/7/2023 - turn the contrasts of interest list-o-lists into a CSV argument that can be read in

contrasts_info <- as.list(read.csv(list_contrasts_file, header=TRUE, sep=",")) 
# each field of the CSV is a list that can then be read into a list-of-lists - that then makes up the contrasts of interest 
contrasts_names <- contrasts_info$cond1_vs_cond2

contrasts_grp <- contrasts_info$group

contrasts_cond1 <- contrasts_info$condition1

contrasts_cond2 <- contrasts_info$condition2

contrasts_list <- list()
# FOR loop for naming the conditions being compared
for(i in 1:length(contrasts_names)){
  contrasts_list[[contrasts_names[i]]] <- contrasts_names[i]
}
# FOR loop for putting 'group' as the first item in the contrasts of interest list
for(i in 1:length(contrasts_grp)){
  contrasts_list[i][1] <- contrasts_grp[i]
}
# FOR loop for putting condition1 as the second item in the contrasts of interest list
for(i in 1:length(contrasts_cond1)){
  contrasts_list[[i]][[2]] <- contrasts_cond1[i]
}
# FOR loop for putting condition1 as the second item in the contrasts of interest list
for(i in 1:length(contrasts_cond2)){
  contrasts_list[[i]][[3]] <- contrasts_cond2[i]
}

# This code is creating a list of contrasts. The list contains as many elements as conditions an individual would want to compare. Each element contains a vector of column names associated with the contrast.




# ===============================================================================


# dds object creation _----------------------------------
dds <- DESeqDataSetFromMatrix(counts_data_colfilter, ColData_rowfilter, design = ~ group) 

#Takes the counts dataframe, the paired ColData dataframe, and the design argument. The 'design' refers to the column in the ColData dataframe will be used as the design matrix for the object. (what samples from the second column that will be taken and compared within the DEseq algorithm)
dds <- DESeq(dds)
k = lapply(contrasts_list,function(cntrst){ # the original list of lists for comparison groups is swapped out for the one that can be read in via CSV
  resMFType <- results(dds,
                       format = "DataFrame",
                       contrast = cntrst)
  return(resMFType)})
#  code is creating a list of results objects from the DESeq object. The lapply function is used to loop over the contrasts in the contrasts_ofint list, and the results function is used to extract the results for each contrast. The format argument is used to specify that the results should be returned as a DataFrame.
k


# ======== Space to work code to include the means of the two conditions being compared into the output CSV ==================== 

# Including each group basemean would have to come from individually calculating the rowMeans(counts(dds, normalized=TRUE)) - per each condition and then added in at the end
dds_counts_df <- data.frame(counts(dds, normalized=TRUE))

# a FOR loop can go through the contrasts of interest list and subset each condition 1 groups into a df
# # a FOR loop can go through the contrasts of interest list and subset each condition 2 groups into a df

condition1_grps <- data.frame()
condition2_grps <- data.frame()
for(i in 1:length(contrasts_list)){
  condition1_grps <- rbind(condition1_grps, data.frame(contrasts_list[[i]][2]))
  condition2_grps <- rbind(condition2_grps, data.frame(contrasts_list[[i]][3]))
}
colnames(condition1_grps) <- c("condition1")
colnames(condition2_grps) <- c("condition2")
contrast_grps_df <- data.frame(condition1_grps,condition2_grps)
# with the contrasts list conditions each in their dataframe, the next step is to match the values in the contrast condition df with column names in the normalized dds counts dataframe - but the column names in dds df do NOT match the names in the contrast_grps_df

# edit so ColData matches column names in dds norm counts
ColData_filter_dds_colnames <- ColData[ColData$sample %in% colnames(dds_counts_df),]

all(colnames(dds_counts_df) %in% ColData_filter_dds_colnames$sample)

dds_counts_df = dds_counts_df[,order(colnames(dds_counts_df))]

ColData_filter_dds_colnames = ColData_filter_dds_colnames[order(ColData_filter_dds_colnames$sample),]

ifelse(colnames(dds_counts_df) == ColData_filter_dds_colnames$sample, colnames(dds_counts_df) <- ColData_filter_dds_colnames$group, colnames(dds_counts_df))

# dds_counts_df have the groups as column names now, meaning it can be sub-selected using the contrast_grps_df

# for condition 1 a list of dataframes is made that is a subset from the larger normalized counts dds_counts_df
basemean_condition1_dfs <- list()

basemean_condition1_dfs <- lapply(contrast_grps_df$condition1, function(row) {
  dds_counts_df_temp <- dds_counts_df[colnames(dds_counts_df) %in% row]
  basemean_condition1_dfs[[row]] <- dds_counts_df_temp
})
names(basemean_condition1_dfs) = lapply(contrasts_list, function(x) x[[2]])

# for condition 2 a list of dataframes is made that is a subset from the larger normalized counts dds_counts_df
#USE lapply() function to avoid weird FOR loop behaviors
basemean_condition2_dfs <- list()

basemean_condition2_dfs <- lapply(contrast_grps_df$condition2, function(row) {
  dds_counts_df_temp <- dds_counts_df[colnames(dds_counts_df) %in% row]
  basemean_condition2_dfs[[row]] <- dds_counts_df_temp
})
names(basemean_condition2_dfs) = lapply(contrasts_list, function(x) x[[3]])

# at this point, both conditions from EACH contrast of interest have the normalized counts from the dds object for each replicate in a list of dataframes

# the next step would be to get the mean off all the replicates for each dataframe in each condition list - and THEN input the mean values into the results table with all the other info

basemean_condition1_dfs <- lapply(basemean_condition1_dfs, rowMeans)

basemean_condition2_dfs <- lapply(basemean_condition2_dfs, rowMeans)
# both condition contrasts have their normalized counts row means calculated and are ready to be inserted into the final results CSV - line 285



# ======================= volcano plot automation ======================


# Create an empty list to store the results of the loop
list_results <- list()

# Create a for loop that swaps ensembleID for the gene symbol - so when it is plotted on as a Volcano Plot, it has human readable points labelled 
for(i in 1:length(k)){
  # Extract the matrix from the k list
  matrix <- k[[i]]
  
  # Create a data frame from the matrix
  df <- data.frame(matrix)
  
  # Add the Gene Symbol column to the data frame
  df$Gene_Symbol <- counts_copy$Gene_Symbol
  
  # Reorder the columns so the Gene Symbol is the first column
  df <- df %>% relocate(Gene_Symbol) 
  
  # Remove duplicate rows
  df <- df %>% distinct(Gene_Symbol, .keep_all = TRUE)
  
  # Set the row names as the Gene Symbol
  rownames(df) <- df$Gene_Symbol
  
  # Remove the Gene Symbol column
  df <- df %>% select(-Gene_Symbol)
  
  # Add the data frame to the list
  list_results[[i]] <- df
}
# list_results contains all the results tables in the same order as they were when ran through the lapply() function 

## use the names from the previous list and apply to the results list
names(list_results)=names(contrasts_list)



# testing the volcano plot using FOR loop for created DEseq datasets
for (i in seq_along(list_results)) {
  EnhancedVolcano(list_results[[i]],
                  lab = rownames(list_results[[i]]),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = names(list_results[i]),
                  subtitle = paste0("Reference Condition: ",contrasts_list[[i]][3]), # volcano plot now lets others know what the values are being compared to 2/7/2023
                  titleLabSize = 10,
                  pointSize = 1.5,
                  labSize = 3)
  ggsave(paste0(outs_directory,names(list_results[i]), "_Volcano_Plot.pdf"), width = 22, height = 15, units = "cm")
}

# a FOR loop that keeps the ensemble gene ID and the gene symbol 

# Create an empty list to store the results of the loop ===========================================
csv_results <- list()

# Create a for loop that adds in a column for the gene symbol - so when it is plotted on as a Volcano Plot, it has human readable points labelled 
for(i in 1:length(k)){
  # Extract the matrix from the k list
  matrix <- k[[i]]
  
  # Create a data frame from the matrix
  df <- data.frame(matrix)
  
  # Add the Gene Symbol column to the data frame
  df$Gene_Symbol <- counts_copy$Gene_Symbol
  
  # Reorder the columns so the Gene Symbol is the first column
  df <- df %>% relocate(Gene_Symbol) 
  
  # 2/7/2023 - to added reference level info to the csv output table
  # Add the data frame to the list
  csv_results[[i]] <- df
  colnames(csv_results[[i]])[which(colnames(csv_results[[i]]) == "log2FoldChange")] <- paste0("log2FoldChange","(",contrasts_list[[i]][2],"/",contrasts_list[[i]][3],")")
}
# name the list items
names(csv_results)=names(contrasts_list)
# this is where I should insert the condition CONTRASTS means from the normalized counts - initially coming from the dds object =============================================
# lapply function should be used here to turn the list of condtion double values into a single dataframe 
basemean_condition1_dfs <- lapply(basemean_condition1_dfs, function(x) data.frame(x))

basemean_condition2_dfs <- lapply(basemean_condition2_dfs, function(x) data.frame(x))


# each condition has normalized means in a list of dataframes and this next part names each column by the corresponding list object

for(i in 1:length(basemean_condition1_dfs)){
  colnames(basemean_condition1_dfs[[i]]) <- paste0(names(basemean_condition1_dfs[i]), " Mean")
}


for(i in 1:length(basemean_condition2_dfs)){
  colnames(basemean_condition2_dfs[[i]]) <- paste0(names(basemean_condition2_dfs[i]), " Mean")
}


combined_conditions <- lapply(seq_along(basemean_condition1_dfs), function(x) cbind(basemean_condition1_dfs[[x]], basemean_condition2_dfs[[x]]))
names(combined_conditions)=names(contrasts_list)
# all the contrasts of interest have their individual sample group means computed
# NEXT STEP: insert sequentially into the csv_results list 
combined_csv_results <- lapply(seq_along(csv_results), function(x) cbind(csv_results[[x]], combined_conditions[[x]]))
names(combined_csv_results)=names(contrasts_list)

# reorganize the columns so that the sample means are next to the baseMean column
combined_csv_results <- map(combined_csv_results, ~ .x[,c(1,2,8,9,3,4,5,6,7)])

# == write all comparison results to a CSV - combine with the steps above. format results, name the df, write the csv ===========================

lapply(1:length(combined_csv_results), function(i) write.csv(combined_csv_results[[i]], 
                                                             file = paste0(outs_directory,names(combined_csv_results[i]), ".csv"),
                                                             row.names = TRUE))
