# Script for Transcriptomics-Informed Stoichiometric Modelling and Network Analysis (TISMAN)
# Claudio Tomi-Andrino (2021)
# other files can be found in https://github.com/CTA-code/TISMAN

# Gather data ---------------------------------
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library(dplyr)
library(cluster)
library(factoextra)
library(pracma)
library(stringr)

# set working directory as file where the script is located
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load gene list from the GSM
gene_list_human_gem <- read.delim(file.choose()) # for RStudio, to manually select the file

GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-GBM")
query_TCGA = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

gbm_res = getResults(query_TCGA) # make results as table
head(gbm_res) # data of the first 6 patients.
colnames(gbm_res) # columns present in the table

GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)
dim(tcga_data)

table(tcga_data@colData$vital_status)
head(assay(tcga_data)[,1:10]) 
head(rowData(tcga_data))

# extract FPKM values ready for rFASTCORMICS ============
# pipeline without clustering
query_TCGA_fast = query_TCGA
df_fastcormics = query_TCGA_fast[[1]][[1]]
df_fastcormics_2 <- filter(df_fastcormics, df_fastcormics$sample_type == 'Primary Tumor') # do not keep normal tissue
query_TCGA_fast[[1]][[1]] = df_fastcormics_2 # overwrite to use the proper object type
tcga_data_fast = GDCprepare(query_TCGA_fast)

dge_fast = DGEList( # creating a DGEList object
  counts=assay(tcga_data_fast),
  samples=colData(tcga_data_fast),
  genes=as.data.frame(rowData(tcga_data_fast)))

write.table(dge_fast$counts, file = "fpkm_fast.txt", sep = "\t",
            row.names = TRUE)

# Change headers so that matlab likes it (no minus sign). Select "fpkm_fast.txt"
re_fpkm_fast = as.data.frame(dge_fast$counts)
headers_fpkm = colnames(re_fpkm_fast)
headers_fpkm_corr = gsub("-", "", fast_headers, fixed = TRUE)
colnames(re_fpkm_fast) <- headers_fpkm_corr
write.table(headers_fpkm_corr, file = "headers_fpkm_corr.txt", sep = "\t",
            row.names = FALSE)

# Differential expression =================================
# https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html
# Tutorial by Tiago Maié & Fabio Ticconi, 21 Nov 2019
# define function
limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}

# call function
limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal" # this is the control group, by default the other one is the 'treated'
)

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)

# optional here
clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
head(design)

dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)
num_genes = dim(fit$genes)[1]
topGenes = topTable(fit, coef=1, number=num_genes, sort.by="logFC")
# now we filter genes by considering the following cut-offs
pCutoff = 10e-16
FCcutoff = 1.5
topGenes_interest_pre = topGenes[1,] # check it fulfills the criteria, use as a primer
rm(i)
rm(j)
i = 1 
j = 1

for(i in 1:dim(topGenes)[1]) {
    tmp_logFC <- topGenes[i,4]
    tmp_AdjPVal <- topGenes[i,8]

    if((tmp_logFC>=FCcutoff)&(tmp_AdjPVal<=pCutoff)) {
      topGenes_interest_pre[j,] = topGenes[i,]
      j = j+1
    }
  rm(tmp_logFC)
  rm(tmp_AdjPVal)
}

# now keep only those which are also defined in the model

topGenes_interest_def = topGenes_interest_pre[FALSE,]

for(i in 1:dim(gene_list_human_gem)[1]) {
  for(j in 1:dim(topGenes_interest_pre)[1]) { 
    tmp_gem <- gene_list_human_gem[i,1]
    tmp_dge <- topGenes_interest_pre[j,1]
    tmp_cmp = strcmp(tmp_gem,tmp_dge)
    if(tmp_cmp == TRUE) {
      topGenes_interest_def[i,] = topGenes_interest_pre[j,]
    }
  }
}

rm(tmp_gem)

write.table(topGenes_interest_def, file = "topGenes_interest_def.txt", sep = "\t",
            row.names = FALSE)

## plot PCA
plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}
res_pca = plot_PCA(limma_res$voomObj, "definition")

### Clustering of tumoural #############################
# Discarded 'solid tissue normal
df_export_1 = query_TCGA[[1]][[1]]
df_export_2 <- filter(df_export_1, df_export_1$sample_type == 'Primary Tumor') # do not keep normal tissue
query_TCGA[[1]][[1]] = df_export_2 # overwrite to use the proper object type
tcga_data_mod = GDCprepare(query_TCGA)

# in order to normalise and remove genes with counts with less than 10 reads, need to convert into DGEList object
df_export_3 = DGEList( # creating a DGEList object
  counts=assay(tcga_data_mod),
  samples=colData(tcga_data_mod),
  genes=as.data.frame(rowData(tcga_data_mod)))

keep = filterByExpr(df_export_3,design=NULL) # this function keeps "genes have sufficiently large counts to be retained in a statistical analysis."

df_export_4 = df_export_3[keep,,keep.lib.sizes=FALSE] 
df_export_norm = calcNormFactors(df_export_4,method="TMM") # this should be suffice, the voom one does not seem to be necessary here
conv_export = as.matrix(df_export_norm) # need to apply filters to discretise (making it binary)
global_threshold = quantile(conv_export, c(.25)) # values below overall 25% are assumed to be not expressed
local_threshold = rowMeans(conv_export) # average expression, gene-specific
  # get a 0ed matrix with the same dimensions and headers
conv_export_bin = conv_export[]
conv_export_bin[] <- 0L
  # binary matrix
for(i in 1:dim(conv_export)[1]) {
  for(j in 1:dim(conv_export)[2]) {
    if( (conv_export[i,j]>global_threshold)&(conv_export[i,j]>local_threshold[i]) ){
      conv_export_bin[i,j] = 1
    }
  }
}
conv_export_bin_t = t(conv_export_bin) # transpose it so that it will cluster samples instead of genes
conv_export_bin_t_df = as.data.frame( conv_export_bin_t[], drop=false) 
final_mat = as.matrix(conv_export_bin_t_df)

local_threshold_export = as.data.frame(local_threshold)
tmp_rownames = rownames(conv_export)
rownames(local_threshold_export) <- tmp_rownames
col_headings <- c('gene_id','FPKM')
names(local_threshold_export) <- col_headings
write.table(local_threshold_export, file = "local_thresholds.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
 
# Estimate number of k clusters using the silhouette method
# manually assess the best k value
fviz_nbclust(final_mat, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# K-medoids
pam.res <- pam(final_mat, 2, metric = "manhattan") # more robust than euclidean, NECESSARY for binary data
print(pam.res)
medoids = pam.res$medoids

# import gene list from HUMAN-GEM: generate a list for each medoid, keeping the same gene order as in the model
t_medoids = t(medoids)
t_medoid_1 = as.data.frame(t_medoids[,1])
t_medoid_2 = as.data.frame(t_medoids[,2])
gene_list_tcga = as.data.frame(row.names(t_medoids)) # indexing of genes: content and number following the order in t_medoids
human_gem_t_medoid_1 = gene_list_human_gem
human_gem_t_medoid_2 = gene_list_human_gem

# for loop with 'if' statements to fill the new vector. Whenever a gene is in the GSM but not in the RNA-seq -> NA
rm(i)
rm(j)
i = 1 
j = 1

for(i in 1:dim(gene_list_human_gem)[1]) {
  for(j in 1:dim(gene_list_tcga)[1]) { # common for both t_medoids
    tmp_gem <- gene_list_human_gem[i,1]
    tmp_tcga <- gene_list_tcga[j,1]
    tmp_cmp = strcmp(tmp_gem,tmp_tcga)
    if(tmp_cmp == TRUE) {
      human_gem_t_medoid_1[i,2] = t_medoid_1[j,1]
      human_gem_t_medoid_2[i,2] = t_medoid_2[j,1]
    }
  }
}

col_headings <- c('gene_id','binary')
names(human_gem_t_medoid_1) <- col_headings
names(human_gem_t_medoid_2) <- col_headings

# export so it's ready for matlab
write.table(human_gem_t_medoid_1, file = "human_gem_t_medoid_1.txt", sep = "\t",
            row.names = FALSE)
write.table(human_gem_t_medoid_2, file = "human_gem_t_medoid_2.txt", sep = "\t",
            row.names = FALSE)

# extract FPKM values ready for rFASTCORMICS ============
# pipeline with clustering

# extract clusters with normalised data
clusters = as.data.frame(pam.res$clustering)
df_conv_export = as.data.frame(conv_export)
df_conv_export_t = t(df_conv_export)
num_medoid_1 = 0
num_medoid_2 = 0
for(i in 1:dim(clusters)[1]) {
  tmp_cl <- clusters[i,1]
    if(tmp_cl == 1) {
      num_medoid_1 = num_medoid_1 + 1
    } else if ( tmp_cl == 2) {
      num_medoid_2 = num_medoid_2 + 1
    }
}

num_genes_exp = dim(df_conv_export_t)[2]
df_medoid_1 = as.data.frame(matrix(numeric(),nrow = num_medoid_1, ncol = num_genes_exp))
df_medoid_2 = as.data.frame(matrix(numeric(),nrow = num_medoid_2, ncol = num_genes_exp))
row_names = row.names(df_conv_export_t)
row_names_df_medoid_1 <- character()[1:num_medoid_1]
row_names_df_medoid_2 <- character()[1:num_medoid_2]

j = 1
for(i in 1:dim(clusters)[1]) {
  tmp_cl <- clusters[i,1]
  if(tmp_cl == 1) {
    df_medoid_1[j,1:dim(df_medoid_1)[2]] = df_conv_export_t[i,1:dim(df_conv_export_t)[2]]
    row_names_df_medoid_1[j] = row_names[i]
    j = j+1
    rm(tmp_cl)
  }
}
row.names(df_medoid_1) <- row_names_df_medoid_1

j = 1
for(i in 1:dim(clusters)[1]) {
  tmp_cl <- clusters[i,1]
  if(tmp_cl == 2) {
    df_medoid_2[j,1:dim(df_medoid_2)[2]] = df_conv_export_t[i,1:dim(df_conv_export_t)[2]]
    row_names_df_medoid_2[j] = row_names[i]
    j = j+1
    rm(tmp_cl)
  }
}
row.names(df_medoid_2) <- row_names_df_medoid_2
col_names = colnames(df_conv_export_t)
colnames(df_medoid_1) <- col_names
colnames(df_medoid_2) <- col_names

# since FASTCORMICS has its own normalisation protocol, we need to export untreated FPKM values
row_names_fast = rownames(re_fpkm_fast)

re_fpkm_fast_cluster_1 = as.data.frame(matrix(numeric(),nrow = length(row_names_fast) , ncol = num_medoid_1))
rownames(re_fpkm_fast_cluster_1) <- row_names_fast
headers_fpkm_fast_1 <- row_names_df_medoid_1
headers_fpkm_fast_1_corr = gsub("-", "", headers_fpkm_fast_1, fixed = TRUE)
colnames(re_fpkm_fast_cluster_1) <- headers_fpkm_fast_1_corr

re_fpkm_fast_cluster_2 = as.data.frame(matrix(numeric(),nrow = length(row_names_fast) , ncol = num_medoid_2))
rownames(re_fpkm_fast_cluster_2) <- row_names_fast
headers_fpkm_fast_2 <- row_names_df_medoid_2
headers_fpkm_fast_2_corr = gsub("-", "", headers_fpkm_fast_2, fixed = TRUE)
colnames(re_fpkm_fast_cluster_2) <- headers_fpkm_fast_2_corr

j = 1
k = 1
for(i in 1:dim(clusters)[1]) {
  tmp_cl <- clusters[i,1]
  if(tmp_cl == 1) {
    re_fpkm_fast_cluster_1[1:dim(re_fpkm_fast_cluster_1)[1],j] = re_fpkm_fast[1:dim(re_fpkm_fast)[1],i]
    j = j+1
  } else if ( tmp_cl == 2) {
    re_fpkm_fast_cluster_2[1:dim(re_fpkm_fast_cluster_2)[1],k] = re_fpkm_fast[1:dim(re_fpkm_fast)[1],i]
    k = k+1
  }
}
write.table(re_fpkm_fast_cluster_1, file = "re_fpkm_fast_cluster_1.txt", sep = "\t",
            row.names = TRUE)
write.table(headers_fpkm_fast_1_corr, file = "headers_fpkm_fast_1_corr.txt", sep = "\t",
            row.names = FALSE)
write.table(re_fpkm_fast_cluster_2, file = "re_fpkm_fast_cluster_2.txt", sep = "\t",
            row.names = TRUE)
write.table(headers_fpkm_fast_2_corr, file = "headers_fpkm_fast_2_corr.txt", sep = "\t",
            row.names = FALSE)
 

# rm(list = ls()) # to clear the workspace








