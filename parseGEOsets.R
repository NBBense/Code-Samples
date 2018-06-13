################################################################################
# Title: parseGEOsets
# Author: Nicholas Bense
# Date: 4/10/18
# Description: R script to parse and separate Gene Expression Omnibus (GEO) expression 
# values and metadata information then subset into disease and control partitions.
################################################################################

# Install/load the required libraries
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
install.packages("Biobase")
library(GEOquery)
library(Biobase)

# Convert GSE file into parsable R object
GSE11138_file<-"GSE11138_series_matrix.txt.gz"
GSE11138<-getGEO(filename=GSE11138_file)

# Obtain gene symbols
GSE11138_gene_symbols<-fData(GSE11138)$Gene_Symbol
GSE11138_gene_symbols<-as.character(GSE11138_gene_symbols)

# Obtain expression data
expression_Data<-as.data.frame(exprs(GSE11138))
expression_Data<-cbind(GSE11138_gene_symbols,expression_Data)

# Examine disease phenotype (also useful to reference supplementary patient file)
pData(GSE11138)$pathology

# Write data into disease (Dz) dataframe
colnames(expression_Data)[1] <- "Gene_Symbol"
GSE11138.Dz.df<-expression_Data

# Convert GSE file into parsable R object
GSE2143_file<-"GSE2143_series_matrix.txt.gz"
GSE2143<-getGEO(filename=GSE2143_file)

# Obtain gene symbols
gene_symbols<-fData(GSE2143)$GENE_NAME
gene_symbols<-as.character(gene_symbols)

# Read in table of known gene symbols
known_gene_symbols_file<-"gene_symbols.txt"
known_gene_symbols<- read.table(known_gene_symbols_file, sep = '\t',header = F, skip = 1,quote='', comment='')
known_gene_symbols<- as.character(known_gene_symbols[,2])

# Extract gene symbols from unstructured metadata
new_gene_symbols = c()
for (symbol in gene_symbols){
  for (known_symbol in known_gene_symbols){
    if (grepl(paste('\\<',known_symbol,'\\>', sep = ""),symbol) == TRUE | grepl(paste('(\\b',known_symbol,'\\b)',sep = ""),symbol)){
      new_gene_symbol<-known_symbol
      break
    } else{
      new_gene_symbol<-symbol
    }
  }
  new_gene_symbols<-c(new_gene_symbols,new_gene_symbol)
}

# Obtain expression data
expression_Data<-as.data.frame(exprs(GSE2143))
expression_Data<-cbind(new_gene_symbols,expression_Data)

# Write data into unknown .CSV
colnames(expression_Data)[1] <- "Gene_Symbol"
write.csv(expression_Data, file = "GEO_Unknown_Type_Expression.csv",row.names=FALSE)
write.table(expression_Data, "GEO_Unknown_Type_Expression.txt", sep="\t",row.names=FALSE)

# Convert GDS file into parsable R object
GDS5083_file<-"GDS5083.soft.gz"
GDS5083<-getGEO(filename=GDS5083_file)

# Convert GDS to ExpressionSet object
GDS5083_eset<-GDS2eSet(GDS5083, do.log2=TRUE)

# Obtain gene symbols
gene_symbols<-fData(GDS5083_eset)$`Gene symbol`
gene_symbols<-as.character(gene_symbols)

# Standardize gene alias lists with single gene symbol for facilitated analysis
new_gene_symbols = c()
for (symbol in gene_symbols){
  for (known_symbol in known_gene_symbols){
    if (grepl(paste('///\\b',known_symbol,'\\b///', sep = ""),symbol) == TRUE | grepl(paste('///\\b',known_symbol,'\\b',sep = ""),symbol) | grepl(paste('\\b',known_symbol,'\\b///',sep = ""),symbol)){
      new_gene_symbol<-known_symbol
      break
    } else{
      new_gene_symbol<-symbol
    }
  }
  new_gene_symbols<-c(new_gene_symbols,new_gene_symbol)
}

# Obtain expression data
expression_Data<-as.data.frame(exprs(GDS5083_eset),stringsAsFactors = FALSE)
expression_Data<-cbind(new_gene_symbols,expression_Data)

# Examine disease phenotype 
pData(GDS5083_eset)$specimen

# Split data into disease and control dataframes and write control file
colnames(expression_Data)[1] <- "Gene_Symbol"
GDS5083.Dz.df<-expression_Data[,1:32]
GDS5083.Control.df<-cbind(new_gene_symbols,expression_Data[,32:64])
colnames(GDS5083.Control.df)[1] <- "Gene_Symbol"
write.csv(GDS5083.Control.df, file = "GEO_Control_Type_Expression.csv",row.names=FALSE)
write.table(GDS5083.Control.df, "GEO_Control_Type_Expression.txt", sep="\t",row.names=FALSE)

# Combine results and write files
final.Dz.df<-merge(GSE11138.Dz.df,GDS5083.Dz.df, all = TRUE)
write.csv(final.Dz.df, file = "GEO_Dz_Type_Expression.csv",row.names=FALSE)
write.table(final.Dz.df, "GEO_Dz_Type_Expression.txt", sep="\t",row.names=FALSE)
