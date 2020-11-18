#!/software/R-3.4.0/bin/R
# Dafni Glinos
# 3 April 2018
#load libraries
suppressMessages(library("foreach"))
suppressMessages(library("doParallel"))
suppressMessages(library("dplyr"))
suppressMessages(library("optparse"))
suppressMessages(library("stringr"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("reshape2"))
suppressMessages(library("IRanges"))
suppressMessages(library("binr"))
suppressMessages(library("annotables"))

#Parse command-line options
option_list <- list(
  make_option(c("-t", "--tablepath"), type="character", default=NULL,
              help="Table of all genes assayed. It should contain the columns chr,
              start, end, length, ensgene, mean_expression and DE_status (0 or 1)", metavar = "type"),
  make_option(c("-i", "--iterations"), type="integer", default=1000,
              help="Number of iterations", metavar = "type")
  make_option(c("-d", "--disease"), type="character", default=NULL,
              help="Bed file of the LD loci of your tested trait", metavar = "type")
  )
opt <- parse_args(OptionParser(option_list=option_list))

#Extract parameters for CMD options
path_table = opt$t
iterations = as.numeric(opt$i)
pathdisease = opt$d

cl<-makeCluster(8)
registerDoParallel(cl)

# Read data inputs
MY_EXPR_DF = read.csv(file = path_table,header=TRUE)
MY_DIS <- read.table(file = path_disease,header=TRUE)

# Overlap the datasets
overlap_control <- IRanges::findOverlapPairs(snps,ld_genes)
overlap_control_sub <- overlap_control@second
overlap_control_sub <- unique(overlap_control_sub)

MY_EXPR_DF$inLD <- ifelse(MY_EXPR_DF$ensgene %in% overlap_control_sub$ensgene, 1, 0)
MY_EXPR_DF <- MY_EXPR_DF[ !duplicated(MY_EXPR_DF$ensgene), ]

print(paste0(length(overlap_control_sub)," genes are in LD loci"),quote=FALSE,file=stderr())
print(paste0(table(overlap_control_sub@elementMetadata$status)[1]," DE genes are in LD loci"),quote=FALSE,file=stderr())

###############
## add step to skip empty bins
###############

MY_EXPR_DF_DIS <- subset(MY_EXPR_DF,inLD == 1)

# set coordinates_length of lengths to samples from
min_size = round(dim(MY_EXPR_DF_DIS)[1]/3)
length_bins <- binr::bins(log2(MY_EXPR_DF_DIS$length), minpts = min_size, 
                 target.bins=round(length(log2(MY_EXPR_DF_DIS$length))/min_size, digits = 0))
coordinates_length <- names(length_bins$binct)
coordinates_length <- as.data.frame(gsub("\\[|\\]", "", coordinates_length))
colnames(coordinates_length) <- c("dummy")
coordinates_length <- stringr::str_split_fixed(coordinates_length$dummy, ",", 2)
coordinates_length <- as.data.frame(coordinates_length)
colnames(coordinates_length) <- c("start","end") 
coordinates_length$start <- as.numeric(as.character(coordinates_length$start))
coordinates_length$end <- as.numeric(as.character(coordinates_length$end))
# set smallest length as the one from the whole dataframe
coordinates_length[1,1] <- min(hist(log2(MY_EXPR_DF$length),plot=FALSE)$breaks)

breaks_length <- c(round(coordinates_length$start,digits=1),max(hist(log2(MY_EXPR_DF$length),plot=FALSE)$breaks))
counts_length <- hist(log2(MY_EXPR_DF_DIS$length),plot=FALSE, breaks = breaks_length)$counts

minG = round(dim(MY_EXPR_DF_DIS)[1]/10)
ls <-
  foreach::foreach(icount(iterations), .combine='c', .errorhandling='pass') %dopar% {
    gene_pool_lengths <- list()
    for (i in 1:length(counts_length)){
      gene_pool_expression <- list()
      tmp1 <- MY_EXPR_DF[log2(MY_EXPR_DF$length) > breaks_length[i] & log2(MY_EXPR_DF$length) <= breaks_length[i+1],]
      mdf <- MY_EXPR_DF_DIS[log2(MY_EXPR_DF_DIS$length) > breaks_length[i] & log2(MY_EXPR_DF_DIS$length) <= breaks_length[i+1],]
      expr_bins <- binr::bins(mdf$baseMean, minpts = minG,target.bins=round(length(mdf$baseMean)/minG, digits = 0))
      coordinates_expr <- names(expr_bins$binct)
      coordinates_expr <- as.data.frame(gsub("\\[|\\]", "", coordinates_expr))
      colnames(coordinates_expr) <- c("dummy")
      coordinates_expr <- stringr::str_split_fixed(coordinates_expr$dummy, ",", 2)
      coordinates_expr <- as.data.frame(coordinates_expr)
      colnames(coordinates_expr) <- c("start","end") 
      coordinates_expr$start <- as.numeric(as.character(coordinates_expr$start))
      coordinates_expr$end <- as.numeric(as.character(coordinates_expr$end))
      coordinates_expr[1,1] <- min(hist(tmp1$baseMean,plot=FALSE)$breaks)
      breakse <- c(round(coordinates_expr$start,digits=1),max(hist(tmp1$baseMean,plot=FALSE)$breaks))
      countse <- hist(mdf$baseMean,plot=FALSE,breaks = breakse)$counts
      for (j in 1:length(countse)){
        tmp2 <- dplyr::sample_n(tmp1[tmp1$baseMean > breakse[j] & tmp1$baseMean <= breakse[j+1],],countse[j],replace=FALSE)
        gene_pool_expression[[j]] <- tmp2
      }
      tmp3 <- do.call("rbind", gene_pool_expression)
      gene_pool_lengths[[i]] <- tmp3
    }
    gene_pool_counts <- do.call("rbind", gene_pool_lengths)
    mysum <- sum(gene_pool_counts$status)
    myobs <- ifelse(nrow(MY_EXPR_DF_DIS[MY_EXPR_DF_DIS$status == 1,]) >= mysum,0,1)
  }

pvalue <- (sum(ls) + 1)/iterations
print(paste0(pvalue," p-value"),quote=FALSE,file=stderr())