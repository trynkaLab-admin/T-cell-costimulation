library(DESeq2)
library(tidyverse)
library(ggrepel)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Calculate all the DEG using condition 7 as baseline
res_naive_1 <- results(dds, contrast=c("group", "hTCR_hCD28Naive", "restingNaive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_naive_1_q <- subset(res_naive_1,padj <= 0.05)
res_naive_1_log <- subset(res_naive_1_q,abs(log2FoldChange) >= 1)
res_naive_1$symbol <- grch38$symbol[match(rownames(res_naive_1),grch38$ensgene)]
write.csv(res_naive_1,"hTCRhCD28_rest_naive_comp.csv",quote=FALSE)

res_naive_2 <- results(dds, contrast=c("group", "lTCR_hCD28Naive", "restingNaive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_naive_2_q <- subset(res_naive_2,padj <= 0.05)
res_naive_2_log <- subset(res_naive_2_q,abs(log2FoldChange) >= 1)
res_naive_2$symbol <- grch38$symbol[match(rownames(res_naive_2),grch38$ensgene)]
write.csv(res_naive_2,"lTCRhCD28_rest_naive_comp.csv",quote=FALSE)

res_naive_3 <- results(dds, contrast=c("group", "hTCR_lCD28Naive", "restingNaive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_naive_3_q <- subset(res_naive_3,padj <= 0.05)
res_naive_3_log <- subset(res_naive_3_q,abs(log2FoldChange) >= 1)
res_naive_3$symbol <- grch38$symbol[match(rownames(res_naive_3),grch38$ensgene)]
write.csv(res_naive_3,"hTCRlCD28_rest_naive_comp.csv",quote=FALSE)

res_naive_4 <- results(dds, contrast=c("group", "lTCR_lCD28Naive", "restingNaive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_naive_4_q <- subset(res_naive_4,padj <= 0.05)
res_naive_4_log <- subset(res_naive_4_q,abs(log2FoldChange) >= 1)
res_naive_4$symbol <- grch38$symbol[match(rownames(res_naive_4),grch38$ensgene)]
write.csv(res_naive_4,"lTCRlCD28_rest_naive_comp.csv",quote=FALSE)

res_naive_5 <- results(dds, contrast=c("group", "hTCRNaive", "restingNaive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_naive_5_q <- subset(res_naive_5,padj <= 0.05)
res_naive_5_log <- subset(res_naive_5_q,abs(log2FoldChange) >= 1)
res_naive_5$symbol <- grch38$symbol[match(rownames(res_naive_5),grch38$ensgene)]
write.csv(res_naive_5,"hTCR_rest_naive_comp.csv",quote=FALSE)

res_naive_6 <- results(dds, contrast=c("group", "hCD28Naive", "restingNaive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_naive_6_q <- subset(res_naive_6,padj <= 0.05)
res_naive_6_log <- subset(res_naive_6_q,abs(log2FoldChange) >= 1)
res_naive_6$symbol <- grch38$symbol[match(rownames(res_naive_6),grch38$ensgene)]
write.csv(res_naive_6,"hCD28_rest_naive_comp.csv",quote=FALSE)

n_vec_pos <- unique(c(rownames(res_naive_1_log[res_naive_1_log$log2FoldChange > 0,]),rownames(res_naive_2_log[res_naive_2_log$log2FoldChange > 0,]),
                  rownames(res_naive_3_log[res_naive_3_log$log2FoldChange > 0,]),rownames(res_naive_4_log[res_naive_4_log$log2FoldChange > 0,]),
                  rownames(res_naive_5_log[res_naive_5_log$log2FoldChange > 0,]),rownames(res_naive_6_log[res_naive_6_log$log2FoldChange > 0,])))

n_vec <- unique(c(rownames(res_naive_1_log),rownames(res_naive_2_log),
                      rownames(res_naive_3_log),rownames(res_naive_4_log),
                      rownames(res_naive_5_log),rownames(res_naive_6_log)))

n_vec_symbol <- grch38$symbol[match(n_vec_pos,grch38$ensgene)]

res_memory_1 <- results(dds, contrast=c("group", "hTCR_hCD28Memory", "restingMemory"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_memory_1_q <- subset(res_memory_1,padj <= 0.05)
res_memory_1_log <- subset(res_memory_1_q,abs(log2FoldChange) >= 1)
res_memory_1$symbol <- grch38$symbol[match(rownames(res_memory_1),grch38$ensgene)]
write.csv(res_memory_1,"hTCRhCD28_rest_memory_comp.csv",quote=FALSE)

res_memory_2 <- results(dds, contrast=c("group", "lTCR_hCD28Memory", "restingMemory"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_memory_2_q <- subset(res_memory_2,padj <= 0.05)
res_memory_2_log <- subset(res_memory_2_q,abs(log2FoldChange) >= 1)
res_memory_2$symbol <- grch38$symbol[match(rownames(res_memory_2),grch38$ensgene)]
write.csv(res_memory_2,"lTCRhCD28_rest_memory_comp.csv",quote=FALSE)

res_memory_3 <- results(dds, contrast=c("group", "hTCR_lCD28Memory", "restingMemory"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_memory_3_q <- subset(res_memory_3,padj <= 0.05)
res_memory_3_log <- subset(res_memory_3_q,abs(log2FoldChange) >= 1)
res_memory_3$symbol <- grch38$symbol[match(rownames(res_memory_3),grch38$ensgene)]
write.csv(res_memory_3,"hTCRlCD28_rest_memory_comp.csv",quote=FALSE)

res_memory_4 <- results(dds, contrast=c("group", "lTCR_lCD28Memory", "restingMemory"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_memory_4_q <- subset(res_memory_4,padj <= 0.05)
res_memory_4_log <- subset(res_memory_4_q,abs(log2FoldChange) >= 1)
res_memory_4$symbol <- grch38$symbol[match(rownames(res_memory_4),grch38$ensgene)]
write.csv(res_memory_4,"lTCRlCD28_rest_memory_comp.csv",quote=FALSE)

res_memory_5 <- results(dds, contrast=c("group", "hTCRMemory", "restingMemory"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_memory_5_q <- subset(res_memory_5,padj <= 0.05)
res_memory_5_log <- subset(res_memory_5_q,abs(log2FoldChange) >= 1)
res_memory_5$symbol <- grch38$symbol[match(rownames(res_memory_5),grch38$ensgene)]
write.csv(res_memory_5,"hTCR_rest_memory_comp.csv",quote=FALSE)

res_memory_6 <- results(dds, contrast=c("group", "hCD28Memory", "restingMemory"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res_memory_6_q <- subset(res_memory_6,padj <= 0.05)
res_memory_6_log <- subset(res_memory_6_q,abs(log2FoldChange) >= 1)
res_memory_6$symbol <- grch38$symbol[match(rownames(res_memory_6),grch38$ensgene)]
write.csv(res_memory_6,"hCD28_rest_memory_comp.csv",quote=FALSE)

m_vec_pos <- unique(c(rownames(res_memory_1_log[res_memory_1_log$log2FoldChange > 0,]),rownames(res_memory_2_log[res_memory_2_log$log2FoldChange > 0,]),
                  rownames(res_memory_3_log[res_memory_3_log$log2FoldChange > 0,]),rownames(res_memory_4_log[res_memory_4_log$log2FoldChange > 0,]),
                  rownames(res_memory_5_log[res_memory_5_log$log2FoldChange > 0,]),rownames(res_memory_6_log[res_memory_6_log$log2FoldChange > 0,])))
m_vec <- unique(c(rownames(res_memory_1_log),rownames(res_memory_2_log),
                  rownames(res_memory_3_log),rownames(res_memory_4_log),
                  rownames(res_memory_5_log),rownames(res_memory_6_log)))
m_vec_symbol <- grch38$symbol[match(m_vec_pos,grch38$ensgene)]

res_naive_1_logP <- res_naive_1_log[res_naive_1_log$log2FoldChange > 0,]
res_naive_2_logP <- res_naive_2_log[res_naive_2_log$log2FoldChange > 0,]
res_naive_3_logP <- res_naive_3_log[res_naive_3_log$log2FoldChange > 0,]
res_naive_4_logP <- res_naive_4_log[res_naive_4_log$log2FoldChange > 0,]
res_naive_5_logP <- res_naive_5_log[res_naive_5_log$log2FoldChange > 0,]
res_naive_6_logP <- res_naive_6_log[res_naive_6_log$log2FoldChange > 0,]
res_memory_1_logP <- res_memory_1_log[res_memory_1_log$log2FoldChange > 0,]
res_memory_2_logP <- res_memory_2_log[res_memory_2_log$log2FoldChange > 0,]
res_memory_3_logP <- res_memory_3_log[res_memory_3_log$log2FoldChange > 0,]
res_memory_4_logP <- res_memory_4_log[res_memory_4_log$log2FoldChange > 0,]
res_memory_5_logP <- res_memory_5_log[res_memory_5_log$log2FoldChange > 0,]
res_memory_6_logP <- res_memory_6_log[res_memory_6_log$log2FoldChange > 0,]

numb_up <- data.frame(Shared = c(dim(res_naive_1_logP[rownames(res_naive_1_logP) %in% rownames(res_memory_1_logP),])[1],
                          dim(res_naive_2_logP[rownames(res_naive_2_logP) %in% rownames(res_memory_2_logP),])[1],
                          dim(res_naive_3_logP[rownames(res_naive_3_logP) %in% rownames(res_memory_3_logP),])[1],
                          dim(res_naive_4_logP[rownames(res_naive_4_logP) %in% rownames(res_memory_4_logP),])[1],
                          dim(res_naive_5_logP[rownames(res_naive_5_logP) %in% rownames(res_memory_5_logP),])[1],
                          dim(res_naive_6_logP[rownames(res_naive_6_logP) %in% rownames(res_memory_6_logP),])[1]),
                    Naive = c(dim(res_naive_1_logP[!rownames(res_naive_1_logP) %in% rownames(res_memory_1_logP),])[1],
                               dim(res_naive_2_logP[!rownames(res_naive_2_logP) %in% rownames(res_memory_2_logP),])[1],
                               dim(res_naive_3_logP[!rownames(res_naive_3_logP) %in% rownames(res_memory_3_logP),])[1],
                               dim(res_naive_4_logP[!rownames(res_naive_4_logP) %in% rownames(res_memory_4_logP),])[1],
                               dim(res_naive_5_logP[!rownames(res_naive_5_logP) %in% rownames(res_memory_5_logP),])[1],
                               dim(res_naive_6_logP[!rownames(res_naive_6_logP) %in% rownames(res_memory_6_logP),])[1]),
                    Memory = c(dim(res_memory_1_logP[!rownames(res_memory_1_logP) %in% rownames(res_naive_1_logP),])[1],
                          dim(res_memory_2_logP[!rownames(res_memory_2_logP) %in% rownames(res_naive_2_logP),])[1],
                          dim(res_memory_3_logP[!rownames(res_memory_3_logP) %in% rownames(res_naive_3_logP),])[1],
                          dim(res_memory_4_logP[!rownames(res_memory_4_logP) %in% rownames(res_naive_4_logP),])[1],
                          dim(res_memory_5_logP[!rownames(res_memory_5_logP) %in% rownames(res_naive_5_logP),])[1],
                          dim(res_memory_6_logP[!rownames(res_memory_6_logP) %in% rownames(res_naive_6_logP),])[1]),
                    s = c("hTCR hCD28","lTCR hCD28","hTCR lCD28","lTCR lCD28","hTCR","hCD28"))
numb_up <- reshape2::melt(numb_up)
ggplot(numb_up, aes(y=value,x=s,fill=variable)) +
  geom_bar(stat = "identity", color="black") +
  theme_classic(base_size=20) +
  scale_fill_manual(name="cell type",values = c("#BFC2C1","#98cdac","#e5b7d6")) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(size=16)) +
  coord_flip() +
  #facet_grid(model~.,scales = "free_y",space = "free_y") +
  xlab("") + ylab("# upregulated genes") 

numb_test <- data.frame(n = c(dim(res_naive_1_log[res_naive_1_log$log2FoldChange > 0,])[1], dim(res_naive_2_log[res_naive_2_log$log2FoldChange > 0,])[1], dim(res_naive_3_log[res_naive_3_log$log2FoldChange > 0,])[1],
                          dim(res_naive_4_log[res_naive_4_log$log2FoldChange > 0,])[1], dim(res_naive_5_log[res_naive_5_log$log2FoldChange > 0,])[1], dim(res_naive_6_log[res_naive_6_log$log2FoldChange > 0,])[1]),
                    m = c(dim(res_memory_1_log[res_memory_1_log$log2FoldChange > 0,])[1], dim(res_memory_2_log[res_memory_2_log$log2FoldChange > 0,])[1], dim(res_memory_3_log[res_memory_3_log$log2FoldChange > 0,])[1],
                          dim(res_memory_4_log[res_memory_4_log$log2FoldChange > 0,])[1], dim(res_memory_5_log[res_memory_5_log$log2FoldChange > 0,])[1], dim(res_memory_6_log[res_memory_6_log$log2FoldChange > 0,])[1]),
                    s = c("1","2","3","4","5","6"))
t.test(numb_test$n,
       numb_test$m,
       paired=TRUE)

data1 <- plotCounts(dds, gene="ENSG00000109471", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data1$rlog <- combat_edata["ENSG00000109471",]
p1 <- ggplot(data1, aes(x=stimulation, y=rlog,fill=stimulation)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("IL2") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  facet_grid(cell_type~., scales = "free") +
  scale_fill_manual(values = c("resting"="#66a61e",
                               "lTCR_lCD28"="#a6761d",
                               "lTCR_hCD28"="#d95f02",
                               "hTCR_lCD28"="#e7298a",
                               "hTCR"="#7570b3",
                               "hCD28"="#e6ab02",
                               "hTCR_hCD28"="#56B4E9"),
                    guide = guide_legend(title = "Stimulation"))
  
data2 <- plotCounts(dds, gene="ENSG00000163599", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data2$rlog <- combat_edata["ENSG00000163599",]
p2 <- ggplot(data2, aes(x=stimulation, y=rlog,fill=stimulation)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("CTLA4") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  facet_grid(cell_type~., scales = "free") +
  scale_fill_manual(values = c("resting"="#66a61e",
                               "lTCR_lCD28"="#a6761d",
                               "lTCR_hCD28"="#d95f02",
                               "hTCR_lCD28"="#e7298a",
                               "hTCR"="#7570b3",
                               "hCD28"="#e6ab02",
                               "hTCR_hCD28"="#56B4E9"),
                    guide = guide_legend(title = "Stimulation"))
data3 <- plotCounts(dds, gene="ENSG00000178562", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data3$rlog <- combat_edata["ENSG00000178562",]
p3 <- ggplot(data3, aes(x=stimulation, y=rlog,fill=stimulation)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("CD28") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  facet_grid(cell_type~., scales = "free") +
  scale_fill_manual(values = c("resting"="#66a61e",
                               "lTCR_lCD28"="#a6761d",
                               "lTCR_hCD28"="#d95f02",
                               "hTCR_lCD28"="#e7298a",
                               "hTCR"="#7570b3",
                               "hCD28"="#e6ab02",
                               "hTCR_hCD28"="#56B4E9"),
                    guide = guide_legend(title = "Stimulation"))
data4 <- plotCounts(dds, gene="ENSG00000114013", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data4$rlog <- combat_edata["ENSG00000114013",]
p4 <- ggplot(data4, aes(x=stimulation, y=rlog,fill=stimulation)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("CD86") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  facet_grid(cell_type~., scales = "free") +
  scale_fill_manual(values = c("resting"="#66a61e",
                               "lTCR_lCD28"="#a6761d",
                               "lTCR_hCD28"="#d95f02",
                               "hTCR_lCD28"="#e7298a",
                               "hTCR"="#7570b3",
                               "hCD28"="#e6ab02",
                               "hTCR_hCD28"="#56B4E9"),
                    guide = guide_legend(title = "Stimulation"))
data5 <- plotCounts(dds, gene="ENSG00000163600", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data5$rlog <- combat_edata["ENSG00000163600",]
p5 <- ggplot(data5, aes(x=stimulation, y=rlog,fill=stimulation)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("ICOS") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  facet_grid(cell_type~., scales = "free") +
  scale_fill_manual(values = c("resting"="#66a61e",
                               "lTCR_lCD28"="#a6761d",
                               "lTCR_hCD28"="#d95f02",
                               "hTCR_lCD28"="#e7298a",
                               "hTCR"="#7570b3",
                               "hCD28"="#e6ab02",
                               "hTCR_hCD28"="#56B4E9"),
                    guide = guide_legend(title = "Stimulation"))

multiplot(p1,p2,p3,p5,p4, cols = 5)