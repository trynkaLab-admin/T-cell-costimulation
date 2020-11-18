# Calculate all the DEG between naive and memory cells

####### LOAD PACKAGES
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(reshape2)
library(ggrepel)
library(annotables)

####### FUNCTIONS
## Volcano plot with "significant" genes labeled
volcano <- function(df,lg2fc,pv,genes_label,xlim1,xlim2,ylim1,ylim2) {
  df$hgnc_symbol <- grch38$symbol[match(rownames(df),grch38$ensgene)]
  df <- as.data.frame(df)
  df$threshold = as.factor(abs(df$log2FoldChange) > lg2fc & df$padj < pv/nrow(df))
  df$legend <- as.factor(df$hgnc_symbol %in% genes_label) 
  ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), label=hgnc_symbol)) +
    geom_point(alpha=0.7, size=1.75, aes(color = threshold)) +
    scale_color_manual(values = c("grey","#89A8D9")) +
    theme_minimal(base_size=20) +
    xlim(xlim1,xlim2) + ylim(ylim1,ylim2) +
    theme(legend.title=element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.position = "none",
          panel.border = element_rect(colour = "grey50", fill=NA, size=1)) +
    xlab("LOG2 fold change") + ylab("-LOG10 p-value") +
    geom_text_repel(aes(label=ifelse(legend==TRUE,as.character(hgnc_symbol),'')),size=6, color="black") +
    geom_vline(xintercept = 0, color="grey50",linetype=2)
}

####### CODE
# Calculate all the DEG between naive and memory cells
res1 <- results(dds, contrast=c("group", "hTCR_hCD28Memory", "hTCR_hCD28Naive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res1_log <- subset(res1,abs(log2FoldChange) >= 1)
res1_q <- subset(res1_log,padj <= 0.05)
res1$symbol <- grch38$symbol[match(rownames(res1),grch38$ensgene)]
write.csv(res1,"naive_memory_hTCRhCD28_comp.csv",quote = FALSE)

res7 <- results(dds, contrast=c("group", "restingMemory", "restingNaive"), alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")
res7_log <- subset(res7,abs(log2FoldChange) >= 1)
res7_q <- subset(res7_log,padj <= 0.05)
res7$symbol <- grch38$symbol[match(rownames(res7),grch38$ensgene)]
write.csv(res7,"naive_memory_rest_comp.csv",quote = FALSE)

plot1 <- results(dds, contrast=c("group", "restingMemory", "restingNaive"))
plot2 <- results(dds, contrast=c("group", "hTCR_hCD28Memory", "hTCR_hCD28Naive"))
genes_label1 <- c("CD58","CXCR3","CXCR5","CCR2","CCR4","CCR5","CCR6","ITGB1",
                  "ITGB5","MAF","IL2RB","IL12RB2","IL18RAP","MCOLN2","SYT11",
                  "AIM2","ITPRIPL1","GPR15","MCF2L2","ADAM19","STOM","NPDC1",
                  "ST8SIA1","LGALS3","THBS1","CPNE7","CYB561","DUSP4","PTGER2",
                  "DUSP5","NOD2","CFH","GZMK","FAS","CDKN1A","TBX21","PECAM1",
                  "CACHD1","LRRN1","NREP","GNAI1","DNTT","SLC15A3","FLT1")
genes_label <- c("CD58","CXCR3","CXCR5","CCR2","CCR5","CCR6","ITGB1",
                 "MAF","MCOLN2","SYT11","CACHD1","PECAM1","GNAI1","SLC15A3",
                 "FLT1","AIM2","GPR15","MCF2L2","ADAM19","NPDC1",
                 "ST8SIA1","LGALS3","THBS1","PTGER2","NOD2","CFH",
                 "GZMK","CDKN1A")

volcano(plot1,1,0.05,genes_label,-3.5,6.5,0,95)
ggsave("volcanoPlot_memory_naive_rest.pdf", plot = last_plot(), width = 14, height = 16, units = "cm", dpi = 576)
volcano(plot2,1,0.05,genes_label,-3.5,6.5,0,95)
ggsave("volcanoPlot_memory_naive_stim.pdf", plot = last_plot(), width = 14, height = 16, units = "cm", dpi = 576)

## ALL GENES
rest_var_N <- apply(combat_edata[ , grepl( "naive_7" , colnames( combat_edata ) ) ], 1, var)
rest_var_M <- apply(combat_edata[ , grepl( "memory_7" , colnames( combat_edata ) ) ], 1, var)
rest_mean_N <- apply(combat_edata[ , grepl( "naive_7" , colnames( combat_edata ) ) ], 1, mean)
rest_mean_M <- apply(combat_edata[ , grepl( "memory_7" , colnames( combat_edata ) ) ], 1, mean)

rest_N <- as.data.frame(cbind(rest_var_N,rest_mean_N))
rest_N$cell <- "naive"
rest_N$stimulation <- "resting"
colnames(rest_N) <- c("var","meanExpr","cell","stimulation")
rest_N$gene <- rownames(rest_N)
rest_M <- as.data.frame(cbind(rest_var_M,rest_mean_M))
rest_M$cell <- "memory"
rest_M$stimulation <- "resting"
colnames(rest_M) <- c("var","meanExpr","cell","stimulation")
rest_M$gene <- rownames(rest_M)

rest_melt_var <- rbind(rest_N,rest_M)

wt <- wilcox.test(var ~ cell, data=rest_melt_var)

rest_melt_var <- cbind(rest_N$var,rest_M$var)
rownames(rest_melt_var) <- rest_M$gene
colnames(rest_melt_var) <- c("naive","memory")
rest_melt_var <- as.data.frame(rest_melt_var)
rest_melt_var$symbol <- grch38$symbol[match(rownames(rest_melt_var),grch38$ensgene)]

ggplot(rest_melt_var,aes(x=naive, y=memory)) +
  geom_point(alpha=0.1, size=3) +
  theme_bw(base_size=18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20)) +
  geom_text_repel(data = subset(rest_melt_var, naive > 1 & memory < 1 | memory > 1 & naive < 1),
                  aes(label = symbol), color="#89A8D9") +
  ylab("memory variance") + xlab("naive variance") +
  annotate("text", label = paste0("Wilcoxon rank sum test p-value: ",round(wt$p.value,digits=3)), x = 2.2, y = 2.7, size = 3.5, colour = "black") 
ggsave("variance_resting_memory_naive.pdf", plot = last_plot(), width = 10, height = 10, units = "cm", dpi = 576)
