library(fgsea)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)

hs <- org.Hs.eg.db
hs_df <- select(hs, keys = rownames(resTCR_M),
       columns = c("ENTREZID", "SYMBOL", "ENSEMBL"),
       keytype = "ENSEMBL")

resfc <- data.frame(resTCR_M_path = as.data.frame(resTCR_M)[,2],
                    resTCR_N_path = as.data.frame(resTCR_N)[,2],
                    resCD28_M_path = as.data.frame(resCD28_M)[,2],
                    resCD28_N_path = as.data.frame(resCD28_N)[,2],
                    row.names = rownames(as.data.frame(resTCR_M)))
resfc$threshold_naive <- ifelse(rownames(resfc) %in% resTCR_N_final, "TCR", 
                                ifelse(rownames(resfc) %in% resCD28_N_final, "CD28", "neither"))
resfc$threshold_memory <- ifelse(rownames(resfc) %in% resTCR_M_final, "TCR",
                                 ifelse(rownames(resfc) %in% resCD28_M_final, "CD28", "neither"))
resfc$Entrez <- hs_df$ENTREZID[match(rownames(resfc), hs_df$ENSEMBL)]

load("human_H_v5p2.rdata")
pathwaysH <- Hs.H
#Conduct analysis
resfc_temp <- na.omit(resfc)

resCD28_N_path <- resfc_temp$resCD28_N_path
names(resCD28_N_path) <- resfc_temp$Entrez
resCD28_N_path <- sort(resCD28_N_path, decreasing = T)
fgseaRes_CD28_N <- fgsea(pathwaysH, resCD28_N_path, minSize=5, maxSize = 500, nperm=1000000)
fgseaRes_CD28_N$Stimulation <- "CD28 sensitive"
fgseaRes_CD28_N$Cell <- "Naive"

resCD28_M_path <- resfc_temp$resCD28_M_path
names(resCD28_M_path) <- resfc_temp$Entrez
resCD28_M_path <- sort(resCD28_M_path, decreasing = T)
fgseaRes_CD28_M <- fgsea(pathwaysH, resCD28_M_path, minSize=5, maxSize = 500, nperm=1000000)
fgseaRes_CD28_M$Stimulation <- "CD28 sensitive"
fgseaRes_CD28_M$Cell <- "Memory"

resTCR_N_path <- resfc_temp$resTCR_N_path
names(resTCR_N_path) <- resfc_temp$Entrez
resTCR_N_path <- sort(resTCR_N_path, decreasing = T)
fgseaRes_TCR_N <- fgsea(pathwaysH, resTCR_N_path, minSize=5, maxSize = 500, nperm=1000000)
fgseaRes_TCR_N$Stimulation <- "TCR sensitive"
fgseaRes_TCR_N$Cell <- "Naive"

resTCR_M_path <- resfc_temp$resTCR_M_path
names(resTCR_M_path) <- resfc_temp$Entrez
resTCR_M_path <- sort(resTCR_M_path, decreasing = T)
fgseaRes_TCR_M <- fgsea(pathwaysH, resTCR_M_path, minSize=5, maxSize = 500, nperm=1000000)
fgseaRes_TCR_M$Stimulation <- "TCR sensitive"
fgseaRes_TCR_M$Cell <- "Memory"

fgseaRes <- rbind(fgseaRes_CD28_N,fgseaRes_CD28_M,fgseaRes_TCR_N,fgseaRes_TCR_M)
pathways <- unique(fgseaRes[fgseaRes$pval <= 0.01 & fgseaRes$ES > 0.5,]$pathway)
fgseaRes_sub <- fgseaRes[fgseaRes$pathway %in% pathways,]
fgseaRes_sub <- fgseaRes_sub[!fgseaRes_sub$pathway=="HALLMARK_MYC_TARGETS_V2",]
fgseaRes_sub$pathway <- gsub("HALLMARK_","",fgseaRes_sub$pathway)
ggplot(fgseaRes_sub, aes(y=ES,x=pathway,fill=Cell)) +
  geom_bar(stat="identity", position="dodge", color ="#1D1D1B") +
  geom_hline(yintercept=0.5, color="red") +
  coord_flip() +
  theme_classic(base_size=20) +
  facet_grid(.~Stimulation,scales = "free",space = "free") +
  ylab("Enrichment score") + xlab("") +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#e5b7d6","#98cdac"))

ggplot(fgseaRes_sub, aes(y=-log10(pval),x=pathway,fill=Cell)) +
  geom_bar(stat="identity", position="dodge", color ="#1D1D1B") +
  geom_hline(yintercept=2, color="red") +
  coord_flip() +
  theme_classic(base_size=20) +
  facet_grid(.~Stimulation,scales = "free",space = "free") +
  ylab("-LOG10 p-value") + xlab("") +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = c("#e5b7d6","#98cdac"))

altplotEnrichment <- function (pathway, stats, gseaParam = 1, mycolor) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) +
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = mycolor) +
    geom_line(color = mycolor) +
    theme_classic(base_size=16) + 
    ylim(-0.2,0.75) +
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = 0.2, color = mycolor) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "Rank", y = "Enrichment score")
  return(g)
}

altplotEnrichment(pathwaysH[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], resCD28_M_path, mycolor="#e5b8d7")
altplotEnrichment(pathwaysH[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], resCD28_N_path, mycolor="#99cdad")
altplotEnrichment(pathwaysH[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], resTCR_M_path, mycolor="#e5b8d7")
altplotEnrichment(pathwaysH[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], resTCR_N_path, mycolor="#99cdad")

data1 <- plotCounts(dds, gene="ENSG00000065328", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data1$rlog <- combat_edata["ENSG00000065328",]
data1 <- data1[data1$stimulation %in% c("hCD28","hTCR", "resting"),]
p1 <- ggplot(data1, aes(x=stimulation, y=rlog,fill=cell_type)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("MCM10") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  scale_fill_manual(values = c("Naive" = "#98cdac", "Memory" = "#e5b7d6"))

data2 <- plotCounts(dds, gene="ENSG00000072274", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data2$rlog <- combat_edata["ENSG00000072274",]
data2 <- data2[data2$stimulation %in% c("hCD28","hTCR", "resting"),]
p2 <- ggplot(data2, aes(x=stimulation, y=rlog,fill=cell_type)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("TFRC") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  scale_fill_manual(values = c("Naive" = "#98cdac", "Memory" = "#e5b7d6"))

data3 <- plotCounts(dds, gene="ENSG00000149554", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data3$rlog <- combat_edata["ENSG00000149554",]
data3 <- data3[data3$stimulation %in% c("hCD28","hTCR", "resting"),]
p3 <- ggplot(data3, aes(x=stimulation, y=rlog,fill=cell_type)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("CHEK1") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  scale_fill_manual(values = c("Naive" = "#98cdac", "Memory" = "#e5b7d6"))

data4 <- plotCounts(dds, gene="ENSG00000094804", intgroup=c("stimulation","cell_type"), returnData=TRUE)
data4$rlog <- combat_edata["ENSG00000094804",]
data4 <- data4[data4$stimulation %in% c("hCD28","hTCR", "resting"),]
p4 <- ggplot(data4, aes(x=stimulation, y=rlog,fill=cell_type)) +
  geom_violin(draw_quantiles = c(0.5)) +
  theme_classic(base_size = 20) +
  ggtitle("CDC6") +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 20)) +
  ylab("log2 counts") + xlab("") +
  scale_fill_manual(values = c("Naive" = "#98cdac", "Memory" = "#e5b7d6"))

multiplot(p1,p2,p3,p4, cols = 1)