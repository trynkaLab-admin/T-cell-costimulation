library(DESeq2)

naive_cells <- cost_table[ , grepl( "naive" , colnames( cost_table ) ) ]

sample.info <- read.table("costim.samples_info.txt",header=T, sep="\t")
rownames(sample.info) <- sample.info$name
sample.info <- sample.info [,-c(1)]
sample.info <- sample.info[row.names(sample.info) %in% colnames(naive_cells),]

# Form a data frame analogous to expression data that will hold the clinical traits.
traitRows = match(colnames(naive_cells), rownames(sample.info));
sample.info = sample.info[traitRows,];
sample.info$sample <- as.factor(sample.info$sample)

# Generate the DESeqDataSet
DESeq.dsN <- DESeqDataSetFromMatrix(countData = naive_cells,
                                    colData = sample.info,
                                    design = ~ sample + stimulation)
dds_N <- DESeq(DESeq.dsN)

design(dds_N) <- formula(~  sample + CD28_level)
dds_N <- DESeq(dds_N, test="LRT", reduced = ~ sample)
resCD28_N <- results(dds_N)
resCD28_N$symbol <- grch38$symbol[match(rownames(resCD28_N),grch38$ensgene)]
write.csv(resCD28_N,"resCD28_N_linear.csv")
resCD28_N_q <- subset(resCD28_N,padj <= 0.05)
resCD28_N_log <- subset(resCD28_N_q,log2FoldChange >= 0.584)

design(dds_N) <- formula(~  sample + TCR_level)
dds_N <- DESeq(dds_N, test="LRT", reduced = ~ sample)
resTCR_N <- results(dds_N)
resTCR_N$symbol <- grch38$symbol[match(rownames(resTCR_N),grch38$ensgene)]
write.csv(resTCR_N,"resTCR_N_linear.csv")
resTCR_N_q <- subset(resTCR_N,padj <= 0.05)
resTCR_N_log <- subset(resTCR_N_q,log2FoldChange >= 0.584)

design(dds_N) <- formula(~sample + CD28_presence)
dds_N <- DESeq(dds_N, test="LRT", reduced = ~ sample)
resCD28_N_pres <- results(dds_N)
resCD28_N_pres$symbol <- grch38$symbol[match(rownames(resCD28_N_pres),grch38$ensgene)]
write.csv(resCD28_N_pres,"resCD28_N_switch.csv")
resCD28_N_pres_q <- subset(resCD28_N_pres,padj <= 0.05)
resCD28_N_pres_log <- subset(resCD28_N_pres_q,log2FoldChange >= 1)

resCD28_N_co <- resCD28_N_pres_log[!rownames(resCD28_N_pres_log) %in% rownames(resCD28_N_log),]

design(dds_N) <- formula(~sample + TCR_presence)
dds_N <- DESeq(dds_N, test="LRT", reduced = ~ sample)
resTCR_N_pres <- results(dds_N)
resTCR_N_pres$symbol <- grch38$symbol[match(rownames(resTCR_N_pres),grch38$ensgene)]
write.csv(resTCR_N_pres,"resTCR_N_switch.csv")
resTCR_N_pres_q <- subset(resTCR_N_pres,padj <= 0.05)
resTCR_N_pres_log <- subset(resTCR_N_pres_q,log2FoldChange >= 1)

resTCR_N_co <- resTCR_N_pres_log[!rownames(resTCR_N_pres_log) %in% rownames(resTCR_N_log),]

resTCR_N_all <- c(rownames(resTCR_N_co),rownames(resTCR_N_log))
resCD28_N_all <- c(rownames(resCD28_N_co),rownames(resCD28_N_log))

resTCR_N_final <- resTCR_N_all[!resTCR_N_all %in% resCD28_N_all]
resCD28_N_final <- resCD28_N_all[!resCD28_N_all %in% resTCR_N_all]

memory_cells <- cost_table[ , grepl( "memory" , colnames( cost_table ) ) ]

sample.info <- read.table("costim.samples_info.txt",header=T, sep="\t")
rownames(sample.info) <- sample.info$name
sample.info <- sample.info [,-c(1)]

# Form a data frame analogous to expression data that will hold the clinical traits.
traitRows = match(colnames(memory_cells), rownames(sample.info));
sample.info = sample.info[traitRows,];
sample.info$stimulation <- as.factor(sample.info$stimulation)
sample.info$sample <- as.factor(sample.info$sample)
sample.info$age <- as.factor(sample.info$age)

# Generate the DESeqDataSet
DESeq.dsM <- DESeqDataSetFromMatrix(countData = memory_cells,
                                    colData = sample.info,
                                    design = ~ sample + stimulation)

dds_M <- DESeq(DESeq.dsM)
design(dds_M) <- formula(~  sample + CD28_level)
dds_M <- DESeq(dds_M, test="LRT", reduced = ~ sample)
resCD28_M <- results(dds_M)
resCD28_M$symbol <- grch38$symbol[match(rownames(resCD28_M),grch38$ensgene)]
write.csv(resCD28_M,"resCD28_M_linear.csv")
resCD28_M_q <- subset(resCD28_M,padj <= 0.05)
resCD28_M_log <- subset(resCD28_M_q,log2FoldChange >= 0.584)

design(dds_M) <- formula(~  sample + TCR_level)
dds_M <- DESeq(dds_M, test="LRT", reduced = ~ sample)
resTCR_M <- results(dds_M)
resTCR_M$symbol <- grch38$symbol[match(rownames(resTCR_M),grch38$ensgene)]
write.csv(resTCR_M,"resTCR_M_linear.csv")
resTCR_M_q <- subset(resTCR_M,padj <= 0.05)
resTCR_M_log <- subset(resTCR_M_q,log2FoldChange >= 0.584)

design(dds_M) <- formula(~sample + CD28_presence)
dds_M <- DESeq(dds_M, test="LRT", reduced = ~ sample)
resCD28_M_pres <- results(dds_M)
resCD28_M_pres$symbol <- grch38$symbol[match(rownames(resCD28_M_pres),grch38$ensgene)]
write.csv(resCD28_M_pres,"resCD28_M_switch.csv")
resCD28_M_pres_q <- subset(resCD28_M_pres,padj <= 0.05)
resCD28_M_pres_log <- subset(resCD28_M_pres_q,log2FoldChange >= 1)

resCD28_M_co <- resCD28_M_pres_log[!rownames(resCD28_M_pres_log) %in% rownames(resCD28_M_log),]

design(dds_M) <- formula(~sample + TCR_presence)
dds_M <- DESeq(dds_M, test="LRT", reduced = ~ sample)
resTCR_M_pres <- results(dds_M)
resTCR_M_pres$symbol <- grch38$symbol[match(rownames(resTCR_M_pres),grch38$ensgene)]
write.csv(resTCR_M_pres,"resTCR_M_switch.csv")
resTCR_M_pres_q <- subset(resTCR_M_pres,padj <= 0.05)
resTCR_M_pres_log <- subset(resTCR_M_pres_q,log2FoldChange >= 1)

resTCR_M_co <- resTCR_M_pres_log[!rownames(resTCR_M_pres_log) %in% rownames(resTCR_M_log),]

resTCR_M_all <- c(rownames(resTCR_M_co),rownames(resTCR_M_log))
resCD28_M_all <- c(rownames(resCD28_M_co),rownames(resCD28_M_log))

resTCR_M_final <- resTCR_M_all[!resTCR_M_all %in% resCD28_M_all]
resCD28_M_final <- resCD28_M_all[!resCD28_M_all %in% resTCR_M_all]

## CD28-dependant
data <- plotCounts(dds, gene="ENSG00000169245",intgroup=c("TCR_level","CD28_level","cell_type"), returnData=TRUE)
data$logcount <- combat_edata["ENSG00000169245",]
data$TCR_level <- as.factor(data$TCR_level)
data$CD28_level <- as.factor(data$CD28_level)
levels(data$TCR_level) = c("no TCR", "low TCR", "high TCR")
levels(data$CD28_level) = c("0", "0.04", "0.4")
#data <- data[data$cell_type == "Memory",]
p1 <- ggplot(data,aes(x=as.factor(CD28_level), y=logcount)) + 
  geom_violin(aes(fill = cell_type),draw_quantiles = 0.5) +
  theme_classic(base_size = 20) +
  facet_grid(cell_type~., scales = "free") +
  ggtitle("CXCL10") +
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("CD28 (cell:CHO ratio)") + ylab("log2 counts") +
  scale_fill_manual(values = c("#e5b7d6","#98cdac"))+
  scale_color_manual(values = c("#e5b7d6","#98cdac"))

data <- plotCounts(dds, gene="ENSG00000169194",intgroup=c("TCR_level","CD28_level","cell_type"), returnData=TRUE)
data$logcount <- combat_edata["ENSG00000169194",]
data$TCR_level <- as.factor(data$TCR_level)
data$CD28_level <- as.factor(data$CD28_level)
levels(data$TCR_level) = c("no TCR", "low TCR", "high TCR")
levels(data$CD28_level) = c("0", "0.04", "0.4")
#data <- data[data$cell_type == "Memory",]
p2 <- ggplot(data,aes(x=as.factor(CD28_level), y=logcount)) + 
  geom_violin(aes(fill = cell_type),draw_quantiles = 0.5) +
  theme_classic(base_size = 20) +
  facet_grid(cell_type~., scales = "free") +
  ggtitle("IL13") +
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("CD28 (cell:CHO ratio)") + ylab("log2 counts") +
  scale_fill_manual(values = c("#e5b7d6","#98cdac"))+
  scale_color_manual(values = c("#e5b7d6","#98cdac"))

data <- plotCounts(dds, gene="ENSG00000112116",intgroup=c("TCR_level","CD28_level","cell_type"), returnData=TRUE)
data$logcount <- combat_edata["ENSG00000112116",]
data$TCR_level <- as.factor(data$TCR_level)
data$CD28_level <- as.factor(data$CD28_level)
levels(data$TCR_level) = c("no TCR", "low TCR", "high TCR")
levels(data$CD28_level) = c("0", "0.04", "0.4")
#data <- data[data$cell_type == "Memory",]
p3 <- ggplot(data,aes(x=as.factor(CD28_level), y=logcount)) + 
  geom_violin(aes(fill = cell_type),draw_quantiles = 0.5) +
  theme_classic(base_size = 20) +
  ggtitle("IL17F") +
  facet_grid(cell_type~., scales = "free") +
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("CD28 (cell:CHO ratio)") + ylab("log2 counts") +
  scale_fill_manual(values = c("#e5b7d6","#98cdac"))+
  scale_color_manual(values = c("#e5b7d6","#98cdac"))

data <- plotCounts(dds, gene="ENSG00000136244",intgroup=c("TCR_level","CD28_level","cell_type"), returnData=TRUE)
data$logcount <- combat_edata["ENSG00000136244",]
data$TCR_level <- as.factor(data$TCR_level)
data$CD28_level <- as.factor(data$CD28_level)
levels(data$TCR_level) = c("no TCR", "low TCR", "high TCR")
levels(data$CD28_level) = c("0", "0.04", "0.4")
#data <- data[data$cell_type == "Memory",]
p4 <- ggplot(data,aes(x=as.factor(CD28_level), y=logcount)) + 
  geom_violin(aes(fill = cell_type),draw_quantiles = 0.5) +
  theme_classic(base_size = 20) +
  ggtitle("IL6") +
  facet_grid(cell_type~., scales = "free") +
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("CD28 (cell:CHO ratio)") + ylab("log2 counts") +
  scale_fill_manual(values = c("#e5b7d6","#98cdac"))+
  scale_color_manual(values = c("#e5b7d6","#98cdac"))

data <- plotCounts(dds, gene="ENSG00000120217",intgroup=c("TCR_level","CD28_level","cell_type"), returnData=TRUE)
data$logcount <- combat_edata["ENSG00000120217",]
data$TCR_level <- as.factor(data$TCR_level)
data$CD28_level <- as.factor(data$CD28_level)
levels(data$TCR_level) = c("no TCR", "low TCR", "high TCR")
levels(data$CD28_level) = c("0", "0.04", "0.4")
#data <- data[data$cell_type == "Memory",]
p5 <- ggplot(data,aes(x=as.factor(CD28_level), y=logcount)) + 
  geom_violin(aes(fill = cell_type),draw_quantiles = 0.5) +
  theme_classic(base_size = 20) +
  ggtitle("CD274") +
  facet_grid(cell_type~., scales = "free") +
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("CD28 (cell:CHO ratio)") + ylab("log2 counts") +
  scale_fill_manual(values = c("#e5b7d6","#98cdac"))+
  scale_color_manual(values = c("#e5b7d6","#98cdac"))

multiplot(p1,p4,p5,p2,p3, cols = 5)