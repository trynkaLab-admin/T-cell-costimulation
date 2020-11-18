library(annotables)
library(DESeq2)
library(ggplot2)

## load a count table a process it to appropriate format
read.counts <- read.table("costim.RNA_featureCounts.txt", header = TRUE )
row.names(read.counts) <- read.counts$Geneid                         # the gene IDs should be stored as row. names
read.counts$Chr <- gsub("\\;.*","",read.counts$Chr) # remove contig info from name
read.counts <- read.counts[!grepl("G", read.counts$Chr),]
read.counts <- read.counts[!grepl("K", read.counts$Chr),]
read.counts <- read.counts[!grepl("J", read.counts$Chr),]
read.counts <- subset(read.counts, !(Chr=="chrY"))      # Remove Y chromosome
row.names(read.counts) <- gsub("\\..*","",rownames(read.counts)) # remove contig info from name
read.counts <- read.counts [,-c(1:6)]                            # exclude all columns that do not contain read counts (chr, start, end, strand, length)
rc.keep <- rowSums(read.counts > 20) >= 3                           # remove genes for which only few reads could be mapped
cost_table <- read.counts [ rc.keep, ]
cost_table$biotype <- grch38$biotype[match(row.names(cost_table),grch38$ensgene)]
cost_table <- subset(cost_table,biotype=="protein_coding")
cost_table <- cost_table[,-c(75)]                            # exclude all columns that do not contain read counts (chr, start, end, strand, length)

###SCRIPT
sample.info <- read.table("costim.samples_info.txt",header=T, sep="\t")
rownames(sample.info) <- sample.info$name
sample.info <- sample.info [,-c(1)]

# Form a data frame analogous to expression data that will hold the clinical traits.
traitRows = match(colnames(cost_table), rownames(sample.info));
sample.info = sample.info[traitRows,];
sample.info$sample <- as.factor(sample.info$sample)

# Generate the DESeqDataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = cost_table,
                                   colData = sample.info,
                                   design = ~ cell_type)
dds <- DESeq(DESeq.ds)
dds$group <- factor(paste0(dds$stimulation, dds$cell_type))
design(dds) <- formula(~ sample + group)
dds = DESeq(dds)

# Regularized log transformation
rld <- rlogTransformation(dds, blind=TRUE)
rlog.norm.counts <- assay(rld)

# there are two types of variables that are being considered: (1) adjustment variables and (2) variables of interest
modcombat = model.matrix(~1, data=sample.info)
batch = sample.info$batch

combat_edata = sva::ComBat(dat=rlog.norm.counts, batch=batch, mod=modcombat)
temporary_file = c("costim.FeatureCounts.rlog_combat.rds")
saveRDS(combat_edata, temporary_file)
combat_edata = readRDS("costim.FeatureCounts.rlog_combat.rds")

PCA <- prcomp(t(combat_edata), scale = T)                      # use the function prcomp to compute the principal components
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
barplot(percentVar, main="variation explained", ylab="% of variation explained", xlab="principal components")
dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], PC3 = PCA$x[,3], PC4 = PCA$x[,4], PC5 = PCA$x[,5],
                    PC6 = PCA$x[,6], PC7 = PCA$x[,7], PC8 = PCA$x[,8], PC9 = PCA$x[,9], PC10 = PCA$x[,10],
                    PC11 = PCA$x[,11], PC12 = PCA$x[,12], PC13 = PCA$x[,13], PC14 = PCA$x[,14], PC15 = PCA$x[,15],
                    sampleNO = colData(rld))

p1 <- qplot(PC1, PC2, data = dataGG, color =  sampleNO.stimulation, shape = sampleNO.cell_type,
           main = "PC1 vs PC2", size = I(5)) +
  labs(x = paste0("PC1: ", round(percentVar[1],4),"%"),
       y = paste0("PC2: ", round(percentVar[2],4),"%")) +
  theme(strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=18),
        panel.border = element_rect(colour = "#d6d6d6", fill=NA, size=2)) +
  scale_color_manual(values = c("resting"="#66a61e",
                               "lTCR_lCD28"="#a6761d",
                               "lTCR_hCD28"="#d95f02",
                               "hTCR_lCD28"="#e7298a",
                               "hTCR"="#7570b3",
                               "hCD28"="#e6ab02",
                               "hTCR_hCD28"="#56B4E9"),
                    guide = guide_legend(title = "Stimulation")) +
  scale_shape(guide = guide_legend(title = "Cell type"))

ggsave("PCA_PC1_PC2.pdf", plot = p1, width = 16, height = 16, units = "cm", dpi = 576)

p2 <- qplot(PC3, PC4, data = dataGG, color =  sampleNO.stimulation, shape = sampleNO.cell_type,
      main = "PC3 vs PC4", size = I(5)) +
  labs(x = paste0("PC3: ", round(percentVar[3],4),"%"),
       y = paste0("PC4: ", round(percentVar[4],4),"%")) +
  theme(strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=18),
        panel.border = element_rect(colour = "#d6d6d6", fill=NA, size=2)) +
  scale_color_manual(values = c("resting"="#66a61e",
                                "lTCR_lCD28"="#a6761d",
                                "lTCR_hCD28"="#d95f02",
                                "hTCR_lCD28"="#e7298a",
                                "hTCR"="#7570b3",
                                "hCD28"="#e6ab02",
                                "hTCR_hCD28"="#56B4E9"),
                     guide = guide_legend(title = "Stimulation")) +
  scale_shape(guide = guide_legend(title = "Cell type"))
ggsave("PCA_PC3_PC4.pdf", plot = p2, width = 16, height = 16, units = "cm", dpi = 576)