#all upregulated genes in naive cells pairwise
all_naive <- unique(c(rownames(res_naive_1_log[res_naive_1_log$log2FoldChange > 0,]),
                      rownames(res_naive_2_log[res_naive_2_log$log2FoldChange > 0,]),
                      rownames(res_naive_3_log[res_naive_3_log$log2FoldChange > 0,]),
                      rownames(res_naive_4_log[res_naive_4_log$log2FoldChange > 0,]),
                      rownames(res_naive_5_log[res_naive_5_log$log2FoldChange > 0,]),
                      rownames(res_naive_6_log[res_naive_6_log$log2FoldChange > 0,])))
#all downregulated genes in naive cells pairwise
all_naive_dn <- unique(c(rownames(res_naive_1_log[res_naive_1_log$log2FoldChange < 0,]),
                         rownames(res_naive_2_log[res_naive_2_log$log2FoldChange < 0,]),
                         rownames(res_naive_3_log[res_naive_3_log$log2FoldChange < 0,]),
                         rownames(res_naive_4_log[res_naive_4_log$log2FoldChange < 0,]),
                         rownames(res_naive_5_log[res_naive_5_log$log2FoldChange < 0,]),
                         rownames(res_naive_6_log[res_naive_6_log$log2FoldChange < 0,])))

#all upregulated genes in naive cells that do not fit linear or switch models for either stimulus
ns_naive <- all_naive[!all_naive %in% resTCR_N_final]
ns_naive <- ns_naive[!ns_naive %in% resCD28_N_final]

#all upregulated genes in naive cells that do not fit linear or switch models for either stimulus and require both stimuli
both_naive <- ns_naive[ns_naive %in% unique(c(rownames(res_naive_1_log[res_naive_1_log$log2FoldChange > 0,]),
                                              rownames(res_naive_2_log[res_naive_2_log$log2FoldChange > 0,]),
                                              rownames(res_naive_3_log[res_naive_3_log$log2FoldChange > 0,]),
                                              rownames(res_naive_4_log[res_naive_4_log$log2FoldChange > 0,])))]

#all upregulated genes in naive cells that do not fit linear or switch models for either stimulus and require both stimuli
# and either stimulus alone can lead to upregulation
either_naive <- both_naive[both_naive %in% unique(c(rownames(res_naive_6_log[res_naive_6_log$log2FoldChange > 0,]),
                                                    rownames(res_naive_5_log[res_naive_5_log$log2FoldChange > 0,])))]

#upregulated genes in naive cells that only go up with hCD28 alone
cd28_naive <- ns_naive[ns_naive %in% rownames(res_naive_6_log[res_naive_6_log$log2FoldChange > 0,])]
cd28_naive <- cd28_naive[!cd28_naive %in% both_naive]
#upregulated genes in naive cells that only go up with hTCR alone
tcr_naive <- ns_naive[ns_naive %in% rownames(res_naive_5_log[res_naive_5_log$log2FoldChange > 0,])]
tcr_naive <- tcr_naive[!tcr_naive %in% both_naive]

cd28_naive1 <- cd28_naive[!cd28_naive %in% tcr_naive]
tcr_naive1 <- tcr_naive[!tcr_naive %in% cd28_naive]

either_alone_naive <- tcr_naive[tcr_naive %in% cd28_naive]
cd28_naive <- cd28_naive1
tcr_naive <- tcr_naive1

both_naive <- both_naive[!both_naive %in% either_naive]

not_de_naive <- rownames(res_naive_1)[!rownames(res_naive_1) %in% all_naive]
not_de_naive <- not_de_naive[!not_de_naive %in% all_naive_dn]
not_de_naive <- not_de_naive[!not_de_naive %in% resTCR_N_final]
not_de_naive <- not_de_naive[!not_de_naive %in% resCD28_N_final]

write.table(not_de_naive,"not_de_naive_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

all_memory <- unique(c(rownames(res_memory_1_log[res_memory_1_log$log2FoldChange > 0,]),
                           rownames(res_memory_2_log[res_memory_2_log$log2FoldChange > 0,]),
                           rownames(res_memory_3_log[res_memory_3_log$log2FoldChange > 0,]),
                           rownames(res_memory_4_log[res_memory_4_log$log2FoldChange > 0,]),
                           rownames(res_memory_5_log[res_memory_5_log$log2FoldChange > 0,]),
                           rownames(res_memory_6_log[res_memory_6_log$log2FoldChange > 0,])))
all_memory_dn <- unique(c(rownames(res_memory_1_log[res_memory_1_log$log2FoldChange < 0,]),
                           rownames(res_memory_2_log[res_memory_2_log$log2FoldChange < 0,]),
                           rownames(res_memory_3_log[res_memory_3_log$log2FoldChange < 0,]),
                           rownames(res_memory_4_log[res_memory_4_log$log2FoldChange < 0,]),
                           rownames(res_memory_5_log[res_memory_5_log$log2FoldChange < 0,]),
                           rownames(res_memory_6_log[res_memory_6_log$log2FoldChange < 0,])))

ns_memory <- all_memory[!all_memory %in% resTCR_M_final]
ns_memory <- ns_memory[!ns_memory %in% resCD28_M_final]
both_memory <- ns_memory[ns_memory %in% unique(c(rownames(res_memory_1_log[res_memory_1_log$log2FoldChange > 0,]),
                                                     rownames(res_memory_2_log[res_memory_2_log$log2FoldChange > 0,]),
                                                     rownames(res_memory_3_log[res_memory_3_log$log2FoldChange > 0,]),
                                                     rownames(res_memory_4_log[res_memory_4_log$log2FoldChange > 0,])))]
either_memory <- both_memory[both_memory %in% unique(c(rownames(res_memory_6_log[res_memory_6_log$log2FoldChange > 0,]),
                                                           rownames(res_memory_5_log[res_memory_5_log$log2FoldChange > 0,])))]
    
cd28_memory <- ns_memory[ns_memory %in% rownames(res_memory_6_log[res_memory_6_log$log2FoldChange > 0,])]
cd28_memory <- cd28_memory[!cd28_memory %in% both_memory]
tcr_memory <- ns_memory[ns_memory %in% rownames(res_memory_5_log[res_memory_5_log$log2FoldChange > 0,])]
tcr_memory <- tcr_memory[!tcr_memory %in% both_memory]

cd28_memory1 <- cd28_memory[!cd28_memory %in% tcr_memory]
tcr_memory1 <- tcr_memory[!tcr_memory %in% cd28_memory]

either_alone_memory <- tcr_memory[tcr_memory %in% cd28_memory]
cd28_memory <- cd28_memory1
tcr_memory <- tcr_memory1
both_memory <- both_memory[!both_memory %in% either_memory]

not_de_memory <- rownames(res_memory_1)[!rownames(res_memory_1) %in% all_memory]
not_de_memory <- not_de_memory[!not_de_memory %in% all_memory_dn]
not_de_memory <- not_de_memory[!not_de_memory %in% resTCR_M_final]
not_de_memory <- not_de_memory[!not_de_memory %in% resCD28_M_final]

write.table(not_de_memory,"not_de_memory_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

sensitivity_summary <- data.frame(Memory = c(length(resTCR_M_final), length(resCD28_M_final),
                                             length(both_memory), length(either_memory),
                                             length(tcr_memory), length(cd28_memory), length(either_alone_memory)),
                                  Naive = c(length(resTCR_N_final), length(resCD28_N_final),
                                             length(both_naive), length(either_naive),
                                             length(tcr_naive), length(cd28_naive), length(either_alone_naive)),
                                  Group = c("TCR-sensitive","CD28-sensitive","requires both stimuli",
                                            "either CD28 or TCR","hTCR alone","hCD28 alone","either alone"))

chi1 <- chisq.test(sensitivity_summary[,-3])

CD28_fish <-
  matrix(c(sum(sensitivity_summary[sensitivity_summary$Group=="CD28-sensitive",]$Memory),
           sum(sensitivity_summary[sensitivity_summary$Group=="CD28-sensitive",]$Naive),
           sum(sensitivity_summary[sensitivity_summary$Group!="CD28-sensitive",]$Memory),
           sum(sensitivity_summary[sensitivity_summary$Group!="CD28-sensitive",]$Naive)),
         nrow = 2,
         dimnames = list(Cell = c("Memory", "Naive"),
                         Stimulus = c("CD28", "Other")))
fish1 <- fisher.test(CD28_fish, alternative = "two.sided")

TCR_fish <-
  matrix(c(sum(sensitivity_summary[sensitivity_summary$Group=="TCR-sensitive",]$Memory),
           sum(sensitivity_summary[sensitivity_summary$Group=="TCR-sensitive",]$Naive),
           sum(sensitivity_summary[sensitivity_summary$Group!="TCR-sensitive",]$Memory),
           sum(sensitivity_summary[sensitivity_summary$Group!="TCR-sensitive",]$Naive)),
         nrow = 2,
         dimnames = list(Cell = c("Memory", "Naive"),
                         Stimulus = c("TCR", "Other")))
fish2 <- fisher.test(TCR_fish, alternative = "two.sided")

sensitivity_summary <- melt(sensitivity_summary)
ggplot(sensitivity_summary,aes(y=value,x=variable,fill=Group)) +
  geom_bar(stat = "identity", color="black", position="dodge") +
  theme_minimal(base_size=18) +
  #scale_fill_manual(name="cell type",values = c("#e5b7d6","#98cdac","#BFC2C1")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(size=16)) +
  #coord_flip() +
  xlab("") + ylab("# upregulated genes") +
  scale_fill_brewer(palette = "Set1")