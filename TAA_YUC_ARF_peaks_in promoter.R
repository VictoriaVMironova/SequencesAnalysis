library(ggplot2)
library(dplyr)

df <- read.csv("TAA_YUC_ARFpeaks_summary.csv") 

df$Normalized_Start <- df$Promoter_start + 2500
df$Normalized_End <- df$Promoter_end + 2500

# Split data based on TF
df_arf2 <- df %>% filter(TF == "ARF2")
df_arf5 <- df %>% filter(TF == "ARF5/MP")

genes <- unique(df$Gene)

df_arf2 <- df_arf2 %>%
  mutate(Label_Pos = gene_mapping[Gene])  # Adjust this offset as needed

df_arf5 <- df_arf5 %>%
  mutate(Label_Pos = gene_mapping[Gene] + 0.4)  # Adjust this offset as needed

# Base plot with ARF2 data
p <- ggplot() +
  geom_segment(data = df_arf2, aes(x = Normalized_Start, xend = Normalized_End, y = gene_mapping[Gene]+ 0.1, yend = gene_mapping[Gene]+ 0.1, color = "ARF2"), size = 5) +
  geom_segment(data = df_arf5, aes(x = Normalized_Start, xend = Normalized_End, y = gene_mapping[Gene] + 0.5, yend = gene_mapping[Gene] + 0.5, color = "ARF5/MP"), size = 5) +
  geom_hline(data = df_arf2, aes(yintercept = gene_mapping[Gene]), color = "grey")+
  geom_hline(data = df_arf5, aes(yintercept = gene_mapping[Gene]), color = "grey")


# Add labels for ARF2
p <- p + geom_text(data = df_arf2, aes(x = (Normalized_Start + Normalized_End) / 2, y = Label_Pos, label = "ARF2"), color = "black", size = 3, vjust = -0.5)

# Add labels for ARF5
p <- p + geom_text(data = df_arf5, aes(x = (Normalized_Start + Normalized_End) / 2, y = Label_Pos, label = "ARF5/MP"), color = "black", size = 3, vjust = -0.5)

# Final touches
p <- p + theme_minimal() +
  labs(x = "Position relative to TSS (-2500 to +100)", y = "Gene") +
  theme(axis.text.x = element_text(hjust = 1)) +
  scale_x_continuous(name = "Position relative to TSS (-2500 to +100)",
                     breaks = seq(0, 2500, by = 500),
                     labels = seq(-2400, 100, by = 500)) +
  coord_cartesian(xlim = c(0, 2500)) +
  scale_color_manual(values = c("ARF2" = "lightgreen", "ARF5/MP" = "aquamarine3")) +
  scale_y_continuous(name = "Gene", breaks = seq_along(genes), labels = genes)

# p<-p+ geom_gene_arrow(IR8_gene_pos,
#                      aes(xmin = start, xmax = end,forward = direction, y = gene_mapping[GeneName], fill = repeats, col = repeats))

p

ggsave(filename = "TAA_YUC_ARF_DAP-Seq_mapping.jpg", p)
