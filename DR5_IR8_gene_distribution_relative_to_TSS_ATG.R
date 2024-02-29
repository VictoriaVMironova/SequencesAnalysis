library("biomaRt")
library("dplyr")


auxin_genes <- data.frame(
  ID = c("AT1G70560", "AT1G23320", "AT4G24670", "AT4G32540", "AT4G13260", "AT1G04610", "AT5G11320", "AT5G43890", "AT5G25620",
  "AT2G33230", "AT4G28720", "AT1G04180", "AT1G48910", "AT1G21430"),
  GeneName = c("TAA1", "TAR1", "TAR2", "YUC1", "YUC2", "YUC3", "YUC4", "YUC5", "YUC6", "YUC7", "YUC8", "YUC9", "YUC10", "YUC11")
)


# 5'UTR data processing (not necessary) -----------------------------------

mart <- useMart('plants_mart', host="plants.ensembl.org", dataset='athaliana_eg_gene')

IR8_5utr <- getBM(
  attributes=c("ensembl_gene_id", "5utr", "5_utr_start", "5_utr_end", "chromosome_name", 
               "strand", "start_position", "end_position"),
  filters=c("ensembl_gene_id"),
  values=auxin_genes$ID,
  mart=mart,
  checkFilters=FALSE,
  bmHeader=TRUE, uniqueRows = T
)

IR8_5utr_plus <- IR8_5utr[IR8_5utr$`5' UTR start` %in% IR8_5utr$`Gene start (bp)`, ]
IR8_5utr_minus <- IR8_5utr[IR8_5utr$`5' UTR end` %in% IR8_5utr$`Gene end (bp)`, ]

setdiff_ir8 <- setdiff(setdiff(auxin_genes, IR8_5utr_minus$`Gene stable ID`), IR8_5utr_plus$`Gene stable ID`)
setdiff_ir8 <- IR8_5utr[IR8_5utr$`Gene stable ID` %in% setdiff_ir8, ]
setdiff_ir8 <- setdiff_ir8[c(1:3), ]

IR8_5utr <- rbind(IR8_5utr_plus, IR8_5utr_minus, setdiff_ir8)
IR8_5utr <- arrange(IR8_5utr, IR8_5utr$`GeneName`)


# loading promoter sequences ----------------------------------------------
load("E:/Datasets/upstream_seq_A.thaliana.Rdata")
#IR8
IR8_seq2000 <- right_join(seq2000, auxin_genes, by = c(`Gene stable ID` = "ID"))

# IR8_seq2000 <- seq2000[seq2000$`Gene stable ID` %in% auxin_genes$ID, ]
# IR8_seq2000 <- left_join(IR8_seq2000, auxin_genes, by = c(`Gene stable ID`, "ID"))

IR8_seq2000 <- arrange(IR8_seq2000, IR8_seq2000$`GeneName`)
IR8_seq2000 <- split(IR8_seq2000, IR8_seq2000$`Gene stable ID`)

# IR8_5utr$Length <- nchar(IR8_5utr$`5' UTR`)
# IR8_5utr <- arrange(IR8_5utr, IR8_5utr$`Gene stable ID`)
# ups1500_5utr_IR8 <- data.frame(seq = paste(IR8_seq1500$`Flank (Gene)`, IR8_5utr$`5' UTR`), 
#                            ID = IR8_5utr$`Gene stable ID`)
# ups1500_5utr_IR8 <- arrange(ups1500_5utr_IR8, ID)


#ups1500_5utr_DR5 <- split(ups1500_5utr_DR5, ups1500_5utr_DR5$ID)

#find IR8 for 21 genes index 1 and find DR5 for 18 genes index 2
#index 1

tgtcnn_pos1 <- lapply(IR8_seq2000, function(x) unlist(gregexpr("TGTC..", x$`Flank (Gene)`)))
tgtcnn_pos1 <- lapply(tgtcnn_pos1, function(x) data.frame(x, direction = c("TRUE")))
nngaca_pos1 <- lapply(IR8_seq2000, function(x) unlist(gregexpr("..GACA", x$`Flank (Gene)`)))
nngaca_pos1<- lapply(nngaca_pos1, function(x) data.frame(x, direction = c("FALSE")))

gene_pos1 <- Map(tgtcnn_pos1, nngaca_pos1, f = bind_rows)
gene_pos1 <- lapply(gene_pos1, data.frame)

for(i in seq_along(gene_pos1)){
   names(gene_pos1[[i]])[1] <- names(gene_pos1)[i]
}
gene_pos1 <- lapply(gene_pos1, function(x) data.frame(x, gene=rep(names(x)[1], each=length(x$direction))))

for(i in seq_along(gene_pos1)){
  colnames(gene_pos1[[i]])[1] <- "start"
}

gene_pos1 <- lapply(gene_pos1, function(x) arrange(x, start))
gene_pos1 <- do.call(rbind, gene_pos1)

add_vIR8 <- c(gene_pos1$start[2:257] - gene_pos1$start[1:256], 0) #find position when IR8 hide 
add_vIR8 <- c(which(add_vIR8 == 14), which(add_vIR8 == 14)+1)
add_vIR8 <- add_vIR8[order(add_vIR8)]
gene_pos1$repeats <- "all"
gene_pos1[add_vIR8, 4] <- "IR8"

row.names(gene_pos1) <- NULL


add_vDR5forward <- gene_pos1[gene_pos1$direction == "TRUE", ]
add_DR5v1<- c(add_vDR5forward$start[2:110] - add_vDR5forward$start[1:109], 0)
add_DR5v1 <- c(which(add_DR5v1== 11), which(add_DR5v1 == 11)+1)
add_vDR5forward$repeats <- "all"
add_vDR5forward[add_DR5v1,4] <- "DR5"

add_vDR5reverse <- gene_pos1[gene_pos1$direction == "FALSE", ]
add_DR5v2 <- c(add_vDR5reverse$start[2:110] - add_vDR5reverse$start[1:109], 0)
add_DR5v2 <- c(which(add_DR5v2 == 11), which(add_DR5v2 == 11)+1)
add_vDR5reverse$repeats <- "all"
add_vDR5reverse[add_DR5v2,4] <- "DR5"

gene_pos1 <- rbind(add_vDR5forward, add_vDR5reverse)

#find DR5 position 
#atg_dr5 <- DR5_5utr$Length+1500
#atg_ir8 <- IR8_5utr$Length+1500

gene_pos1$end <- gene_pos1$start+6
gene_pos2$end <- gene_pos2$start+6

IR8_gene_pos <- gene_pos1
IR8_gene_pos <- left_join(gene_pos1, auxin_genes, by = c("gene" = "ID"))
IR8_order <- levels(IR8_gene_pos$GeneName)

#gene_pos <- rbind(DR5_gene_pos, IR8_gene_pos)
#gene_pos$gene <- factor(gene_pos$gene, levels = c(DR5_order, IR8_order))

atg_repeats <- c(atg_dr5, atg_ir8)
high <- seq(0.95, 38.95, by = 1)

library("ggplot2")
library("gggenes")
library("svglite")
ggplot(IR8_gene_pos,
             aes(xmin = start, xmax = end,forward = direction, y = GeneName, fill = repeats, col = repeats))+
             scale_color_manual(values = c("black", "lightcoral", "darkorchid", "gray40"))+
             scale_fill_manual(values = c("black", "lightcoral", "darkorchid", "gray40"))+
             scale_y_discrete(limits = rev(levels(IR8_gene_pos$GeneName)))+
             geom_gene_arrow()+
             theme_genes()+
            theme(legend.position = "none")+
            scale_x_continuous(name = "Position relative to TSS (-1500 to TSS)",
                     breaks = seq(0, 2000, by = 500),
                     labels = seq(-2000, 0, by = 500)) +
            coord_cartesian(xlim = c(0, 2000))

ggsave("TGTC_TAA_YUC_genes2.jpg", height = 8, width = 11)

+
  
  geom_hline(yintercept = seq(1.5, 39.5, by = 1))
plot <-  plot + theme(legend.position = "none")+
         annotate("text", x = atg_repeats, y = high, label = rep("ATG", 39), size = 3.5)+
         geom_hline(yintercept = seq(1.5, 39.5, by = 1))

plot

