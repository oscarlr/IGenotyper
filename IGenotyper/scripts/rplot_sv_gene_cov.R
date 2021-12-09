#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

infn = args[1]
outfn = args[2]

x = read.table(infn)

#x <- melt(x, id=c("V1","V2","V3","V4","V5")) 

genes = c("IGHV7-4-1",
          "IGHV3-64D","IGHV5-10-1",
          "IGHV3-23","IGHV3-23D",
          "IGHV4-28","IGHV3-30","IGHV4-30-2","IGHV3-30-3","IGHV4-30-4","IGHV3-30-5","IGHV4-31",
          "IGHV4-38-2","IGHV3-43D","IGHV3-38-3","IGHV1-38-4",
          "IGHV1-69","IGHV2-70D","IGHV1-69-2","IGHV1-69D","IGHV2-70",
          "IGHV1-8","IGHV3-9")

sv = c("1",
       "2","2",
       "3","3",
       "4","4","4","4","4","4","4",
       "5","5","5","5",
       "6","6","6","6","6",
       "7","7")

x = x[x$V4 %in% genes,]
x = x[order(x$V2),]
x$sv = sv

x <- melt(x, id=c("V1","V2","V3","V4","V5","sv")) 

x$V6 = substr(x$V4,4,100)

ggplot(x,aes(reorder(V6,V2),value,fill=variable)) +
  geom_bar(stat="identity") +
  facet_grid(. ~ sv,scale = "free_x",space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
  ylab("Coverage") +
  xlab("Genes") +
  scale_fill_discrete(name = "hap",labels = c("unphased","1","2"))
  
ggsave(outfn, width = 16, height = 10,units = "cm",dpi=1000)
