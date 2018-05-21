library(ggplot2)
library(scales)

gef <- data.frame(read.delim("D:/schatz_project/comparisons/gtex_exact_filtered/compare.tsv"))
geu <- data.frame(read.delim("D:/schatz_project/comparisons/gtex_exact_unfiltered/compare.tsv"))
sef <- data.frame(read.delim("D:/schatz_project/comparisons/short_exact_filtered/compare.tsv"))
seu <- data.frame(read.delim("D:/schatz_project/comparisons/short_exact_unfiltered/compare.tsv"))
reu <- data.frame(read.delim("D:/schatz_project/comparisons/refseq_exact_unfiltered/compare.tsv"))

ggplot(gef, aes(lr_min_reads)) + 
  geom_line(aes(y=lr_recall, color="gtex_lr"), linetype=2, size=2) + geom_line(aes(y=target_recall, color="gtex_target"), linetype=2, size=2) +
  geom_line(aes(y=geu$lr_recall, color="gtex_lr"), size=2) + geom_line(aes(y=geu$target_recall, color="gtex_target"), size=2) +
  geom_line(aes(y=sef$lr_recall, color="short_lr"), linetype=2, size=2) + geom_line(aes(y=sef$target_recall, color="short_target"), linetype=2, size=2) +
  geom_line(aes(y=seu$lr_recall, color="short_lr"), size=2) + geom_line(aes(y=seu$target_recall, color="short_target"), size=2) +
  geom_line(aes(y=reu$lr_recall, color="refseq_lr"), size=2) + geom_line(aes(y=reu$target_recall, color="refseq_target"), size=2) +
  scale_x_log10(breaks = df1$lr_min_reads) +scale_y_continuous(breaks=c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00)) +
  xlab("Minimum # of reads supporting a junction") + ylab("Recall") + theme_bw() + labs(color="Color/Linetype:\ndashed=filtered\nsolid=unfiltered")
 