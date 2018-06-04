library(ggplot2)
library(scales)

novel <- data.frame(read.delim("D:/schatz_project/comparisons/original_plus_4_wiggle.juncs.sorted.full.novel"))
g<-ggplot(novel, aes(x=reads_required, y=novel_junctions))
g+geom_line()+scale_y_log10(breaks=novel$novel_junctions) + scale_x_log10(breaks=novel$reads_required) +
  xlab("Minimum # of reads required to support a novel junction (log10)") + ylab("# of novel junctions (log10)") + theme_bw()

novel20_lens <- data.frame(read.delim("D:/schatz_project/comparisons/novel_20.unique_rads.ids.lengths"))


gef <- data.frame(read.delim("D:/schatz_project/comparisons/gtex_exact_filtered/compare.tsv"))
geu <- data.frame(read.delim("D:/schatz_project/comparisons/gtex_exact_unfiltered/compare.tsv"))
sef <- data.frame(read.delim("D:/schatz_project/comparisons/short_exact_filtered/compare.tsv"))
seu <- data.frame(read.delim("D:/schatz_project/comparisons/short_exact_unfiltered/compare.tsv"))
reu <- data.frame(read.delim("D:/schatz_project/comparisons/refseq_exact_unfiltered/compare.tsv"))

#sra
aef <- data.frame(read.delim("D:/schatz_project/comparisons/sra_exact_filtered/compare.tsv"))
aeu <- data.frame(read.delim("D:/schatz_project/comparisons/sra_exact_unfiltered/compare.tsv"))

#pacbio
pef <- data.frame(read.delim("D:/schatz_project/comparisons/pacbio_exact_filtered/compare.tsv"))
peu <- data.frame(read.delim("D:/schatz_project/comparisons/pacbio_exact_unfiltered/compare.tsv"))


#filtered
ggplot(gef, aes(lr_min_reads)) + 
  geom_line(aes(y=lr_recall, color="gtex_lr"), linetype=2, size=2) + geom_line(aes(y=target_recall, color="gtex_target"), linetype=2, size=2) +
  geom_line(aes(y=sef$lr_recall, color="short_lr"), linetype=2, size=2) + geom_line(aes(y=sef$target_recall, color="short_target"), linetype=2, size=2) +
  geom_line(aes(y=aef$lr_recall, color="sra_lr"), linetype=2, size=2) + geom_line(aes(y=aef$target_recall, color="sra_target"), linetype=2, size=2) +
  geom_line(aes(y=pef$lr_recall, color="pacbio_lr"), linetype=2, size=2) + geom_line(aes(y=pef$target_recall, color="pacbio_target"), linetype=2, size=2) +
  geom_line(aes(y=pef$lr_recall, color="pacbio_lr"), linetype=2, size=2) + geom_line(aes(y=pef$target_recall, color="pacbio_target"), linetype=2, size=2) +
  scale_x_log10(breaks = gef$lr_min_reads) +scale_y_continuous(breaks=c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00)) +
  xlab("Minimum # of reads supporting a junction") + ylab("Recall") + ggtitle("Target Filtered") +theme_bw()

#unfiltered
ggplot(gef, aes(lr_min_reads)) + 
  geom_line(aes(y=geu$lr_recall, color="gtex_lr"), size=2) + geom_line(aes(y=geu$target_recall, color="gtex_target"), size=2) +
  geom_line(aes(y=seu$lr_recall, color="short_lr"), size=2) + geom_line(aes(y=seu$target_recall, color="short_target"), size=2) +
  geom_line(aes(y=aeu$lr_recall, color="sra_lr"), size=2) + geom_line(aes(y=aeu$target_recall, color="sra_target"), size=2) +
  geom_line(aes(y=peu$lr_recall, color="pacbio_lr"), size=2) + geom_line(aes(y=peu$target_recall, color="pacbio_target"), size=2) +
  geom_line(aes(y=reu$lr_recall, color="refseq_lr"), size=2) + geom_line(aes(y=reu$target_recall, color="refseq_target"), size=2) +
  scale_x_log10(breaks = gef$lr_min_reads) +scale_y_continuous(breaks=c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00)) +
  xlab("Minimum # of reads supporting a junction") + ylab("Recall") + ggtitle("Target Unfiltered") + theme_bw()



#both
ggplot(gef, aes(lr_min_reads)) + 
  geom_line(aes(y=lr_recall, color="gtex_lr"), linetype=2, size=2) + geom_line(aes(y=target_recall, color="gtex_target"), linetype=2, size=2) +
  geom_line(aes(y=geu$lr_recall, color="gtex_lr"), size=2) + geom_line(aes(y=geu$target_recall, color="gtex_target"), size=2) +
  geom_line(aes(y=sef$lr_recall, color="short_lr"), linetype=2, size=2) + geom_line(aes(y=sef$target_recall, color="short_target"), linetype=2, size=2) +
  geom_line(aes(y=seu$lr_recall, color="short_lr"), size=2) + geom_line(aes(y=seu$target_recall, color="short_target"), size=2) +
  geom_line(aes(y=aef$lr_recall, color="sra_lr"), linetype=2, size=2) + geom_line(aes(y=aef$target_recall, color="sra_target"), linetype=2, size=2) +
  geom_line(aes(y=aeu$lr_recall, color="sra_lr"), size=2) + geom_line(aes(y=aeu$target_recall, color="sra_target"), size=2) +
  geom_line(aes(y=pef$lr_recall, color="pacbio_lr"), linetype=2, size=2) + geom_line(aes(y=pef$target_recall, color="pacbio_target"), linetype=2, size=2) +
  geom_line(aes(y=peu$lr_recall, color="pacbio_lr"), size=2) + geom_line(aes(y=peu$target_recall, color="pacbio_target"), size=2) +
  geom_line(aes(y=reu$lr_recall, color="refseq_lr"), size=2) + geom_line(aes(y=reu$target_recall, color="refseq_target"), size=2) +
  scale_x_log10(breaks = gef$lr_min_reads) +scale_y_continuous(breaks=c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00)) +
  xlab("Minimum # of reads supporting a junction") + ylab("Recall") + theme_bw() + labs(color="Color/Linetype:\ndashed=target filtered\nsolid=target unfiltered")
 


#wiggle compare

gef <- data.frame(read.delim("D:/schatz_project/comparisons/gtex_wiggle_filtered_compare.tsv"))
geu <- data.frame(read.delim("D:/schatz_project/comparisons/gtex_wiggle_unfiltered_compare.tsv"))
sef <- data.frame(read.delim("D:/schatz_project/comparisons/short_wiggle_filtered_compare.tsv"))
seu <- data.frame(read.delim("D:/schatz_project/comparisons/short_wiggle_unfiltered_compare.tsv"))
reu <- data.frame(read.delim("D:/schatz_project/comparisons/refseq_wiggle_unfiltered_compare.tsv"))

#sra
aef <- data.frame(read.delim("D:/schatz_project/comparisons/sra_wiggle_filtered_compare.tsv"))
aeu <- data.frame(read.delim("D:/schatz_project/comparisons/sra_wiggle_unfiltered_compare.tsv"))

#pacbio
pef <- data.frame(read.delim("D:/schatz_project/comparisons/pacbio_wiggle_filtered_compare.tsv"))
peu <- data.frame(read.delim("D:/schatz_project/comparisons/pacbio_wiggle_unfiltered_compare.tsv"))



#filtered
ggplot(gef, aes(lr_min_reads)) + 
  geom_line(aes(y=lr_recall, color="gtex_lr"), linetype=2, size=2) + geom_line(aes(y=target_recall, color="gtex_target"), linetype=2, size=2) +
  geom_line(aes(y=sef$lr_recall, color="short_lr"), linetype=2, size=2) + geom_line(aes(y=sef$target_recall, color="short_target"), linetype=2, size=2) +
  geom_line(aes(y=aef$lr_recall, color="sra_lr"), linetype=2, size=2) + geom_line(aes(y=aef$target_recall, color="sra_target"), linetype=2, size=2) +
  geom_line(aes(y=pef$lr_recall, color="pacbio_lr"), linetype=2, size=2) + geom_line(aes(y=pef$target_recall, color="pacbio_target"), linetype=2, size=2) +
  geom_line(aes(y=pef$lr_recall, color="pacbio_lr"), linetype=2, size=2) + geom_line(aes(y=pef$target_recall, color="pacbio_target"), linetype=2, size=2) +
  scale_x_log10(breaks = gef$lr_min_reads) +scale_y_continuous(breaks=c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00)) +
  xlab("Minimum # of reads supporting a junction (wiggle of 9 bp)") + ylab("Recall") + ggtitle("Target Filtered") +theme_bw()

#unfiltered
ggplot(gef, aes(lr_min_reads)) + 
  geom_line(aes(y=geu$lr_recall, color="gtex_lr"), size=2) + geom_line(aes(y=geu$target_recall, color="gtex_target"), size=2) +
  geom_line(aes(y=seu$lr_recall, color="short_lr"), size=2) + geom_line(aes(y=seu$target_recall, color="short_target"), size=2) +
  geom_line(aes(y=aeu$lr_recall, color="sra_lr"), size=2) + geom_line(aes(y=aeu$target_recall, color="sra_target"), size=2) +
  geom_line(aes(y=peu$lr_recall, color="pacbio_lr"), size=2) + geom_line(aes(y=peu$target_recall, color="pacbio_target"), size=2) +
  geom_line(aes(y=reu$lr_recall, color="refseq_lr"), size=2) + geom_line(aes(y=reu$target_recall, color="refseq_target"), size=2) +
  scale_x_log10(breaks = gef$lr_min_reads) +scale_y_continuous(breaks=c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00)) +
  xlab("Minimum # of reads supporting a junction (wiggle of 9 bp)") + ylab("Recall") + ggtitle("Target Unfiltered") + theme_bw()


#both
ggplot(gef, aes(lr_min_reads)) + 
  geom_line(aes(y=lr_recall, color="gtex_lr"), linetype=2, size=2) + geom_line(aes(y=target_recall, color="gtex_target"), linetype=2, size=2) +
  geom_line(aes(y=geu$lr_recall, color="gtex_lr"), size=2) + geom_line(aes(y=geu$target_recall, color="gtex_target"), size=2) +
  geom_line(aes(y=sef$lr_recall, color="short_lr"), linetype=2, size=2) + geom_line(aes(y=sef$target_recall, color="short_target"), linetype=2, size=2) +
  geom_line(aes(y=seu$lr_recall, color="short_lr"), size=2) + geom_line(aes(y=seu$target_recall, color="short_target"), size=2) +
  geom_line(aes(y=aef$lr_recall, color="sra_lr"), linetype=2, size=2) + geom_line(aes(y=aef$target_recall, color="sra_target"), linetype=2, size=2) +
  geom_line(aes(y=aeu$lr_recall, color="sra_lr"), size=2) + geom_line(aes(y=aeu$target_recall, color="sra_target"), size=2) +
  geom_line(aes(y=pef$lr_recall, color="pacbio_lr"), linetype=2, size=2) + geom_line(aes(y=pef$target_recall, color="pacbio_target"), linetype=2, size=2) +
  geom_line(aes(y=peu$lr_recall, color="pacbio_lr"), size=2) + geom_line(aes(y=peu$target_recall, color="pacbio_target"), size=2) +
  geom_line(aes(y=reu$lr_recall, color="refseq_lr"), size=2) + geom_line(aes(y=reu$target_recall, color="refseq_target"), size=2) +
  scale_x_log10(breaks = gef$lr_min_reads) +scale_y_continuous(breaks=c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00)) +
  xlab("Minimum # of reads supporting a junction (wiggle of 9 bp)") + ylab("Recall") + theme_bw() + labs(color="Color/Linetype:\ndashed=target filtered\nsolid=target unfiltered")
