#in the paper
#----------------------------------------
library(matrixStats)
cwd <- getwd()

update_cwd = paste(cwd, "AP_clustering_output" ,sep = "/")
print(update_cwd)
setwd(update_cwd)

print(getwd())

data = read.table("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_ER_100",sep=",")
labels_ = read.table("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_ER_100_labels", head=TRUE, sep=",")[,2]
pdf("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_ER_100.pdf")
#----------------------------------------


LL <- list()
for(i in seq(1, max(labels_)))
{
LL<-c(LL,sum(labels_==i))
}

sorted_by_size<-order(unlist(LL), decreasing = TRUE)

j = 0


x <- seq(1,dim(data)[2],1)


old_max <- 0
old_min <- 0
for(i in sorted_by_size) {
d<-data.matrix(data[labels_==i,], rownames.force = NA)
m<-colMeans(d, na.rm=TRUE)
std<-rowSds(t(d))
if (min(m - std) < old_min) {old_min = min(m - std)}
if (max(m + std) > old_max) {old_max = max(m + std)}
}



iter <- 0

for(i in sorted_by_size) {
if (j %% 12 == 0) {par(mfrow=c(3,4), mar = c(1.5,2.8,0.3,1.3), oma=c(0,0,0,0))}


d<-data.matrix(data[labels_==i,], rownames.force = NA)
m<-colMeans(d, na.rm=TRUE)
std<-rowSds(t(d))

#matplot(t(d), type="l",lty=1 ,lwd=1, col = 1, main = unlist(LL)[i], xlab="",ylab="") # plots individual plots

matplot(m-std, type="l", lty=1, lwd = 2, xlim=c(1, dim(data)[2]), ylim=c(old_min, old_max), col=rgb(0, 0, 1, 0.2), xlab="",ylab="", axes=F)
title(unlist(LL)[i], line = -1.75, cex.main = 3.0)
matlines(m+std, type="l", lty=1, lwd = 2, col=rgb(0, 0, 1, 0.2), xlab="",ylab="")
matlines(m, type="l", lty=1, lwd = 1, col="black")
polygon(c(x,rev(x)),c(m+std, rev(m-std)), col=rgb(0, 0, 1, 0.4),  border = NA)


if (j %% 4 == 0) {axis(2, at = c(-2, 0, 2), cex.axis = 2.5)}
#if (j %%7 %in% (1:4)) {axis(1, at = 1:8, labels = c("0", "5", "10", "20", "40", "80", "160", "320"), cex.axis = 1)}
if (j %in% (c(8:11, 9:12 + 11, 9:12+23))) {axis(1, at = 1:8, labels = c("0", "5", "10", "20", "40", "80", "160", "320"), cex.axis = 2.5, line = -1.55)}#
#box(lty=1, col="black")
 
j<-j+1



}
dev.off()

proc.time()


