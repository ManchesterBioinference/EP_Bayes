library(matrixStats)
cwd <- getwd()

update_cwd = paste(cwd, "AP_clustering_output" ,sep = "/")
print(update_cwd)
setwd(update_cwd)

print(getwd())

data = read.table("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_PolII_2012-03_PolII_30_cor_0.2",sep=",")
labels_ = read.table("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_PolII_2012-03_PolII_30_cor_0.2_labels", head=TRUE, sep=",")[,2]
pdf("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_PolII_2012-03_PolII_30_cor_0.2.pdf")


LL <- list()
for(i in seq(1, max(labels_)))
{
LL<-c(LL,sum(labels_==i))
}
 
sorted_by_size<-order(unlist(LL), decreasing = TRUE)
 
j = 0
 
x <- seq(1,dim(data)[2]/2,1)
 
#pdf("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_unfiltered_count_concat_2012-03_PolII_PolII_100_0.2_clusters_.pdf")
 
 
 
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
 
#matplot(t(d), type="l",lty=1 ,lwd=1, col = 1, main = paste( j, ", " , unlist(LL)[i],sep = ""), xlab="",ylab="") # plots individual plots
 
matplot(m[1:8]-std[1:8], type="l", lty=1, lwd = 2, xlim=c(1, dim(data)[2]/2), ylim=c(old_min, old_max), col=rgb(0, 0, 1, 0.2), xlab="",ylab="",  axes=F) # main = unlist(LL)[i],
title(unlist(LL)[i], line = -1.75, cex.main = 3.0)
matlines(m[1:8]+std[1:8], type="l", lty=1, lwd = 2, col=rgb(0, 0, 1, 0.2), xlab="",ylab="")
matlines(m[1:8], type="l", lty=1, lwd = 1, col="deepskyblue")
polygon(c(x,rev(x)),c(m[1:8]+std[1:8], rev(m[1:8]-std[1:8])), col=rgb(0, 0, 1, 0.4),  border = NA)
 
matlines((m[9:16]-std[9:16]), type="l", lty=1, lwd = 2, col= rgb(1, 0, 0, 0.2))#, xlim=c(1, dim(data)[2]/2), ylim=c(min(m[9:16]-std[9:16]), max(m[9:16]std[9:16])), col="pink", xlab="",ylab="")
matlines(m[9:16]+std[9:16], type="l", lty=1, lwd = 2,  col= rgb(1, 0, 0, 0.2), xlab="",ylab="")
matlines(m[9:16], type="l", lty=1, lwd = 1,col="magenta")
 
polygon(c(x,rev(x)),c(m[9:16]+std[9:16], rev(m[9:16]-std[9:16])), col=rgb(1, 0, 0, 0.4),  border = NA)
 
if (j %% 4 == 0) {axis(2, at = c(-2, 0, 2), cex.axis = 2.5)}
#if (j %%7 %in% (1:4)) {axis(1, at = 1:8, labels = c("0", "5", "10", "20", "40", "80", "160", "320"), cex.axis = 1)}
if (j %in% (c(8:11, 9:12 + 11, 9:12+23))) {axis(1, at = 1:8, labels = c("0", "5", "10", "20", "40", "80", "160", "320"), cex.axis = 2.5, line = -1.55)}#
#box(lty=1, col="black")
 
j<-j+1
}
 
dev.off()

