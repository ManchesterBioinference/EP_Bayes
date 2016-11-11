#------------------------------------------------------------------------------------------------------------------------------------------
#--------------
#in the paper
#Joined Pol2/ER

library(matrixStats)
cwd <- getwd()

update_cwd = paste(cwd, "AP_clustering_output" ,sep = "/")
print(update_cwd)
setwd(update_cwd)

print(getwd())

data = read.table("common_region_peaks_extended_less_time_points_corrected_0_indexed_unfiltered_count_concat_PolII_ER_200",sep=",")
labels_ = read.table("common_region_peaks_extended_less_time_points_corrected_0_indexed_unfiltered_count_concat_PolII_ER_200_labels", head=TRUE, sep=",")[,2]
pdf("common_region_peaks_extended_less_time_points_corrected_0_indexed_unfiltered_count_concat_PolII_ER_200_clusters.pdf")

LL <- list()
for(i in seq(1, max(labels_)))
{
LL<-c(LL,sum(labels_==i))
}

sorted_by_size<-order(unlist(LL), decreasing = TRUE)

j = 0


x <- c(seq(1,dim(data)[2]/2,1), seq(1,dim(data)[2]/2,1))

x2 <- c(seq(1,dim(data)[2],1))

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
if (j %% 9 == 0) {par(mfrow=c(3,3), mar = c(2.,2.5,0.0,1.25), oma=c(0,0,0,0))}


d<-data.matrix(data[labels_==i,], rownames.force = NA)
m<-colMeans(d, na.rm=TRUE)
std<-rowSds(t(d))

#matplot(t(d), type="l",lty=1 ,lwd=1, col = 1, main = unlist(LL)[i], xlab="",ylab="") # plots individual plots


matplot((m-std)[1:8], type="l", lty=1, lwd = 2, xlim=c(1, dim(data)[2]), ylim=c(old_min, old_max), col=rgb(255,0,255, 50, maxColorValue=255), xlab="",ylab="", axes = FALSE)
title(unlist(LL)[i], line = -1.75)
matlines((m+std)[1:8], type="l", lty=1, lwd = 2, col=rgb(255,0,255, 50, maxColorValue=255), xlab="",ylab="")

matlines(x2[9:16] , (m-std)[9:16], type="l", lty=1, lwd = 2, xlim=c(1, dim(data)[2]), ylim=c(min(m-std), max(m+std)), col=rgb(0, 0, 1, 0.2), main = paste(unlist(LL)[i]), xlab="",ylab="")
matlines(x2[9:16], (m+std)[9:16], type="l", lty=1, lwd = 2, col=rgb(0, 0, 1, 0.2), xlab="",ylab="")


matlines(m[1:8], type="l", lty=1, lwd = 1, col="black")
matlines(9:16, m[9:16], type="l", lty=1, lwd = 1, col="black")
#polygon(c(x2[1:8],rev(x2[1:8])), c(m+std, rev(m-std)), col=rgb(0, 0, 1, 0.2),  border = NA, axes = FALSE)

polygon(c(x2[9:16],rev(x2[9:16])), c((m+std)[9:16], rev((m-std)[9:16])), col=rgb(0, 0, 1, 0.2),  border = NA)
polygon(c(x2[1:8],rev(x2[1:8])), c((m+std)[1:8], rev((m-std)[1:8])), col=rgb(255,0,255, 50, maxColorValue=255),  border = NA)

#abline(v=8.5,col=3, lty=3)
lines(x=c(8.5,8.5), y=c(-2,2), lty=3, col = c("#7570B3"))
#box()

if (j %% 3 == 0) {axis(2, at = c(-2, 0, 2), cex.axis = 1.8)} 
#if (j %%7 %in% (1:4)) {axis(1, at = 1:8, labels = c("0", "5", "10", "20", "40", "80", "160", "320"), cex.axis = 1)}
if (j %in% (c(6:8, 15:17, 24:26, 33:35))) {axis(1, at = 1:16, labels = c("0", "5", "10", "20", "40", "80", "160", "320", "0", "5", "10", "20", "40", "80", "160", "320"), cex.axis = 2.25, line = -1.25)#


j<-j+1

}
dev.off()

proc.time()



