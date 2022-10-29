data_file=read.table("GSE156922_read_counts.txt",row.names  = 1,header = TRUE )
head(data_file)

#Create a count per matrix
cpmatrix=data_file
for(i in 1:ncol(data_file)){
  cpmatrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
}

#Calculate a log of cpm
logcpm=log2(cpmatrix+1)
saveRDS(logcpm,file="logCPM.rds")
summary(logcpm)

#Calculate a z score
library(matrixStats)
z_score = (logcpm - rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]
z_score

#Calculate variance using log 
variance = apply(logcpm, 1, var)
variance = sort(variance,decreasing = T)
top50 = variance[1:50]
pmat = z_score[names(top50),]

#Create a heatmap
library(ComplexHeatmap)
Heatmap(pmat)
