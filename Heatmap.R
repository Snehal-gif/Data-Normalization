read.csv("C:/Users/Admin/OneDrive/Documents/cancer genomics/GSE183913_Karanjin.csv")
cancer<-matrix(rnorm(100,0,5),nrow = 10,ncol = 10)

colnames(cancer)<-paste0("gene",1:10)
rownames(cancer)<-paste0("species",1:10)

heatmap(cancer)
heatmap(cancer,xlab = "species",ylab="gene")

cpm = apply(cancer, 2, function(x) (x/sum(x))*1000000)
print(cpm)

vargenes = apply(x_stand1,1,var)
print(vargenes)

vargenes = sort(vargenes,decreasing = T)
top50 = vargenes[1:50]

mat = matrix(NA, ncol= 4, nrow= nrow(cancer))
rownames(mat)= rownames(cancer)
colnames(mat)= c('grp1','grp2','pval','log2FC')
print(mat)

for(i in 1:nrow(cancer)){
  vec1 = as.numeric(cancer[i,1:4])
  vec2 = as.numeric(cancer[i,5:7])
  res = t.test(vec1, vec2, paired = F, alternative = 'two.sided')
  mat[i,1]= res$estimate[[1]]
  mat[i,2]= res$estimate[[2]]
  mat[i,3]= res$p.value
  mat[1,4]=mat[i,1]-mat[i,2]
  
}

mat= as.data.frame(mat)
num = which(is.nan(mat$pval))
mat[num,'pval']=1

library("EnhancedVolcano")

EnhancedVolcano(mat,lab = rownames(mat),x='log2FC',y= 'pval')
