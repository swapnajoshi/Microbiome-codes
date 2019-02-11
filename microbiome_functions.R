# Functions for microbiome analysis

################################################################################
# Acess element of vector separated by semicolon
################################################################################
extract.name.level = function(x, level){
  a=c(unlist(strsplit(x,';')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}
################################################################################
# get taxonomy for the OTU number
################################################################################
otu2taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa))
  {taxa = colnames(x)}
  if(length(taxa)!=dim(x)[2])
  {print("ERROR: taxonomy should have the same length
         as the number of columns in OTU table")
    return;}
  level.names = sapply(as.character(taxa),
                       function(x)
                         extract.name.level(x,level=level))
  t(apply(x, 1, function(y) tapply(y,level.names,sum)))
}
################################################################################
# make ggplots scatter plot from 2 vectors and a color variable in a data frame 
################################################################################
plotAllScatter <- function(z) {
  P<-ggplot(z, aes(z[,x],z[,y])) + geom_point(size=2.5) 
  q <- P + geom_point(aes(color=as.factor(z$Group))) + xlab(colnames(z)[x]) + 
  ylab(colnames(z)[y])
  ggsave(q, filename = paste(colnames(z)[x],colnames(z)[y],".png"))
}
# usage
# for(x in 1:2) for(y in 3:4) {
#   plotAllScatter(data3)
# }
################################################################################
# make ggplots box plot from 2 vectors, (one numeric, one factor) in a dataframe 
################################################################################
plotAllBox <- function(z) {
  P <- ggplot(z, aes(as.factor(z[,x]), z[,y], fill=as.factor(z[,x]))) 
  q <- P + geom_boxplot() + xlab(colnames(z)[x]) +  ylab(colnames(z)[y])
  ggsave(q, filename = paste(colnames(z)[x],colnames(z)[y],".png"))
}
# usage
# for(y in 3:4) {
#   x = 'Group'
#   plotAllBox(data3)
# }
################################################################################
# Scatterplots with regression line with geom_smooth function
################################################################################
plotAllScatterRegr <- function(z) {
  P<-ggplot(z, aes(z[,x],z[,y])) + geom_point(size=2.5) 
  q <- P + geom_point(aes()) + xlab(colnames(z)[x]) + 
  ylab(colnames(z)[y]) + geom_smooth(aes(z[,x],z[,y]), method=lm, se=FALSE)
  ggsave(q, filename = paste(colnames(z)[x],colnames(z)[y],".png"))
}
# usage
# for(x in 1:2) for(y in 3:4) {
#   plotAllScatterRegr(data3)
# }
################################################################################
# Function for correlation between two dataframes
################################################################################
cor2dfs <- function(x,y) {
  corList <- rcorr(as.matrix(x), type="spearman")
  corVal  <- corList$r
  adjP <- apply(corList$P,2,p.adjust)
  highCor <- ifelse(adjP <= 0.1, 1,0)
  corInt <- highCor[1:y, y+1:dim(x)[2]]
  corInt$sum1<-apply(corInt,1,sum)
  corInt1<-corInt[corInt$sum1>0,]
  corInt2<-apply(corInt1,2,sum)
  df2Int<-names(corInt2[corInt2>0])
  SigCor <- corInt2[,df2Int]
  return(SigCor)
}
# Usage
# df1Sel <- t(microbMat2.ibs[row.names(microbMat2.ibs)%in%
#                                  row.names(cor2dfs),]); dim(cor.microb)
# df2Sel <- t(mirna1.ibs[row.names(mirna1.ibs) %in% colnames(cor2dfs),]); 
# dim(cor.mirna)
# 
# data2<-as.data.frame(cbind(df1Sel,df2Sel));dim(data2)

