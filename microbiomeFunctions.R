clust1 = cutreeStatic(sampleTree, cutHeight = 4, minSize = 1)
table(clust1)
lab<-map2[sampleTree$labels,]$Dx
lab1<-map2[sampleTree$labels,]$NDPNum
df2<-as.data.frame(cbind(clust1,as.character(lab)))
head(df2)
table(df2$clust1, df2$V2) #   42% vs 30% hcs in cluster 1 vs cluster 2

# set.seed(123)
# rainbow(n=6)
# [1] "#FF0000FF" "#FFFF00FF" "#00FF00FF" "#00FFFFFF" "#0000FFFF" "#FF00FFFF"

df2$color1<-NA
for ( i in 1: 85)
{if (df2[i,1]=="1")
{df2[i,3]<-"red"}
  else if (df2[i,1]=="2")
  {df2[i,3]<-"cyan"}
  else if (df2[i,1]=="3")
  {df2[i,3]<-"orange"}
}

df2$colorDx<-NA
for ( i in 1: 85)
{if (df2[i,2]=="IBS")
{df2[i,4]<-"cyan"}
  else if (df2[i,2]=="HC")
  {df2[i,4]<-"red"}
}
df2<-cbind(df2, as.character(lab1))

library(dendextend)
png("BrayCurtiswardHierarchical.png", res = 200, height = 2000, width = 2500)
plot(sampleTree, main = "Heirachical Clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, labels=labels1)
abline(h = 4, col = "red")
colored_bars(colors= df2[,c(3,4)], dend = sampleTree, y_shift = -1.0, rowLabels = c("Cluster", "Dx"),sort_by_labels_order 	
             = TRUE)
dev.off()

#plot
plot(sampleTree, main = "Heirachical Clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, labels=labels1)
abline(h = 4, col = "red")
colored_bars(colors= df2[,c(3,4)], dend = sampleTree, y_shift = -1.0, rowLabels = c("Cluster", "Dx"),sort_by_labels_order 	
             = TRUE)

# cluster memberships
map2$IBS_HC_bc_ward_clu_memb <- df2[,1]