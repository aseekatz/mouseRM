# Heatmap (with clustering) created in R

### Anna Seekatz
### 12.20.14

> load all of the necessary packages
Load your packages:
```{r, error=FALSE, message=FALSE, warning=FALSE}
library(vegan)
library(gplots)
library(RColorBrewer)
library(ggplot2)

```

Change into your working directory: 
```{r}
setwd("~/Desktop/umich/projects/Relapse_model2/analysis/mothur_analysis")
list.files(path="~/Desktop/umich/projects/Relapse_model2/analysis/mothur_analysis")

```

Files used:
	+ mouseRM_otucounts.txt (this is modified to be a matrix from the mothur-produced file, mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared)
	+ mouse.var.txt

To make heatmap:

```{r}

# read in your otutable 
otu<-read.table(file="mouseRM_otucounts.txt", header=TRUE, row.names=1)
otuf<-otu/rowSums(otu)*100												#this calculates relative abundance
maxab <- apply(otuf, 2, max)											
n2 <- names(which(maxab < 2.0))											#eliminated OTUs present at less than 2%
otu2 <- otuf[, -which(names(otuf) %in% n2)]								
otuf.dist <- vegdist(otu, method = "morisita")							#calculates a distance matrix (you can customize to be whatever distance you want)	
row.clus <- hclust(otuf.dist, "aver")									#clusters based on the pairwise distance matrix (used in your heatmap to cluster)

# now, you can add some metadata as a column against your heatmap--you just need a supplemental file for that
#make sure that your variable file has the samples listed in the same order as your initial otuf/otu2! otherwise, you will not get matched clustering
mouse.var<-read.table(file="mouse.var.txt", header=TRUE)
mouse.var<- mouse.var[order(mouse.var$sampleID),]
var1<-mouse.var[,c("sampleID", "group_time_col")]						# I had previously added these colors into the metadata, but you could also substitute your variables with colors if necessary
cbind(row.names(otu2), var1)

#can use group instead:
var2<-mouse.var[,c("sampleID", "group_col")]
cbind(row.names(otu2), var2)

# define your color breaks according to the relative abundance you want shown per color group (always have 1 more break than color)
col_breaks = c(seq(-1,0,length=1),
  seq(0,1,length=1),
  seq(1,5,length=1),
  seq(5,10,length=1),
  seq(10,50,length=1),
  seq(50,80,length=1),
  seq(80,100,length=1)
  )
  
#define your color groups:
my_palette <- colorRampPalette(c("azure", "skyblue1", "steelblue1", "dodgerblue", "royalblue3", "blue", "darkblue"))(n = 6)

# then, plot the 
heatmap.2(as.matrix(otu2), Rowv = as.dendrogram(row.clus), dendrogram="row", Colv=NA, col = my_palette, breaks=col_breaks, margins = c(10, 26), RowSideColors = as.character(var1$group_time_col), lhei = c(2, 6), density.info = "none", trace="none", symbreaks=FALSE, symkey=FALSE)
legend("bottomright",legend=c("pre-treatment", "ABX only: during cef", "ABX only: post-cef", "ABX only: during vanco", "ABX only: post-vanco", "CDI: pre-vanco", "CDI: post-vanco treatment", "CDI: post-IP clinda", "CDI: post-FMT", "CDI: post-FMT, post-clinda", "donor FMT"), col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"), cex=0.8, pch=15)
legend("topright",legend=c("not detected", "<1%", "1-5%", "5-10%", "10-50%", "50-80%", "80-100%"), col=c("azure", "skyblue1", "steelblue1", "dodgerblue", "royalblue3", "blue", "darkblue"), cex=0.8, pch=15)

# The figure was edited in Adobe Illustrator (margins in R are annoying...)

That's it!

```