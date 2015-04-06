# Making barplots of relative abundance of top genera (per sample)

### Anna Seekatz
### 12.17.14

Change into your working directory: 
```{r}
setwd("~/Desktop/umich/projects/Relapse_model2/analysis/mothur_analysis")
list.files(path="~/Desktop/umich/projects/Relapse_model2/analysis/mothur_analysis")
```

Please note: this code is not completely reproducible in R!
	+ As indicated in the comments, I also used excel to organize and sort the data to my liking
	+ I have tried to include examples of the output of my excel modifications
	+ The methods I used can probably be adapted to R, but would require a defined taxonomic file to your liking (since I organize by phyla, and handpicked my colors)

The files I used for this analysis:
	+ mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary (this will eventually be sorted)
	+ mouseRM_genbar.txt (modified as described--this is just to show an example of what I filtered)
	+ mouseRM_genfrac2p.txt (further filtered to only include the top 98% of genera; also converted values into relative abundance)
	+ mouseRM_phylofrac2p_ordered.txt (same file, only ordered by phylum and relative abundance within each phylum)
	
	
To create the files, I used the following methods:

```{r}
#all data in directory: /Users/aseekatz/Desktop/umich/projects/KPC_rush/rush_analysis/rush2
#took files:
	mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
	#saved as: mouseRM_phylotype.xlsx
	#organized sheet to include genus only
	#deleted genera with <100 sequences
	#added unclassified info
	#added phyla
	#reordered samples by day, treatment, cage, mouse
	#saved as: mouseRM_genbar.txt
```

Then in R:

```{r}

#in R:

#For genus barplot:
#wd: /Users/aseekatz/Desktop/umich/projects/Relapse_model2/analysis/mothur_analysis
genus<-read.table(file="mouseRM_genbar.txt", header=TRUE, row.names=1)
genus<-as.data.frame(t(genus))
genus<-genus/rowSums(genus)*100
genus <- genus[, order(-unlist(lapply(genus, median)))]
genus<-t(genus)
write.table(genus, file="mouseRM_genfrac.txt", quote=FALSE, sep="\t")
genus1<- genus[rowSums(genus>=1)>=1,] #will list all genera above 1%
genus2<- genus[rowSums(genus>=2)>=2,] #will list all genera above 2%
#nrow(genus2) #this will give you the number of taxa -->using 1p, since only 52 taxa (82 in 1p)
write.table(genus2, 'mouseRM_genfrac2p.txt',quote=FALSE,sep="\t")

#took 2% file, and sorted by phyla (using VLOOKUP), then abundance; added colors
	#=VLOOKUP(A2,genus_only!C:D,2,0)
#to collapse some of the smaller percentages, within each phyla, any phylotypes less than 0.2% average abundance were combined into 'other_phylaname'
#saved as mouseRM_phylofrac2p_ordered.txt

##For color palette:
#I picked these colors to distinguish between different phylum-level classification
col.rm_g<-c("darkgreen","green1","greenyellow","green3","palegreen2","limegreen","midnightblue","blue3","blue","dodgerblue4","deepskyblue3","dodgerblue2","cornflowerblue","royalblue1","deepskyblue","cadetblue3","lightskyblue1","lightskyblue3","cyan","darkcyan","orchid4","purple4","yellow2","gold","orchid1","tomato4","hotpink","maroon4","red","grey47")

#Plotting:
rm_g<-read.table(file="mouseRM_phylofrac2p_ordered.txt", header=TRUE, row.names=1)
rm_g<-as.matrix(rm_g)
par(mar=c(5,4,2,5))
par(xpd=T)
barplot(rm_g, las=2, main="Mouse Relapse Model (>2% genera)", ylab="Relative abundance-genera (%)", cex.names=0.4, ylim=c(0,100), col=col.rm_g)
legend(350,100,legend=rownames(rm_g),col=col.rm_g,fill=col.rm_g,cex=0.8)

#if you only wanted to make a bargraph of pre-abx mice, you could filter all of the 'dn7' samples:
filtered.g <- rm_g[,grepl("dn7",colnames(rm_g))]
test<-filtered.g[which(rowSums(filtered.g)>0),]
filtered.g<-as.matrix(filtered.g)

# you can order the samples as you please
# the graph was further edited to fit the page in Adobe illustrator

```
