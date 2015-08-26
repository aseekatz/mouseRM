# Calculating median (or mean) relative abundance of top genera within each treatment group
### Anna Seekatz
### 2.11.15

> load all of the necessary packages
Load your packages:
```{r, error=FALSE, message=FALSE, warning=FALSE}
library(ddply)
```

Change into your working directory: 
```{r}
setwd("/Users/aseekatz/Desktop/umich/projects/Relapse_model2/mouseRM_analysis/mothur_analysis")
list.files(path="/Users/aseekatz/Desktop/umich/projects/Relapse_model2/mouseRM_analysis/mothur_analysis")
```

> These are the files I used:
	+ mouseRM_phylofrac2p_for.barplot.txt (this is a modified version of all of the genera above 2% of the full dataset, created when making previous barplots of each sample)
	+ mouse.var.txt

```{r}
##Making bargraphs comparing FMT and no FMT:

# these are the explanations for the metadata in mouse.var.txt (used in this graph):
	--bargroup:
		-pre=donor and pre-abx samples
		-cdi=after cef, before vanco (days 0-4 for all mice)
		-post'n'=post-vanco for FMT or noFMT mice (days 12-35)
		-postclinda'n' for either FMT or noFMT mice (days 36-42)
	--barfmt: FMT, noFMT, or 'other' group

# added some extra rows to mouseRM_phylofrac2p_ordered.txt (saved as mouseRM_phylofrac2p_for.barplot.txt):
	# added 'other_Bacteroidetes (sum of all Bact. other than uncl_porph, bacteroides, and barnesiella)
	# added 'other_Firmicutes' (sum of all except Lacto, Turicibacter, uncl_lachno, Clostridium_XI)
	
# merged previous 'mouseRM_phylofrac2p_ordered.txt' and mouse.var.txt:
mouse.var<-read.table(file="mouse.var.txt", header=TRUE)
rm_g<-read.table(file="mouseRM_phylofrac2p_for.barplot.txt", header=TRUE, row.names=1)
barg<-as.data.frame(t(rm_g))
barg<- barg[, order(-unlist(lapply(barg, mean)))]
barg$sampleID<-rownames(barg)
bar<-merge(mouse.var, barg, by.x=c("sampleID"), by.y=c("sampleID"))
#write.table(bar, file="mouseRM_bar.graphs.txt", quote=FALSE, sep="\t", col.names=NA)

# first, get the median and min/max for your groups:
#colnames(bar)			#this will show you the correct columns to choose (I am choosing the top 9 since I ordered by mean earlier)
median<- sapply(levels(bar$bargroup), function(class) sapply(bar[bar$bargroup == class, c("unclassified_Porphyromonadaceae", "Bacteroides", "Barnesiella", "other_bacteroidetes", "Lactobacillus", "Turicibacter", "unclassified_Lachnospiraceae", "Clostridium_XI", "other_Firmicutes", "unclassified_Enterobacteriaceae", "Akkermansia")], median))
min<- sapply(levels(bar$bargroup), function(class) sapply(bar[bar$bargroup == class, c("unclassified_Porphyromonadaceae", "Bacteroides", "Barnesiella", "other_bacteroidetes", "Lactobacillus", "Turicibacter", "unclassified_Lachnospiraceae", "Clostridium_XI", "other_Firmicutes", "unclassified_Enterobacteriaceae", "Akkermansia")], function(x) quantile(x)[2]))
max<- sapply(levels(bar$bargroup), function(class) sapply(bar[bar$bargroup == class, c("unclassified_Porphyromonadaceae", "Bacteroides", "Barnesiella", "other_bacteroidetes", "Lactobacillus", "Turicibacter", "unclassified_Lachnospiraceae", "Clostridium_XI", "other_Firmicutes", "unclassified_Enterobacteriaceae", "Akkermansia")], function(x) quantile(x)[4]))

# then select and subset the columns you want:
	#first row: pre, cdi, post
median1<-median[,c("pre", "cdi", "post10")]  #this picks the specific columns
min1<-min[,c("pre", "cdi", "post10")]
max1<-max[,c("pre", "cdi", "post10")]
	#second row: post27_noFMT, postclinda36_noFMT, postclinda42_noFMT
median2<-median[,c("post27_noFMT", "post_clinda36_noFMT", "post_clinda42_noFMT")]
min2<-min[,c("post27_noFMT", "post_clinda36_noFMT", "post_clinda42_noFMT")]
max2<-max[,c("post27_noFMT", "post_clinda36_noFMT", "post_clinda42_noFMT")]
	#third row: post27_FMT, post_clinda36_FMT, post_clinda42_FMT
median3<-median[,c("post27_FMT", "post_clinda36_FMT", "post_clinda42_FMT")]
min3<-min[,c("post27_FMT", "post_clinda36_FMT", "post_clinda42_FMT")]
max3<-max[,c("post27_FMT", "post_clinda36_FMT", "post_clinda42_FMT")]
	
# then plot:
bar.col=c("darkgreen", "green1", "greenyellow", "limegreen", "midnightblue", "blue3", "blue", "dodgerblue4", "cyan", "yellow2", "hotpink")
par(mfrow=c(3,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
par(xpd=TRUE)
mp1<-barplot(median1, beside=TRUE,col=bar.col,las=1, ylim=c(0,100), ylab="Relative Abundance")
segments(mp1, max1, mp1, min1)
legend(36, 100, c("Porphyromonadaceae fam", "Bacteroides", "Barnesiella", "other Bacteroidetes", "Lactobacillus", "Turicibacter", "Lachnospiraceae fam.", "Clostridium XI", "other Firmicutes", "Enterobacteriaceae fam.", "Akkermansia"), col=bar.col, pch=15)
mp2<-barplot(median2, beside=TRUE,col=bar.col,las=1, ylim=c(0,100), ylab="Relative Abundance")
segments(mp2, max2, mp2, min2)
mp3<-barplot(median3, beside=TRUE,col=bar.col,las=1, ylim=c(0,100), xlab="Day pi", ylab="Relative Abundance")
segments(mp3, max3, mp3, min3)

```

For stats on the bargraphs:

```{r}

# Stats for bargraphs:
#read in merged table that was created previously:
barplots<-read.table(file="mouseRM_bar.graphs.txt", header=TRUE)

levels(barplots$bargroup)				#gives the classes I had defined before

# wilcoxon test:

#pre-abx compared to post-cef:

#subset to include only pre-abx and d1-4 post-cef:
pre.cef<-subset(barplots, (subset=bargroup %in% c("pre", "cdi")))
wilcox.test(unclassified_Porphyromonadaceae~bargroup, data=pre.cef)
	W = 0, p-value = 6.026e-07
# could do this individually for each phylotype...however, can create a file (although nonvector) for each of the desired columns (in this case, columns 14:24)
pre.v.cef<-apply(pre.cef[14:24],2,function(x) wilcox.test(x~pre.cef$bargroup))
#if you want this as a list, you can use the following:


## could also do pairwise wilcoxon test, since there are three groups within figure A that we want to look at:
#subset to include pre-abx, post-cef, and post-vanco:
pre.cef.vanco<-subset(barplots, (subset=bargroup %in% c("pre", "cdi", "post10")))
pre.v.cef.v.vanco<-apply(pre.cef.vanco[14:24],2,function(x) pairwise.wilcox.test(x, pre.cef.vanco$bargroup, exact=FALSE, p.adj="bonferroni"))
	#this is better!
	
# this command will present your p-values in a nicer way (as a list):
porph<-melt(pre.v.cef.v.vanco$unclassified_Porphyromonadaceae[[3]])
porph$phylotype<- c("uncl_Porphyromonadaceae")
bact<-melt(pre.v.cef.v.vanco$Bacteroides[[3]])
bact$phylotype<- c("Bacteroides")
barn<-melt(pre.v.cef.v.vanco$Barnesiella[[3]])
barn$phylotype<- c("Barnesiella")
lacto<-melt(pre.v.cef.v.vanco$Lactobacillus[[3]])
lacto$phylotype<- c("Lactobacillus")
turi<-melt(pre.v.cef.v.vanco$Turicibacter[[3]])
turi$phylotype<- c("Turicibacter")
lachno<-melt(pre.v.cef.v.vanco$unclassified_Lachnospiraceae[[3]])
lachno$phylotype<- c("uncl_Lachnospiraceae")
cdiff<-melt(pre.v.cef.v.vanco$Clostridium_XI[[3]])
cdiff$phylotype<- c("Clostridium_XI")
akk<-melt(pre.v.cef.v.vanco$Akkermansia[[3]])
akk$phylotype<- c("Akkermansia")
entero<-melt(pre.v.cef.v.vanco$unclassified_Enterobacteriaceae[[3]])
entero$phylotype<- c("uncl_Enterobacteriaceae")
other_b<-melt(pre.v.cef.v.vanco$other_bacteroidetes[[3]])
other_b$phylotype<- c("other_bacteroidetes")
other_f<-melt(pre.v.cef.v.vanco$other_Firmicutes[[3]])
other_f$phylotype<- c("other_Firmicutes")

figA_stats<-rbind(porph,bact,barn,lacto,turi,lachno,cdiff,akk,entero,other_b,other_f)
figA_stats<-na.omit(figA_stats)
figA_stats<-figA_stats[order(figA_stats$Var1, figA_stats$Var2, -figA_stats$value) ,]
figA_stats_filtered<-subset(figA_stats, value<=0.05)
write.table(figA_stats, file="barplots_stats/figA_stats.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(figA_stats_filtered, file="barplots_stats/figA_stats_filtered.txt", quote=FALSE, sep="\t", col.names=NA)

# Fig. B stats over time:
nofmt.bar<-subset(barplots, (subset=bargroup %in% c("post27_noFMT", "post_clinda36_noFMT", "post_clinda42_noFMT")))
nofmt.v.time<-apply(nofmt.bar[14:24],2,function(x) pairwise.wilcox.test(x, nofmt.bar$bargroup, exact=FALSE, p.adj="bonferroni"))
porph<-melt(nofmt.v.time$unclassified_Porphyromonadaceae[[3]])
porph$phylotype<- c("uncl_Porphyromonadaceae")
bact<-melt(nofmt.v.time$Bacteroides[[3]])
bact$phylotype<- c("Bacteroides")
barn<-melt(nofmt.v.time$Barnesiella[[3]])
barn$phylotype<- c("Barnesiella")
lacto<-melt(nofmt.v.time$Lactobacillus[[3]])
lacto$phylotype<- c("Lactobacillus")
turi<-melt(nofmt.v.time$Turicibacter[[3]])
turi$phylotype<- c("Turicibacter")
lachno<-melt(nofmt.v.time$unclassified_Lachnospiraceae[[3]])
lachno$phylotype<- c("uncl_Lachnospiraceae")
cdiff<-melt(nofmt.v.time$Clostridium_XI[[3]])
cdiff$phylotype<- c("Clostridium_XI")
akk<-melt(nofmt.v.time$Akkermansia[[3]])
akk$phylotype<- c("Akkermansia")
entero<-melt(nofmt.v.time$unclassified_Enterobacteriaceae[[3]])
entero$phylotype<- c("uncl_Enterobacteriaceae")
other_b<-melt(nofmt.v.time$other_bacteroidetes[[3]])
other_b$phylotype<- c("other_bacteroidetes")
other_f<-melt(nofmt.v.time$other_Firmicutes[[3]])
other_f$phylotype<- c("other_Firmicutes")

figB_stats<-rbind(porph,bact,barn,lacto,turi,lachno,cdiff,akk,entero,other_b,other_f)
figB_stats<-na.omit(figB_stats)
figB_stats<-figB_stats[order(figB_stats$Var1, figB_stats$Var2, -figB_stats$value) ,]
figB_stats_filtered<-subset(figB_stats, value<=0.05)
write.table(figB_stats, file="barplots_stats/figB_stats.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(figB_stats_filtered, file="barplots_stats/figB_stats_filtered.txt", quote=FALSE, sep="\t", col.names=NA)

# Fig. C stats over time:
fmt.bar<-subset(barplots, (subset=bargroup %in% c("post27_FMT", "post_clinda36_FMT", "post_clinda42_FMT")))
fmt.v.time<-apply(fmt.bar[14:24],2,function(x) pairwise.wilcox.test(x, fmt.bar$bargroup, exact=FALSE, p.adj="bonferroni"))
porph<-melt(fmt.v.time$unclassified_Porphyromonadaceae[[3]])
porph$phylotype<- c("uncl_Porphyromonadaceae")
bact<-melt(fmt.v.time$Bacteroides[[3]])
bact$phylotype<- c("Bacteroides")
barn<-melt(fmt.v.time$Barnesiella[[3]])
barn$phylotype<- c("Barnesiella")
lacto<-melt(fmt.v.time$Lactobacillus[[3]])
lacto$phylotype<- c("Lactobacillus")
turi<-melt(fmt.v.time$Turicibacter[[3]])
turi$phylotype<- c("Turicibacter")
lachno<-melt(fmt.v.time$unclassified_Lachnospiraceae[[3]])
lachno$phylotype<- c("uncl_Lachnospiraceae")
cdiff<-melt(fmt.v.time$Clostridium_XI[[3]])
cdiff$phylotype<- c("Clostridium_XI")
akk<-melt(fmt.v.time$Akkermansia[[3]])
akk$phylotype<- c("Akkermansia")
entero<-melt(fmt.v.time$unclassified_Enterobacteriaceae[[3]])
entero$phylotype<- c("uncl_Enterobacteriaceae")
other_b<-melt(fmt.v.time$other_bacteroidetes[[3]])
other_b$phylotype<- c("other_bacteroidetes")
other_f<-melt(fmt.v.time$other_Firmicutes[[3]])
other_f$phylotype<- c("other_Firmicutes")

figC_stats<-rbind(porph,bact,barn,lacto,turi,lachno,cdiff,akk,entero,other_b,other_f)
figC_stats<-na.omit(figC_stats)
figC_stats<-figC_stats[order(figC_stats$Var1, figC_stats$Var2, -figC_stats$value) ,]
figC_stats_filtered<-subset(figC_stats, value<=0.05)
write.table(figC_stats, file="barplots_stats/figC_stats.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(figC_stats_filtered, file="barplots_stats/figC_stats_filtered.txt", quote=FALSE, sep="\t", col.names=NA)


```


# to get number of time points (n=) per groupings in bargraph, as per reviewer comments:
length<- sapply(levels(bar$bargroup), function(class) sapply(bar[bar$bargroup == class, c("unclassified_Porphyromonadaceae", "Bacteroides")], length))


pre<-bar[bar$bargroup==c("pre"), ]
