# Using R to calculate diversity measures

### Anna Seekatz
### 12.17.14

> load all of the necessary packages
Load your packages:
```{r, error=FALSE, message=FALSE, warning=FALSE}
library(vegan)
library(gplots)
library(RColorBrewer)
library(Heatplus)
library(gtools)
library(GMD)
library(devtools)
library(plyr)
library(ggplot2)
```

Change into your working directory: 
```{r}
setwd("~/Desktop/umich/projects/Relapse_model2/analysis/mothur_analysis")
list.files(path="~/Desktop/umich/projects/Relapse_model2/analysis/mothur_analysis")
```

> These are the files I used:
+ group summary file from mothur: mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.groups.summary
+ metadata file created previously: mouse.var.txt
+ (optional) pcoa axes added from mothur file: mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.pcoa.axes
+ (optional) nmds axes added from mothur file: mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.nmds.axes

First, you must read in your files and merge them together with the metadata:

```{r}
#read in relevant files:

mouse.var<-read.table(file="mouse.var.txt", header=TRUE)
head(mouse.var)

mouse.sum<-read.table(file="mouseRM_mothur.files/mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.groups.summary", header=TRUE)
#head(mouse.sum)
colnames(mouse.sum)[which(names(mouse.sum) == "group")] <- "sample_ID"		#changed the column name since already have a column named 'group' in mouse.var

mouse.pcoa<-read.table(file="mouseRM_mothur.files/mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.pcoa.axes", header=TRUE)
colnames(mouse.pcoa)[which(names(mouse.pcoa) == "group")] <- "sample_ID"

mouse.nmds<-read.table(file="mouseRM_mothur.files/mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.nmds.axes", header=TRUE)
colnames(mouse.nmds)[which(names(mouse.nmds) == "group")] <- "sample_ID"

#combine files for a 'master' file:
#add sample information for 'group1':
v<-merge(mouse.var, mouse.sum, by.x=c("sampleID"), by.y=c("sample_ID")) #this merges the data based on the sampleID match
v1<-v[,c("sampleID", "cage", "day", "mouse", "exp", "treatment", "treatment2", "group_time", "group", "group_col", "group_time_col", "invsimpson", "shannon", "npshannon", "sobs", "nseqs")]  #this picks the specific columns
v1<-rename(v1, c("treatment2"="treat")) #another way to rename columns

v2<-merge(v1, mouse.pcoa, by.x=c("sampleID"), by.y=c("sample_ID"))
v2<-v2[,c("sampleID", "cage", "day", "mouse", "exp", "treatment", "treat", "group_time", "group", "group_col", "group_time_col", "invsimpson", "shannon", "npshannon", "sobs", "nseqs", "axis1", "axis2", "axis3", "axis4")] 
v2<-rename(v2, c("axis1"="pcoa_axis1", "axis2"="pcoa_axis2", "axis3"="pcoa_axis3", "axis4"="pcoa_axis4"))

v3<-merge(v2, mouse.nmds, by.x=c("sampleID"), by.y=c("sample_ID"))
v3<-v3[,c("sampleID", "cage", "day", "mouse", "exp", "treatment", "treat", "group_time", "group", "group_col", "group_time_col", "invsimpson", "shannon", "npshannon", "sobs", "nseqs", "pcoa_axis1", "pcoa_axis2", "pcoa_axis3", "pcoa_axis4", "axis1", "axis2", "axis3", "axis4")] 
mouse.div<-rename(v3, c("axis1"="nmds_axis1", "axis2"="nmds_axis2", "axis3"="nmds_axis3", "axis4"="nmds_axis4"))

head(mouse.div)
	#this file should have all of your metadata, diversity measures, and pcoa/nmds axes
write.table(mouse.div, file="mouseRM_div.thetayc.txt", quote=FALSE, sep="\t", col.names=NA)

```

> **Plotting your data:**

Now that you have your file ready, you can plot your data. You can do this either by individual points, summary statistics, or both. I've included some options below:

```{R}	
##plotting invsimpson over time (individually):

#separate individual points by group:
#first split into treatment groups:
all<-split(mouse.div, mouse.div$group)
abx.all<-all$'abx'
f1.all<-all$'CDI_FMT1'
f2.all<-all$'CDI_FMT2'
no.all<-all$'CDI_noFMT'
don.all<-all$'donor'

#plot the points by group, individually:
#par(mfrow=c(2,1))
plot(jitter(f1.all$day, amount=0.2), f1.all$invsimpson, 
	pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="chartreuse4", 
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
points(jitter(f2.all$day, amount=0.2), f2.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="yellowgreen")
points(jitter(don.all$day, amount=0.2), don.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="black")
points(jitter(abx.all$day, amount=0.2), abx.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="violetred1")
points(jitter(no.all$day, amount=0.2), no.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="blue4")
legend("bottomleft", c("FMT-d11", "FMT-d12","no FMT", "Donor", "Abx only"), pch=21, col="black", pt.bg=c("chartreuse4", "yellowgreen", "blue4", "black", "violetred1"), cex=0.6)
axis(1, at=f1.all$day, labels=f1.all$day)
#this looks really busy, so might want to stick to plotting the median or something
	#vs
plot(jitter(f1.all$day, amount=0.2), f1.all$invsimpson, 
	pch=19, xlim=c(-10,45), xaxt='n', col="chartreuse4", 
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
points(jitter(f2.all$day, amount=0.2), f2.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="yellowgreen")
points(jitter(don.all$day, amount=0.2), don.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="black")
points(jitter(abx.all$day, amount=0.2), abx.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="violetred1")
points(jitter(no.all$day, amount=0.2), no.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="blue4")
legend("bottomleft", c("FMT-d11", "FMT-d12","no FMT", "Donor", "Abx only"), col=c("chartreuse4", "yellowgreen", "blue4", "black", "violetred1"), pch=19, cex=0.6)
axis(1, at=f1.all$day, labels=f1.all$day)
```
You could also plot the average (or median, etc):

```{R}
#plotting with the average:

#calculate averages, sd:

#to get averages of multiple measures, you can use this function:
mouse.avg <- sapply(levels(mouse.div$group), function(x) sapply(mouse.div[mouse.div$group == x, c("shannon", "sobs")], mean)) 	#you can calculate anything you want with this, such as sd
mouse.avg
#               abx   CDI_FMT1   CDI_FMT2 CDI_noFMT      donor
#shannon  0.8984109   2.206616   2.056727  1.188316   3.397824
#sobs    48.1851852 118.524590 116.564706 41.240000 164.750000
#you could do this with any summary statistic, which will produce one summary statistic for many diversity measures.

#you can also produce one summary table of multiple summary statistics for one diversity measure. This is what I did:
#use another function to get average, sd, se in the same table:
mdata <- ddply(mouse.div, c("group", "day"), summarise,
               N    = length(invsimpson),
               mean = mean(invsimpson),
               sd   = sd(invsimpson),
               se   = sd / sqrt(N) )
head(mdata)
#       group day  N      mean          sd           se
#1        abx  -7  3 19.013596 2.242470002 1.2946906592
#2        abx  -2  3  5.338352 2.439503777 1.4084481623
#3        abx   0  3  1.003559 0.001164831 0.0006725153

#separate into separate groups:
abx<-subset(mdata, (subset=group %in% c("abx")))				#subset data by group
cdi.f1<-subset(mdata, (subset=group %in% c("CDI_FMT1")))
cdi.f2<-subset(mdata, (subset=group %in% c("CDI_FMT2")))
cdi<-subset(mdata, (subset=group %in% c("CDI_noFMT")))
don<-subset(mdata, (subset=group %in% c("donor")))

#now, you can graph the lines on each of the individual graphs:

#par(mfrow=c(2,1))
plot(jitter(f1.all$day, amount=0.2), f1.all$invsimpson, 
	pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="chartreuse4", cex=0.8, 
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
lines(cdi.f1$day, cdi.f1$mean, col="chartreuse4")
lines(cdi.f2$day, cdi.f2$mean, col="yellowgreen")
lines(cdi$day, cdi$mean, col="blue4")
lines(abx$day, abx$mean, col="violetred1")
points(jitter(f2.all$day, amount=0.2), f2.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="yellowgreen", cex=0.8)
points(jitter(don.all$day, amount=0.2), don.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="black", cex=0.8)
points(jitter(abx.all$day, amount=0.2), abx.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="violetred1", cex=0.8)
points(jitter(no.all$day, amount=0.2), no.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="blue4", cex=0.8)
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=2)
#legend("bottomleft", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=2, merge=FALSE, trace=TRUE)
	#note: this legend is not perfect (I cannot separate the circles and the lines), and I cannot figure out how to add lines using the info given by trace and using the 'arrows' command
axis(1, at=f1.all$day, labels=f1.all$day)
	#OR:
#this looks really busy, so might want to stick to plotting the median or something
	#vs
plot(jitter(f1.all$day, amount=0.2), f1.all$invsimpson, 
	pch=19, xlim=c(-10,45), xaxt='n', col="chartreuse4", cex=0.8,
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
lines(cdi.f1$day, cdi.f1$mean, col="chartreuse4")
lines(cdi.f2$day, cdi.f2$mean, col="yellowgreen")
lines(cdi$day, cdi$mean, col="blue4")
lines(abx$day, abx$mean, col="violetred1")
points(jitter(f2.all$day, amount=0.2), f2.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="yellowgreen", cex=0.8)
points(jitter(don.all$day, amount=0.2), don.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="black", cex=0.8)
points(jitter(abx.all$day, amount=0.2), abx.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="violetred1", cex=0.8)
points(jitter(no.all$day, amount=0.2), no.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="blue4", cex=0.8)
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=2)
axis(1, at=f1.all$day, labels=f1.all$day)

#You could also only plot the the means, with the sd:
plot(cdi.f1$day, cdi.f1$mean, 
	pch=19, xlim=c(-10,45), ylim=c(-1,25), xaxt='n', col="chartreuse4", cex=1, 
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
points(cdi.f2$day, cdi.f2$mean, pch=19, xlim=c(-10,45), xaxt='n', col="yellowgreen", cex=1)
points(don$day, don$mean, pch=19, xlim=c(-10,45), xaxt='n', col="black", cex=1)
points(abx$day, abx$mean, pch=19, xlim=c(-10,45), xaxt='n', col="violetred1", cex=1)
points(cdi$day, cdi$mean, pch=19, xlim=c(-10,45), xaxt='n', col="blue4", cex=1)
arrows(cdi.f1$day, cdi.f1$mean-cdi.f1$sd, cdi.f1$day, cdi.f1$mean+cdi.f1$sd, length=0.05, angle=90, code=3, col="chartreuse4")
arrows(cdi.f2$day, cdi.f2$mean-cdi.f2$sd, cdi.f2$day, cdi.f2$mean+cdi.f2$sd, length=0.05, angle=90, code=3, col="yellowgreen")
arrows(cdi$day, cdi$mean-cdi$sd, cdi$day, cdi$mean+cdi$sd, length=0.05, angle=90, code=3, col="blue4")
arrows(abx$day, abx$mean-abx$sd, abx$day, abx$mean+abx$sd, length=0.05, angle=90, code=3, col="violetred1")
arrows(don$day, don$mean-don$sd, don$day, don$mean+don$sd, length=0.05, angle=90, code=3, col="black")
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=2)
axis(1, at=f1.all$day, labels=f1.all$day)
		#optional lines:
#lines(cdi.f1$day, cdi.f1$mean, col="chartreuse4")
#lines(cdi.f2$day, cdi.f2$mean, col="yellowgreen")
#lines(cdi$day, cdi$mean, col="blue4")
#lines(abx$day, abx$mean, col="violetred1")
```
As you can see, there are multiple ways of presenting your data. Below, I have some other options:

```{R}
##Different ways of presenting the inverse Simpson, altogether:

par(mfrow=c(2,2))
plot(jitter(f1.all$day, amount=0.2), f1.all$invsimpson, 
	pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="chartreuse4", cex=0.8, 
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
lines(cdi.f1$day, cdi.f1$mean, col="chartreuse4")
lines(cdi.f2$day, cdi.f2$mean, col="yellowgreen")
lines(cdi$day, cdi$mean, col="blue4")
lines(abx$day, abx$mean, col="violetred1")
points(jitter(f2.all$day, amount=0.2), f2.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="yellowgreen", cex=0.8)
points(jitter(don.all$day, amount=0.2), don.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="black", cex=0.8)
points(jitter(abx.all$day, amount=0.2), abx.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="violetred1", cex=0.8)
points(jitter(no.all$day, amount=0.2), no.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="blue4", cex=0.8)
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=1)
axis(1, at=f1.all$day, labels=f1.all$day)

plot(jitter(f1.all$day, amount=0.2), f1.all$invsimpson, 
	pch=19, xlim=c(-10,45), xaxt='n', col="chartreuse4", cex=0.8,
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
lines(cdi.f1$day, cdi.f1$mean, col="chartreuse4")
lines(cdi.f2$day, cdi.f2$mean, col="yellowgreen")
lines(cdi$day, cdi$mean, col="blue4")
lines(abx$day, abx$mean, col="violetred1")
points(jitter(f2.all$day, amount=0.2), f2.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="yellowgreen", cex=0.8)
points(jitter(don.all$day, amount=0.2), don.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="black", cex=0.8)
points(jitter(abx.all$day, amount=0.2), abx.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="violetred1", cex=0.8)
points(jitter(no.all$day, amount=0.2), no.all$invsimpson, pch=19, xlim=c(-10,45), xaxt='n', col="blue4", cex=0.8)
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=1)
axis(1, at=f1.all$day, labels=f1.all$day)

plot(cdi.f1$day, cdi.f1$mean, 
	pch=19, xlim=c(-10,45), ylim=c(-1,35), xaxt='n', col="chartreuse4", cex=1, 
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
points(cdi.f2$day, cdi.f2$mean, pch=19, xlim=c(-10,45), xaxt='n', col="yellowgreen", cex=1)
points(don$day, don$mean, pch=19, xlim=c(-10,45), xaxt='n', col="black", cex=1)
points(abx$day, abx$mean, pch=19, xlim=c(-10,45), xaxt='n', col="violetred1", cex=1)
points(cdi$day, cdi$mean, pch=19, xlim=c(-10,45), xaxt='n', col="blue4", cex=1)
arrows(cdi.f1$day, cdi.f1$mean-cdi.f1$sd, cdi.f1$day, cdi.f1$mean+cdi.f1$sd, length=0.05, angle=90, code=3, col="chartreuse4")
arrows(cdi.f2$day, cdi.f2$mean-cdi.f2$sd, cdi.f2$day, cdi.f2$mean+cdi.f2$sd, length=0.05, angle=90, code=3, col="yellowgreen")
arrows(cdi$day, cdi$mean-cdi$sd, cdi$day, cdi$mean+cdi$sd, length=0.05, angle=90, code=3, col="blue4")
arrows(abx$day, abx$mean-abx$sd, abx$day, abx$mean+abx$sd, length=0.05, angle=90, code=3, col="violetred1")
arrows(don$day, don$mean-don$sd, don$day, don$mean+don$sd, length=0.05, angle=90, code=3, col="black")
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=1)
axis(1, at=f1.all$day, labels=f1.all$day)

plot(cdi.f1$day, cdi.f1$mean, 
	pch=21, xlim=c(-10,45), ylim=c(-1,35), xaxt='n', col="black", bg="chartreuse4", cex=1, 
	main="InvSimpson overtime (individual mice)", ylab="Inverse Simpson diversity index", xlab="Day")
arrows(cdi.f1$day, cdi.f1$mean-cdi.f1$sd, cdi.f1$day, cdi.f1$mean+cdi.f1$sd, length=0.05, angle=90, code=3, col="chartreuse4")
arrows(cdi.f2$day, cdi.f2$mean-cdi.f2$sd, cdi.f2$day, cdi.f2$mean+cdi.f2$sd, length=0.05, angle=90, code=3, col="yellowgreen")
arrows(cdi$day, cdi$mean-cdi$sd, cdi$day, cdi$mean+cdi$sd, length=0.05, angle=90, code=3, col="blue4")
arrows(abx$day, abx$mean-abx$sd, abx$day, abx$mean+abx$sd, length=0.05, angle=90, code=3, col="violetred1")
arrows(don$day, don$mean-don$sd, don$day, don$mean+don$sd, length=0.05, angle=90, code=3, col="black")
lines(cdi.f1$day, cdi.f1$mean, col="chartreuse4")
lines(cdi.f2$day, cdi.f2$mean, col="yellowgreen")
lines(cdi$day, cdi$mean, col="blue4")
lines(abx$day, abx$mean, col="violetred1")
points(cdi.f2$day, cdi.f2$mean, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="yellowgreen", cex=1)
points(don$day, don$mean, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="black", cex=1)
points(abx$day, abx$mean, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="violetred1", cex=1)
points(cdi$day, cdi$mean, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="blue4", cex=1)
points(cdi.f1$day, cdi.f1$mean, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="chartreuse4", cex=1) #to get the points on top of the graph, again
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=0.6, lty=1)
axis(1, at=f1.all$day, labels=f1.all$day)
legend(-10:45, -1:35)
```

> To do stats between FMT/no_FMT groups:

```{r}
# statistics for simpson/change over time (Fig. 5):

#previously created table from above for thetayc pcoa/nmdd/diversity:
mouse.div<-read.table(file="mouseRM_div.thetayc.txt", header=TRUE)

# let's remove the 'abx_only' group for the time being, since we are interested in the FMT vs. no_FMT mice:
sub.mouse.div<-subset(mouse.div, (subset=treat %in% c("FMT", "no_FMT")))

# split by days:
div.split<-split(sub.mouse.div, sub.mouse.div$day)
div.d10<-div.split$'10'
div.d12<-div.split$'12'
div.d16<-div.split$'16'
div.d19<-div.split$'19'
div.d22<-div.split$'22'
div.d27<-div.split$'27'
div.d33<-div.split$'33'
div.d37<-div.split$'37'
div.d36<-div.split$'36'
div.d42<-div.split$'42'

# do stats:
wilcox.test(invsimpson~treat, data=div.d12, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d10, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d16, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d19, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d22, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d27, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d33, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d36, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d37, paired=FALSE)
wilcox.test(invsimpson~treat, data=div.d42, paired=FALSE)

# Could also do Kruskal-Wallis to ensure that differences are there between all 3 groups:
kruskal.test(invsimpson~treatment, data=div.d12)
kruskal.test(invsimpson~treatment, data=div.d10)
kruskal.test(invsimpson~treatment, data=div.d16)
kruskal.test(invsimpson~treatment, data=div.d19)
kruskal.test(invsimpson~treatment, data=div.d22)
kruskal.test(invsimpson~treatment, data=div.d27)
kruskal.test(invsimpson~treatment, data=div.d33)
kruskal.test(invsimpson~treatment, data=div.d36)
kruskal.test(invsimpson~treatment, data=div.d37)
kruskal.test(invsimpson~treatment, data=div.d42)

# some of these values look weird, so investigated post-hoc tests for multiple groups
# found that there is a 'pairwise.wilcox.test' that can be applied:

pairwise.wilcox.test(div.d12$invsimpson, div.d12$treatment, p.adj="bonferroni")
pairwise.wilcox.test(div.d16$invsimpson, div.d16$treatment, p.adj="bonferroni")
pairwise.wilcox.test(div.d19$invsimpson, div.d19$treatment, p.adj="bonferroni")

```

############

Final code used for Fig. 5A:

```{r}
#simpson diversity:
plot(jitter(f1.all$day, amount=0.2), f1.all$invsimpson, 
	pch=21, xlim=c(-8,45), xaxt='n', col="black", bg="chartreuse4", cex=1, 
	ylab="Inverse Simpson diversity index")
lines(cdi.f1$day, cdi.f1$mean, col="chartreuse4")
lines(cdi.f2$day, cdi.f2$mean, col="yellowgreen")
lines(cdi$day, cdi$mean, col="blue4")
lines(abx$day, abx$mean, col="violetred1")
points(jitter(f2.all$day, amount=0.2), f2.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="yellowgreen", cex=1)
points(jitter(don.all$day, amount=0.2), don.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="black", cex=1)
points(jitter(abx.all$day, amount=0.2), abx.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="violetred1", cex=1)
points(jitter(no.all$day, amount=0.2), no.all$invsimpson, pch=21, xlim=c(-10,45), xaxt='n', col="black", bg="blue4", cex=1)
legend("topright", c("FMT-d11", "FMT-d12","no FMT", "Abx only", "Donor"), pch=21, col=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), pt.bg=c("chartreuse4", "yellowgreen", "blue4", "violetred1", "black"), cex=1, lty=1)
axis(1, at=f1.all$day, labels=f1.all$day)

```

# 6.8.15
# to test for statistical differences between FMTs on different dates (reported in IAI reviewer comments):
```{r}

#previously created table from above for thetayc pcoa/nmdd/diversity:
mouse.div<-read.table(file="mouseRM_div.thetayc.txt", header=TRUE)

# let's remove the 'abx_only' group for the time being, since we are interested in the FMT vs. no_FMT mice:
sub.mouse.div<-subset(mouse.div, (subset=treat %in% c("FMT", "no_FMT")))

div.split<-split(sub.mouse.div, sub.mouse.div$day)
div.d16<-div.split$'16'
div.d19<-div.split$'19'

d16<-split(div.d16, div.d16$treatment)
fmt2<-d16$'FMT2'
fmt3<-d16$'FMT3'
fmt4<-d16$'FMT4'
all<-rbind(fmt2, fmt3, fmt4)

d19<-split(div.d19, div.d19$treatment)
fmt2<-d19$'FMT2'
fmt3<-d19$'FMT3'
fmt4<-d19$'FMT4'
all19<-rbind(fmt2, fmt3, fmt4)

d22<-split(div.d22, div.d22$treatment)
fmt2<-d22$'FMT2'
fmt3<-d22$'FMT3'
fmt4<-d22$'FMT4'
all22<-rbind(fmt2, fmt3, fmt4)

d22<-split(div.d22, div.d22$treatment)
fmt2<-d22$'FMT2'
fmt3<-d22$'FMT3'
fmt4<-d22$'FMT4'
all22<-rbind(fmt2, fmt3, fmt4)

d27<-split(div.d27, div.d27$treatment)
fmt2<-d27$'FMT2'
fmt3<-d27$'FMT3'
fmt4<-d27$'FMT4'
all27<-rbind(fmt2, fmt3, fmt4)

kruskal.test(invsimpson~treatment, data=all)			#d16 is still different
kruskal.test(invsimpson~treatment, data=all19)			#however, by the time C. diff is 'cleared', no statistical difference between different FMT days
kruskal.test(invsimpson~treatment, data=all22)
kruskal.test(invsimpson~treatment, data=all27)