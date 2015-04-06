# Different thetayc calculations for longitudinal graphs
### Anna Seekatz
### 1.20.15

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
setwd("/Users/aseekatz/Desktop/umich/projects/Relapse_model2/mouseRM_analysis/mouseRM_non.mb.data")
list.files(path="/Users/aseekatz/Desktop/umich/projects/Relapse_model2/mouseRM_analysis/mouseRM_non.mb.data")
```

> These are the files I used:
+ tox1: RM1.2.FMT_col.tox.txt  			#this contains all toxin information from first 2 repeats of the relapse model
+ tox2: RM3.FMT_col.tox.txt				#this contains all toxin info from the 3rd (FMT) experiment
+ all.var: all.samples_mouse.var.txt	#this contains all the metadata necessary to categorize and group the data

> Merge your information with the metadata:

```{r}
tox1<-read.table(file="RM1.2.FMT_col.tox.txt", header=TRUE)
tox2<-read.table(file="RM3.FMT_col.tox.txt", header=TRUE)
all.var<-read.table(file="all.samples_mouse.var.txt", header=TRUE)

#head(tox1)			#to give you an idea of what the data looks like
#sampleID group treatment day cage mouse exp colonization  type weight log_colonization toxin
#1 d1_901_1   ABX   control   1  901     1   1          100 stool  103.4         2.000000    NA
#2 d1_901_2   ABX   control   1  901     2   1          100 stool   98.5         2.000000    NA
#3 d1_901_3   ABX   control   1  901     3   1          100 stool  103.5         2.000000    NA

#head(all.var)
#sampleID Day cage_mouse exp cage mouse treatment treatment2 group_time group
#1 dn7_901_1  -7      901_1   1  901     1       abx   abx_only        pre   abx
#2 dn7_901_2  -7      901_2   1  901     2       abx   abx_only        pre   abx
#3 dn7_901_3  -7      901_3   1  901     3       abx   abx_only        pre   abx

#combine to get the same variables that you want as before:

#merge to combine the metadata for each experimental group:
t1<-merge(tox1, all.var, by.x=c("sampleID"), by.y=c("sampleID")) #this merges the data based on sampleID (note: the columns do not have to match)
t1.all<-t1[,c("sampleID", "day", "cage.x", "mouse.x", "treatment.y", "treatment2", "group_time", "group.y", "type", "log_colonization", "colonization", "toxin")]  #this picks the specific columns you want to incorporate

t2<-merge(tox2, all.var, by.x=c("sampleID"), by.y=c("sampleID")) #this merges the data based on sampleID (note: the columns do not have to match)
t2.all<-t2[,c("sampleID", "day", "cage.x", "mouse.y", "treatment.y", "treatment2", "group_time", "group", "type", "log_colonization", "colonization", "toxin")]
head(t2.all)
#sampleID day cage.x mouse.y treatment.y treatment2 group_time     group type log_colonization colonization toxin
#1 d0_1006_1   0   1006       1         abx     no_FMT   post_cef CDI_noFMT <NA>               NA           NA    NA
#2 d0_1006_2   0   1006       2         abx     no_FMT   post_cef CDI_noFMT <NA>               NA           NA    NA
	# you now have the appropriate metadata for your analysis
```
> You are now ready to plot both the toxin and cfu counts.

> Toxin:

```{r}
#########
# For toxin:
######### 

#within each toxin set, calculate the median and interquartile range:

#Now we calculate the mean and sd for whatever variable we want to include (by treatment, etc) within each day:

#tox1:
t1.sum <- ddply(t1.all, c("treatment2", "day"), na.rm=TRUE, summarise,
               N    = length(toxin),
               median = median(toxin, na.rm=TRUE),
               iqr   = IQR(toxin, na.rm=TRUE),
               lowq	 = quantile(toxin, 0.25, na.rm=TRUE),
               highq  = quantile(toxin, 0.75, na.rm=TRUE)
               )
#note: you must use na.rm=TRUE within the functions, since there is technically missing data from the original N (the mice were sacrificed at some points, so there is no data available)

head(t1.sum)
#   treatment2 day na.rm N median iqr
#1   abx_only   1  TRUE 7     NA  NA
#2   abx_only   4  TRUE 6     NA  NA
#3   abx_only   7  TRUE 2     NA  NA

#split up summary by group:
t1.split.sum<-split(t1.sum, t1.sum$treatment2)
t1.fmt<-t1.split.sum$'FMT'
t1.nofmt<-t1.split.sum$'no_FMT'
t1.novanco<-t1.split.sum$'no_vanco'

#also split up individual points:
t1.split.all<-split(t1.all, t1.all$treatment2)
t1.ind.fmt<-t1.split.all$'FMT'
t1.ind.nofmt<-t1.split.all$'no_FMT'
t1.ind.novanco<-t1.split.all$'no_vanco'

#plot the groups:

#if separating FMT from no_fmt:
plot(0.25+jitter(t1.ind.fmt$day, amount=0.5), t1.ind.fmt$toxin, 
	pch=19, xaxt='n', col="chartreuse4", cex=0.8, ylim=c(0,6),
	main="toxin (cecal)", ylab="log10 (reciprocal dilution toxin per g of sample)", xlab="Day")
axis(1, at=t1.ind.fmt$day, labels=t1.ind.fmt$day)
arrows(0.25+t1.fmt$day, t1.fmt$lowq, 0.25+t1.fmt$day, t1.fmt$highq, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
arrows(0.25+t1.fmt$day-0.05, t1.fmt$median, 0.25+t1.fmt$day+0.05, t1.fmt$median, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
points(jitter(t1.ind.nofmt$day, amount=0.5)-0.25, t1.ind.nofmt$toxin, pch=19, col="blue4", cex=0.8)
arrows(t1.nofmt$day-0.25, t1.nofmt$lowq, t1.nofmt$day-0.25, t1.nofmt$highq, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(t1.nofmt$day-0.05-0.25, t1.nofmt$median, t1.nofmt$day+0.05-0.25, t1.nofmt$median, angle=180, code=3, col="blue4", lwd=4, length=.10)
points(jitter(t1.ind.novanco$day, amount=0.5), t1.ind.novanco$toxin, pch=19, col="orange", cex=0.8)
arrows(t1.novanco$day, t1.novanco$lowq, t1.novanco$day, t1.novanco$highq, length=0.10, angle=90, code=3, col="orange", lwd=2)
arrows(t1.novanco$day-0.05, t1.novanco$median, t1.novanco$day+0.05, t1.novanco$median, angle=180, code=3, col="orange", lwd=4, length=.10) 
legend("bottomleft", c("Relapsed (no FMT)", "Relapsed (FMT)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "chartreuse4", "orange"), pt.bg=c("blue4", "chartreuse4", "orange"), cex=0.8, lty=1)

#tox2:
t2.sum <- ddply(t2.all, c("treatment2", "day"), na.rm=TRUE, summarise,
               N    = length(toxin),
               median = median(toxin, na.rm=TRUE),
               iqr   = IQR(toxin, na.rm=TRUE),
               lowq	 = quantile(toxin, 0.25, na.rm=TRUE),
               highq  = quantile(toxin, 0.75, na.rm=TRUE)
               )

#split up summary by group:
t2.split.sum<-split(t2.sum, t2.sum$treatment2)
t2.fmt<-t2.split.sum$'FMT'
t2.nofmt<-t2.split.sum$'no_FMT'

#also split up individual points:
t2.split.all<-split(t2.all, t2.all$treatment2)
t2.ind.fmt<-t2.split.all$'FMT'
t2.ind.nofmt<-t2.split.all$'no_FMT'

#plot the groups:
plot(0.25+jitter(t2.ind.fmt$day, amount=0.5), t2.ind.fmt$toxin, 
	pch=19, xaxt='n', col="chartreuse4", cex=0.8, ylim=c(0,7),
	main="toxin (cecal)", ylab="log10 (reciprocal dilution toxin per g of sample)", xlab="Day")
axis(1, at=t2.ind.fmt$day, labels=t2.ind.fmt$day)
arrows(0.25+t2.fmt$day, t2.fmt$lowq, 0.25+t2.fmt$day, t2.fmt$highq,, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
arrows(0.25+t2.fmt$day-0.05, t2.fmt$median, 0.25+t2.fmt$day+0.05, t2.fmt$median, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
points(jitter(t2.ind.nofmt$day, amount=0.5)-0.25, t2.ind.nofmt$toxin, pch=19, col="blue4", cex=0.8)
arrows(t2.nofmt$day-0.25, t2.nofmt$lowq, t2.nofmt$day-0.25, t2.nofmt$highq, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(t2.nofmt$day-0.05-0.25, t2.nofmt$median, t2.nofmt$day+0.05-0.25, t2.nofmt$median, angle=180, code=3, col="blue4", lwd=4, length=.10)
legend("bottomleft", c("FMT", "no FMT"), pch=21, col=c("chartreuse4", "blue4"), pt.bg=c("chartreuse4", "blue4"), cex=0.8, lty=1)

--
##plotting both experiments:

par(mfrow=c(2,1))
plot(0.25+jitter(t1.ind.fmt$day, amount=0.5), t1.ind.fmt$toxin, 
	pch=19, xaxt='n', col="chartreuse4", cex=0.8, ylim=c(0,6),
	main="toxin (cecal)", ylab="log10 (reciprocal dilution toxin per g of sample)", xlab="Day")
axis(1, at=t1.ind.fmt$day, labels=t1.ind.fmt$day)
arrows(0.25+t1.fmt$day, t1.fmt$lowq, 0.25+t1.fmt$day, t1.fmt$highq, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
arrows(0.25+t1.fmt$day-0.05, t1.fmt$median, 0.25+t1.fmt$day+0.05, t1.fmt$median, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
points(jitter(t1.ind.nofmt$day, amount=0.5)-0.25, t1.ind.nofmt$toxin, pch=19, col="blue4", cex=0.8)
arrows(t1.nofmt$day-0.25, t1.nofmt$lowq, t1.nofmt$day-0.25, t1.nofmt$highq, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(t1.nofmt$day-0.05-0.25, t1.nofmt$median, t1.nofmt$day+0.05-0.25, t1.nofmt$median, angle=180, code=3, col="blue4", lwd=4, length=.10)
points(jitter(t1.ind.novanco$day, amount=0.5), t1.ind.novanco$toxin, pch=19, col="orange", cex=0.8)
arrows(t1.novanco$day, t1.novanco$lowq, t1.novanco$day, t1.novanco$highq, length=0.10, angle=90, code=3, col="orange", lwd=2)
arrows(t1.novanco$day-0.05, t1.novanco$median, t1.novanco$day+0.05, t1.novanco$median, angle=180, code=3, col="orange", lwd=4, length=.10) 
legend("bottomleft", c("Relapsed (no FMT)", "Relapsed (FMT)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "chartreuse4", "orange"), pt.bg=c("blue4", "chartreuse4", "orange"), cex=0.8, lty=1)

plot(0.25+jitter(t2.ind.fmt$day, amount=0.5), t2.ind.fmt$toxin, 
	pch=19, xaxt='n', col="chartreuse4", cex=0.8, ylim=c(0,7),
	main="toxin (cecal)", ylab="log10 (reciprocal dilution toxin per g of sample)", xlab="Day")
axis(1, at=t2.ind.fmt$day, labels=t2.ind.fmt$day)
arrows(0.25+t2.fmt$day, t2.fmt$lowq, 0.25+t2.fmt$day, t2.fmt$highq,, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
arrows(0.25+t2.fmt$day-0.05, t2.fmt$median, 0.25+t2.fmt$day+0.05, t2.fmt$median, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
points(jitter(t2.ind.nofmt$day, amount=0.5)-0.25, t2.ind.nofmt$toxin, pch=19, col="blue4", cex=0.8)
arrows(t2.nofmt$day-0.25, t2.nofmt$lowq, t2.nofmt$day-0.25, t2.nofmt$highq, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(t2.nofmt$day-0.05-0.25, t2.nofmt$median, t2.nofmt$day+0.05-0.25, t2.nofmt$median, angle=180, code=3, col="blue4", lwd=4, length=.10)
legend("bottomleft", c("FMT", "no FMT"), pch=21, col=c("chartreuse4", "blue4"), pt.bg=c("chartreuse4", "blue4"), cex=0.8, lty=1)

```
> CFU:

Change into your working directory: 
```{r}
#########
# For CFU:
######### 

#You can use the same sheet from above (t1.all and t2.all)

#first, t1 data set:

#calculate summary statistics for colonization instead of toxin--except we use mean and sd again:
               
cfu1.sum <- ddply(t1.all, c("treatment2", "day"), summarise,
               N    = length(colonization),
               mean = mean(log10(colonization), na.rm=TRUE),
               sd   = sd(log10(colonization), na.rm=TRUE),
               se   = sd / sqrt(N) )
               
#separate each group, again:
#split up summary by group:
cfu1.split.sum<-split(cfu1.sum, cfu1.sum$treatment2)
cfu1.fmt<-cfu1.split.sum$'FMT'
cfu1.nofmt<-cfu1.split.sum$'no_FMT'
cfu1.novanco<-cfu1.split.sum$'no_vanco'

#also split up individual points:
cfu1.split.all<-split(t1.all, t1.all$treatment2)
cfu1.ind.fmt<-cfu1.split.all$'FMT'
cfu1.ind.nofmt<-cfu1.split.all$'no_FMT'
cfu1.ind.novanco<-cfu1.split.all$'no_vanco'

#plot using log scale:

#all points:
plot(0.3+jitter(cfu1.ind.fmt$day, amount=0.2), log10(cfu1.ind.fmt$colonization), 
	pch=19, xaxt='n', col="chartreuse4", cex=0.8, log="y",
	main="C. difficile load", ylab="CFU/g sample", xlab="Day", yaxt="n", xlim=c(0.5,16.5))
axis(1, at=cfu1.ind.fmt$day, labels=cfu1.ind.fmt$day)
labels <- sapply(aty,function(i)
            as.expression(bquote(10^ .(i)))
          )															#this labels your logscale axis nicely
axis(2,at=aty,labels=labels)
arrows(0.3+cfu1.fmt$day, cfu1.fmt$mean-cfu1.fmt$sd, 0.3+cfu1.fmt$day, cfu1.fmt$mean+cfu1.fmt$sd, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
arrows(0.3+cfu1.fmt$day-0.05, cfu1.fmt$mean, 0.3+cfu1.fmt$day+0.05, cfu1.fmt$mean, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
points(jitter(cfu1.ind.nofmt$day, amount=0.2), log10(cfu1.ind.nofmt$colonization), 
	pch=19, col="blue4", cex=0.8,)
arrows(cfu1.nofmt$day, cfu1.nofmt$mean-cfu1.nofmt$sd, cfu1.nofmt$day, cfu1.nofmt$mean+cfu1.nofmt$sd, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(cfu1.nofmt$day-0.05, cfu1.nofmt$mean, cfu1.nofmt$day+0.05, cfu1.nofmt$mean, angle=180, code=3, col="blue4", lwd=4, length=.10)
points(jitter(cfu1.ind.novanco$day, amount=0.2)-0.3, log10(cfu1.ind.novanco$colonization), 
	pch=19, col="orange", cex=0.8,)
arrows(cfu1.novanco$day-0.3, cfu1.novanco$mean-cfu1.novanco$sd, cfu1.novanco$day-0.3, cfu1.novanco$mean+cfu1.novanco$sd, length=0.10, angle=90, code=3, col="orange", lwd=2)
arrows(cfu1.novanco$day-0.05-0.3, cfu1.novanco$mean, cfu1.novanco$day+0.05-0.3, cfu1.novanco$mean, angle=180, code=3, col="orange", lwd=4, length=.10)
legend("bottomleft", c("Relapsed (no FMT)", "Relapsed (FMT)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "chartreuse4", "orange"), pt.bg=c("blue4", "chartreuse4", "orange"), cex=0.8, lty=1)

#vs.
#only the mean plus sd:

plot(0.3+cfu1.fmt$day, cfu1.fmt$mean, 
	pch=21, xaxt='n', bg="chartreuse4", cex=1, log="y", col="black",
	main="C. difficile load", ylab="CFU/g sample", xlab="Day", xlim=c(0.5,16.5))
axis(1, at=cfu1.ind.fmt$day, labels=cfu1.ind.fmt$day)
arrows(0.3+cfu1.fmt$day, cfu1.fmt$mean-cfu1.fmt$sd, 0.3+cfu1.fmt$day, cfu1.fmt$mean+cfu1.fmt$sd, length=0.10, angle=90, code=3, col="chartreuse4", lwd=1)
arrows(cfu1.nofmt$day, cfu1.nofmt$mean-cfu1.nofmt$sd, cfu1.nofmt$day, cfu1.nofmt$mean+cfu1.nofmt$sd, length=0.10, angle=90, code=3, col="blue4", lwd=1)
arrows(cfu1.novanco$day-0.3, cfu1.novanco$mean-cfu1.novanco$sd, cfu1.novanco$day-0.3, cfu1.novanco$mean+cfu1.novanco$sd, length=0.10, angle=90, code=3, col="orange", lwd=1)
lines(0.3+cfu1.fmt$day, cfu1.fmt$mean, col="chartreuse4")
lines(cfu1.nofmt$day, cfu1.nofmt$mean, col="blue4", na.rm=TRUE)
lines(cfu1.novanco$day-0.3, cfu1.novanco$mean, col="orange")
points(0.3+cfu1.fmt$day, cfu1.fmt$mean, pch=21, bg="chartreuse4", cex=1, col="black")
points(cfu1.nofmt$day, cfu1.nofmt$mean, pch=21, bg="blue4", cex=1, col="black")
points(cfu1.novanco$day-0.3, cfu1.novanco$mean, pch=21, bg="orange", cex=1, col="black")
legend("bottomleft", c("Relapsed (no FMT)", "Relapsed (FMT)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "chartreuse4", "orange"), pt.bg=c("blue4", "chartreuse4", "orange"), cex=0.8, lty=1)

#for exp2:

cfu2.sum <- ddply(t2.all, c("treatment2", "day"), summarise,
               N    = length(colonization),
               mean = mean(log10(colonization), na.rm=TRUE),
               sd   = sd(log10(colonization), na.rm=TRUE),
               se   = sd / sqrt(N) )
               
#separate each group, again:
#split up summary by group:
cfu2.split.sum<-split(cfu2.sum, cfu2.sum$treatment2)
cfu2.fmt<-cfu2.split.sum$'FMT'
cfu2.nofmt<-cfu2.split.sum$'no_FMT'

#also split up individual points:
cfu2.split.all<-split(t2.all, t2.all$treatment2)
cfu2.ind.fmt<-cfu2.split.all$'FMT'
cfu2.ind.nofmt<-cfu2.split.all$'no_FMT'

#have some missing data, so need to eraser the rows where there is no data (days 2, 14):
cfu2.fmt <- cfu2.fmt[-c(1,7), ]
cfu2.nofmt <- cfu2.nofmt[-c(1,7), ]

#plot using log scale:

#all points:
plot(0.3+jitter(cfu2.ind.fmt$day, amount=0.2), log10(cfu2.ind.fmt$colonization), 
	pch=19, xaxt='n', col="chartreuse4", cex=0.8, log="y", ylim=range(c(cfu2.fmt$mean-cfu2.fmt$sd, cfu2.fmt$mean+cfu2.fmt$sd)),
	main="C. difficile load", ylab="CFU/g sample", xlab="Day", yaxt="n", xlim=c(0.5,42.5))
axis(1, at=cfu2.ind.fmt$day, labels=cfu2.ind.fmt$day)
aty <- axTicks(2)
labels <- sapply(aty,function(i)
            as.expression(bquote(10^ .(i)))
          )
axis(2,at=aty,labels=labels)														#this labels your logscale axis nicely
arrows(0.3+cfu2.fmt$day, cfu2.fmt$mean-cfu2.fmt$sd, 0.3+cfu2.fmt$day, cfu2.fmt$mean+cfu2.fmt$sd, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
arrows(0.3+cfu2.fmt$day-0.05, cfu2.fmt$mean, 0.3+cfu2.fmt$day+0.05, cfu2.fmt$mean, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
points(jitter(cfu2.ind.nofmt$day, amount=0.2), log10(cfu2.ind.nofmt$colonization), 
	pch=19, col="blue4", cex=0.8,)
arrows(cfu2.nofmt$day, cfu2.nofmt$mean-cfu2.nofmt$sd, cfu2.nofmt$day, cfu2.nofmt$mean+cfu2.nofmt$sd, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(cfu2.nofmt$day-0.05, cfu2.nofmt$mean, cfu2.nofmt$day+0.05, cfu2.nofmt$mean, angle=180, code=3, col="blue4", lwd=4, length=.10)
legend("bottomleft", c("FMT", "no FMT"), pch=21, col=c("chartreuse4", "blue4"), pt.bg=c("chartreuse4", "blue4"), cex=0.8, lty=1)

#vs.
#only the mean plus sd:

plot(0.3+cfu2.fmt$day, cfu2.fmt$mean, yaxt='n', 
	pch=21, xaxt='n', bg="chartreuse4", cex=1, log="y", col="black",
	main="C. difficile load", ylab="CFU/g sample", xlab="Day", xlim=c(0.5,42.5), ylim=range(c(cfu2.fmt$mean-cfu2.fmt$sd, cfu2.fmt$mean+cfu2.fmt$sd)))
axis(1, at=cfu2.ind.fmt$day, labels=cfu2.ind.fmt$day)
labels <- sapply(aty,function(i)
            as.expression(bquote(10^ .(i)))
          )	
axis(2,at=aty,labels=labels)
arrows(0.3+cfu2.fmt$day, cfu2.fmt$mean-cfu2.fmt$sd, 0.3+cfu2.fmt$day, cfu2.fmt$mean+cfu2.fmt$sd, length=0.10, angle=90, code=3, col="chartreuse4", lwd=1)
arrows(cfu2.nofmt$day, cfu2.nofmt$mean-cfu2.nofmt$sd, cfu2.nofmt$day, cfu2.nofmt$mean+cfu2.nofmt$sd, length=0.10, angle=90, code=3, col="blue4", lwd=1)
lines(0.3+cfu2.fmt$day, cfu2.fmt$mean, col="chartreuse4")
lines(cfu2.nofmt$day, cfu2.nofmt$mean, col="blue4", na.rm=TRUE)
points(0.3+cfu2.fmt$day, cfu2.fmt$mean, pch=21, bg="chartreuse4", cex=1, col="black")
points(cfu2.nofmt$day, cfu2.nofmt$mean, pch=21, bg="blue4", cex=1, col="black")
legend("bottomleft", c("FMT", "no FMT"), pch=21, col=c("chartreuse4", "blue4"), pt.bg=c("chartreuse4", "blue4"), cex=0.8, lty=1)

```

Great!

> Now, graph the CFU and toxin together, per each experiment:

Change into your working directory: 
```{r}
#########
# For CFU:
#########

##for experiment 1:

#exp 1 (Fig 1):
par(mfrow=c(2,1))

# CFU:
#note: excluding the FMT points for now
plot(0.3+cfu1.nofmt$day, cfu1.nofmt$mean, 
	pch=21, xaxt='n', yaxt="n", bg="blue4", cex=1, log="y", col="black",
	main="C. difficile load", ylab="CFU/g sample", xlab=NA, xlim=c(0.5,16.5), ylim=range(c(cfu1.nofmt$mean-cfu1.nofmt$sd, cfu1.nofmt$mean+cfu1.nofmt$sd)))
axis(1, at=cfu1.ind.nofmt$day, labels=cfu1.ind.nofmt$day, cex.axis=0.8)
aty <- axTicks(2)
labels <- sapply(aty,function(i)
            as.expression(bquote(10^ .(i)))
          )	
axis(2,at=aty,labels=labels, cex.axis=0.8)
#arrows(0.3+cfu1.fmt$day, cfu1.fmt$mean-cfu1.fmt$sd, 0.3+cfu1.fmt$day, cfu1.fmt$mean+cfu1.fmt$sd, length=0.10, angle=90, code=3, col="chartreuse4", lwd=1)
arrows(cfu1.nofmt$day, cfu1.nofmt$mean-cfu1.nofmt$sd, cfu1.nofmt$day, cfu1.nofmt$mean+cfu1.nofmt$sd, length=0.10, angle=90, code=3, col="blue4", lwd=1)
arrows(cfu1.novanco$day-0.3, cfu1.novanco$mean-cfu1.novanco$sd, cfu1.novanco$day-0.3, cfu1.novanco$mean+cfu1.novanco$sd, length=0.10, angle=90, code=3, col="orange", lwd=1)
#lines(0.3+cfu1.fmt$day, cfu1.fmt$mean, col="chartreuse4")
lines(cfu1.nofmt$day, cfu1.nofmt$mean, col="blue4", na.rm=TRUE)
lines(cfu1.novanco$day-0.3, cfu1.novanco$mean, col="orange")
#points(0.3+cfu1.fmt$day, cfu1.fmt$mean, pch=21, bg="chartreuse4", cex=1, col="black")
points(cfu1.nofmt$day, cfu1.nofmt$mean, pch=21, bg="blue4", cex=1, col="black")
points(cfu1.novanco$day-0.3, cfu1.novanco$mean, pch=21, bg="orange", cex=1, col="black")
#legend("bottomleft", c("Relapsed (no FMT)", "Relapsed (FMT)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "chartreuse4", "orange"), pt.bg=c("blue4", "chartreuse4", "orange"), cex=0.8, lty=1)
legend("bottomleft", c("Relapsed (vanco)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "orange"), pt.bg=c("blue4", "orange"), cex=0.8, lty=1)


# Toxin (note: changed the points so that they are a different type from above):
#note: excluding the FMT points for now
plot(0.25+jitter(t1.ind.nofmt$day, amount=0.5), t1.ind.nofmt$toxin, 
	pch=21, xaxt='n', col="black", bg="blue4", cex=0.75, ylim=c(0,6), cex.axis=0.8,
	main="toxin (cecal)", ylab="log10 (reciprocal dilution toxin per g of sample)", xlab=NA)
axis(1, at=t1.ind.fmt$day, labels=t1.ind.fmt$day, cex.axis=0.8)
#arrows(t1.fmt$day-0.25, t1.fmt$lowq, t1.fmt$day-0.25, t1.fmt$highq, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
#arrows(t1.fmt$day-0.05-0.25, t1.fmt$median, t1.fmt$day+0.05-0.25, t1.fmt$median, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
#points(jitter(t1.ind.fmt$day, amount=0.5)-0.25, t1.ind.fmt$toxin, pch=21, col="black", bg="chartreuse4", cex=0.75)
arrows(t1.nofmt$day+0.25, t1.nofmt$lowq, t1.nofmt$day+0.25, t1.nofmt$highq, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(t1.nofmt$day-0.05+0.25, t1.nofmt$median, t1.nofmt$day+0.05+0.25, t1.nofmt$median, angle=180, code=3, col="blue4", lwd=4, length=.10)
points(jitter(t1.ind.novanco$day, amount=0.5), t1.ind.novanco$toxin, pch=21, col="black", bg="orange", cex=0.75)
arrows(t1.novanco$day, t1.novanco$lowq, t1.novanco$day, t1.novanco$highq, length=0.10, angle=90, code=3, col="orange", lwd=2)
arrows(t1.novanco$day-0.05, t1.novanco$median, t1.novanco$day+0.05, t1.novanco$median, angle=180, code=3, col="orange", lwd=4, length=.10) 
#legend("bottomleft", c("Relapsed (no FMT)", "Relapsed (FMT)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "chartreuse4", "orange"), pt.bg=c("blue4", "chartreuse4", "orange"), cex=0.75, lty=1)
legend("bottomleft", c("Relapsed (vanco)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "orange"), pt.bg=c("blue4", "orange"), cex=0.75, lty=1)


#exp 2 (Fig 2):

# CFU:
par(mfrow=c(2,1))
plot(0.3+cfu2.fmt$day, cfu2.fmt$mean, yaxt='n', 
	pch=21, xaxt='n', bg="chartreuse4", cex=1, log="y", col="black",
	main="C. difficile load", ylab="CFU/g sample", xlab="Day", xlim=c(0.5,42.5), ylim=range(c(cfu2.fmt$mean-cfu2.fmt$sd, cfu2.fmt$mean+cfu2.fmt$sd)))
axis(1, at=cfu2.ind.fmt$day, labels=cfu2.ind.fmt$day, cex.axis=0.8)
labels <- sapply(aty,function(i)
            as.expression(bquote(10^ .(i)))
          )	
axis(2,at=aty,labels=labels, cex.axis=0.8)
arrows(0.3+cfu2.fmt$day, cfu2.fmt$mean-cfu2.fmt$sd, 0.3+cfu2.fmt$day, cfu2.fmt$mean+cfu2.fmt$sd, length=0.10, angle=90, code=3, col="chartreuse4", lwd=1)
arrows(cfu2.nofmt$day, cfu2.nofmt$mean-cfu2.nofmt$sd, cfu2.nofmt$day, cfu2.nofmt$mean+cfu2.nofmt$sd, length=0.10, angle=90, code=3, col="blue4", lwd=1)
lines(0.3+cfu2.fmt$day, cfu2.fmt$mean, col="chartreuse4")
lines(cfu2.nofmt$day, cfu2.nofmt$mean, col="blue4", na.rm=TRUE)
points(0.3+cfu2.fmt$day, cfu2.fmt$mean, pch=21, bg="chartreuse4", cex=1, col="black")
points(cfu2.nofmt$day, cfu2.nofmt$mean, pch=21, bg="blue4", cex=1, col="black")
legend("bottomleft", c("FMT", "no FMT"), pch=21, col=c("chartreuse4", "blue4"), pt.bg=c("chartreuse4", "blue4"), cex=0.7, lty=1)

# Toxin:
plot(0.25+jitter(t2.ind.fmt$day, amount=0.5), t2.ind.fmt$toxin, 
	pch=21, xaxt='n', col="black", bg="chartreuse4", cex=0.75, ylim=c(0,6),
	main="toxin (cecal)", ylab="log10 (reciprocal dilution toxin per g of sample)", xlab="Day", cex.axis=0.8)
axis(1, at=t2.ind.fmt$day, labels=t2.ind.fmt$day, cex.axis=0.8)
arrows(0.25+t2.fmt$day, t2.fmt$lowq, 0.25+t2.fmt$day, t2.fmt$highq,, length=0.10, angle=90, code=3, col="chartreuse4", lwd=2)
arrows(0.25+t2.fmt$day-0.05, t2.fmt$median, 0.25+t2.fmt$day+0.05, t2.fmt$median, angle=180, code=3, col="chartreuse4", lwd=4, length=.10)
points(jitter(t2.ind.nofmt$day, amount=0.5)-0.25, t2.ind.nofmt$toxin, pch=21, col="black", bg="blue4", cex=0.75)
arrows(t2.nofmt$day-0.25, t2.nofmt$lowq, t2.nofmt$day-0.25, t2.nofmt$highq, length=0.10, angle=90, code=3, col="blue4", lwd=2)
arrows(t2.nofmt$day-0.05-0.25, t2.nofmt$median, t2.nofmt$day+0.05-0.25, t2.nofmt$median, angle=180, code=3, col="blue4", lwd=4, length=.10)
legend("bottomleft", c("FMT", "no FMT"), pch=21, col=c("chartreuse4", "blue4"), pt.bg=c("chartreuse4", "blue4"), cex=0.7, lty=1)


```

# stats for toxin and CFU:

```{r}

## doing stats on toxin/cfu:
--since technically the toxin levels are ordinal (or continuous), you can use the Wilcoxon rank sum test (Mann-Whitney test, for two independent variables)

# exp 1:

# subset to only include groups 'no_FMT' and 'vanco_only':
exp1<-subset(t1.all, (subset=treatment2 %in% c("no_FMT", "no_vanco")))

# split data groups up into days:
exp1.split<-split(exp1, exp1$day)
exp1.d7<-exp1.split$'7'
exp1.d9<-exp1.split$'9'
exp1.d10<-exp1.split$'10'
exp1.d11<-exp1.split$'11'
exp1.d12<-exp1.split$'12'
exp1.d16<-exp1.split$'16'

# for toxin on days 9, 16:

wilcox.test(toxin~treatment2, data=exp1.d9, paired=FALSE)
#	Wilcoxon rank sum test with continuity correction
#data:  toxin by treatment2
#W = 0, p-value = 0.0002034
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(toxin~treatment2, data=exp1.d16, paired=FALSE)
#	Wilcoxon rank sum test with continuity correction\
#data:  toxin by treatment2
#W = 80, p-value = 0.0003687
#alternative hypothesis: true location shift is not equal to 0

# for colonization on days 7, 9, 10, 11, 12:
wilcox.test(log_colonization~treatment2, data=exp1.d7, paired=FALSE)
	#W = 0, p-value = 0.0004198
wilcox.test(log_colonization~treatment2, data=exp1.d9, paired=FALSE)
	#W = 0, p-value = 0.0001463
wilcox.test(log_colonization~treatment2, data=exp1.d10, paired=FALSE)
	#W = 0, p-value = 0.008911
wilcox.test(log_colonization~treatment2, data=exp1.d11, paired=FALSE)
	#W = 6, p-value = 0.01259
wilcox.test(log_colonization~treatment2, data=exp1.d12, paired=FALSE)
	#W = 14, p-value = 0.7619

##
# exp2:

# subset to only include groups 'no_FMT' and 'vanco_only':
exp2<-subset(t2.all, (subset=treatment2 %in% c("no_FMT", "FMT")))

# split data groups up into days:
exp2.split<-split(exp2, exp2$day)
exp2.d7<-exp2.split$'7'
exp2.d10<-exp2.split$'10'
exp2.d12<-exp2.split$'12'
exp2.d16<-exp2.split$'16'
exp2.d19<-exp2.split$'19'
exp2.d22<-exp2.split$'22'
exp2.d27<-exp2.split$'27'
exp2.d33<-exp2.split$'33'
exp2.d36<-exp2.split$'36'
exp2.d37<-exp2.split$'37'
exp2.d42<-exp2.split$'42'

# toxin on days 12, 15, 19, 22, 33, 36, 42:
wilcox.test(toxin~treatment2, data=exp2.d12, paired=FALSE)
	#W = 48.5, p-value = 0.5219
wilcox.test(toxin~treatment2, data=exp2.d16, paired=FALSE)
	#W = 2.5, p-value = 0.003726
wilcox.test(toxin~treatment2, data=exp2.d19, paired=FALSE)
	#W = 0, p-value = 2.424e-06
wilcox.test(toxin~treatment2, data=exp2.d22, paired=FALSE)	
	#W = 0, p-value = 6.726e-06
wilcox.test(toxin~treatment2, data=exp2.d33, paired=FALSE)	
	#W = 0, p-value = 2.424e-06
wilcox.test(toxin~treatment2, data=exp2.d37, paired=FALSE)	
	#W = 0, p-value = 1.091e-05
etc.

# colonization on days 12, 15, 19, 22, 33, 36, 42:
wilcox.test(log_colonization~treatment2, data=exp2.d12, paired=FALSE)
	#W = 81, p-value = 0.1337
wilcox.test(log_colonization~treatment2, data=exp2.d16, paired=FALSE)
	#W = 7, p-value = 0.01291
wilcox.test(log_colonization~treatment2, data=exp2.d19, paired=FALSE)
	#W = 0, p-value = 0.0065
wilcox.test(log_colonization~treatment2, data=exp2.d22, paired=FALSE)	
	#W = 0, p-value = 6.726e-06
wilcox.test(log_colonization~treatment2, data=exp2.d33, paired=FALSE)	
	#W = 2, p-value = 2.595e-05
wilcox.test(log_colonization~treatment2, data=exp2.d36, paired=FALSE)	
	#W = 0, p-value = 2.548e-06
wilcox.test(log_colonization~treatment2, data=exp2.d37, paired=FALSE)	
	#W = 0, p-value = 0.0002712
wilcox.test(log_colonization~treatment2, data=exp2.d36, paired=FALSE)	
	#W = 0, p-value = 3.194e-05

```