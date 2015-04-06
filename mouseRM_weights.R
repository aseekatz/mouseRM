# Plotting weights of mice by experiment

### Anna Seekatz
### 12.19.14

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
This tutorial will cover how to graph your weights by your chosen variable. 

> Reshaping the data:

I have collated the following data sheets to reflect all of the toxin, cfu, and weight data:

### Data sheets with all of the information for days that cfu's were counted:
+ mouseRM_toxin.col.weight.xlsx (this has several sheets)
+ RM1.2.FMT_col.tox.txt
+ RM3.FMT_col.tox.txt
+ all.samples_mouse.var.xlsx and all.samples_mouse.var.txt

### Data sheets of weight matrices:
+ RM1_weight.matrix.txt
+ RM2_weight.matrix.txt
+ RM3.FMT_weight.matrix.txt

First, we need to merge the weights together and add the metadata:

```{r}
#First, I need to reshape the data so that the weights are in a list form, not a matrix:

r1<-as.matrix(read.table(file="RM1_weight.matrix.txt", header=TRUE))	
head(r1)																#this is how the list looked before (as a matrix)
r1.list <- reshape(as.data.frame(r1), idvar="Day", timevar="mouse", 
        	varying = -1, direction = "long", sep = "", ids = 2:30(r1))
#r1.list <- reshape(as.data.frame(r1), idvar="Day", timevar="mouse", 
        	varying = -1, direction = "long", sep = "", ids = 2:30(r1), v.names="weight")
        #for whatever reason, this code will add the name of the variable (weight), but it will not use the name of the mouse--instead, it just numbers it
        #thus, I just rename the colum ID as before:
r1.list<-rename(r1.list, c("X"="weight"))
r1.list$experiment<- c(1)
head(r1.list)															#this is how the list looks now (as a list)

r2<-as.matrix(read.table(file="RM2_weight.matrix.txt", header=TRUE))
r2.list <- reshape(as.data.frame(r2), idvar="DAY", timevar="mouse", 
        	varying = -1, direction = "long", sep = "", ids = 2:37(r2))
r2.list<-rename(r2.list, c("X"="weight"))
r2.list<-rename(r2.list, c("DAY"="Day"))
r2.list$experiment<- c(2)

r3<-as.matrix(read.table(file="RM3.FMT_weight.matrix.txt", header=TRUE))
r3.list <- reshape(as.data.frame(r3), idvar="Day", timevar="mouse", 
        	varying = -1, direction = "long", sep = "", ids = 2:31(r3))
r3.list<-rename(r3.list, c("X"="weight"))
r3.list$experiment<- c(3)

#now you can add them all together (we still need to fix the sampleID name to match those before):
all.weights<-rbind(r1.list, r2.list, r3.list)

#I will now save this as a table, then use the column add function to rename each of the mice in the 'd*_cage_mouse' format
write.table(all.weights, file="all_weights.txt", quote=FALSE, sep="\t", col.names=TRUE)

#after modifying the file, we should be able to add some metadata using a full all.samples_mouse.var.txt file that has metadata per sample:

all.var<-read.table(file="all.samples_mouse.var.txt", header=TRUE)
weight<-read.table(file="all_weights.txt", header=TRUE)
head(all.var)															#checking what the files look like (head displays the first 6 rows of a file)
head(weight)

#merge to combine the metadata:
w<-merge(weight, all.var, by.x=c("sampleID"), by.y=c("sampleID")) #this merges the data based on sampleID (note: the columns do not have to match)
w.all<-w[,c("sampleID", "weight", "Day.x", "cage_mouse", "mouse.y", "exp", "treatment", "treatment2", "group_time", "group")]  #this picks the specific columns you want to incorporate
head(w.all)
	#looks good
```

> Merging the data:

I have three repeats of the same experiment. I would like to show a graph of the weights per experiment.
To do this, I will split up the data by experiment, and then calculate the mean/sd/se for the weights by treatment group and day:

```{r}
#now we can subset by experiment (I want to graph each repeat separately):
w1<-subset(w.all, (subset=exp %in% c(1)))
w2<-subset(w.all, (subset=exp %in% c(2)))
w3<-subset(w.all, (subset=exp %in% c(3)))

#another way of splitting these up:
all<-split(w.all, w.all$exp)
w1<-all$'1'
w2<-all$'2'
w3<-all$'3'


#Now we calculate the mean and sd for whatever variable we want to include (by treatment, etc) within each day:
#for exp. 1:
wd1 <- ddply(w1, c("treatment2", "Day.x"), na.rm=TRUE, summarise,
               N    = length(weight),
               mean = mean(weight, na.rm=TRUE),
               sd   = sd(weight, na.rm=TRUE),
               se   = sd / sqrt(N) )
#note: you must use na.rm=TRUE within the functions, since there is technically missing data from the original N (the mice were sacrificed at some points, so there is no data available)

head(wd1)
#   treatment2 Day.x  N      mean       sd        se
#1    abx_only    -7  4  93.77283 1.934594 0.9672970
#2    abx_only    -2  4  98.59980 3.262564 1.6312820
#3    abx_only     0  4 100.00000 0.000000 0.0000000

#split up means, sd by group:
ave<-split(wd1, wd1$treatment2)
abx1<-ave$'abx_only'
vanco1<-ave$'no_FMT'
nov1<-ave$'no_vanco'

#for exp 2:
wd2 <- ddply(w2, c("treatment2", "Day.x"), summarise,
               N    = length(weight),
               mean = mean(weight, na.rm=TRUE),
               sd   = sd(weight, na.rm=TRUE),
               se   = sd / sqrt(N) )

#split up means, sd by group:
ave<-split(wd2, wd2$treatment2)
abx2<-ave$'abx_only'
vanco2<-ave$'no_FMT'
nov2<-ave$'no_vanco'
fmt2<-ave$'FMT'

#for exp 3:
wd3 <- ddply(w3, c("treatment2", "Day.x"), na.rm=TRUE, summarise,
               N    = length(weight),
               mean = mean(weight, na.rm=TRUE),
               sd   = sd(weight, na.rm=TRUE),
               se   = sd / sqrt(N) )

#split up means, sd by group:
ave<-split(wd3, wd3$treatment2)
fmt3<-ave$'FMT'
nof3<-ave$'no_FMT'
```

> Plotting the data:

Now, we are ready to plot!

```{r}
par(mfrow=c(3,1))
plot(abx1$Day.x, abx1$mean, 
	pch=21, xlim=c(-8,18), ylim=c(80, 110), xaxt='n', col="black", bg="violetred1", cex=1, 
	main="Relapse model (exp 1)", ylab="mean weight change from d0 (%)", xlab="Day")
axis(1, at=abx1$Day.x, labels=abx1$Day.x)
lines(abx1$Day.x, abx1$mean, col="violetred1")
lines(vanco1$Day.x, vanco1$mean, col="blue4")
lines(nov1$Day.x, nov1$mean, col="orange")
arrows(abx1$Day.x, abx1$mean-abx1$sd, abx1$Day.x, abx1$mean+abx1$sd, length=0.05, angle=90, code=3, col="violetred1")
arrows(vanco1$Day.x, vanco1$mean-vanco1$sd, vanco1$Day.x, vanco1$mean+vanco1$sd, length=0.05, angle=90, code=3, col="blue4")
arrows(nov1$Day.x, nov1$mean-nov1$sd, nov1$Day.x, nov1$mean+nov1$sd, length=0.05, angle=90, code=3, col="orange")
points(vanco1$Day.x, vanco1$mean, pch=21, col="black", bg="blue4", xlim=c(-8,18), xaxt='n')
points(nov1$Day.x, nov1$mean, pch=21, col="black", bg="orange", xlim=c(-8,18), xaxt='n')
points(abx1$Day.x, abx1$mean, pch=21, col="black", bg="violetred1", xlim=c(-8,18), xaxt='n')		#replot the points so that they are over the lines
legend("bottomright", c("abx only", "Vanco (relapse)", "No vanco"), pch=21, col=c("black"), pt.bg=c("violetred1", "blue4", "orange"), cex=1)

plot(abx2$Day.x, abx2$mean, 
	pch=21, xlim=c(-8,18), ylim=c(80, 110), xaxt='n', col="black", bg="violetred1", cex=1, 
	main="Relapse model (exp 2)", ylab="mean weight change from d0 (%)", xlab="Day")
axis(1, at=abx2$Day.x, labels=abx2$Day.x)
lines(abx2$Day.x, abx2$mean, col="violetred1")
lines(vanco2$Day.x, vanco2$mean, col="blue4")
lines(nov2$Day.x, nov2$mean, col="orange")
lines(fmt2$Day.x, fmt2$mean, col="chartreuse4")
arrows(abx2$Day.x, abx2$mean-abx2$sd, abx2$Day.x, abx2$mean+abx2$sd, length=0.05, angle=90, code=3, col="violetred1")
arrows(vanco2$Day.x, vanco2$mean-vanco2$sd, vanco2$Day.x, vanco2$mean+vanco2$sd, length=0.05, angle=90, code=3, col="blue4")
arrows(nov2$Day.x, nov2$mean-nov2$sd, nov2$Day.x, nov2$mean+nov2$sd, length=0.05, angle=90, code=3, col="orange")
arrows(fmt2$Day.x, fmt2$mean-nov2$sd, nov2$Day.x, nov2$mean+nov2$sd, length=0.05, angle=90, code=3, col="chartreuse4")
points(vanco2$Day.x, vanco2$mean, pch=21, col="black", bg="blue4", xlim=c(-8,18), xaxt='n')
points(nov2$Day.x, nov2$mean, pch=21, col="black", bg="orange", xlim=c(-8,18), xaxt='n')
points(fmt2$Day.x, fmt2$mean, pch=21, col="black", bg="chartreuse4", xlim=c(-8,18), xaxt='n')
points(abx2$Day.x, abx2$mean, pch=21, col="black", bg="violetred1", xlim=c(-8,18), xaxt='n')		#replot the points so that they are over the lines
legend("bottomright", c("abx only", "Vanco (relapse)", "No vanco", "FMT"), pch=21, col=c("black"), pt.bg=c("violetred1", "blue4", "orange", "chartreuse4"), cex=1)

plot(fmt3$Day.x, fmt3$mean, 
	pch=21, xlim=c(-8,45), ylim=c(80, 120), xaxt='n', col="black", bg="chartreuse4", cex=1, 
	main="Relapse model 1 (exp 3)", ylab="mean weight change from d0 (%)", xlab="Day")
axis(1, at=fmt3$Day.x, labels=fmt3$Day.x)
lines(fmt3$Day.x, fmt3$mean, col="chartreuse4")
lines(nof3$Day.x, nof3$mean, col="blue4")
arrows(fmt3$Day.x, fmt3$mean-fmt3$sd, fmt3$Day.x, fmt3$mean+fmt3$sd, length=0.05, angle=90, code=3, col="chartreuse4")
arrows(nof3$Day.x, nof3$mean-nof3$sd, nof3$Day.x, nof3$mean+nof3$sd, length=0.05, angle=90, code=3, col="blue4")
points(fmt3$Day.x, fmt3$mean, pch=21, col="black", bg="chartreuse4", xlim=c(-8,45), xaxt='n')
points(nof3$Day.x, nof3$mean, pch=21, col="black", bg="blue4", xlim=c(-8,45), xaxt='n')
legend("bottomright", c("FMT", "no FMT"), pch=21, col=c("black"), pt.bg=c("chartreuse4", "blue4"), cex=1)

```