# Different thetayc calculations for longitudinal graphs
### Anna Seekatz
### 12.15.14

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
+ mouseRM_shared.dist.txt: has pairwise distances (originally was mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.summary)
+ mouse.var.txt: this is all of the metadata by sample (note: must be on same order as the dist.txt file)

> Create a master file from which you can annotate data:
```{r}
mouse.var<-read.table(file="mouse.var.txt", header=TRUE)
mdist<-read.table(file="mouseRM_shared.dist.txt", header=TRUE)
head(mouse.var)															#checking what the files look like (head displays the first 6 rows of a file)
head(mdist)

#add sample information for 'group1':
m<-merge(mouse.var, mdist, by.x=c("sampleID"), by.y=c("group1")) #this merges the data based on the sampleID/group1 match
m1<-m[,c("sampleID", "group2", "cage", "day", "mouse", "exp", "treatment", "treatment2", "group_time", "sharedsobs", "braycurtis", "spearman", "thetayc")]  #this picks the specific columns you want to incorporate
m1<-rename(m1, c("sampleID"="sample1", "cage"="cage1", "day"="day1", "mouse"="mouse1", "exp"="exp1", "treatment"="treatment1", "treatment2"="treat1", "group_time"="group_time1")) #this renames the columns so you know which metadata matches group1 (or, renamed, sample1)

#add sample information for 2nd sample: ('group2'):
m2<-merge(mouse.var, m1, by.x=c("sampleID"), by.y=c("group2"))
m3<-m2[,c("sample1", "sampleID", "cage1", "day1", "mouse1", "exp1", "treatment1", "treat1", "group_time1", "cage", "day", "mouse", "exp", "treatment", "treatment2", "group_time", "sharedsobs", "braycurtis", "spearman", "thetayc")]
mshare<-rename(m3, c("sampleID"="sample2", "cage"="cage2", "day"="day2", "mouse"="mouse2", "exp"="exp2", "treatment"="treatment2", "treatment2"="treat2", "group_time"="group_time2")) #this renames the columns so you know which metadata matches group1 (or, renamed, sample1)
mshare[1:20, 1:20] #check to see if it looks correct
#write.table(mshare, file="mouse.var_shared.dist.txt", quote=FALSE, sep="\t")
head(mshare)
```
The file you created (mshare) now has all of the information necessary to parse out specific comparisons between the two samples compared pairwise

> Parsing out your information 

> **Example1--comparing changes within treatment groups over time (day-to-day comparison):**
```{r}
#first graph, between days, over time, within each treatment group:

#Between FMT and pre (day-7):
don1<-mshare[(mshare$day1==-8) & (mshare$day2==-7), ]
don2<-mshare[(mshare$day1==-7) & (mshare$day2==-8), ]					#empty--but for future reference
a1<-rbind(don1, don2)
a1$DAY<-c(-1)
a1$GROUP<-c("Donor")

#FMT mice:
#check which day comparisons need to be made:
fmt<-subset(mshare, (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("FMT")))
levels(as.factor(fmt$day2))
 	#[1] "-7" "1"  "4"  "10" "12" "16" "19" "22" "27" "33" "36" "37" "42"
levels(as.factor(fmt$day1))
 	#[1] "-7" "1"  "4"  "10" "12" "16" "19" "22" "27" "33" "36" "37" "42"
#they match--good

	#day -7 to 1:
sub1<-subset(fmt, (subset=day1 %in% c(1)) & (subset=day2 %in% c(-7)))	#there may also be comparisons between the two on opposite sides
sub2<-subset(fmt, (subset=day1 %in% c(-7)) & (subset=day2 %in% c(1)))
d1<-rbind(sub1, sub2)													#this joins the two comparisons together
d1$DAY<- c(1)
	#day 1 to 4:
sub1<-subset(fmt, (subset=day1 %in% c(1)) & (subset=day2 %in% c(4)))
sub2<-subset(fmt, (subset=day1 %in% c(4)) & (subset=day2 %in% c(1)))
d2<-rbind(sub1, sub2)
d2$DAY<- c(4)
	#day 4 to 10:
sub1<-subset(fmt, (subset=day1 %in% c(10)) & (subset=day2 %in% c(4)))
sub2<-subset(fmt, (subset=day1 %in% c(4)) & (subset=day2 %in% c(10)))
d3<-rbind(sub1, sub2)
d3$DAY<- c(10)
	#day 10 to 12:
sub1<-subset(fmt, (subset=day1 %in% c(10)) & (subset=day2 %in% c(12)))
sub2<-subset(fmt, (subset=day1 %in% c(12)) & (subset=day2 %in% c(10)))
d4<-rbind(sub1, sub2)
d4$DAY<- c(12)
	#day 12 to 16:
sub1<-subset(fmt, (subset=day1 %in% c(16)) & (subset=day2 %in% c(12)))
sub2<-subset(fmt, (subset=day1 %in% c(12)) & (subset=day2 %in% c(16)))
d5<-rbind(sub1, sub2)
d5$DAY<- c(16)	
	#day 16 to 19:
sub1<-subset(fmt, (subset=day1 %in% c(16)) & (subset=day2 %in% c(19)))
sub2<-subset(fmt, (subset=day1 %in% c(19)) & (subset=day2 %in% c(16)))
d6<-rbind(sub1, sub2)
d6$DAY<- c(19)
	#day 19 to 22:
sub1<-subset(fmt, (subset=day1 %in% c(22)) & (subset=day2 %in% c(19)))
sub2<-subset(fmt, (subset=day1 %in% c(19)) & (subset=day2 %in% c(22)))
d7<-rbind(sub1, sub2)
d7$DAY<- c(22)
	#day 22 to 27:
sub1<-subset(fmt, (subset=day1 %in% c(22)) & (subset=day2 %in% c(27)))
sub2<-subset(fmt, (subset=day1 %in% c(27)) & (subset=day2 %in% c(22)))
d8<-rbind(sub1, sub2)
d8$DAY<- c(22)
	#day 27 to 33:
sub1<-subset(fmt, (subset=day1 %in% c(33)) & (subset=day2 %in% c(27)))
sub2<-subset(fmt, (subset=day1 %in% c(27)) & (subset=day2 %in% c(33)))
d9<-rbind(sub1, sub2)
d9$DAY<- c(33)
	#day 33 to 36:
sub1<-subset(fmt, (subset=day1 %in% c(33)) & (subset=day2 %in% c(36)))
sub2<-subset(fmt, (subset=day1 %in% c(36)) & (subset=day2 %in% c(33)))
d10<-rbind(sub1, sub2)
d10$DAY<- c(36)
	#day 36 to 37:
sub1<-subset(fmt, (subset=day1 %in% c(37)) & (subset=day2 %in% c(36)))
sub2<-subset(fmt, (subset=day1 %in% c(36)) & (subset=day2 %in% c(37)))
d11<-rbind(sub1, sub2)
d11$DAY<- c(37)
	#day 37 to 42:
sub1<-subset(fmt, (subset=day1 %in% c(37)) & (subset=day2 %in% c(42)))
sub2<-subset(fmt, (subset=day1 %in% c(42)) & (subset=day2 %in% c(37)))
d12<-rbind(sub1, sub2)
d12$DAY<- c(42)
	#finally, combine these together:
fmt.all<-rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12)
fmt.all$GROUP<-c("FMT")

#no_FMT mice:
#check which day comparisons need to be made:
no.fmt<-subset(mshare, (subset=treat1 %in% c("no_FMT")) & (subset=treat2 %in% c("no_FMT")))

#day -7 to 1:
sub1<-subset(no.fmt, (subset=day1 %in% c(1)) & (subset=day2 %in% c(-7)))
sub2<-subset(no.fmt, (subset=day1 %in% c(-7)) & (subset=day2 %in% c(1)))
n1<-rbind(sub1, sub2)
n1$DAY<- c(1)
	#day 1 to 4:
sub1<-subset(no.fmt, (subset=day1 %in% c(1)) & (subset=day2 %in% c(4)))
sub2<-subset(no.fmt, (subset=day1 %in% c(4)) & (subset=day2 %in% c(1)))
n2<-rbind(sub1, sub2)
n2$DAY<- c(4)
	#day 4 to 10:
sub1<-subset(no.fmt, (subset=day1 %in% c(10)) & (subset=day2 %in% c(4)))
sub2<-subset(no.fmt, (subset=day1 %in% c(4)) & (subset=day2 %in% c(10)))
n3<-rbind(sub1, sub2)
n3$DAY<- c(10)
	#day 10 to 12:
sub1<-subset(no.fmt, (subset=day1 %in% c(10)) & (subset=day2 %in% c(12)))
sub2<-subset(no.fmt, (subset=day1 %in% c(12)) & (subset=day2 %in% c(10)))
n4<-rbind(sub1, sub2)
n4$DAY<- c(12)
	#day 12 to 16:
sub1<-subset(no.fmt, (subset=day1 %in% c(16)) & (subset=day2 %in% c(12)))
sub2<-subset(no.fmt, (subset=day1 %in% c(12)) & (subset=day2 %in% c(16)))
n5<-rbind(sub1, sub2)
n5$DAY<- c(16)	
	#day 16 to 19:
sub1<-subset(no.fmt, (subset=day1 %in% c(16)) & (subset=day2 %in% c(19)))
sub2<-subset(no.fmt, (subset=day1 %in% c(19)) & (subset=day2 %in% c(16)))
n6<-rbind(sub1, sub2)
n6$DAY<- c(19)
	#day 19 to 22:
sub1<-subset(no.fmt, (subset=day1 %in% c(22)) & (subset=day2 %in% c(19)))
sub2<-subset(no.fmt, (subset=day1 %in% c(19)) & (subset=day2 %in% c(22)))
n7<-rbind(sub1, sub2)
n7$DAY<- c(22)
	#day 22 to 27:
sub1<-subset(no.fmt, (subset=day1 %in% c(22)) & (subset=day2 %in% c(27)))
sub2<-subset(no.fmt, (subset=day1 %in% c(27)) & (subset=day2 %in% c(22)))
n8<-rbind(sub1, sub2)
n8$DAY<- c(22)
	#day 27 to 33:
sub1<-subset(no.fmt, (subset=day1 %in% c(33)) & (subset=day2 %in% c(27)))
sub2<-subset(no.fmt, (subset=day1 %in% c(27)) & (subset=day2 %in% c(33)))
n9<-rbind(sub1, sub2)
n9$DAY<- c(33)
	#day 33 to 36:
sub1<-subset(no.fmt, (subset=day1 %in% c(33)) & (subset=day2 %in% c(36)))
sub2<-subset(no.fmt, (subset=day1 %in% c(36)) & (subset=day2 %in% c(33)))
n10<-rbind(sub1, sub2)
n10$DAY<- c(36)
	#day 36 to 37:
sub1<-subset(no.fmt, (subset=day1 %in% c(37)) & (subset=day2 %in% c(36)))
sub2<-subset(no.fmt, (subset=day1 %in% c(36)) & (subset=day2 %in% c(37)))
n11<-rbind(sub1, sub2)
n11$DAY<- c(37)
	#day 37 to 42:
sub1<-subset(no.fmt, (subset=day1 %in% c(37)) & (subset=day2 %in% c(42)))
sub2<-subset(no.fmt, (subset=day1 %in% c(42)) & (subset=day2 %in% c(37)))
n12<-rbind(sub1, sub2)
n12$DAY<- c(42)
	#finally, combine these together:
no.fmt.all<-rbind(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12)
no.fmt.all$GROUP<-c("no_FMT")
```
You have now created a sheet that only includes the distance measured ONLY between day 1 and -7, within each of the groups.

> Creating summary statistics for graphing:
```{r}

#combine fmt, no.fmt datasets:
theta.overtime<-rbind(fmt.all, no.fmt.all)

#calculate summary of theta.overtime:
tdata <- ddply(theta.overtime, c("DAY", "GROUP"), summarise,
               N    = length(thetayc),
               mean = mean(thetayc),
               sd   = sd(thetayc),
               se   = sd / sqrt(N) )
head(tdata)									#example of what your data looks like
               
#option 2 for getting these results:
#must definte summarySE function before hand
dfc <- summarySE(theta.overtime, measurevar="thetayc", groupvars=c("DAY","GROUP"))
head(dfc)
```

> Graphing the data:

There are many ways to graph the data. The following code gives some instruction on how to use:
+ ggplot2 package (although I only do one plot)
+ interaction plots (a quick way to view your data)
+ plotting overtime (no packages required!)

```{r}
#graph the data:
ggplot(dfc, aes(x=DAY, y=thetayc, colour=GROUP)) + 
    geom_errorbar(aes(ymin=thetayc-se, ymax=thetayc+se), width=.1) +
    geom_line() +
    geom_point()
#looks ok, but I don't necessarily like the style


#interaction plot:
meantheta = aggregate(theta.overtime$thetayc, list(theta.overtime$GROUP, theta.overtime$DAY), 
    mean)
interaction.plot(meantheta$Group.2, meantheta$Group.1, meantheta$x, type = "b", col = c(1:2), pch = c(18, 24))
 	#this worked, and looks nice   
	#attempting to add standard error bars:
sd.theta = aggregate(theta.overtime$thetayc, list(theta.overtime$GROUP, theta.overtime$DAY), sd)
	#OR
#sd.up = meantheta$x+sd.theta$x
#sd.down = meantheta$x-sd.theta$x
	#then:
#p1<-interaction.plot(meantheta$Group.2, meantheta$Group.1, meantheta$x, type = "b", col = c(1:5), pch = 19,)
#arrows(p1,sd.down,p1,sd.up,code=3,length=0.2,angle=90,col='red')
	#this did not work--however, the interaction plot might be a nice quick way to look at things over time (just use the first 2 lines at the top)
	

#another way of doing this (using the means, sd defined earlier):
#group1:
sub.td<-subset(tdata, (subset=GROUP %in% c("FMT")))				#subset data by group
avg<-as.numeric(sub.td$mean)
sdev<-as.numeric(sub.td$sd)
#group2:
sub2.td<-subset(tdata, (subset=GROUP %in% c("no_FMT")))
avg2<-as.numeric(sub2.td$mean)
sdev2<-as.numeric(sub2.td$sd)

x<-sub.td$DAY														
	#this defines where your graph is (both mean and sd are of length n, or how many samples you have (if you want categorical, define x<1:11, or the number of observations you are plotting)
plot(x, avg,
    ylim=range(c(avg2-sdev2, avg2+sdev2)),
    pch=19, xlab="Day sampled", ylab="Mean distance to prior sampling",
    main="change over time", xaxt='n')
arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)	#adds sd to the graph
axis(1, at=x, labels=x)											#adds labels to the graph
lines(x, avg)													#adds lines to the graph

#now add the other data group (no_FMT):
points(x, avg2, col="blue", pch=19)
lines(x, avg2, col="blue", pch=19)
arrows(x, avg2-sdev2, x, avg2+sdev2, length=0.05, angle=90, code=3, col="blue")
legend("bottomleft", c("FMT", "no FMT"), col=c("black", "blue"), pch=19, cex=0.8)
```

> **Example2--comparing changes within treatment groups over time compared to PRE status:**

This has all of the code required to create a graph 

```{r}
mshare<-read.table(file="mouse.var_shared.dist.txt", header=TRUE)

#check which day comparisons need to be made:
fmt<-subset(mshare, (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("FMT")))
no.fmt<-subset(mshare, (subset=treat1 %in% c("no_FMT")) & (subset=treat2 %in% c("no_FMT")))

	#day 1:
sub1<-subset(mshare, (subset=day1 %in% c(1)) & (subset=day2 %in% c(-7, -8)))
sub2<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(1)))
f1<-rbind(sub1, sub2)
f1$DAY<- c(1)
	#day 4:
sub1<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(4)))
sub2<-subset(mshare, (subset=day1 %in% c(4)) & (subset=day2 %in% c(-7, -8)))
f2<-rbind(sub1, sub2)
f2$DAY<- c(4)
	#day 10:
sub1<-subset(mshare, (subset=day1 %in% c(10)) & (subset=day2 %in% c(-7, -8)))
sub2<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(10)))
f3<-rbind(sub1, sub2)
f3$DAY<- c(10)
	#day 12:
sub1<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(12)))
sub2<-subset(mshare, (subset=day1 %in% c(12)) & (subset=day2 %in% c(-7, -8)))
f4<-rbind(sub1, sub2)
f4$DAY<- c(12)
	#day 16:
sub1<-subset(mshare, (subset=day1 %in% c(16)) & (subset=day2 %in% c(-7, -8)))
sub2<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(16)))
f5<-rbind(sub1, sub2)
f5$DAY<- c(16)	
	#day 19:
sub1<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(19)))
sub2<-subset(mshare, (subset=day1 %in% c(19)) & (subset=day2 %in% c(-7, -8)))
f6<-rbind(sub1, sub2)
f6$DAY<- c(19)
	#day 22:
sub1<-subset(mshare, (subset=day1 %in% c(22)) & (subset=day2 %in% c(-7, -8)))
sub2<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(22)))
f7<-rbind(sub1, sub2)
f7$DAY<- c(22)
	#day 27:
sub1<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(27)))
sub2<-subset(mshare, (subset=day1 %in% c(27)) & (subset=day2 %in% c(-7, -8)))
f8<-rbind(sub1, sub2)
f8$DAY<- c(27)
	#day 33:
sub1<-subset(mshare, (subset=day1 %in% c(33)) & (subset=day2 %in% c(-7, -8)))
sub2<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(33)))
f9<-rbind(sub1, sub2)
f9$DAY<- c(33)
	#day 36:
sub1<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(36)))
sub2<-subset(mshare, (subset=day1 %in% c(36)) & (subset=day2 %in% c(-7, -8)))
f10<-rbind(sub1, sub2)
f10$DAY<- c(36)
	#day 37:
sub1<-subset(mshare, (subset=day1 %in% c(37)) & (subset=day2 %in% c(-7, -8)))
sub2<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(37)))
f11<-rbind(sub1, sub2)
f11$DAY<- c(37)
	#day 42:
sub1<-subset(mshare, (subset=day1 %in% c(-7, -8)) & (subset=day2 %in% c(42)))
sub2<-subset(mshare, (subset=day1 %in% c(42)) & (subset=day2 %in% c(-7, -8)))
f12<-rbind(sub1, sub2)
f12$DAY<- c(42)
	
	#finally, combine these together:
pre.comp.all<-rbind(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12)
	#you can compare these across cages:
levels(as.factor(fmt$cage1))
	#"982"  "1007" "1008" "1009" "1010" "1015"
levels(as.factor(no.fmt$cage1))
	#[1] "1006"
levels(as.factor(pre.comp.all$cage1))
 	# [1] "972"  "982"  "986"  "1001" "1002" "1006" "1007" "1008" "1009" "1010" "1015"
#use these cage assignments to group the mice into the groups to be compared (FMT or no_FMT against pre/abx time points)

#to make the FMT group:	
presub1<-subset(pre.comp.all, (subset=cage1 %in% c(982, 1007, 1008, 1009, 1010, 1015)) & (subset=cage2 %in% c(1001, 1002, 972, 986)))
presub2<-subset(pre.comp.all, (subset=cage2 %in% c(982, 1007, 1008, 1009, 1010, 1015)) & (subset=cage1 %in% c(1001, 1002, 972, 986)))
pre.fmt<-rbind(presub1, presub2)
pre.fmt$GROUP<-c("FMT")

#to make the no_FMT group:
presub1<-subset(pre.comp.all, (subset=cage1 %in% c(1006)) & (subset=cage2 %in% c(1001, 1002, 972, 986)))
presub2<-subset(pre.comp.all, (subset=cage2 %in% c(1006)) & (subset=cage1 %in% c(1001, 1002, 972, 986)))
pre.nofmt<-rbind(presub1, presub2)
pre.nofmt$GROUP<-c("no_FMT")

#combine them together to make a large group:
pre.all<-rbind(pre.fmt, pre.nofmt)

#write.table(pre.all, file="mouseRM_dist.to.preabx.txt", quote=FALSE, sep="\t", col.names=NA)

#make summary files:

#calculate summary of pre.all:
pdata <- ddply(pre.all, c("DAY", "GROUP"), summarise,
               N    = length(thetayc),
               mean = mean(thetayc),
               sd   = sd(thetayc),
               se   = sd / sqrt(N) )

#subset data by group:
#group1:
sub.pd<-subset(pdata, (subset=GROUP %in% c("FMT")))				
p.ave<-as.numeric(sub.pd$mean)
p.sd<-as.numeric(sub.pd$sd)
#group2:
sub2.pd<-subset(pdata, (subset=GROUP %in% c("no_FMT")))
p.ave2<-as.numeric(sub2.pd$mean)
p.sd2<-as.numeric(sub2.pd$sd)

x<-sub.pd$DAY														
plot(x, p.ave,
    ylim=range(c(p.ave-p.sd, p.ave+p.sd)),
    pch=19, xlab="Day sampled", ylab="Mean distance to pre-abx community",
    main="change over time (compared to pre-abx status)", xaxt='n')
arrows(x, p.ave-p.sd, x, p.ave+p.sd, length=0.05, angle=90, code=3)	#adds sd to the graph
axis(1, at=x, labels=x)												#adds labels to the graph
lines(x, p.ave)														#adds lines to the graph

#now add the other data group (no_FMT):
points(sub2.pd$DAY, p.ave2, col="blue", pch=19)						#note: since this group has no day 4 comparison, you use the actual variable (sub2.pd) instead of x
lines(sub2.pd$DAY, p.ave2, col="blue", pch=19)
arrows(sub2.pd$DAY, p.ave2-p.sd2, sub2.pd$DAY, p.ave2+p.sd2, length=0.05, angle=90, code=3, col="blue")
legend("bottomleft", c("FMT", "no FMT"), col=c("black", "blue"), pch=19, cex=0.8)
```

> **Example3--comparing treatment (FMT vs. no_FMT) overtime:**

This has all of the code required to create a graph 

```{r}
##example 3--comparing differences directly between FMT and no_FMT over time:

mshare<-read.table(file="mouse.var_shared.dist.txt", header=TRUE)

#for this comparison, you will be comparing FMT to no_FMT on each day
#thus, you can just choose all of the FMT to no_FMT comparisons on a specific day

#check the days by 

#day -7:
sub1<-subset(mshare, (subset=day1 %in% c(-7)) & (subset=day2 %in% c(-7)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(-7)) & (subset=day2 %in% c(-7)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a1<-rbind(sub1, sub2)
a1$DAY<- c(-7)
#day 1:
sub1<-subset(mshare, (subset=day1 %in% c(1)) & (subset=day2 %in% c(1)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(1)) & (subset=day2 %in% c(1)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a2<-rbind(sub1, sub2)
a2$DAY<- c(1)
#day 4:
sub1<-subset(mshare, (subset=day1 %in% c(4)) & (subset=day2 %in% c(4)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(4)) & (subset=day2 %in% c(4)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a3<-rbind(sub1, sub2)
a3$DAY<- c(4)
#day 10:
sub1<-subset(mshare, (subset=day1 %in% c(10)) & (subset=day2 %in% c(10)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(10)) & (subset=day2 %in% c(10)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a4<-rbind(sub1, sub2)
a4$DAY<- c(10)
#day 12:
sub1<-subset(mshare, (subset=day1 %in% c(12)) & (subset=day2 %in% c(12)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(12)) & (subset=day2 %in% c(12)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a5<-rbind(sub1, sub2)
a5$DAY<- c(12)
#day 16:
sub1<-subset(mshare, (subset=day1 %in% c(16)) & (subset=day2 %in% c(16)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(16)) & (subset=day2 %in% c(16)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a6<-rbind(sub1, sub2)
a6$DAY<- c(16)
#day 19:
sub1<-subset(mshare, (subset=day1 %in% c(19)) & (subset=day2 %in% c(19)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(19)) & (subset=day2 %in% c(19)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a7<-rbind(sub1, sub2)
a7$DAY<- c(19)
#day 22:
sub1<-subset(mshare, (subset=day1 %in% c(22)) & (subset=day2 %in% c(22)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(22)) & (subset=day2 %in% c(22)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a8<-rbind(sub1, sub2)
a8$DAY<- c(22)
#day 27:
sub1<-subset(mshare, (subset=day1 %in% c(27)) & (subset=day2 %in% c(27)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(27)) & (subset=day2 %in% c(27)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a9<-rbind(sub1, sub2)
a9$DAY<- c(27)
#day 33:
sub1<-subset(mshare, (subset=day1 %in% c(33)) & (subset=day2 %in% c(33)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(33)) & (subset=day2 %in% c(33)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a10<-rbind(sub1, sub2)
a10$DAY<- c(33)
#day 36:
sub1<-subset(mshare, (subset=day1 %in% c(36)) & (subset=day2 %in% c(36)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(36)) & (subset=day2 %in% c(36)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a11<-rbind(sub1, sub2)
a11$DAY<- c(36)
#day 37:
sub1<-subset(mshare, (subset=day1 %in% c(37)) & (subset=day2 %in% c(37)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(37)) & (subset=day2 %in% c(37)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a12<-rbind(sub1, sub2)
a12$DAY<- c(37)
#day 42:
sub1<-subset(mshare, (subset=day1 %in% c(42)) & (subset=day2 %in% c(42)) & (subset=treat1 %in% c("FMT")) & (subset=treat2 %in% c("no_FMT")),)
sub2<-subset(mshare, (subset=day1 %in% c(42)) & (subset=day2 %in% c(42)) & (subset=treat2 %in% c("FMT")) & (subset=treat1 %in% c("no_FMT")),)
a13<-rbind(sub1, sub2)
a13$DAY<- c(42)

#merge the datasets together:
f.nof.comp<-rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13)

#calculate summary of f.nof.comp:
fdata <- ddply(f.nof.comp, c("DAY"), summarise,
               N    = length(thetayc),
               mean = mean(thetayc),
               sd   = sd(thetayc),
               se   = sd / sqrt(N) )

f.ave<-as.numeric(fdata$mean)
f.sd<-as.numeric(fdata$sd)

x<-fdata$DAY														
plot(x, f.ave,
    ylim=range(c(f.ave-f.sd, f.ave+f.sd)),
    pch=19, xlab="Day sampled", ylab="Mean distance between FMT and no FMT groups",
    main="change over time (FMT vs. no FMT)", xaxt='n')
arrows(x, f.ave-f.sd, x, f.ave+f.sd, length=0.05, angle=90, code=3)	#adds sd to the graph
axis(1, at=x, labels=x)												#adds labels to the graph
#lines(x, f.ave)
```

> **all graphs together:**

```{r}

##adding all graphs together:

#pdf("mnouseRM_tyc_time.pdf") #if you want this printed in a pdf

par(mfrow=c(3,1))
x<-sub.td$DAY														
plot(x, avg,
    ylim=range(c(avg2-sdev2, avg2+sdev2)),
    pch=19, xlab="Day sampled", ylab="Mean distance to prior sampling",
    main="change over time (consecutive day comparison)", xaxt='n')
arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)	#adds sd to the graph
axis(1, at=x, labels=x)											#adds labels to the graph
lines(x, avg)													#adds lines to the graph
points(x, avg2, col="blue", pch=19)
lines(x, avg2, col="blue", pch=19)
arrows(x, avg2-sdev2, x, avg2+sdev2, length=0.05, angle=90, code=3, col="blue")
legend("bottomleft", c("FMT", "no FMT"), col=c("black", "blue"), pch=19, cex=0.8)

t<-sub.pd$DAY														
plot(t, p.ave,
    ylim=range(c(p.ave-p.sd, p.ave+p.sd)),
    pch=19, xlab="Day sampled", ylab="Mean distance to pre-abx community",
    main="change over time (compared to pre-abx status)", xaxt='n')
arrows(t, p.ave-p.sd, t, p.ave+p.sd, length=0.05, angle=90, code=3)	#adds sd to the graph
axis(1, at=t, labels=t)												#adds labels to the graph
lines(t, p.ave)														#adds lines to the graph
points(sub2.pd$DAY, p.ave2, col="blue", pch=19)						#note: since this group has no day 4 comparison, you use the actual variable (sub2.pd) instead of x
lines(sub2.pd$DAY, p.ave2, col="blue", pch=19)
arrows(sub2.pd$DAY, p.ave2-p.sd2, sub2.pd$DAY, p.ave2+p.sd2, length=0.05, angle=90, code=3, col="blue")
legend("bottomleft", c("FMT", "no FMT"), col=c("black", "blue"), pch=19, cex=0.8)

f<-fdata$DAY														
plot(f, f.ave,
    ylim=range(c(f.ave-f.sd, f.ave+f.sd)),
    pch=19, xlab="Day sampled", ylab="Mean distance between FMT and no FMT groups",
    main="change over time (FMT vs. no FMT)", xaxt='n')
arrows(f, f.ave-f.sd, f, f.ave+f.sd, length=0.05, angle=90, code=3)	#adds sd to the graph
axis(1, at=f, labels=f)												#adds labels to the graph
#lines(f, f.ave)

#dev.off()
```

> Stats for the thetayc differences:

```{r}

## stats for change over time (Fig. 5B):
# this only shows stats for the distance between pre-abx/donor vs FMT or no_FMT groups:

# previously, had created a file called 'pre.all' that has thetayc pairwise distance ONLY for FMT/pre-abx and no_FMT/pre-abx comparisons (B)
preabx.comp<-read.table(file="mouseRM_dist.to.preabx.txt", header=TRUE)

# split by DAY:
pre.split<-split(preabx.comp, preabx.comp$DAY)
pre.d1<-pre.split$'1'
pre.d10<-pre.split$'10'
pre.d12<-pre.split$'12'
pre.d16<-pre.split$'16'
pre.d19<-pre.split$'19'
pre.d22<-pre.split$'22'
pre.d27<-pre.split$'27'
pre.d33<-pre.split$'33'
pre.d37<-pre.split$'37'
pre.d36<-pre.split$'36'
pre.d42<-pre.split$'42'

# do stats:
wilcox.test(thetayc~GROUP, data=pre.d10)
wilcox.test(thetayc~GROUP, data=pre.d12)
wilcox.test(thetayc~GROUP, data=pre.d16)
wilcox.test(thetayc~GROUP, data=pre.d19)
wilcox.test(thetayc~GROUP, data=pre.d22, paired=FALSE)
wilcox.test(thetayc~GROUP, data=pre.d27, paired=FALSE)
wilcox.test(thetayc~GROUP, data=pre.d33, paired=FALSE)
wilcox.test(thetayc~GROUP, data=pre.d36, paired=FALSE)
wilcox.test(thetayc~GROUP, data=pre.d37, paired=FALSE)
wilcox.test(thetayc~GROUP, data=pre.d42, paired=FALSE)	
wilcox.test(thetayc~GROUP, data=pre.d1, paired=FALSE)

```

# code used for figure 5B,C:


```{r}
#distance to pre-abx community:
x<-sub.pd$DAY														
plot(x, p.ave,
    ylim=range(c(p.ave-p.sd, p.ave+p.sd)), 
    pch=19, xlab="Day sampled", ylab="Mean distance to pre-abx community", cex=1.2,
    xaxt='n', col="chartreuse4", xlim=c(-8,45))
arrows(x, p.ave-p.sd, x, p.ave+p.sd, length=0.05, angle=90, code=3, col="chartreuse4")	#adds sd to the graph
axis(1, at=x, labels=x)												#adds labels to the graph
lines(x, p.ave, col="chartreuse4")														#adds lines to the graph
points(sub2.pd$DAY, p.ave2, col="blue4", pch=19, cex=1.2)						#note: since this group has no day 4 comparison, you use the actual variable (sub2.pd) instead of x
lines(sub2.pd$DAY, p.ave2, col="blue4", pch=19)
arrows(sub2.pd$DAY, p.ave2-p.sd2, sub2.pd$DAY, p.ave2+p.sd2, length=0.05, angle=90, code=3, col="blue4")
legend("bottomleft", c("FMT", "no FMT"), col=c("chartreuse4", "blue4"), pch=19, cex=1)													#adds lines to the graph

#difference/distance between FMT/no_FMT groups:
f<-fdata$DAY														
plot(f, f.ave,
    ylim=range(c(f.ave-f.sd, f.ave+f.sd)),
    pch=19, xlab="Day sampled", ylab="Mean distance between FMT and no FMT groups",
    xaxt='n', cex=1.2, xlim=c(-8,45))
arrows(f, f.ave-f.sd, f, f.ave+f.sd, length=0.05, angle=90, code=3)	#adds sd to the graph
axis(1, at=f, labels=f)

```