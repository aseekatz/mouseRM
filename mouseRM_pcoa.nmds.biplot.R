# PCOA or NMDS oordination graphs, with biplot coordinates of significant OTUs
### Anna Seekatz
### 1.20.15

> load all of the necessary packages

```{r, error=FALSE, message=FALSE, warning=FALSE}

library(shape)

```

Change into your working directory: 

```{r}
# this is only applicable to me, so use your own path...
setwd("/Users/aseekatz/Desktop/umich/projects/Relapse_model2/mouseRM_analysis/mouseRM_non.mb.data")
list.files(path="/Users/aseekatz/Desktop/umich/projects/Relapse_model2/mouseRM_analysis/mouseRM_non.mb.data")

```

The files you will need for this include:
	+ mouseRM_div.thetayc.txt (see below on how to create this)
		#OR# if you do not have this file:
	+ group summary file from mothur: mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.groups.summary
	+ metadata file created previously: mouse.var.txt
	+ pcoa axes added from mothur file: mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.pcoa.axes
	+ nmds axes added from mothur file: mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.nmds.axes

To add biplot coordinates, you will need the following files generated from mothur's corr.axes function:
	+ mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.nmds.spearman.corr.axes
	+ mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.pcoa.spearman.corr.axes


```{r}

##note: I had previously created this file for alpha diversity measures (Fig. 5A), but to recreate this file, you can follow this:
#to create the 'mouseRM_div.thetayc.txt' file:

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

# you know have a file to proceed

```

To graph PCOA oordination of your data:

```{r}

#pcoa:

#in R:
pcoa<-read.table(file="mouseRM_div.thetayc.txt", header=TRUE)
par(mfrow=c(1,3))
#class(pcoa$group_time)  #if integer data type, must change to character
#pcoa$group_time<-as.character(pcoa$group_time)
plot(pcoa$pcoa_axis2~pcoa$pcoa_axis1, pch=19, col=time.col(pcoa$group_time))
plot(pcoa$pcoa_axis2~pcoa$pcoa_axis3, pch=19, col=time.col(pcoa$group_time))
plot(pcoa$pcoa_axis3~pcoa$pcoa_axis1, pch=19, col=time.col(pcoa$group_time))
legend("bottomright",legend=c("pre-treatment", "ABX only: during cef", "ABX only: post-cef", "ABX only: during vanco", "ABX only: post-vanco", "CDI: pre-vanco", "CDI: post-vanco treatment", "CDI: post-IP clinda", "CDI: post-FMT", "CDI: post-FMT, post-clinda", "donor FMT"), col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"), cex=0.8, pch=19)

#1 graph at a time:
plot(pcoa$pcoa_axis3~pcoa$pcoa_axis1, pch=19, col=time.col(pcoa$group_time) )
legend("bottomright",legend=c("pre-treatment", "ABX only: during cef", "ABX only: post-cef", "ABX only: during vanco", "ABX only: post-vanco", "CDI: pre-vanco", "CDI: post-vanco treatment", "CDI: post-IP clinda", "CDI: post-FMT", "CDI: post-FMT, post-clinda", "donor FMT"), col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"), cex=0.8, pch=19)
#note: you can add the % variance explained per axis by looking at the .loadings file produced alongside your .axes file in mothur

##color palettes:
#Color palette by treatment group over time:
time.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "pre" ) {
colorvec[i] = "green4"
}
if (n[i] == "fmt" ) {
colorvec[i] = "green3"
}
if( n[i] == "cef") {
colorvec[i] = "hotpink"
}
if( n[i] == "post_cef") {
colorvec[i] = "magenta1"
}
if ( n[i] == "vanco" ) {
colorvec[i] = "magenta4"
}
if (n[i] == "post_vanco" ) {
colorvec[i] = "red4"
}
if( n[i] == "post_CDI") {
colorvec[i] = "cyan"
}
if( n[i] == "post_CDI_vanco") {
colorvec[i] = "deepskyblue1"
}
if (n[i] == "post_CDI_vanco_clinda" ) {
colorvec[i] = "blue1"
}
if( n[i] == "post_FMT") {
colorvec[i] = "olivedrab2"
}
if( n[i] == "post_FMT_clinda") {
colorvec[i] = "gold1"
}
}
c(colorvec)
}

#some diversity indeces:
pcoa<-read.table(file="mouseRM1_pcoa.txt", header=TRUE)
pcoa$group<-factor(pcoa$group_time, levels=c("pre", "cef", "post_cef", "vanco", "post_vanco", "post_CDI", "post_CDI_vanco", "post_CDI_vanco_clinda", "post_FMT", "post_FMT_clinda", "fmt"))
	#this orders the axis
plot(pcoa$shannon~pcoa$group_time, pch=19, las=2, ylab="Shannon diversity index", col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"))

--
#pcoa overtime
#this was not published, but I wanted to get an idea of what the first axis looks like over time

#used file: mouseRM1_pcoa.txt
	#note: changed day within FMT to be -6

pcoa<-read.table(file="mouseRM1_pcoa.txt", header=TRUE)
#par(mfrow=c(1,3))
#class(pcoa$group)  #if integer data type, must change to character
#pcoa$group<-as.character(pcoa$group)
plot(pcoa$pcoa_axis1~pcoa$day, pch=19, col=group.col(pcoa$group), ylab="PCOA axis 1 (10.6%)", xlab="Time (days)")
legend("bottomleft",legend=c("FMT donor", "ABX only", "CDI: no FMT", "CDI: FMT (d11)", "CDI: FMT (d12)"), col=c("green4", "magenta1", "deepskyblue1", "olivedrab2", "green3"), cex=0.8, pch=19)

#color palette for group:
group.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "abx" ) {
colorvec[i] = "magenta1"
}
if (n[i] == "donor" ) {
colorvec[i] = "green4"
}
if( n[i] == "CDI_noFMT") {
colorvec[i] = "deepskyblue1"
}
if( n[i] == "CDI_FMT1") {
colorvec[i] = "olivedrab2"
}
if ( n[i] == "CDI_FMT2" ) {
colorvec[i] = "green3"
}
}
c(colorvec)
}

```

For NMDS:

```{r}

pcoa<-read.table(file="mouseRM_div.thetayc.txt", header=TRUE)
par(mfrow=c(1,3))
#class(pcoa$group_time)  #if integer data type, must change to character
#pcoa$group_time<-as.character(pcoa$group_time)
plot(pcoa$nmds_axis2~pcoa$nmds_axis1, pch=19, col=time.col(pcoa$group_time))
plot(pcoa$nmds_axis2~pcoa$nmds_axis3, pch=19, col=time.col(pcoa$group_time))
plot(pcoa$nmds_axis3~pcoa$nmds_axis1, pch=19, col=time.col(pcoa$group_time))
legend("bottomright",legend=c("pre-treatment", "ABX only: during cef", "ABX only: post-cef", "ABX only: during vanco", "ABX only: post-vanco", "CDI: pre-vanco", "CDI: post-vanco treatment", "CDI: post-IP clinda", "CDI: post-FMT", "CDI: post-FMT, post-clinda", "donor FMT"), col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"), cex=0.8, pch=19)

#1 graph at a time:
plot(pcoa$nmds_axis3~pcoa$nmds_axis1, pch=19, col=time.col(pcoa$group_time) )
legend("bottomright",legend=c("pre-treatment", "ABX only: during cef", "ABX only: post-cef", "ABX only: during vanco", "ABX only: post-vanco", "CDI: pre-vanco", "CDI: post-vanco treatment", "CDI: post-IP clinda", "CDI: post-FMT", "CDI: post-FMT, post-clinda", "donor FMT"), col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"), cex=0.8, pch=19)


##color palettes:
#Color palette by treatment group over time:
time.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "pre" ) {
colorvec[i] = "green4"
}
if (n[i] == "fmt" ) {
colorvec[i] = "green3"
}
if( n[i] == "cef") {
colorvec[i] = "hotpink"
}
if( n[i] == "post_cef") {
colorvec[i] = "magenta1"
}
if ( n[i] == "vanco" ) {
colorvec[i] = "magenta4"
}
if (n[i] == "post_vanco" ) {
colorvec[i] = "red4"
}
if( n[i] == "post_CDI") {
colorvec[i] = "cyan"
}
if( n[i] == "post_CDI_vanco") {
colorvec[i] = "deepskyblue1"
}
if (n[i] == "post_CDI_vanco_clinda" ) {
colorvec[i] = "blue1"
}
if( n[i] == "post_FMT") {
colorvec[i] = "olivedrab2"
}
if( n[i] == "post_FMT_clinda") {
colorvec[i] = "gold1"
}
}
c(colorvec)
}

```

I also calculated biplot values using the corr.axes function in mothur. 
	+ This was used to define the OTUs that drove the biplot
	+ files used: corr.axes files (listed at beginning of file)
This was also the code used for Fig. 4 


```{r}

#note: the beginning part of this code required filtering some of the identified correlating OTUs
#this means that I defined my cutoffs by looking at all of the numbers, and then chose the coordinates that fit those cutoffs
#this will be different for each dataset, but I have explained the process of elimination for my plots

# NMDS biplot (Fig. 4):

#opened file: mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.nmds.spearman.corr.axes
	#saved as: mouseRM.corraxes_nmds.xlsx
	#added taxonomy with file mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy
	#filtered by: 
		-axis1_pvalue, axis2_pvalue, then axis3_pvalue (nested sorting): 
			this gave 262 OTUs that were sign. (p<0.001) for axis1 alone (blue), 243 that were significant for axes 1 and 2 (purple), and 193 OTUs that were sign. for all three axes (green)
			(need more sorting)
		-within each category, also sorted by size (number of OTUs) and length
		-For now, from each of these, will choose OTUs that have >100,000 in abundance (size) and >0.5 in length (note: none were isolated from the 'only axis 1' group)
		-this totals 9 OTUs
		-once the graph is plotted, can go back and check if any other OTUs might be interesting additions
		
#graphing onto nmds plot:

#recreate nmds plot:

#Color palette by treatment group over time:
time.col <- function(n) {
colorvec <- vector(mode="character", length=length(n))
for (i in 1:length(n)) {
colorvec[i] = "light grey"
if ( n[i] == "pre" ) {
colorvec[i] = "green4"
}
if (n[i] == "fmt" ) {
colorvec[i] = "green3"
}
if( n[i] == "cef") {
colorvec[i] = "hotpink"
}
if( n[i] == "post_cef") {
colorvec[i] = "magenta1"
}
if ( n[i] == "vanco" ) {
colorvec[i] = "magenta4"
}
if (n[i] == "post_vanco" ) {
colorvec[i] = "red4"
}
if( n[i] == "post_CDI") {
colorvec[i] = "cyan"
}
if( n[i] == "post_CDI_vanco") {
colorvec[i] = "deepskyblue1"
}
if (n[i] == "post_CDI_vanco_clinda" ) {
colorvec[i] = "blue1"
}
if( n[i] == "post_FMT") {
colorvec[i] = "olivedrab2"
}
if( n[i] == "post_FMT_clinda") {
colorvec[i] = "gold1"
}
}
c(colorvec)
}

#plot it:
pcoa<-read.table(file="mouseRM_div.thetayc.txt", header=TRUE)
class(pcoa$group_time)  														#for color function to work
pcoa$group_time<-as.character(pcoa$group_time)

#testing out arrows:
plot(pcoa$nmds_axis2~pcoa$nmds_axis1, pch=19, col=time.col(pcoa$group_time), ylim=c(-0.8, 0.8), xlim=c(-0.8, 0.8))
arrows(0, 0, x1=-0.777181, y1=-0.705047, lty=1, length=0.1)						#this is using default R packages

library(shape)
plot(pcoa$nmds_axis2~pcoa$nmds_axis1, pch=19, col=time.col(pcoa$group_time), ylim=c(-0.8, 0.8), xlim=c(-0.8, 0.8))
Arrows(0, 0, x1=-0.777181, y1=-0.705047, lty=1, arr.length=0.3, arr.type="triangle")					#this is using the package 'shape' (it has more options)
text(-0.777181, -0.705047, label="Lactobacillus", cex=.5, pos=1)

#altogether:
plot(pcoa$nmds_axis2~pcoa$nmds_axis1, pch=19, col=time.col(pcoa$group_time), 
		ylim=c(-0.8, 0.8), xlim=c(-0.8, 0.8), xlab="NMDS axis 1", ylab="NMDS axis 2")
Arrows(0, 0, x1=-0.777181, y1=-0.705047, lty=1, arr.length=0.3, arr.type="triangle")
text(-0.777181, -0.705047, label="Lactobacillus", cex=.8, pos=4)

Arrows(0, 0, x1=0.784767, y1=0.456718, lty=1, arr.length=0.3, arr.type="triangle")
text(0.784767, 0.456718, label="Bacteroides", cex=.8, pos=2)

Arrows(0, 0, x1=0.746129, y1=0.63305, lty=1, arr.length=0.3, arr.type="triangle")						#all of these are Porphyromonadaceae, so just labeled once (although different OTUs)
Arrows(0, 0, x1=0.716142, y1=0.709681, lty=1,  arr.length=0.3, arr.type="triangle")
Arrows(0, 0, x1=0.727474, y1=0.691207, lty=1, arr.length=0.3, arr.type="triangle")
Arrows(0, 0, x1=0.75561, y1=0.621556, lty=1, arr.length=0.3, arr.type="triangle")
Arrows(0, 0, x1=0.721277, y1=0.654462, lty=1, arr.length=0.3, arr.type="triangle")
text(0.716142, 0.709681, label="Porphyromonadaceae (5 OTUs)", cex=.8, pos=2)

Arrows(0, 0, x1=0.569287, y1=-0.227633, lty=1, arr.length=0.3, arr.type="triangle")
text(0.569287, -0.227633, label="Enterobacteriaceae", cex=.8, pos=1)

Arrows(0, 0, x1=-0.426376, y1=-0.303904, lty=1, arr.length=0.3, arr.type="triangle")
text(-0.426376, -0.303904, label="Clostridium XI", cex=.8, pos=2)

legend("topleft",legend=c("pre-treatment", "ABX only: during cef", "ABX only: post-cef", "ABX only: during vanco", "ABX only: post-vanco", "CDI: pre-vanco", "CDI: post-vanco treatment", "CDI: post-IP clinda", "CDI: post-FMT", "CDI: post-FMT, post-clinda", "donor FMT"), col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"), cex=0.8, pch=19)

--
#for pcoa:

#opened up file:
	#copied onto sheet 'corr.axes_pcoa' in same excel file as the nmds data
	#added taxonomy info from same sheet
	#sorted in same way:
		#only axis 1 (blue): 291 OTUs that were sign (p<0.001)
		#axis 1 and 2 (purple): only 23!
		#all 3 axes (green): only 10 OTUs
	#sorted by size and length, again (<100,000 in size, and <0.5 in length):
		#this gave 3 OTUs that were significant across all axes, and 8 that were significant ONLY for axis 1 (which explains most of the difference)

#plotted using the pcoa (same file as with nmds):

plot(pcoa$pcoa_axis2~pcoa$pcoa_axis1, pch=19, col=time.col(pcoa$group_time), 
		ylim=c(-1, 1), xlim=c(-1, 1), xlab="PCOA 1 (15.5%)", ylab="PCOA 2 (10.7%)")
Arrows(0, 0, x1=-0.367244, y1=0.912414, lty=1, arr.length=0.3, arr.type="triangle")
text(-0.367244, 0.912414, label="Akkermansia", cex=.8, pos=4)

Arrows(0, 0, x1=-0.723751, y1=-0.262051, lty=1, arr.length=0.3, arr.type="triangle")
text(-0.723751, -0.262051, label="Lactobacillus", cex=.8, pos=4)

Arrows(0, 0, x1=0.338274, y1=-0.385487, lty=1, arr.length=0.3, arr.type="triangle")
text(0.338274, -0.385487, label="Enterobactericeae", cex=.8, pos=4)

Arrows(0, 0, x1=0.766286, y1=-0.008976, lty=1, arr.length=0.3, arr.type="triangle")
text(0.766286, -0.008976, label="Bacteroides", cex=.8, pos=1)

Arrows(0, 0, x1=0.19629, y1=0.023317, lty=1, arr.length=0.3, arr.type="triangle")
text(0.19629, 0.023317, label="Turicibacter", cex=.8, pos=1)

Arrows(0, 0, x1=0.735185, y1=0.090545, lty=1, arr.length=0.3, arr.type="triangle")
Arrows(0, 0, x1=0.805379, y1=0.082323, lty=1, arr.length=0.3, arr.type="triangle")
Arrows(0, 0, x1=0.807154, y1=0.055562, lty=1, arr.length=0.3, arr.type="triangle")
Arrows(0, 0, x1=0.774691, y1=0.134099, lty=1, arr.length=0.3, arr.type="triangle")
text(0.774691, 0.134099, label="Porphyromonadaceae (4 OTUs)", cex=.8, pos=3)

Arrows(0, 0, x1=0.770544, y1=0.081607, lty=1, arr.length=0.3, arr.type="triangle")
text(0.770544, 0.081607, label="Barnesiella", cex=.8, pos=4)

Arrows(0, 0, x1=-0.480743, y1=0.178618, lty=1, arr.length=0.3, arr.type="triangle")
text(-0.480743, 0.178618, label="Clostridium XI", cex=.8, pos=4)

legend("topright",legend=c("pre-treatment", "ABX only: during cef", "ABX only: post-cef", "ABX only: during vanco", "ABX only: post-vanco", "CDI: pre-vanco", "CDI: post-vanco treatment", "CDI: post-IP clinda", "CDI: post-FMT", "CDI: post-FMT, post-clinda", "donor FMT"), col=c("green4", "hotpink", "magenta1", "magenta4", "red4", "cyan", "deepskyblue1", "blue1", "olivedrab2", "gold1", "green3"), cex=0.8, pch=19)

```


	
