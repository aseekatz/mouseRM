##Plotting histology results using R
# A. Seekatz
# 1.19.15

##
--------------
##Plotting histology results:

##made histo.txt, which contains all histo results (both cecum and colon)

# merged histo.txt with mouse.var.txt to combine with metadata:

mouse.var<-read.table(file="all.samples_mouse.var.txt", header=TRUE)
histo<-read.table(file="histo.txt", header=TRUE)
h<-merge(histo, mouse.var, by.x=c("sampleID"), by.y=c("sampleID")) #this merges the data based on sampleID (note: the column names do not have to match)

# subset into cecum and colon results:

#split up summary by group:
h.split<-split(h, h$tissue3)
cecum<-h.split$'cecum'
colon<-h.split$'colon'

--
# plotting summary_score cecum individually, adding a median line:

# levels(cecum$treatment2)			#lets you know how many treatment groups you have
#for whatever reason, also shows all levels from mouse.var...(this includes FMT, which is NOT in this data set)

#get the median summary_score:

cecum.h.sum <- ddply(cecum, c("treatment2", "Day"), summarise,
               N    = length(summary_score),
               median = median(summary_score, na.rm=TRUE),
               iqr   = IQR(summary_score, na.rm=TRUE),
               lowq	 = quantile(summary_score, 0.25, na.rm=TRUE),
               highq  = quantile(summary_score, 0.75, na.rm=TRUE)
               )

head(cecum.h.sum)				#take a look at your table
#  treatment2 Day N median  iqr lowq highq
#1   abx_only  16 4      2 0.75  1.5  2.25
#2     no_FMT   4 3      7 1.00  6.5  7.50
#3     no_FMT   9 3      1 1.00  0.5  1.50

#split up summary by group:
split.cecum.h.sum<-split(cecum.h.sum, cecum.h.sum$treatment2)
cecum.abx.sum<-split.cecum.h.sum$'abx_only'
cecum.vanco.sum<-split.cecum.h.sum$'no_FMT'
cecum.novanco.sum<-split.cecum.h.sum$'no_vanco'

#also split up individual points:
cecum.split.all<-split(cecum, cecum$treatment2)
cecum.abx<-cecum.split.all$'abx_only'
cecum.vanco<-cecum.split.all$'no_FMT'
cecum.novanco<-cecum.split.all$'no_vanco'

#plot the groups:
plot(0.25+jitter(cecum.vanco$Day, amount=0.25), cecum.vanco$summary_score, 
	pch=21, xaxt='n', col="black", bg="blue4", cex=0.8,
	ylab="Summary score (cecum)", xlab="Day", xlim=c(0,16), ylim=c(0,10))
axis(1, at=cecum$Day, labels=cecum$Day)
arrows(0.25+cecum.vanco.sum$Day-0.05, cecum.vanco.sum$median, 0.25+cecum.vanco.sum$Day+0.05, cecum.vanco.sum$median, angle=180, code=3, col="blue4", lwd=2, length=.10)
points(jitter(cecum.novanco$Day, amount=0.5)-0.25, cecum.novanco$summary_score, pch=21, col="black", bg="orange", cex=0.8)
arrows(cecum.novanco.sum$Day-0.05-0.25, cecum.novanco.sum$median, cecum.novanco.sum$Day+0.05-0.25, cecum.novanco.sum$median, angle=180, code=3, col="orange", lwd=2, length=.10)
legend("bottomleft", c("Relapsed (vanco)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "orange"), pt.bg=c("blue4", "orange"), cex=0.8, lty=1)

# do the same for the colon:

colon.h.sum <- ddply(colon, c("treatment2", "Day"), summarise,
               N    = length(summary_score),
               median = median(summary_score, na.rm=TRUE),
               iqr   = IQR(summary_score, na.rm=TRUE),
               lowq	 = quantile(summary_score, 0.25, na.rm=TRUE),
               highq  = quantile(summary_score, 0.75, na.rm=TRUE)
               )

#split up summary by group:
split.colon.h.sum<-split(colon.h.sum, colon.h.sum$treatment2)
colon.abx.sum<-split.colon.h.sum$'abx_only'
colon.vanco.sum<-split.colon.h.sum$'no_FMT'
colon.novanco.sum<-split.colon.h.sum$'no_vanco'

#also split up individual points:
colon.split.all<-split(colon, colon$treatment2)
colon.abx<-colon.split.all$'abx_only'
colon.vanco<-colon.split.all$'no_FMT'
colon.novanco<-colon.split.all$'no_vanco'

plot(0.25+jitter(colon.vanco$Day, amount=0.25), colon.vanco$summary_score, 
	pch=21, xaxt='n', col="black", bg="blue4", cex=0.8,
	ylab="Summary score (colon)", xlab="Day", xlim=c(0,16), ylim=c(0,10))
axis(1, at=colon$Day, labels=colon$Day)
arrows(0.25+colon.vanco.sum$Day-0.05, colon.vanco.sum$median, 0.25+colon.vanco.sum$Day+0.05, colon.vanco.sum$median, angle=180, code=3, col="blue4", lwd=2, length=.10)
points(jitter(colon.novanco$Day, amount=0.5)-0.25, colon.novanco$summary_score, pch=21, col="black", bg="orange", cex=0.8)
arrows(colon.novanco.sum$Day-0.05-0.25, colon.novanco.sum$median, colon.novanco.sum$Day+0.05-0.25, colon.novanco.sum$median, angle=180, code=3, col="orange", lwd=2, length=.10)
legend("bottomleft", c("Relapsed (vanco)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "orange"), pt.bg=c("blue4", "orange"), cex=0.8, lty=1)

--
# for both tissue summary scores:

par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
plot(0.25+jitter(cecum.vanco$Day, amount=0.25), cecum.vanco$summary_score, 
	pch=21, xaxt='n', col="black", bg="blue4", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Summary score (cecum)", xlab=NA, xlim=c(3,17), ylim=c(0,10))
axis(1, at=cecum$Day, labels=cecum$Day, cex.axis=0.8)
arrows(0.25+cecum.vanco.sum$Day-0.05, cecum.vanco.sum$median, 0.25+cecum.vanco.sum$Day+0.05, cecum.vanco.sum$median, angle=180, code=3, col="blue4", lwd=2, length=.10)
points(jitter(cecum.novanco$Day, amount=0.5)-0.25, cecum.novanco$summary_score, pch=21, col="black", bg="orange", cex=0.8)
arrows(cecum.novanco.sum$Day-0.05-0.25, cecum.novanco.sum$median, cecum.novanco.sum$Day+0.05-0.25, cecum.novanco.sum$median, angle=180, code=3, col="orange", lwd=2, length=.10)
legend("bottomleft", c("Relapsed (vanco)", "Primary CDI (no vanco)"), pch=21, col=c("blue4", "orange"), pt.bg=c("blue4", "orange"), cex=0.6, lty=1)

plot(0.25+jitter(colon.vanco$Day, amount=0.25), colon.vanco$summary_score, 
	pch=21, xaxt='n', col="black", bg="blue4", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Summary score (colon)", xlab="Day", xlim=c(3,17), ylim=c(0,10))
axis(1, at=colon$Day, labels=colon$Day, cex.axis=0.8)
arrows(0.25+colon.vanco.sum$Day-0.05, colon.vanco.sum$median, 0.25+colon.vanco.sum$Day+0.05, colon.vanco.sum$median, angle=180, code=3, col="blue4", lwd=2, length=.10)
points(jitter(colon.novanco$Day, amount=0.5)-0.25, colon.novanco$summary_score, pch=21, col="black", bg="orange", cex=0.8)
arrows(colon.novanco.sum$Day-0.05-0.25, colon.novanco.sum$median, colon.novanco.sum$Day+0.05-0.25, colon.novanco.sum$median, angle=180, code=3, col="orange", lwd=2, length=.10)

###
--
###

# For supplemental figure that shows all scores separately:

#split data into cecum and colon (as above)
#split into treatment
#graph points/median/etc for each score around each day per treatment

#using h data set, the cecal and colonic results were already split into cecum and colon (see above)

# define medians (individually) for each of factors you want to plot within each group:

cecum.sum <- ddply(cecum, c("treatment2", "Day"), summarise,
               N    = length(summary_score),
               med.edema = median(edema, na.rm=TRUE),
               med.inflam   = median(inflammation, na.rm=TRUE),
               med.epi	 = median(epithelial_damage, na.rm=TRUE)
               )

#split into treatment groups:
split.cecum<-split(cecum.sum, cecum.sum$treatment2)
cecum.abx<-split.cecum$'abx_only'
cecum.vanco<-split.cecum$'no_FMT'
cecum.novanco<-split.cecum$'no_vanco'

#split individual data points by treatment groups:
split.cecum.all<-split(cecum, cecum$treatment2)
cecum.abx.al<-split.cecum.all$'abx_only'
cecum.vanco.all<-split.cecum.all$'no_FMT'
cecum.novanco.all<-split.cecum.all$'no_vanco'

# plot:

par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
par(xpd=TRUE)
	# vanco-treated:
plot(0.75+jitter(cecum.vanco.all$Day, amount=0.25), cecum.vanco.all$edema, 
	pch=21, xaxt='n', col="black", bg="darkslategray3", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Histopathological scores (cecum)", ylim=c(0,5), xlim=c(3,17), xlab=NA, main="Relapse CDI (vanco)")
axis(1, at=cecum$Day, labels=cecum$Day, cex.axis=0.8)
arrows(0.75+cecum.vanco$Day-0.05, cecum.vanco$med.edema, 0.75+cecum.vanco$Day+0.05, cecum.vanco$med.edema, angle=180, code=3, col="darkslategray3", lwd=2, length=.10)
points(jitter(cecum.vanco.all$Day, amount=0.25)-0.75, cecum.vanco.all$inflammation, pch=21, col="black", bg="deepskyblue", cex=0.8)
arrows(cecum.vanco$Day-0.05-0.75, cecum.vanco$med.inflam, cecum.vanco$Day+0.05-0.75, cecum.vanco$med.inflam, angle=180, code=3, col="deepskyblue", lwd=2, length=.10)
points(jitter(cecum.vanco.all$Day, amount=0.25), cecum.vanco.all$epithelial_damage, pch=21, col="black", bg="dodgerblue4", cex=0.8)
arrows(cecum.vanco$Day-0.05, cecum.vanco$med.epi, cecum.vanco$Day+0.05, cecum.vanco$med.epi, angle=180, code=3, col="dodgerblue4", lwd=2, length=.10)
#legend(3,-1, c("edema", "inflammation", "epithelial damage"), pch=21, col=c("black"), title="Score:", pt.bg=c("darkslategray3", "deepskyblue", "dodgerblue4"), cex=0.8, lty=1)
	# no-vanco:
plot(0.75+jitter(cecum.novanco.all$Day, amount=0.25), cecum.novanco.all$edema, 
	pch=21, xaxt='n', col="black", bg="darkslategray3", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Histopathological scores (cecum)", ylim=c(0,5), xlim=c(3,17), xlab="Day", main="Primary CDI (no vanco)")
axis(1, at=cecum$Day, labels=cecum$Day, cex.axis=0.8)
arrows(0.75+cecum.novanco$Day-0.05, cecum.novanco$med.edema, 0.75+cecum.novanco$Day+0.05, cecum.novanco$med.edema, angle=180, code=3, col="darkslategray3", lwd=2, length=.10)
points(jitter(cecum.novanco.all$Day, amount=0.25)-0.75, cecum.novanco.all$inflammation, pch=21, col="black", bg="deepskyblue", cex=0.8)
arrows(cecum.novanco$Day-0.05-0.75, cecum.novanco$med.inflam, cecum.novanco$Day+0.05-0.75, cecum.novanco$med.inflam, angle=180, code=3, col="deepskyblue", lwd=2, length=.10)
points(jitter(cecum.novanco.all$Day, amount=0.25), cecum.novanco.all$epithelial_damage, pch=21, col="black", bg="dodgerblue4", cex=0.8)
arrows(cecum.novanco$Day-0.05, cecum.novanco$med.epi, cecum.novanco$Day+0.05, cecum.novanco$med.epi, angle=180, code=3, col="dodgerblue4", lwd=2, length=.10)
legend(3,-1, c("edema", "inflammation", "epithelial damage"), pch=21, col=c("black"), pt.bg=c("darkslategray3", "deepskyblue", "dodgerblue4"), cex=0.8, lty=1)

# do same for colon:

colon.sum <- ddply(colon, c("treatment2", "Day"), summarise,
               N    = length(summary_score),
               med.edema = median(edema, na.rm=TRUE),
               med.inflam   = median(inflammation, na.rm=TRUE),
               med.epi	 = median(epithelial_damage, na.rm=TRUE)
               )

#split into treatment groups:
split.colon<-split(colon.sum, colon.sum$treatment2)
colon.abx<-split.colon$'abx_only'
colon.vanco<-split.colon$'no_FMT'
colon.novanco<-split.colon$'no_vanco'

#split individual data points by treatment groups:
split.colon.all<-split(colon, colon$treatment2)
colon.abx.al<-split.colon.all$'abx_only'
colon.vanco.all<-split.colon.all$'no_FMT'
colon.novanco.all<-split.colon.all$'no_vanco'

# plot all four graphs together:

par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
par(xpd=TRUE)

#cecum:
	# vanco-treated:
plot(0.75+jitter(cecum.vanco.all$Day, amount=0.25), cecum.vanco.all$edema, 
	pch=21, xaxt='n', col="black", bg="darkslategray3", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Histopathological scores (cecum)", ylim=c(0,5), xlim=c(3,17), xlab=NA, main="Relapse CDI (vanco)")
axis(1, at=cecum$Day, labels=cecum$Day, cex.axis=0.8)
arrows(0.75+cecum.vanco$Day-0.05, cecum.vanco$med.edema, 0.75+cecum.vanco$Day+0.05, cecum.vanco$med.edema, angle=180, code=3, col="darkslategray3", lwd=2, length=.10)
points(jitter(cecum.vanco.all$Day, amount=0.25)-0.75, cecum.vanco.all$inflammation, pch=21, col="black", bg="deepskyblue", cex=0.8)
arrows(cecum.vanco$Day-0.05-0.75, cecum.vanco$med.inflam, cecum.vanco$Day+0.05-0.75, cecum.vanco$med.inflam, angle=180, code=3, col="deepskyblue", lwd=2, length=.10)
points(jitter(cecum.vanco.all$Day, amount=0.25), cecum.vanco.all$epithelial_damage, pch=21, col="black", bg="dodgerblue4", cex=0.8)
arrows(cecum.vanco$Day-0.05, cecum.vanco$med.epi, cecum.vanco$Day+0.05, cecum.vanco$med.epi, angle=180, code=3, col="dodgerblue4", lwd=2, length=.10)
	# no-vanco:
plot(0.75+jitter(cecum.novanco.all$Day, amount=0.25), cecum.novanco.all$edema, 
	pch=21, xaxt='n', col="black", bg="darkslategray3", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Histopathological scores (cecum)", ylim=c(0,5), xlim=c(3,17), xlab="Day", main="Primary CDI (no vanco)")
axis(1, at=cecum$Day, labels=cecum$Day, cex.axis=0.8)
arrows(0.75+cecum.novanco$Day-0.05, cecum.novanco$med.edema, 0.75+cecum.novanco$Day+0.05, cecum.novanco$med.edema, angle=180, code=3, col="darkslategray3", lwd=2, length=.10)
points(jitter(cecum.novanco.all$Day, amount=0.25)-0.75, cecum.novanco.all$inflammation, pch=21, col="black", bg="deepskyblue", cex=0.8)
arrows(cecum.novanco$Day-0.05-0.75, cecum.novanco$med.inflam, cecum.novanco$Day+0.05-0.75, cecum.novanco$med.inflam, angle=180, code=3, col="deepskyblue", lwd=2, length=.10)
points(jitter(cecum.novanco.all$Day, amount=0.25), cecum.novanco.all$epithelial_damage, pch=21, col="black", bg="dodgerblue4", cex=0.8)
arrows(cecum.novanco$Day-0.05, cecum.novanco$med.epi, cecum.novanco$Day+0.05, cecum.novanco$med.epi, angle=180, code=3, col="dodgerblue4", lwd=2, length=.10)
legend("topright", c("edema", "inflammation", "epithelial damage"), pch=21, col=c("black"), pt.bg=c("darkslategray3", "deepskyblue", "dodgerblue4"), cex=0.8, lty=1)

#colon:
	# vanco-treated:
plot(0.75+jitter(colon.vanco.all$Day, amount=0.25), colon.vanco.all$edema, 
	pch=21, xaxt='n', col="black", bg="darkslategray3", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Histopathological scores (colon)", ylim=c(0,5), xlim=c(3,17), xlab=NA, main="Relapse CDI (vanco)")
axis(1, at=colon$Day, labels=colon$Day, cex.axis=0.8)
arrows(0.75+colon.vanco$Day-0.05, colon.vanco$med.edema, 0.75+colon.vanco$Day+0.05, colon.vanco$med.edema, angle=180, code=3, col="darkslategray3", lwd=2, length=.10)
points(jitter(colon.vanco.all$Day, amount=0.25)-0.75, colon.vanco.all$inflammation, pch=21, col="black", bg="deepskyblue", cex=0.8)
arrows(colon.vanco$Day-0.05-0.75, colon.vanco$med.inflam, colon.vanco$Day+0.05-0.75, colon.vanco$med.inflam, angle=180, code=3, col="deepskyblue", lwd=2, length=.10)
points(jitter(colon.vanco.all$Day, amount=0.25), colon.vanco.all$epithelial_damage, pch=21, col="black", bg="dodgerblue4", cex=0.8)
arrows(colon.vanco$Day-0.05, colon.vanco$med.epi, colon.vanco$Day+0.05, colon.vanco$med.epi, angle=180, code=3, col="dodgerblue4", lwd=2, length=.10)
	# no-vanco:
plot(0.75+jitter(colon.novanco.all$Day, amount=0.25), colon.novanco.all$edema, 
	pch=21, xaxt='n', col="black", bg="darkslategray3", cex=0.8, cex.lab=0.8, cex.axis=0.8,
	ylab="Histopathological scores (colon)", ylim=c(0,5), xlim=c(3,17), xlab="Day", main="Primary CDI (no vanco)")
axis(1, at=colon$Day, labels=colon$Day, cex.axis=0.8)
arrows(0.75+colon.novanco$Day-0.05, colon.novanco$med.edema, 0.75+colon.novanco$Day+0.05, colon.novanco$med.edema, angle=180, code=3, col="darkslategray3", lwd=2, length=.10)
points(jitter(colon.novanco.all$Day, amount=0.25)-0.75, colon.novanco.all$inflammation, pch=21, col="black", bg="deepskyblue", cex=0.8)
arrows(colon.novanco$Day-0.05-0.75, colon.novanco$med.inflam, colon.novanco$Day+0.05-0.75, colon.novanco$med.inflam, angle=180, code=3, col="deepskyblue", lwd=2, length=.10)
points(jitter(colon.novanco.all$Day, amount=0.25), colon.novanco.all$epithelial_damage, pch=21, col="black", bg="dodgerblue4", cex=0.8)
arrows(colon.novanco$Day-0.05, colon.novanco$med.epi, colon.novanco$Day+0.05, colon.novanco$med.epi, angle=180, code=3, col="dodgerblue4", lwd=2, length=.10)


########
--------
########

# to plot histology results with day 13 collapsed into day 16, as a factor:

# plotting when the time is a factor:
#this is a bit tricky, and you must use stripchart to get the data:

h$day_collapsed<-as.factor(h$day_collapsed)

# For summary score figure:

#split up histo data by tissue group:
h.split<-split(h, h$tissue3)
cecum<-h.split$'cecum'
colon<-h.split$'colon'

#split up by group:
#note: we will be using another way to get the medians in this case
cecum.split.all<-split(cecum, cecum$treatment2)
cecum.abx<-cecum.split.all$'abx_only'
cecum.vanco<-cecum.split.all$'no_FMT'
cecum.novanco<-cecum.split.all$'no_vanco'

colon.split.all<-split(colon, colon$treatment2)
colon.abx<-colon.split.all$'abx_only'
colon.vanco<-colon.split.all$'no_FMT'
colon.novanco<-colon.split.all$'no_vanco'

# plot it:
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))	
par(mfrow=c(1,2))
par(xpd=TRUE)
	#for cecum:
	#plot points:
loc.strip<-1:3					#this gives you a specific location
stripchart(summary_score~day_collapsed, data=cecum.vanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="blue4", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip, cex.axis=0.8, ylab="Summary score (cecum)", cex.lab=0.8, xaxt='n')
stripchart(summary_score~day_collapsed, data=cecum.novanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="orange", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip-0.25, add=TRUE)
stripchart(summary_score~day_collapsed, data=cecum.abx, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="violetred1", cex=0.8, xlim=c(0.5, 3.5), at=3.25, add=TRUE)
axis(1, at=loc.strip, labels=c("primary CDI\n(4 dpi)", "post-vanco\n(9 dpi)", "CDI relapse\n(13-16 dpi)"), cex.axis=0.6)
	#add groups:
cecum.med.sum.vanco<-tapply(cecum.vanco$summary_score, cecum.vanco$day_collapsed, median)
segments(loc.strip-0.15, cecum.med.sum.vanco, loc.strip+0.15, cecum.med.sum.vanco, col="blue4", lwd=2)
cecum.med.sum.novanco<-tapply(cecum.novanco$summary_score, cecum.novanco$day_collapsed, median)
segments((loc.strip-0.25)-0.15, cecum.med.sum.novanco, (loc.strip-0.25)+0.15, cecum.med.sum.novanco, col="orange", lwd=3)
cecum.med.sum.abx<-tapply(cecum.abx$summary_score, cecum.abx$day_collapsed, median)
segments(3.25-0.15, cecum.med.sum.abx, 3.25+0.15, cecum.med.sum.abx, col="violetred1", lwd=3)

	#for colon:
	#plot points:
stripchart(summary_score~day_collapsed, data=colon.vanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="blue4", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip, cex.axis=0.8, ylab="Summary score (colon)", xaxt='n', cex.lab=0.8)
stripchart(summary_score~day_collapsed, data=colon.novanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="orange", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip-0.25, add=TRUE)
stripchart(summary_score~day_collapsed, data=colon.abx, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="violetred1", cex=0.8, xlim=c(0.5, 3.5), at=3.25, add=TRUE)
axis(1, at=loc.strip, labels=c("primary CDI\n(4 dpi)", "post-vanco\n(9 dpi)", "CDI relapse\n(13-16 dpi)"), cex.axis=0.6)
	#add groups:
colon.med.sum.vanco<-tapply(colon.vanco$summary_score, colon.vanco$day_collapsed, median)
segments(loc.strip-0.15, colon.med.sum.vanco, loc.strip+0.15, colon.med.sum.vanco, col="blue4", lwd=2)
colon.med.sum.novanco<-tapply(colon.novanco$summary_score, colon.novanco$day_collapsed, median)
segments((loc.strip-0.25)-0.15, colon.med.sum.novanco, (loc.strip-0.25)+0.15, colon.med.sum.novanco, col="orange", lwd=3)
colon.med.sum.abx<-tapply(colon.abx$summary_score, colon.abx$day_collapsed, median)
segments(3.25-0.15, colon.med.sum.abx, 3.25+0.15, colon.med.sum.abx, col="violetred1", lwd=3)
legend("topleft", c("Relapsed (vanco)", "Primary CDI (no vanco)", "abx control"), pch=21, col=c("blue4", "orange", "violetred1"), pt.bg=c("blue4", "orange", "violetred1"), cex=0.6, lty=1)

###
# For supplemental histology figures:

# since we are using another way of splitting the data to get the median, you do not need to do it beforehand
par(mar=c(4,4,1,1))
par(oma=c(2,0,0,0))	
par(mfrow=c(2,2))
par(xpd=TRUE)
	
# cecum for relapse group:
loc.strip<-1:3					
stripchart(edema~day_collapsed, data=cecum.vanco, vertical=TRUE, ylim=c(0,6), method="jitter", pch=21, col="black", bg="darkslategray3", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip, cex.axis=0.8, ylab="Histopathological score\n(cecum)", cex.lab=0.8, xaxt='n', main="CDI Relapse (vanco)")
stripchart(inflammation~day_collapsed, data=cecum.vanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="deepskyblue", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip-0.25, add=TRUE)
stripchart(epithelial_damage~day_collapsed, data=cecum.vanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="dodgerblue4", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip+0.25, add=TRUE)
axis(1, at=loc.strip, labels=c("primary CDI\n(4 dpi)", "post-vanco\n(9 dpi)", "CDI relapse\n(13-16 dpi)"), cex.axis=0.6)
	#add groups:
cecum.med.ed.vanco<-tapply(cecum.vanco$edema, cecum.vanco$day_collapsed, median)
segments(loc.strip-0.15, cecum.med.ed.vanco, loc.strip+0.15, cecum.med.ed.vanco, col="darkslategray3", lwd=2)
cecum.med.inf.vanco<-tapply(cecum.vanco$inflammation, cecum.vanco$day_collapsed, median)
segments((loc.strip-0.25)-0.15, cecum.med.inf.vanco, (loc.strip-0.25)+0.15, cecum.med.inf.vanco, col="deepskyblue", lwd=3)
cecum.med.epi.vanco<-tapply(cecum.vanco$epithelial_damage, cecum.vanco$day_collapsed, median)
segments((loc.strip+0.25)-0.15, cecum.med.epi.vanco, (loc.strip+0.25)+0.15, cecum.med.epi.vanco, col="dodgerblue4", lwd=3)
	
# cecum for primary CDI group:
stripchart(edema~day_collapsed, data=cecum.novanco, vertical=TRUE, ylim=c(0,6), method="jitter", pch=21, col="black", bg="darkslategray3", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip, cex.axis=0.8, ylab="Histopathological score\n(cecum)", cex.lab=0.8, xaxt='n', main="CDI Relapse (no vanco)")
stripchart(inflammation~day_collapsed, data=cecum.novanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="deepskyblue", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip-0.25, add=TRUE)
stripchart(epithelial_damage~day_collapsed, data=cecum.novanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="dodgerblue4", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip+0.25, add=TRUE)
axis(1, at=loc.strip, labels=c("primary CDI\n(4 dpi)", "post-novanco\n(9 dpi)", "CDI relapse\n(13-16 dpi)"), cex.axis=0.6)
	#add groups:
cecum.med.ed.novanco<-tapply(cecum.novanco$edema, cecum.novanco$day_collapsed, median)
segments(loc.strip-0.15, cecum.med.ed.novanco, loc.strip+0.15, cecum.med.ed.novanco, col="darkslategray3", lwd=2)
cecum.med.inf.novanco<-tapply(cecum.novanco$inflammation, cecum.novanco$day_collapsed, median)
segments((loc.strip-0.25)-0.15, cecum.med.inf.novanco, (loc.strip-0.25)+0.15, cecum.med.inf.novanco, col="deepskyblue", lwd=3)
cecum.med.epi.novanco<-tapply(cecum.novanco$epithelial_damage, cecum.novanco$day_collapsed, median)
segments((loc.strip+0.25)-0.15, cecum.med.epi.novanco, (loc.strip+0.25)+0.15, cecum.med.epi.novanco, col="dodgerblue4", lwd=3)
legend("topright", c("edema", "inflammation", "epithelial damage"), pch=21, col=c("black"), pt.bg=c("darkslategray3", "deepskyblue", "dodgerblue4"), cex=0.8, lty=1)

# colon for relapse group:
loc.strip<-1:3					
stripchart(edema~day_collapsed, data=colon.vanco, vertical=TRUE, ylim=c(0,6), method="jitter", pch=21, col="black", bg="darkslategray3", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip, cex.axis=0.8, ylab="Histopathological score\n(colon)", cex.lab=0.8, xaxt='n', main="CDI Relapse (vanco)")
stripchart(inflammation~day_collapsed, data=colon.vanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="deepskyblue", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip-0.25, add=TRUE)
stripchart(epithelial_damage~day_collapsed, data=colon.vanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="dodgerblue4", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip+0.25, add=TRUE)
axis(1, at=loc.strip, labels=c("primary CDI\n(4 dpi)", "post-vanco\n(9 dpi)", "CDI relapse\n(13-16 dpi)"), cex.axis=0.6)
	#add groups:
colon.med.ed.vanco<-tapply(colon.vanco$edema, colon.vanco$day_collapsed, median)
segments(loc.strip-0.15, colon.med.ed.vanco, loc.strip+0.15, colon.med.ed.vanco, col="darkslategray3", lwd=2)
colon.med.inf.vanco<-tapply(colon.vanco$inflammation, colon.vanco$day_collapsed, median)
segments((loc.strip-0.25)-0.15, colon.med.inf.vanco, (loc.strip-0.25)+0.15, colon.med.inf.vanco, col="deepskyblue", lwd=3)
colon.med.epi.vanco<-tapply(colon.vanco$epithelial_damage, colon.vanco$day_collapsed, median)
segments((loc.strip+0.25)-0.15, colon.med.epi.vanco, (loc.strip+0.25)+0.15, colon.med.epi.vanco, col="dodgerblue4", lwd=3)
	
	# colon for primary CDI group:
stripchart(edema~day_collapsed, data=colon.novanco, vertical=TRUE, ylim=c(0,6), method="jitter", pch=21, col="black", bg="darkslategray3", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip, cex.axis=0.8, ylab="Histopathological score\n(colon)", cex.lab=0.8, xaxt='n', main="CDI Relapse (no vanco)")
stripchart(inflammation~day_collapsed, data=colon.novanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="deepskyblue", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip-0.25, add=TRUE)
stripchart(epithelial_damage~day_collapsed, data=colon.novanco, vertical=TRUE, ylim=c(0,10), method="jitter", pch=21, col="black", bg="dodgerblue4", cex=0.8, xlim=c(0.5, 3.5), at=loc.strip+0.25, add=TRUE)
axis(1, at=loc.strip, labels=c("primary CDI\n(4 dpi)", "post-novanco\n(9 dpi)", "CDI relapse\n(13-16 dpi)"), cex.axis=0.6)
	#add groups:
colon.med.ed.novanco<-tapply(colon.novanco$edema, colon.novanco$day_collapsed, median)
segments(loc.strip-0.15, colon.med.ed.novanco, loc.strip+0.15, colon.med.ed.novanco, col="darkslategray3", lwd=2)
colon.med.inf.novanco<-tapply(colon.novanco$inflammation, colon.novanco$day_collapsed, median)
segments((loc.strip-0.25)-0.15, colon.med.inf.novanco, (loc.strip-0.25)+0.15, colon.med.inf.novanco, col="deepskyblue", lwd=3)
colon.med.epi.novanco<-tapply(colon.novanco$epithelial_damage, colon.novanco$day_collapsed, median)
segments((loc.strip+0.25)-0.15, colon.med.epi.novanco, (loc.strip+0.25)+0.15, colon.med.epi.novanco, col="dodgerblue4", lwd=3)

```

# stats for both Fig. 1D:

```{r}

# statistics for summary score in Fig. 1:

# (in case you need to read in data again):
mouse.var<-read.table(file="all.samples_mouse.var.txt", header=TRUE)
histo<-read.table(file="histo.txt", header=TRUE)
h<-merge(histo, mouse.var, by.x=c("sampleID"), by.y=c("sampleID")) #this merges the data based on sampleID (note: the column names do not have to match)

# same data reading as before:
h$day_collapsed<-h$Day					
h$day_collapsed[h$day_collapsed==13]<-16

#split up histo data by tissue group:
h.split<-split(h, h$tissue3)
cecum<-h.split$'cecum'
colon<-h.split$'colon'

# split up by day:

#cecum:
c.split<-split(cecum, cecum$day_collapsed)
c.d4<-c.split$'4'
c.d9<-c.split$'9'
c.d16<-c.split$'16'

#colon:
co.split<-split(colon, colon$day_collapsed)
co.d4<-co.split$'4'
co.d9<-co.split$'9'
co.d16<-co.split$'16'

#for d16, remove abx:
a.c.d16<-subset(c.d16, (subset=treatment2 %in% c("no_vanco", "no_FMT")))
a.co.d16<-subset(co.d16, (subset=treatment2 %in% c("no_vanco", "no_FMT")))

# do the stats:

wilcox.test(summary_score~treatment2, data=c.d4, paired=FALSE)
	#W = 4.5, p-value = 0.5536
wilcox.test(summary_score~treatment2, data=c.d9, paired=FALSE)
	#W = 0.5, p-value = 0.04083
wilcox.test(summary_score~treatment2, data=a.c.d16, paired=FALSE)
	#W = 23.5, p-value = 0.01794
wilcox.test(summary_score~treatment2, data=co.d4, paired=FALSE)
	#W = 1, p-value = 0.3329
wilcox.test(summary_score~treatment2, data=co.d9, paired=FALSE)
	#W = 3.5, p-value = 0.2308
wilcox.test(summary_score~treatment2, data=a.co.d16, paired=FALSE)
	#W = 24, p-value = 0.01306
	
# out of curiosity, checked what results would be for all three groups on d16 (no_FMT, no_vanco, abx_only):

kruskal.test(summary_score~treatment2, data=c.d16)
	# Kruskal-Wallis chi-squared = 9.8721, df = 2, p-value = 0.007183
kruskal.test(summary_score~treatment2, data=co.d16)
	# Kruskal-Wallis chi-squared = 10.303, df = 2, p-value = 0.005791

```