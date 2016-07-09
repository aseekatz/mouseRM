# mouseRM
Code and files for generation of figures in Mouse Relapse paper

#### Mouse Relapse Model READ ME file#####
3.12.15
Anna M. Seekatz

this file explains what each of the files in this directory are, and what files are necessary for each figure

this depository should include all of the files necessary to reproduce the same graphs (including raw data files)

note: I am not a coder, so please excuse the messy code--if you have the capability of making this code more efficient, please feel free to adjust! (the addition of functions are encouraged)

note2: I tried to comment where necessary, but some things may not be clear

use this code as you want, but please cite the Mouse Relapse Model paper when applicable


###
------
###

##### Sequences were processed as follows:
- mothur v.1.33.3
- RDP database (v9, trainset adapted and downloaded from mothur)
- SILVA alignment (v109, downloaded from mothur)
- batch was ran on the University of Michigan FLUX HPC cluster
- list of files are available in file: mouseRM1.files

#####mothur batch file, or pipeline of sequence processing:

make.contigs(file=mouseRM1.files, processors=10)
summary.seqs(fasta=mouseRM1.trim.contigs.fasta, processors=10)
screen.seqs(fasta=mouseRM1.trim.contigs.fasta, group=mouseRM1.contigs.groups, maxambig=0, maxlength=275, processors=10)
unique.seqs(fasta=mouseRM1.trim.contigs.good.fasta)
count.seqs(name=mouseRM1.trim.contigs.good.names, group=mouseRM1.contigs.good.groups)
summary.seqs(count=mouseRM1.trim.contigs.good.count_table, processors=10)
pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=10)
system(mv silva.bacteria.pcr.fasta silva.v4.fasta)
summary.seqs(fasta=silva.v4.fasta, processors=10)
align.seqs(fasta=mouseRM1.trim.contigs.good.unique.fasta, reference=silva.v4.fasta, processors=10)
summary.seqs(fasta=mouseRM1.trim.contigs.good.unique.align, count=mouseRM1.trim.contigs.good.count_table, processors=10)
screen.seqs(fasta=mouseRM1.trim.contigs.good.unique.align, count=mouseRM1.trim.contigs.good.count_table, summary=mouseRM1.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8, processors=10)
summary.seqs(fasta=current, count=current, processors=10)
filter.seqs(fasta=mouseRM1.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=10)
unique.seqs(fasta=mouseRM1.trim.contigs.good.unique.good.filter.fasta, count=mouseRM1.trim.contigs.good.good.count_table)
pre.cluster(fasta=mouseRM1.trim.contigs.good.unique.good.filter.unique.fasta, count=mouseRM1.trim.contigs.good.unique.good.filter.count_table, diffs=2, processors=10)
chimera.uchime(fasta=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t, processors=10)
remove.seqs(fasta=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.accnos)
summary.seqs(fasta=current, count=current, processors=10)
remove.lineage(fasta=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, taxonomy=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
count.seqs(name=current, group=current)

remove.groups(count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, fasta=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, groups=d10_1008_2)
classify.seqs(fasta=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80)

dist.seqs(fasta=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.15)
cluster(column=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table)
make.shared(list=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, label=0.03)
classify.otu(list=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pds.wang.taxonomy, label=0.03)
phylotype(taxonomy=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pds.wang.taxonomy)
classify.otu(list=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.list, count=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.pick.count_table, taxonomy=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pds.wang.taxonomy, label=1)

######extra commands for mouseRM1.plate in mothur:
dist.shared(shared=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, calc=thetayc-jclass)

######to get thetayc distances:
pcoa(phylip=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.jclass.0.03.lt.dist)
pcoa(phylip=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.dist)
nmds(phylip=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.dist, mindim=4, maxdim=4)
- Output for these=.axes and .loadings for pcoa
        
######for shared distances:
summary.shared(shared=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, calc=sharedsobs-braycurtis-spearman-thetayc)

######for diversity metrics (single):
summary.single(shared=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, calc=simpsoneven-simpson-invsimpson-shannon-npshannon-sobs-chao-nseqs)

amova(phylip=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.dist, design=mouseRM_group.time.design.txt)
metastats(shared=mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, design=mouseRM_group.time.design.txt, processors=10)


-  note: negative extraction/PCR controls and positive PCR (mock) controls were ran alongside sequencing of this run, all of which were on the same plate
- negative controls did not amplify during processing


###
------
###

#####These files were used to analyze the data and create the figures published in Seekatz et al.:

##### Mean CFU (colonization) and median toxin (reciprocal log) overtime (Fig. 1BC, 2BC):
- R code: mouseRM_toxin.cfu.R
- files:
	- tox1: RM1.2.FMT_col.tox.txt (contains all toxin information from first 2 repeats of the relapse model)
	- tox2: RM3.FMT_col.tox.txt	(contains all toxin info from the 3rd (FMT) experiment)
	- all.samples_mouse.var.txt	(contains all the metadata necessary to categorize and group the data for all samples collected)

##### Histology results (Fig. 1D/E, S2):
- R code: mouseRM_histology.R
- files:
	- all.samples_mouse.var.txt	(mouse metadata for all samples)
	- histo.txt (this is a dataframe of the histology scores for each sample)

##### Mouse % weight loss (Fig. S1):
- R code: mouseRM_weights.R
- files:
	- all.samples_mouse.var.txt	(mouse metadata for all samples)
	- RM1_weight.matrix.txt
	- RM2_weight.matrix.txt
	- RM3_weight.matrix.txt
		
##### Sample Barplots of the relative abundance of top genera (above 2%) (Fig. S3):
- R code: 
- files:
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary (this will eventually be sorted)
	- mouseRM_genbar.txt (modified as described--this is just to show an example of what I filtered)
	- mouseRM_genfrac2p.txt (further filtered to only include the top 98% of genera; also converted values into relative abundance)
	- mouseRM_phylofrac2p_ordered.txt (same file, only ordered by phylum and relative abundance within each phylum)
	
##### Median bargraphs of top abundant genera (Fig. 3):
- R code: mouseRM_top.genera.bargraphs.R
- files:
	- mouseRM_phylofrac2p_for.barplot.txt (this is a modified version of all of the genera above 2% of the full dataset, created when making previous barplots of each sample)
	- mouse.var.txt
		
##### PCOA or NMDS biplot oordination (Fig. 4):
- R code: mouseRM_pcoa.nmds.biplot.R
- files: 
	- mouse.var.txt (this file has the necessary mouse metadata for sequenced samples, i.e. cage, treatment, different groupings, etc)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.groups.summary (alpha diversity measures calculated using mothur)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.nmds.axes (nmds axis loadings)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.pcoa.axes (pcoa axis loadings)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.nmds.spearman.corr.axes (for identifying correlated OTUs for a biplot)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.pcoa.spearman.corr.axes

##### Heatmap of top OTUs by sample (with unsupervised clustering) (Fig. S4):
- R code: mouseRM_heatmap.R
- files:
	- mouseRM_otucounts.txt (this is modified to be a matrix from the mothur-produced file, mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared)
	- mouse.var.txt

		
##### Alpha Diversity measures (Fig. 5A):
- R code: mouseRM_div.measures.R
- files:
	- mouse.var.txt (this file has the necessary mouse metadata for sequenced samples, i.e. cage, treatment, different groupings, etc)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.groups.summary (alpha diversity measures calculated using mothur)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.nmds.axes (nmds axis loadings)
	- mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.thetayc.0.03.lt.pcoa.axes (pcoa axis loadings)
		
##### Beta diversity measures, or pairwise comparisons of samples, such as theta yc similarity measure (Fig. 5B, C):
- R code: mouseRM_thetayc.R
- files:
	- mouseRM_shared.dist.txt: has pairwise distances (modified from mouseRM1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.summary, created in mothur)
	- mouse.var.txt (mouse metadata)
		
###
------
###

The End
