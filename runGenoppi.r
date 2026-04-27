##########################################################################################
# R script for processing and analyzing Taplin IP-MS data
# to be run for each IP-MS sample set (bait IP + IgG control replicates)

# STEPS:
# (1) extract relevant columns (accession, gene name, uniq_peps, sum intensity)
# (2) log2 transfomation + median normalization of protein intensities in each sample
# (3) remove proteins with no human accession/gene name, w/ <2 unique peptides,
#     or detected in <2 bait samples
# (4) OUTPUT log of filtered proteins
# (5) impute missing values for the remaining proteins in each sample
# (6) perform two-sample moderated t-test to identify significant interactors
#     (log2 FC > 0 & FDR <= 0.1)
# (7) OUTPUT full results table (with InWeb annotation)
# (8) OUTPUT ggpairs plots of sample correlations (color: significance, imputation)
# (9) OUTPUT volcano plots (color: significance, InWeb overlay)

# Authors: Yu-Han Hsu, Samantha Cassity
# Last updated: 2026-04-24
##########################################################################################

rm(list=ls())

library(fs)
library(readxl)
library(stringr)
library(dplyr)
library(genoppi) # development branch >= v1.1.0
library(GGally)
library(ggrepel)

# ----------------------------------------------------------------------------------------
# input arguments

# Example
if (FALSE) {
msFile <- 'RawData/101772_786dataA_15.xlsx'
msSheet <- 'Sheet2'
msSite <- 'Taplin'
outDir <- '260423_Analysis_YH'
date <- '230207'
bait <- 'TCF4'
cellType <- 'NPC'
samples <- unlist(strsplit('4 TCF4 Rep1,5 TCF4 Rep2,6 TCF4 Rep3',','))
controls <-unlist(strsplit('#1 Igg Rep1,#2 Igg Rep2,#3 Igg Rep3',',')) 
baitInWeb <- 'TCF4'
}

# Read arguments from command line
if (TRUE) {
args <- commandArgs(trailingOnly=T)
msFile <- args[1]
msSheet <- args[2]
msSite <- args[3]
outDir <- args[4]
date <- args[5]
bait <- args[6]
cellType <- args[7]
samples <- unlist(strsplit(args[8],','))
controls <- unlist(strsplit(args[9],','))
baitInWeb <- args[10] 
}

# create output directory if needed
ipName <- paste(date,bait,cellType,msSite,sep='.')
baitDir <- fs::path(outDir,ipName)
fs::dir_create(baitDir) 

cat('*** created output directory:',baitDir,'\n')

# ----------------------------------------------------------------------------------------
# (1) extract relevant columns (accession, gene name, uniq_peps, sum intensity)

# read in data from MS report file
msDf <- read_excel(msFile,sheet=msSheet)

refStrs <- strsplit(msDf$reference,'\\|')
msDf['accession_number'] <- sapply(refStrs,function(x) ifelse(length(x)==3,x[2],NA))
	# rows without 3 "|" delimited substrings in reference column
	# set to NA (usually contaminant/non-human protein entries)

msDf['gene'] <- str_extract(msDf$Annotation,'GN=(\\S+) ',group=1)
	# re-extract gene names to avoid Excel string/date conversion errors 

msDf['uniq_peps'] <- apply(msDf[c(paste('Unique',samples),paste('Unique',controls))],1,max)
	# take maximum of # unique peps found across all samples

# extract protein intensity in each sample or control, set missing values as NA
intDf <- msDf[c(paste('Sum Intensity',samples),paste('Sum Intensity',controls))]
intDf[intDf=='NF'] <- NA # replace NF (Not Found) with NA
intDf <- data.frame(apply(intDf,2,as.numeric))
intDf[intDf==0] <- NA # replace zero with NA
scCols <- c(paste('sample',1:length(samples),sep=''),
	paste('control',1:length(controls),sep=''))
names(intDf) <- scCols

intDf <- cbind(msDf[c('accession_number','gene','uniq_peps')],intDf)

cat('*** read in input data for',nrow(intDf),'proteins\n')
cat('*** extracted relevant columns:',names(intDf),'\n')

# ----------------------------------------------------------------------------------------
# (2) log2 transfomation + median normalization of protein intensities in each sample

intDf[scCols] <- log2(intDf[scCols])
intDf[scCols] <- scale(intDf[scCols],center=apply(intDf[scCols],2,median,na.rm=T),scale=F)
	# subtract column median from values in each column (ignore NA)
	# i.e. center median of each column at ~0

cat('*** performed log2 transformation and median normalization\n')

# ----------------------------------------------------------------------------------------
# (3) remove proteins with no human accession/gene name, w/ <2 unique peptides,
# or detected in <2 bait samples

humanFilterInds <- is.na(intDf$accession_number) | is.na(intDf$gene)
pepFilterInds <- intDf$uniq_peps < 2
sampFilterInds <- rowSums(!is.na(intDf[paste('sample',1:length(samples),sep='')])) < 2

outDf <- intDf[!(humanFilterInds | pepFilterInds | sampFilterInds),]

cat('*** filtered out',sum(humanFilterInds | pepFilterInds | sampFilterInds),'proteins,',
	nrow(outDf),'proteins remaining\n')

# ----------------------------------------------------------------------------------------
# (4) OUTPUT log of filtered proteins

filteredDf <- data.frame(reference=c(msDf$reference[humanFilterInds],
	msDf$reference[pepFilterInds],msDf$reference[sampFilterInds]),
	reason=c(rep('no human accession/gene name',sum(humanFilterInds)),
		rep('<2 unique peptides',sum(pepFilterInds)),
		rep('<2 bait IP samples',sum(sampFilterInds))))

filteredDf <- filteredDf %>% group_by(reference) %>% summarise(reason=toString(reason))

output_filter <- fs::path(baitDir,paste0(ipName,'.Genoppi.FilteredProteins.txt'))
write.table(filteredDf,output_filter,quote=F,sep='\t',row.names=F)

cat('*** output filtered proteins:',paste0(ipName,'.Genoppi.FilteredProteins.txt\n'))

# ----------------------------------------------------------------------------------------
# (5) impute missing values for the remaining proteins in each sample

set.seed(123) # so rerunning gives same imputation results

outDf <- cbind(outDf[,!colnames(outDf) %in% scCols],
	impute_na(outDf[scCols],shift=1.8,width=0.3))

cat('*** performed missing value imputation\n')

# ----------------------------------------------------------------------------------------
# (6) perform two-sample moderated t-test to identify significant interactors

outDf <- calc_mod_ttest(outDf,iter=2000,two_sample=T,eBayes_trend=T)
outDf$significant <- outDf$logFC > 0 & outDf$FDR <= 0.1

cat('*** performed two-sample moderated t-test\n')

# ----------------------------------------------------------------------------------------
# (7) OUTPUT full results table (with InWeb annotation)

# annotate InWeb interactors and calculate overlap enrichment
inwebDf <- get_inweb_list(baitInWeb,type='all')
if (is.null(inwebDf)) {
	inwebInts <- ''
	inwebOverlap <- NULL
} else {
	inwebInts <- subset(inwebDf,significant)$gene
	inwebOverlap <- calc_hyper(outDf,inwebDf,
		data.frame(listName='InWeb',intersectN=T),bait=bait)
}
outDf$InWeb <- outDf$gene %in% inwebInts

# output results table 
output_stats <- fs::path(baitDir,paste0(ipName,'.Genoppi.Stats.txt'))
write.table(outDf,output_stats,quote=F,sep='\t',row.names=F)

cat('*** output Genoppi stats table:',paste0(ipName,'.Genoppi.Stats.txt\n'))

# ----------------------------------------------------------------------------------------
# (8) OUTPUT ggpairs plots of sample correlations (color: significance, imputation)

output_ggpairs <- fs::path(baitDir,paste0(ipName,'.Genoppi.SampleCorrelations.pdf'))
pdf(output_ggpairs,height=6,width=6)

# TRUE if protein has any imputed values across samples and controls
outDf$imp <- apply(outDf[paste(scCols,'_imp',sep='')],1,any)

# color by significance (plot non-sig points first)
ggDf <- outDf
if (sum(ggDf$gene==bait) > 0) {
	ggDf[ggDf$gene==bait,]$significant <- 'BAIT'
}
ggDf$significant <- factor(ggDf$significant,levels=c('FALSE','TRUE','BAIT'))

ggpairs(ggDf[order(ggDf$significant),],columns=scCols,
	lower = list(mapping=aes(color=significant),
		continuous=wrap('points',size=0.6,alpha=0.7)),
	title=paste(ipName,', color = significance',sep=''),
	xlab=bquote(log[2]*" protein intensity"),ylab=bquote(log[2]*" protein intensity")) +
	#scale_color_manual(values=c('BAIT'='red','TRUE'='dodgerblue','FALSE'='grey')) +
	scale_color_manual(values=c('grey','dodgerblue','red')) +
	theme_bw() + theme(axis.text=element_text(size=7))

# color by imputation status (plot non-imputed points first)
if (sum(ggDf$gene==bait) > 0) {
	ggDf[ggDf$gene==bait,]$imp <- 'BAIT'
}
ggDf$imp <- factor(ggDf$imp,levels=c('FALSE','TRUE','BAIT'))

ggpairs(ggDf[order(ggDf$imp),],columns=scCols,
	lower = list(mapping=aes(color=imp),
		continuous=wrap('points',size=0.6,alpha=0.7)),
	title=paste(ipName,', color = imputation',sep=''),
	xlab=bquote(log[2]*" protein intensity"),ylab=bquote(log[2]*" protein intensity")) +
	scale_color_manual(values=c('grey','orange','red')) +
	theme_bw() + theme(axis.text=element_text(size=7))

dev.off()

cat('*** output sample correlation plots:',
	paste0(ipName,'.Genoppi.SampleCorrelations.pdf\n'))

# ----------------------------------------------------------------------------------------
# (9) OUTPUT volcano plots (color: significance, InWeb overlay)

output_volcanos <- fs::path(baitDir,paste0(ipName,'.Genoppi.VolcanoPlots.pdf'))
pdf(output_volcanos,height=4,width=4)

ggplot(outDf,aes(x=logFC,y=-log10(pvalue))) +
xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"(P-value)")) +
geom_hline(yintercept=0,linetype='dashed') +
geom_vline(xintercept=0,linetype='dashed') +

# plot all proteins (dodgerblue = significant, grey = not significant)
geom_point(size=1.2,color=ifelse(outDf$significant,"dodgerblue","grey")) +

# label bait (red = significant, orange = not significant)
geom_point(subset(outDf,gene==bait & significant),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +
geom_point(subset(outDf,gene==bait & !significant),
	mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +
geom_point(subset(outDf,gene==bait),
	mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="black",shape=1) +

# text label for bait only
geom_text_repel(subset(outDf, gene==bait), mapping=aes(label=gene),
	arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.1, "lines"),
	point.padding=unit(0.15, "lines"), color="black", size=2) +

ggtitle(paste(ipName,': ',sum(outDf$significant),' sig/',nrow(outDf),' detected',sep='')) +
theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))

# plot InWeb overlay volcanos if bait found in InWeb
if (!is.null(inwebOverlap)) {
	inwebTitle <- paste('InWeb: ',inwebOverlap$statistic$successInSample_count,' sig/',
		inwebOverlap$statistic$sample_count,' detected, p=',
		signif(inwebOverlap$statistic$pvalue,3),sep='')

	# InWeb volcano without labeling InWeb ints
	p <- ggplot(outDf,aes(x=logFC,y=-log10(pvalue))) +
	xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"(P-value)")) +
	geom_hline(yintercept=0,linetype='dashed') +
	geom_vline(xintercept=0,linetype='dashed') +

	# plot all proteins (dodgerblue = significant, grey = not significant)
	geom_point(size=1.2,color=ifelse(outDf$significant,"dodgerblue","grey")) +

	# label InWeb interactors (yellow = significant, white = not significant)
	geom_point(subset(outDf, InWeb & significant),
		mapping=aes(x=logFC, y=-log10(pvalue)), size=1.2, color="yellow") +
	geom_point(subset(outDf, InWeb & !significant),
		mapping=aes(x=logFC, y=-log10(pvalue)), size=1.2, color="white") +
	geom_point(subset(outDf, InWeb),
		mapping=aes(x=logFC, y=-log10(pvalue)),size=1.2,color="black",shape=1) +

	# label bait (red = significant, orange = not significant)
	geom_point(subset(outDf,gene==bait & significant),
		mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +
	geom_point(subset(outDf,gene==bait & !significant),
		mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +
	geom_point(subset(outDf,gene==bait),
		mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="black",shape=1) +

	# text label for bait only
	geom_text_repel(subset(outDf, gene==bait), mapping=aes(label=gene),
		arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.1, "lines"),
		point.padding=unit(0.15, "lines"), color="black", size=2) +

	# title with InWeb overlap P
	ggtitle(inwebTitle) +
	theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))

	print(p)

	# InWeb volcano with text labels for InWeb ints
	p <- ggplot(outDf,aes(x=logFC,y=-log10(pvalue))) +
	xlab(bquote(log[2]*"(Fold change)")) + ylab(bquote(-log[10]*"(P-value)")) +
	geom_hline(yintercept=0,linetype='dashed') +
	geom_vline(xintercept=0,linetype='dashed') +

	# plot all proteins (dodgerblue = significant, grey = not significant)
	geom_point(size=1.2,color=ifelse(outDf$significant,"dodgerblue","grey")) +

	# label InWeb interactors (yellow = significant, white = not significant)
	geom_point(subset(outDf, InWeb & significant),
		mapping=aes(x=logFC, y=-log10(pvalue)), size=1.2, color="yellow") +
	geom_point(subset(outDf, InWeb & !significant),
		mapping=aes(x=logFC, y=-log10(pvalue)), size=1.2, color="white") +
	geom_point(subset(outDf, InWeb),
		mapping=aes(x=logFC, y=-log10(pvalue)),size=1.2,color="black",shape=1) +

	# label bait (red = significant, orange = not significant)
	geom_point(subset(outDf,gene==bait & significant),
		mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="red") +
	geom_point(subset(outDf,gene==bait & !significant),
		mapping=aes(x=logFC, y=-log10(pvalue)), size=2, color="orange") +
	geom_point(subset(outDf,gene==bait),
		mapping=aes(x=logFC,y=-log10(pvalue)),size=2,color="black",shape=1) +

	# text label for bait + InWeb ints
	geom_text_repel(subset(outDf, gene==bait | InWeb), mapping=aes(label=gene),
		arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.1, "lines"),
		point.padding=unit(0.15, "lines"), color="black", size=2) +

	# title with InWeb overlap P
	ggtitle(inwebTitle) +
	theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))

	print(p)
}

dev.off()

cat('*** output volcano plots:',paste0(ipName,'.Genoppi.VolcanoPlots.pdf\n'))
cat('*** GENOPPI SCRIPT COMPLETED\n')

