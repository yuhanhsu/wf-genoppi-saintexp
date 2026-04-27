##########################################################################################
# R script for processing and analyzing Taplin IP-MS data using SAINTexpress-spc
# to be run for each IP-MS sample set (bait IP + IgG control replicates)

# STEPS:
# (1) extract relevant columns (accession, gene name, uniq_peps, total peptide count)
# (2) remove proteins with no human accession/gene name, w/ <2 unique peptides,
#     or detected in <2 bait samples
# (3) OUTPUT log of filtered proteins
# (4) create intermediate files for input to SAINTexpress
#     Bait file: sample name, bait name, sample type (T/C)
#     Prey file: accession_number, gene
#     Interaction file: sample name, bait name, accession_number, total peptide count
#     (skip entries with NA count)
# (5) run SAINTexpress-spc via command line
# (6) OUTPUT reformatted results table (with InWeb annotation)
# (7) OUTPUT volcano plots (color: significance, InWeb overlay)

# Authors: Yu-Han Hsu, Samantha Cassity
# Last updated: 2026-04-24
##########################################################################################

rm(list=ls())

library(fs)
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(genoppi) # so we can use get_inweb_list

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

# extract total count from each sample or control, set missing values as NA
intDf <- msDf[c(paste('Total',samples),paste('Total',controls))]
intDf[intDf=='NF'] <- NA # replace NF (Not Found) with NA
intDf <- data.frame(apply(intDf,2,as.numeric))
intDf[intDf==0] <- NA # replace zero with NA
scCols <- c(paste0('sample',1:length(samples)),paste0('control',1:length(controls)))
names(intDf) <- scCols

intDf <- cbind(msDf[c('accession_number','gene','uniq_peps')],intDf)

cat('*** read in input data for',nrow(intDf),'proteins\n')
cat('*** extracted relevant columns:',names(intDf),'\n')

# ----------------------------------------------------------------------------------------
# (2) remove proteins with no human accession/gene name, w/ <2 unique peptides,
# or detected in <2 bait samples

humanFilterInds <- is.na(intDf$accession_number) | is.na(intDf$gene)
pepFilterInds <- intDf$uniq_peps < 2
sampFilterInds <- rowSums(!is.na(intDf[paste0('sample',1:length(samples))])) < 2

outDf <- intDf[!(humanFilterInds | pepFilterInds | sampFilterInds),]

cat('*** filtered out',sum(humanFilterInds | pepFilterInds | sampFilterInds),'proteins,',
	nrow(outDf),'proteins remaining\n')

# ----------------------------------------------------------------------------------------
# (3) OUTPUT log of filtered proteins

filteredDf <- data.frame(reference=c(msDf$reference[humanFilterInds],
                                     msDf$reference[pepFilterInds],msDf$reference[sampFilterInds]),
                         reason=c(rep('no human accession/gene name',sum(humanFilterInds)),
                                  rep('<2 unique peptides',sum(pepFilterInds)),
                                  rep('<2 bait IP samples',sum(sampFilterInds))))

filteredDf <- filteredDf %>% group_by(reference) %>% summarise(reason=toString(reason))

output_filter <- fs::path(baitDir,paste0(ipName,'.SAINTexpressSpc.FilteredProteins.txt'))
write.table(filteredDf,output_filter,quote=F,sep='\t',row.names=F)

cat('*** output filtered proteins:',
	paste0(ipName,'.SAINTexpressInt.FilteredProteins.txt\n'))

# ----------------------------------------------------------------------------------------
# (4) create intermediate files for input to SAINTexpress
#     Bait file: sample name, bait name, sample type (T/C)
#     Prey file: accession_number, gene
#     Interaction file: sample name, bait name, accession_number, total peptide count
#     (skip entries with NA count)

baitDf <- data.frame(Sample=scCols,
	# modify bait name so SAINTexpress would not ignore the bait entry in .Int file
	# use different control bait names (IgG_1,IgG_2,...) as recommended
	Bait=c(rep(paste0(bait,'_bait'),length(samples)),
		paste0('IgG_',1:length(controls))),
	Type=c(rep('T',length(samples)),rep('C',length(controls))))

output_bait <- fs::path(baitDir,paste0(ipName,'.SAINTexpressSpc.Bait'))
write.table(baitDf,output_bait,quote=F,sep='\t',row.names=F,col.names=F)

output_prey <- fs::path(baitDir,paste0(ipName,'.SAINTexpressSpc.Prey'))
write.table(outDf[c('accession_number','gene')],output_prey,
	quote=F,sep='\t',row.names=F,col.names=F)

interDf <- NULL
for (samp in scCols) {
	tempDf <- mutate(outDf,sample=samp,bait=baitDf$Bait[baitDf$Sample==samp]) %>%
		select(sample,bait,accession_number,{{samp}}) %>%
		na.omit %>% rename(total={{samp}})
	
	interDf <- rbind(interDf,tempDf)
}

output_int <- fs::path(baitDir,paste0(ipName,'.SAINTexpressSpc.Int'))
write.table(interDf,output_int,quote=F,sep='\t',row.names=F,col.names=F)

cat('*** created SAINTexpress input: *.Bait, *.Prey, *.Int\n')

# ----------------------------------------------------------------------------------------
# (5) run SAINTexpress-spc via command line

system(str_glue(
	'cd {baitDir} && ',
	'SAINTexpress-spc ',
	'{ipName}.SAINTexpressSpc.Int ',
	'{ipName}.SAINTexpressSpc.Prey ',
	'{ipName}.SAINTexpressSpc.Bait ',
	'&& mv list.txt {ipName}.SAINTexpressSpc.Stats.txt'
))

cat('*** ran SAINTexpress-spc\n')

# ----------------------------------------------------------------------------------------
# (6) OUTPUT reformatted results table (with InWeb annotation)

# read in the saintexpress stats file
output_stats <- fs::path(baitDir,paste0(ipName,'.SAINTexpressSpc.Stats.txt'))
SaintStats <- read.delim(output_stats) %>% rename(gene = "PreyGene")

# replace zero BFDR with 1/2 of next smallest BFDR to avoid plotting -log10(0)
min_nonzero <- min(SaintStats$BFDR[SaintStats$BFDR > 0]) # minimum non-zero BFDR

SaintStats$BFDR_plot <- ifelse(
	SaintStats$BFDR == 0, # if BFDR is zero:
	min_nonzero / 2, # BFDR_plot == 1/2 of min(BFDR)
	SaintStats$BFDR # else: BFDR_plot == BFDR 
)

# make a logFC column that takes log2 of the FoldChange column
SaintStats$logFC <- log2(SaintStats$FoldChange)

# make a significance column
SaintStats$significant <- SaintStats$BFDR <= 0.1

# annotate InWeb interactors and calculate overlap enrichment
inwebDf <- get_inweb_list(baitInWeb,type='all')
if (is.null(inwebDf)) {
	inwebInts <- ''
	inwebOverlap <- NULL
} else {
	inwebInts <- subset(inwebDf,significant)$gene
	inwebOverlap <- calc_hyper(SaintStats,inwebDf,
		data.frame(listName='InWeb',intersectN=T),bait=bait)
}

SaintStats$InWeb <- SaintStats$gene %in% inwebInts

# output results table 
write.table(SaintStats,output_stats,quote=F,sep='\t',row.names=F)

cat('*** output SAINTexpress-spc stats table:',
	paste0(ipName,'.SAINTexpressSpc.Stats.txt\n'))

# ----------------------------------------------------------------------------------------
# (7) OUTPUT volcano plots (color: significance, InWeb overlay)

output_volcanos <- fs::path(baitDir,paste0(ipName,'.SAINTexpressSpc.VolcanoPlots.pdf'))
pdf(output_volcanos,height=4,width=4)

p <- ggplot(SaintStats,aes(x=logFC,y=-log10(BFDR_plot))) +
  xlab(bquote(log[2]*"(FoldChange)")) + ylab(bquote(-log[10]*"(BFDR)")) +
  geom_hline(yintercept=0,linetype='dashed') +
  geom_vline(xintercept=0,linetype='dashed') +
  
  # plot all proteins (dodgerblue = significant, grey = not significant)
  geom_point(size=1.2,color=ifelse(SaintStats$significant,"dodgerblue","grey")) +
  
  # label bait (red = significant, orange = not significant)
  geom_point(subset(SaintStats,gene==bait & significant),
             mapping=aes(x=logFC,y=-log10(BFDR_plot)),size=2,color="red") +
  geom_point(subset(SaintStats,gene==bait & !significant),
             mapping=aes(x=logFC, y=-log10(BFDR_plot)), size=2, color="orange") +
  geom_point(subset(SaintStats,gene==bait),
             mapping=aes(x=logFC,y=-log10(BFDR_plot)),size=2,color="black",shape=1) +
  
  # text label for bait only
  geom_text_repel(subset(SaintStats, gene==bait), mapping=aes(label=gene),
                  arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.1, "lines"),
                  point.padding=unit(0.15, "lines"), color="black", size=2) +
  
  ggtitle(paste(ipName,': ',sum(SaintStats$significant),' sig/',nrow(SaintStats),' detected',sep='')) +
  theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))

print(p)

# plot InWeb overlay volcanos if bait found in InWeb
if (!is.null(inwebOverlap)) {
  inwebTitle <- paste('InWeb: ',inwebOverlap$statistic$successInSample_count,' sig/',
                      inwebOverlap$statistic$sample_count,' detected, p=',
                      signif(inwebOverlap$statistic$pvalue,3),sep='')
  
  # InWeb volcano without labeling InWeb ints
  p <- ggplot(SaintStats,aes(x=logFC,y=-log10(BFDR_plot))) +
    xlab(bquote(log[2]*"(FoldChange)")) + ylab(bquote(-log[10]*"(BFDR)")) +
    geom_hline(yintercept=0,linetype='dashed') +
    geom_vline(xintercept=0,linetype='dashed') +
    
    # plot all proteins (dodgerblue = significant, grey = not significant)
    geom_point(size=1.2,color=ifelse(SaintStats$significant,"dodgerblue","grey")) +
    
    # label InWeb interactors (yellow = significant, white = not significant)
    geom_point(subset(SaintStats, InWeb & significant),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)), size=1.2, color="yellow") +
    geom_point(subset(SaintStats, InWeb & !significant),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)), size=1.2, color="white") +
    geom_point(subset(SaintStats, InWeb),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)),size=1.2,color="black",shape=1) +
    
    # label bait (red = significant, orange = not significant)
    geom_point(subset(SaintStats,gene==bait & significant),
               mapping=aes(x=logFC,y=-log10(BFDR_plot)),size=2,color="red") +
    geom_point(subset(SaintStats,gene==bait & !significant),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)), size=2, color="orange") +
    geom_point(subset(SaintStats,gene==bait),
               mapping=aes(x=logFC,y=-log10(BFDR_plot)),size=2,color="black",shape=1) +
    
    # text label for bait only
    geom_text_repel(subset(SaintStats, gene==bait), mapping=aes(label=gene),
                    arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.1, "lines"),
                    point.padding=unit(0.15, "lines"), color="black", size=2) +
    
    # title with InWeb overlap P
    ggtitle(inwebTitle) +
    theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))
  
  print(p)
  
  # InWeb volcano with text labels for InWeb ints
  p <- ggplot(SaintStats,aes(x=logFC,y=-log10(BFDR_plot))) +
    xlab(bquote(log[2]*"(FoldChange)")) + ylab(bquote(-log[10]*"(BFDR)")) +
    geom_hline(yintercept=0,linetype='dashed') +
    geom_vline(xintercept=0,linetype='dashed') +
    
    # plot all proteins (dodgerblue = significant, grey = not significant)
    geom_point(size=1.2,color=ifelse(SaintStats$significant,"dodgerblue","grey")) +
    
    # label InWeb interactors (yellow = significant, white = not significant)
    geom_point(subset(SaintStats, InWeb & significant),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)), size=1.2, color="yellow") +
    geom_point(subset(SaintStats, InWeb & !significant),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)), size=1.2, color="white") +
    geom_point(subset(SaintStats, InWeb),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)),size=1.2,color="black",shape=1) +
    
    # label bait (red = significant, orange = not significant)
    geom_point(subset(SaintStats,gene==bait & significant),
               mapping=aes(x=logFC,y=-log10(BFDR_plot)),size=2,color="red") +
    geom_point(subset(SaintStats,gene==bait & !significant),
               mapping=aes(x=logFC, y=-log10(BFDR_plot)), size=2, color="orange") +
    geom_point(subset(SaintStats,gene==bait),
               mapping=aes(x=logFC,y=-log10(BFDR_plot)),size=2,color="black",shape=1) +
    
    # text label for bait + InWeb ints
    geom_text_repel(subset(SaintStats, gene==bait | InWeb), mapping=aes(label=gene),
                    arrow=arrow(length=unit(0.015, 'npc')), box.padding=unit(0.1, "lines"),
                    point.padding=unit(0.15, "lines"), color="black", size=2) +
    
    # title with InWeb overlap P
    ggtitle(inwebTitle) +
    theme_bw() + theme(axis.line=element_line(color="grey"), plot.title=element_text(size=10))
  
  print(p)
}

dev.off()

cat('*** output volcano plots:',paste0(ipName,'.SAINTexpressSpc.VolcanoPlots.pdf\n'))
cat('*** SAINTEXPRESS-SPC SCRIPT COMPLETED\n')
