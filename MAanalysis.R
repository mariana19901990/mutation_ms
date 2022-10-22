##Analysis of genetic mutations from MA-experiment
library(dplyr)
library(forcats)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(brms)
library(coda)
library(latex2exp)

#Source scripts required for mutation analysis
source("mutationscripts.R")

### * Load the data

### ** Load mutation data

##Load the curated mutation data
aineisto <- read.table("curated_mutations_final.csv", header = T, sep = ",")

### ** Load data about genomic features

### Duplicated regions

##Data about duplicated regions from Wang et al. 2021
dupregions <- read.csv("duplicated.csv", header = T)

#Consolidate duplicated regions, select only those detected by Blast   
dupreg <- consolidate.dupregions(dupregions, method = "Blast")        
sum(dupreg$End - (dupreg$Start - 1)) #Calculate length of all duplicated regions: 6 585 945 bp

dupcheck <- in.dup(aineisto, dupreg) #Check which mutations happened in duplicated regions
##Little over half of the mutations occur in duplicated regions (seems quite a lot)
## What is the number of called bases in duplicated regions?
aineisto$duplicatedreg <- dupcheck

### Centromeric regions

cent <- read.csv("centromeres.csv", header = T)

#Only the seven chromosomes
aineisto <- filter(aineisto, Chromosome != "Supercontig_12.8" & Chromosome != "Supercontig_12.9" & Chromosome != "Supercontig_12.10" & Chromosome != "Supercontig_12.11" & Chromosome != "Supercontig_12.12" & Chromosome != "Supercontig_12.13" & Chromosome != "Supercontig_12.14" & Chromosome != "Supercontig_12.15" & Chromosome != "Supercontig_12.17")

##For each mutation look whether it occurred in a centromeric region
aineisto$centromere <- in.centromeric(aineisto, cent)

sum(aineisto$centromere) #235 mutations happened in centromeric regions. This seems quite a lot


### Regions in H3K27me domains (duplicates removed from chip-seq)
h3k27 <- read.table("2489.H3K27.domains.duprm.bed", header = FALSE, sep = "\t")
colnames(h3k27) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k27[,2] <- h3k27[,2] + 1 #BED files have zero based start coordinate, need to fix this
h3k27 <- filter(h3k27, Chromosome != "mtDNA") #There are reads from mtDNA, but it should not have histones?, but drop it from the dataset
sum(h3k27$End - (h3k27$Start - 1)) #Lenght of all H3K27 domains: 4 538 600 bp

aineisto$H3K27 <- in.dup(aineisto, h3k27) #in.dup function should work for this
sum(aineisto$H3K27) # 128 mutations in H3K27

#H3K27 regions with overlapping H3K9 regions excluded
h3k27exk9 <- read.table("2489.H3K27_exK9.bed", header = FALSE, sep = "\t")
colnames(h3k27exk9) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k27exk9[,2] <- h3k27exk9[,2] + 1 #BED files have zero based start coordinate, need to fix this
sum(h3k27exk9$End - (h3k27exk9$Start - 1)) #Length of H3K27 domains with no K9 3 839 000 bp

aineisto$H3K27exK9 <- in.dup(aineisto, h3k27exk9)
sum(aineisto$H3K27exK9) # 79 mutations in H3K27 without K9


### Regions in H3K9me domtains (duplicates removed from chip-seq)
h3k9 <- read.table("2489.H3K9.domains.duprm.bed", header = FALSE, sep = "\t")
colnames(h3k9) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k9[,2] <- h3k9[,2] + 1 #BED files have zero based start coordinate, need to fix this
h3k9 <- filter(h3k9, Chromosome != "mtDNA") #There are reads from mtDNA, but it should not have histones?, but drop it from the dataset
sum(h3k9$End - (h3k9$Start - 1)) #Length of all H3K9 domains: 7 513 000 bp

aineisto$H3K9 <- in.dup(aineisto, h3k9) #in.dup function should work for this
sum(aineisto$H3K9) # 706 mutations in H3K9 (this is almost all of the mutations, seems a lot)


### Regions in H3K36me domains (duplicates removed from chip-seq)
h3k36 <- read.table("2489.H3K36.domains.duprm.bed", header = FALSE, sep = "\t")
colnames(h3k36) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k36[,2] <- h3k36[,2] + 1 #BED files have zero based start coordinate, need to fix this
h3k36 <- filter(h3k36, Chromosome != "mtDNA") #There are reads from mtDNA, but it should not have histones?, but drop it from the dataset
sum(h3k36$End - (h3k36$Start - 1)) #Length of all H3K36 domains: 32 307 600 bp

aineisto$H3K36 <- in.dup(aineisto, h3k36)
sum(aineisto$H3K36) #500 mutations in H3K36


### Euchromatic regions
euchromatin <- read.table("2489.euchromatin.bed", header = FALSE, sep = "\t")
colnames(euchromatin) <-  c("Chromosome", "Start", "End")
euchromatin[,2] <- euchromatin[,2] + 1 #BED files have zero based start coordinate, need to fix this
euchromatin <- filter(euchromatin, Chromosome != "mtDNA") #filter out mtDNA
sum(euchromatin$End - (euchromatin$Start - 1)) #Length 29 692 084 bp
euchromatin <- euchromatin[1:260,] #Only the seven chromosomes

aineisto$EUCHR <- in.dup(aineisto, euchromatin)
sum(aineisto$EUCHR) #435 mutations in euchromatin


### ** Load data for plotting all of the chromosomes, new mutations, and genomic features

#Neurospora chromosome sizes
chr.sizes <- data.frame(Chromosome = c("Supercontig_12.1", "Supercontig_12.2", "Supercontig_12.3", "Supercontig_12.4", "Supercontig_12.5", "Supercontig_12.6","Supercontig_12.7"), Start = c(1,1,1,1,1,1,1), End = c(9798893, 4478683, 5274802, 6000761, 6436246, 4218384, 4255303) )

labels <- c(Supercontig_12.1 = "Chr I", Supercontig_12.2 = "Chr II", Supercontig_12.3 = "Chr III", Supercontig_12.4 = "Chr IV", Supercontig_12.5 = "Chr V", Supercontig_12.6 = "Chr VI", Supercontig_12.7 = "Chr VII")

h3k27.plot <- filter(h3k27, Chromosome %in% chr.sizes[,1])
h3k9.plot <- filter(h3k9, Chromosome %in% chr.sizes[,1])
h3k36.plot <- filter(h3k36, Chromosome %in% chr.sizes[,1])
h3k27exk9.plot <- filter(h3k27exk9, Chromosome %in% chr.sizes[,1])

#H3K9 domains with centromeric regions excluded
h3k9.ex.cent.plot <- ex.centro(h3k9.plot, cent)[[1]]
ex.cent <- ex.centro(h3k9.plot, cent)[[2]]

#Write these to a BED file
#write.table(h3k9.ex.cent.plot, file= "~/Genomics/Neurospora/mutacc/h3k9.ex.cent.domains.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

aineisto$H3K9.ex.cent <- in.dup(aineisto, h3k9.ex.cent.plot) #in.dup function should work for this
sum(aineisto$H3K9.ex.cent) # 413 mutations H3K9 regions remain

aineisto$ex.cent <- in.dup(aineisto, ex.cent)
sum(aineisto$ex.cent) #290 mutations in excluded regions

### ** Load data for GC-content

###########################################################################
### Load the GC content data, GC data is 200 bp non-overlapping windows ###
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data.RData")
###########################################################################

#GC-content is in file ~/Genomics/Neurospora/mutacc/NC12_200bp_nucstats.txt
#For different window sizes, change the files
#nucstats <- read.table(file = "~/Genomics/Neurospora/mutacc/NC12_1000bp_nucstats.txt", sep = "\t", skip = 1)
#interval <- 200 ### CHANGE THIS!
#nucstats <- read.table(file = "~/Genomics/Neurospora/mutacc/NC12_200bp_nucstats.txt", sep = "\t", skip = 1)
#colnames(nucstats) <- c("Chromosome", "Start", "End", "ATcont", "GCcont", "nA", "nC", "nG", "nT", "nN", "nOther", "seqlen")
#nucstats.plot <- filter(nucstats, Chromosome %in% chr.sizes[,1])[,c(1:3,5)]
#nucstats.plot$Mid <- nucstats.plot$Start + interval/2

#Plotting GC-content across chromosomes, seems to correlate with H3K9 methylation
gc.chr.plot <- ggplot(nucstats.plot, aes()) +
    geom_line(aes(x = Mid, y = GCcont*100)) +
    geom_rect(data = h3k9.plot, aes(xmin = Start, xmax = End, ymax = -30, ymin = -55, fill = "H3K9me")) +
    geom_rect(data = cent, aes(xmin = start, xmax = end, ymin = 0, ymax = -25, fill = "Centromeric")) +
    xlab("Position (bp)") +
    ylab("GC content (%)") +
    scale_y_continuous(breaks = seq(0,100, 50)) +
    scale_fill_manual(breaks = c("Centromeric", "H3K9me"), values = c(Centromeric = "grey", H3K9me = "red")) +    
    facet_grid(Chromosome ~ ., labeller = labeller(Chromosome = labels)) +
    theme(legend.position = "top", legend.justification = c(0.5, 0.5), legend.title = element_blank())

save_plot("./mutac_ms/fig/gc_chr.pdf", gc.chr.plot, base_height = 6, base_width = 9*1.618)

#Distribution of whole genome
gc.genome <- ggplot(nucstats.plot, aes(GCcont*100)) +
    #geom_density() +
    geom_histogram(colour = "black", fill = "white", binwidth = 2.5) +
    xlab("GC content (%)") +
    ylab("") +    
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())   
        
     

##Which intervals overlap with chromatin modification domains?
##These were calculated with bedtools

###

### DO NOT RUN! ###################################################################################
### Save the data, so don't need to load it all the time
### Also can change the files to get different window sizes
h3k27.gc <- read.table("~/Genomics/Neurospora/mutacc/GC.and.H3K27.domains1000.bed", sep = "\t")
colnames(h3k27.gc) <- c("Chromosome", "Start", "End", "ATcont", "GCcont", "nA", "nC", "nG", "nT", "nN", "nOther", "seqlen")
h3k27.gc <- filter(h3k27.gc, Chromosome %in% chr.sizes[,1])[,c(1:3,5)]

h3k27exk9.gc <- read.table("~/Genomics/Neurospora/mutacc/GC.and.H3K27exK9.domains1000.bed", sep = "\t")
colnames(h3k27exk9.gc) <-  c("Chromosome", "Start", "End", "ATcont", "GCcont", "nA", "nC", "nG", "nT", "nN", "nOther", "seqlen")
h3k27exk9.gc <- filter(h3k27exk9.gc, Chromosome %in% chr.sizes[,1])[,c(1:3,5)]

h3k9.gc <- read.table("~/Genomics/Neurospora/mutacc/GC.and.H3K9.domains1000.bed", sep = "\t")
colnames(h3k9.gc) <- c("Chromosome", "Start", "End", "ATcont", "GCcont", "nA", "nC", "nG", "nT", "nN", "nOther", "seqlen")
h3k9.gc <- filter(h3k9.gc, Chromosome %in% chr.sizes[,1])[,c(1:3,5)]

h3k9.ex.cent.gc <- read.table("~/Genomics/Neurospora/mutacc/GC.and.H3K9.ex.cent.domains1000.bed", sep = "\t")
colnames(h3k9.ex.cent.gc) <- c("Chromosome", "Start", "End", "ATcont", "GCcont", "nA", "nC", "nG", "nT", "nN", "nOther", "seqlen")
h3k9.ex.cent.gc <- filter(h3k9.ex.cent.gc, Chromosome %in% chr.sizes[,1])[,c(1:3,5)]

cent.gc <- read.table("~/Genomics/Neurospora/mutacc/GC.and.centromeres1000.bed", sep = "\t")
colnames(cent.gc) <- c("Chromosome", "Start", "End", "ATcont", "GCcont", "nA", "nC", "nG", "nT", "nN", "nOther", "seqlen")
cent.gc <- filter(cent.gc, Chromosome %in% chr.sizes[,1])[,c(1:3,5)]

euchrom.gc <- read.table("~/Genomics/Neurospora/mutacc/GC.and.euchromatic1000.bed", sep = "\t")
colnames(euchrom.gc) <- c("Chromosome", "Start", "End", "ATcont", "GCcont", "nA", "nC", "nG", "nT", "nN", "nOther", "seqlen")
euchrom.gc <- filter(euchrom.gc, Chromosome %in% chr.sizes[,1])[,c(1:3,5)]

gc.data <- data.frame(Domain = c(rep("Euchromatic", nrow(euchrom.gc)), rep("H3K27", nrow(h3k27.gc)), rep("H3K9", nrow(h3k9.gc)), rep("Centromeric", nrow(cent.gc)), rep("H3K9 ex. centromeric", nrow(h3k9.ex.cent.gc)), rep("H3K27exK9", nrow(h3k27exk9.gc))), rbind(euchrom.gc, h3k27.gc, h3k9.gc, cent.gc, h3k9.ex.cent.gc, h3k27exk9.gc))

#save(gc.data, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data.RData")
#save(gc.data, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data400.RData")
#save(gc.data, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data600.RData")
#save(gc.data, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data800.RData")
#save(gc.data, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data1000.RData")
###########################################################################
### Load the GC content data, GC data is 200 bp non-overlapping windows ###
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data.RData")
###########################################################################

#Plotting GC-content in different domains
#Showing H3K27 domains w/ K9 excluded
gc.domains <- ggplot(filter(gc.data, Domain != "H3K27"), aes(x = Domain, y = GCcont*100)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0.1, outlier.shape = NA) +
    scale_x_discrete(labels = c("Centromeric", "Euchromatic",  "H3K27", "H3K9", "H3K9 ex. centromeric")) +
    ylab("GC-content (%)") +
    xlab("") +    
    coord_flip()

#library(ggnewscale)

#First draw rectangles for the chromosomes, centromeric regions, 
mutchr.plot <- ggplot(chr.sizes, aes()) +
     geom_rect( aes(xmin = Start, xmax = End, ymin = 0, ymax = 5, fill = "Chromosome"), alpha = 0.5) +
     geom_rect(data = cent, aes(xmin = start, xmax = end, ymin = 0, ymax = 5, fill = "Centromeric")) +
     geom_rect(data = h3k27.plot, aes(xmin = Start, xmax = End, ymax = -1, ymin = -3, fill = "H3K27me3")) +
     geom_rect(data = h3k9.plot, aes(xmin = Start, xmax = End, ymax = -4, ymin = -6, fill = "H3K9me")) +
     geom_rect(data = dupreg, aes(xmin = Start, xmax = End, ymax = -7, ymin = -9, fill = "Duplicated")) +
     geom_rect(data = h3k36.plot, aes(xmin = Start, xmax = End, ymax = -10, ymin = -13, fill = "H3K36me")) +
     #geom_rect(data = h3k9.ex.cent.plot, aes(xmin = Start, xmax = End, ymax = -10, ymin = -12, fill = "hotpink")) +
     geom_tile(data = aineisto, aes(x = Position, y = 2.5, height = 5, width = 2), colour = "black", alpha = 0.5) +
     xlab("Position (bp)") +
     ylab("") +
     scale_fill_manual(breaks = c("Chromosome", "Centromeric", "H3K27me3", "H3K9me", "Duplicated", "H3K36me"), values = c(Chromosome = "deepskyblue", Centromeric = "grey", H3K27me3 = "blue", H3K9me = "red", Duplicated = "green4", H3K36me = "purple")) +
     guides(fill = guide_legend(nrow = 1)) +
     facet_grid(Chromosome ~ ., labeller = labeller(Chromosome = labels)) +
     theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "top", legend.justification = c(0.5, 0.5), legend.title = element_blank())

save_plot("./mutac_ms/fig/mut_chromosomes.pdf", mutchr.plot, base_height = 6, base_width = 9*1.618)


##Next get the number of called bases in centromeric, duplicated, etc. regions

### ** Load data for called sites

##Load data about called sites

load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/calledsites.RData")

##### No need to run this again #####
calledsites.matA <- read.table(file = "~/Genomics/Neurospora/mutacc/called_sites_MA_matA.txt", header = F, sep = "\t", colClasses = c("character", "integer", "character"))
colnames(calledsites.matA) <- c("Chromosome", "Position", "Base")

#Format the ref base column
calledsites.matA$Base <- substr(calledsites.matA$Base, 1, 1)

calledsites.mata <- read.table(file = "~/Genomics/Neurospora/mutacc/called_sites_MA_mata.txt", header = F, sep = "\t", colClasses = c("character", "integer", "character"))
colnames(calledsites.mata) <- c("Chromosome", "Position", "Base")
#Format the ref base column
calledsites.mata$Base <- substr(calledsites.mata$Base, 1, 1)

called.centromeric.matA <- called.sites(calledsites.matA, cent, type = "cent")
called.duplicated.matA <- called.sites(calledsites.matA, dupreg, type = "dup")
called.H3K27.matA <- called.sites(calledsites.matA, h3k27, type = "dup") #Type "dup" should be OK
called.H3K9.matA <- called.sites(calledsites.matA, h3k9, type = "dup")
called.chr.matA <-  called.sites(calledsites.matA, chr.sizes, type = "cent")
called.H3K9.exc.matA <- called.sites(calledsites.matA, h3k9.ex.cent.plot, type = "dup")
called.exc.matA <- called.sites(calledsites.matA, ex.cent, type = "dup")

called.H3K27.exK9.matA <- called.sites(calledsites.matA, h3k27exk9, type = "dup")
called.euchr.matA <- called.sites(calledsites.matA, euchromatin, type = "dup")

called.centromeric.mata <- called.sites(calledsites.mata, cent, type = "cent")
called.duplicated.mata <- called.sites(calledsites.mata, dupreg, type = "dup")
called.H3K27.mata <- called.sites(calledsites.mata, h3k27, type = "dup") #Type "dup" should be OK
called.H3K9.mata <- called.sites(calledsites.mata, h3k9, type = "dup")
called.chr.mata <- called.sites(calledsites.mata, chr.sizes, type = "cent")
called.H3K9.exc.mata <- called.sites(calledsites.mata, h3k9.ex.cent.plot, type = "dup")
called.exc.mata <- called.sites(calledsites.mata, ex.cent, type = "dup")

called.H3K27.exK9.mata <- called.sites(calledsites.mata, h3k27exk9, type = "dup")
called.euchr.mata <- called.sites(calledsites.mata, euchromatin, type = "dup")

##Since calculating the called sites takes abit long, save the results
##############################################
#save(called.centromeric.matA, called.duplicated.matA, called.H3K27.matA, called.H3K9.matA, called.chr.matA, called.centromeric.mata, called.duplicated.mata, called.H3K27.mata, called.H3K9.mata, called.chr.mata, called.H3K9.exc.matA, called.exc.matA, called.H3K9.exc.mata, called.exc.mata, called.H3K27.exK9.matA, called.euchr.matA, called.H3K27.exK9.mata, called.euchr.mata, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/calledsites.RData")

##############################################


### * Summary statistics

#Summary of mutations that occurred
aineisto.lines <- group_by(aineisto, Line)
mutperline <- summarise(aineisto.lines, nmut = n()) #Mutations per line

med.lines <- median(unlist(mutperline[,2])) #Median of 33 mutations per line

ggplot(mutperline, aes(x = nmut)) +
    geom_histogram(fill = "grey", colour = "black", binwidth = 3) +
    scale_y_continuous(expand = c(0,0), breaks = seq(0,10,1)) +
    scale_x_continuous(breaks = seq(15,105,15)) +
    xlab("Number of mutations") +
    ylab("Number of lines") +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14))     

aineisto.mutypes <- group_by(aineisto, type)
mutpertype <- summarise(aineisto.mutypes, nmut = n()) #Mutations per type

#save(mutpertype, med.lines, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/summary.RData")

#Number of mitoses
transfers <- 40
load("./data/mitoses.RData") #Number of mitoses that happened during the experiment
mitoses.transfer <- quantile(mitoses, probs = c(0.5, 0.025, 0.975)) #25.4 mitoses per transfer
mitoses.exp <- quantile(mitoses*transfers, probs = c(0.5, 0.025, 0.975))

#save(mitoses.transfer, mitoses.exp, file = "./mutac_ms/data/mitosesexp.RData")


### * Draw figure 1
##Check cowplot tutorial and stack overflow

library(magick)

MAexp <- ggdraw() + draw_image(magick::image_read_pdf("~/Documents/tutkijatohtori/epimutation/MA_WGS/MA_exp.pdf", density=600))

Muts <- ggplot(mutperline, aes(x = nmut)) +
    geom_histogram(fill = "grey", colour = "black", binwidth = 3) +
    scale_y_continuous(expand = c(0,0), breaks = seq(0,10,1)) +
    scale_x_continuous(breaks = seq(15,105,15)) +
    xlab("Number of mutations") +
    ylab("Number of lines") +
    theme(text=element_text(size=12), axis.text=element_text(size=12))

#Muts <- ggdraw() + draw_image(magick::image_read_pdf("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutperline.pdf", density=600))

transferplot <- ggdraw() +
draw_image(magick::image_read_pdf("~/Documents/tutkijatohtori/epimutation/MA_WGS/transfer.pdf", density=600))

nucleiplot <- ggdraw() +
draw_image(magick::image_read_pdf("~/Documents/tutkijatohtori/epimutation/MA_WGS/nuclei3.pdf", density=600))


plot_grid(MAexp, Muts, labels = c("A", "B"))

exp.plot <- plot_grid(MAexp, transferplot, labels = c("A", "B"), nrow = 2, rel_heights = c(1, 0.5))

col2plot <- plot_grid(nucleiplot, Muts, labels = c("C", "D"), nrow = 2)

finalplot <- plot_grid(exp.plot, col2plot, rel_widths = c(1, 0.8)) + theme(plot.background = element_rect(fill = "white", colour = "white"))

save_plot(filename = "./mutac_ms/fig/experiment.jpg", base_height = 5, plot = finalplot)



ggplot(linesG40, aes(x = meangr)) +
    geom_histogram(fill = "grey", colour = "black") + #, binwidth = 3) +
    geom_vline(xintercept = c(3.23 - 2*0.0137, 3.23, 3.23 + 2*0.0137), lty = c(2,1,2)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(0,13,1)) + #, breaks = seq(0,10,1)) +
    #scale_x_continuous(breaks = seq(015,105,15)) +
    xlab("Growth rate (mm / h)") +
    ylab("Number of lines") +
    #facet_wrap( ~ mat) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14))    

### * Statistical analysis

### *** All mutations in total

transfers <- 40
load("./mutac_ms/data/mitoses.RData")


model.all <- brm(data = mutperline, family = poisson,
                 nmut ~ 1,
                 prior = c(prior(normal(0, 10), class = Intercept)),
                 iter = 3000, warmup = 1000, chains = 4, cores = 4)

##Extract posterior and calculate mutation rate
post <- posterior_samples(model.all)
mutation.rate.all <- exp(post[,1])/(mitoses*transfers)
rate.all <- c(median(mutation.rate.all), HPDinterval(as.mcmc(mutation.rate.all))) #median, lower, upper

### *** All point mutations

transfers <- 40
load("./mutac_ms/data/mitoses.RData")

#Number of called sites: 40561734 #mat A
#Number of called sites: 40601276 #mat a

called.bases <- (40561734 + 40601276)/2 #Number of sites called in the ancestors

#aineisto.lines <- group_by(aineisto, Line)
#mutperline <- summarise(aineisto.lines, nmut = n()) #Mutations per line

#aineisto.mutations <- group_by(aineisto, Line, mutation, type)
#mutpertype <- summarise(aineisto.mutations, nmut = n()) #Mutations per type

aineisto.point <- filter(aineisto, type == "point")
aineisto.point <- group_by(aineisto.point, Line)
aineisto.point <- summarise(aineisto.point, nmut = n())
#aineisto.point <- filter(aineisto.point, mutation == "point")

model.point <- brm(data = aineisto.point, family = poisson,
                  nmut ~ 1,
                  prior = c(prior(normal(0, 10), class = Intercept)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

##Extract posterior and calculate mutation rate
post <- posterior_samples(model.point)

#Calculate mutation rate
mutation.rate.point <- exp(post[,1])/(mitoses*transfers*called.bases)
#get estimates
rate.point <- c(median(mutation.rate.point), HPDinterval(as.mcmc(mutation.rate.point)))  #median estimate, 2.05*10^-8  #HPD interval [1.92*10^-8, 2.19*10^-8]

### *** Rate of transversions and transitions

#Calculate transition/transversion bias
aineisto.tstv <- filter(aineisto, type == "point")
aineisto.tstv <- group_by(aineisto.tstv, Line, point.type)
aineisto.tstv <- summarise(aineisto.tstv, nmut = n())

#aineisto.transition <- filter(mutpertype, type == "transition")
#aineisto.transversion <- filter(mutpertype, type == "transversion")

model.tstv <- brm(data = aineisto.tstv, family = poisson,
                  nmut ~ -1 + point.type,
                  #prior = c(prior(normal(0, 10), class = Intercept)), #Change prior to b
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.tstv)
transition.rate <- exp(post[,1])/(transfers) #Transition rate per transfer
transversion.rate <- exp(post[,2])/(transfers)

#Calculating transition / transversion bias, kappa
kappa <- transition.rate/transversion.rate
median(kappa) #1.24
HPDinterval(as.mcmc(kappa)) #[1.07, 1.40]
ts.tv.rate <- c(median(kappa), HPDinterval(as.mcmc(kappa)))

#Taking the number of paths into account
koe <- transition.rate/(transversion.rate/2)

ts.rate <- c(median(transition.rate), HPDinterval(as.mcmc(transition.rate)))
tv.rate <- c(median(transversion.rate), HPDinterval(as.mcmc(transversion.rate)))





### *** Each mutation separately

#Mutation spectra

#In mat A
#Number of called A's: 10488800
#Number of called T's: 10491164
#Number of called C's: 9787353
#Number of called G's: 9794417

#In mat a
#Number of called A's: 10498173
#Number of called T's: 10499371
#Number of called C's: 9798851
#Number of called G's: 9804881

As <- (10488800 + 10498173)/2
Ts <- (10491164 + 10499371)/2
Cs <- (9787353 + 9798851)/2
Gs <- (9794417 + 9804881)/2

#As <- 10492539
#Ts <- 10494773
#Cs <- 9791134
#Gs <- 9798146
total <- As + Ts + Cs + Gs
GCcontent <- (Cs + Gs) / total
basefreqs <- c(As, Ts, Cs, Gs)/total
CGs <- Cs + Gs
ATs <- As + Ts

#Total lenght of the reference genome
reftot <- sum(9798893, 4478683, 5274802, 6000761, 6436246, 4218384, 4255303, 192308,142473, 125404, 31696, 19714, 13515, 11565, 9397, 8983, 6701, 6309, 4755, 1646, 64840, 6548)

aineisto.spectra <- filter(aineisto, type == "point")

##Assuming that all mutations are equally likely (1/6)
#possible mutations are:
#transversions
#A:T -> T:A
#C:G -> G:C
#C:G -> A:T
#A:T -> C:G
#transitions
#A:T -> G:C
#C:G -> T:A
mutspectra <- data.frame(mutation = factor(c("A:T → T:A", "C:G → G:C", "C:G → A:T", "A:T → C:G", "A:T → G:C", "C:G → T:A"), levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A")) , type = factor(c(rep("transversion", 4), rep("transition", 2))))
expected <- c((1-GCcontent)*(1/3), GCcontent*(1/3), GCcontent*(1/3), (1-GCcontent)*(1/3), (1-GCcontent)*(1/3), GCcontent*(1/3))
mutspectra$expected <- expected
observed <- rep(0, 6)
observed[1] <- sum( (aineisto.spectra$anc.base == "A" & aineisto.spectra$sample.base == "T") | (aineisto.spectra$anc.base == "T" & aineisto.spectra$sample.base == "A") )
observed[2] <- sum( (aineisto.spectra$anc.base == "C" & aineisto.spectra$sample.base == "G") | (aineisto.spectra$anc.base == "G" & aineisto.spectra$sample.base == "C") )
observed[3] <- sum( (aineisto.spectra$anc.base == "C" & aineisto.spectra$sample.base == "A") | (aineisto.spectra$anc.base == "G" & aineisto.spectra$sample.base == "T") )
observed[4] <- sum( (aineisto.spectra$anc.base == "A" & aineisto.spectra$sample.base == "C") | (aineisto.spectra$anc.base == "T" & aineisto.spectra$sample.base == "G") )
observed[5] <- sum( (aineisto.spectra$anc.base == "A" & aineisto.spectra$sample.base == "G") | (aineisto.spectra$anc.base == "T" & aineisto.spectra$sample.base == "C") )
observed[6] <- sum( (aineisto.spectra$anc.base == "C" & aineisto.spectra$sample.base == "T") | (aineisto.spectra$anc.base == "G" & aineisto.spectra$sample.base == "A") )
mutspectra$observed <- observed
num.exp <- expected*sum(observed)
mutspectra$rel.obs <- observed/num.exp

#Making datasets for each possible point mutations
#AT -> TA
aineisto.indmut <- filter(aineisto.spectra, (anc.base == "A" & sample.base == "T") | (anc.base == "T" & sample.base == "A") )
aineisto.indmut <- group_by(aineisto.indmut, Line, .drop = F) #Need to keep lines with 0 mutations
aineisto.AT.TA <- summarise(aineisto.indmut, nmut = n())

#
aineisto.indmut <- filter(aineisto.spectra, (anc.base == "C" & sample.base == "G") | (anc.base == "G" & sample.base == "C") )
aineisto.indmut <- group_by(aineisto.indmut, Line, .drop = F)
aineisto.CG.GC <- summarise(aineisto.indmut, nmut = n())

aineisto.indmut <- filter(aineisto.spectra, (anc.base == "C" & sample.base == "A") | (anc.base == "G" & sample.base == "T") )
aineisto.indmut <- group_by(aineisto.indmut, Line, .drop = F)
aineisto.CG.AT <- summarise(aineisto.indmut, nmut = n())

aineisto.indmut <- filter(aineisto.spectra, (anc.base == "A" & sample.base == "C") | (anc.base == "T" & sample.base == "G") )
aineisto.indmut <- group_by(aineisto.indmut, Line, .drop = F)
aineisto.AT.CG <- summarise(aineisto.indmut, nmut = n())

aineisto.indmut <- filter(aineisto.spectra, (anc.base == "A" & sample.base == "G") | (anc.base == "T" & sample.base == "C") )
aineisto.indmut <- group_by(aineisto.indmut, Line, .drop = F)
aineisto.AT.GC <- summarise(aineisto.indmut, nmut = n())

aineisto.indmut <- filter(aineisto.spectra, (anc.base == "C" & sample.base == "T") | (anc.base == "G" & sample.base == "A") )
aineisto.indmut <- group_by(aineisto.indmut, Line, .drop = F)
aineisto.CG.TA <- summarise(aineisto.indmut, nmut = n())

aineisto.spectra.mutations <- rbind(aineisto.AT.TA, aineisto.AT.CG, aineisto.CG.GC, aineisto.CG.AT, aineisto.AT.GC, aineisto.CG.TA)
aineisto.spectra.mutations$mutation <- factor(c(rep("A:T → T:A", 39), rep("A:T → C:G", 39), rep("C:G → G:C", 39), rep("C:G → A:T", 39),  rep("A:T → G:C", 39), rep("C:G → T:A", 39)), levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"))

model.spectra <- brm(data = aineisto.spectra.mutations, family = poisson,
                  nmut ~ -1 + mutation,
                  #prior = c(prior(normal(0, 10), class = Intercept)), #Change prior to b
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)
post <- posterior_samples(model.spectra)
rate.AT.TA <- exp(post[,1])/ATs
rate.AT.CG <- exp(post[,2])/ATs
rate.CG.GC <- exp(post[,3])/CGs
rate.CG.AT <- exp(post[,4])/CGs
rate.AT.GC <- exp(post[,5])/ATs
rate.CG.TA <- exp(post[,6])/CGs

### Compare some mutation rates
comp1 <- rate.AT.GC/rate.CG.TA

rate.comps <- data.frame(matrix(rep(0, 4*1), ncol = 4))
colnames(rate.comps) <- c("comparison", "median", "lower", "upper")
rate.comps[1,1] <- "AT->GC / CG->TA"

rate.comps[1,2:4] <- c(median(comp1), HPDinterval(as.mcmc(comp1)))

#Calculating the mutation spectra
mutspectra$rate.med <- c(median(rate.AT.TA), median(rate.CG.GC), median(rate.CG.AT), median(rate.AT.CG), median(rate.AT.GC), median(rate.CG.TA))

#expected rate (i.e. weighted average of all rates)
expected.rate <- (rate.AT.TA*ATs + rate.AT.CG*ATs + rate.CG.GC*CGs + rate.CG.AT*CGs + rate.AT.GC*ATs + rate.CG.TA*CGs)/(ATs*3 + CGs*3)

#Calculate relative rates
rel.rate.AT.TA <- rate.AT.TA / expected.rate
rel.rate.AT.CG <- rate.AT.CG / expected.rate
rel.rate.CG.GC <- rate.CG.GC / expected.rate
rel.rate.CG.AT <- rate.CG.AT / expected.rate
rel.rate.AT.GC <- rate.AT.GC / expected.rate
rel.rate.CG.TA <- rate.CG.TA / expected.rate

mutspectra$rel.rate.med <- c(median(rel.rate.AT.TA), median(rel.rate.CG.GC), median(rel.rate.CG.AT), median(rel.rate.AT.CG), median(rel.rate.AT.GC), median(rel.rate.CG.TA))
mutspectra$rel.rate.low <- c(HPDinterval(as.mcmc(rel.rate.AT.TA))[1], HPDinterval(as.mcmc(rel.rate.CG.GC))[1], HPDinterval(as.mcmc(rel.rate.CG.AT))[1], HPDinterval(as.mcmc(rel.rate.AT.CG))[1], HPDinterval(as.mcmc(rel.rate.AT.GC))[1], HPDinterval(as.mcmc(rel.rate.CG.TA))[1])
mutspectra$rel.rate.high <- c(HPDinterval(as.mcmc(rel.rate.AT.TA))[2], HPDinterval(as.mcmc(rel.rate.CG.GC))[2], HPDinterval(as.mcmc(rel.rate.CG.AT))[2], HPDinterval(as.mcmc(rel.rate.AT.CG))[2], HPDinterval(as.mcmc(rel.rate.AT.GC))[2], HPDinterval(as.mcmc(rel.rate.CG.TA))[2]) 

### Making plots for single nucleotide results ##############################################3

library(latex2exp)
mylabels <- c(TeX("A:T $\\rightarrow$ T:A"), TeX("A:T $\\rightarrow$ C:G"), TeX("C:G $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ A:T"), TeX("A:T $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ T:A") )

#library(extrafont)
#loadfonts()

#Note! Need to use cairo pdf for plotting unicode characters
grDevices::cairo_pdf("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/fig/mutspectra.pdf")
ggplot(mutspectra, aes(y = rel.rate.med, ymin = rel.rate.low, ymax = rel.rate.high, x = mutation, fill = type)) +
    geom_bar(stat = "identity", colour = "black") +
    geom_errorbar(width = 0.1) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2), breaks = c(seq(0,2, by = 0.25))) +
    scale_x_discrete(labels = mylabels) +
    ylab("Relative mutation rate") +
    xlab("Mutation") +
    geom_hline(yintercept = 1, lty = "dashed") +
    theme(legend.position = "none")
dev.off()

#Note! Need to use cairo pdf for plotting unicode characters
grDevices::cairo_pdf("mutspectraFI.pdf")
ggplot(mutspectra, aes(y = rel.obs, x = mutation, fill = type)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_y_continuous(expand = c(0,0), limits = c(0,2), breaks = c(seq(0,2, by = 0.25))) +
    scale_x_discrete(labels = mylabels) +
    ylab("Havaittu / Odotettu") +
    xlab("Mutaatio") +
    geom_hline(yintercept = 1, lty = "dashed") +
    theme(legend.position = "none")
dev.off()

### Save point mutation results #################################################################
save(rate.all, rate.point, ts.tv.rate, rate.comps, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/mutation_rates.RData")
#################################################################################################

####################################################################################################

### *** Effect of local base composition

###Are mutations equally common in all possible triplets?
#The possible (32) triplets are:
# if mutation is A:T base pair
# AAA, TAA, TAT, AAT, CAA, AAC, CAC, GAA, AAG, GAG, TAC, CAT, TAG, GAT, CAG, GAC
# and their complementary codons
# TTT, ATT, ATA, TTA, GTT, TTG, GTG, CTT, TTC, CTC, ATG, GTA, ATC, CTA, GTC, CTG
# These codons are pooled

# if mutation is in C:G base pair
# CCC, GCC, GCG, CCG, ACC, CCA, ACA, TCC, CCT, TCT, GCA, ACG, TCA, ACT, TCG, GCT
# and their complementary codons
# GGG, CGG, CGC, GGC, TGG, GGT, TGT, AGG, GGA, AGA, CGT, TGC, AGT, TGA, AGC, CGA
# These codons are pooled

##First count the codons where mutations occurred. Then use observed (called) base frequencies to
##determine whether mutations occur as they are expected

#Check bases of coordinates -1 and +1 where mutations occurred

#Retrieving the triplets
#Using the aineisto.spectra dataset from previous step

### This step requires that there is access to the reference genome
### DO NOT RUN, but load the data below instead ##################
triplets <- rep(0, nrow(aineisto.spectra))
for(i in 1:nrow(aineisto.spectra)) {
    #Get triplet
    triplets[i] <- retrieve.triplet(aineisto.spectra$Chromosome[i], aineisto.spectra$Position[i])
    #
    ## Perform some QC ##
    if( substr(triplets[i], 2, 2) != aineisto.spectra$anc.base[i] ) {
        print(paste("Warning, check mutation at", aineisto.spectra$Chromosome[i], aineisto.spectra$Position[i], "ancestor base and reference don't match!", sep = " "))
    }
}

aineisto.spectra$triplet <- triplets
################################################################

#Since triplet retrieving requires that the large genome files are accessible, it is good idea to save, so that can continue even if on a different computer
#save(aineisto.spectra, file = "./mutac_ms/data/mutspectra.RData")
load("./mutac_ms/data/mutspectra.RData") #Load the data about point mutations, triplets included
#triplets <- aineisto.spectra$triplet
#save(aineisto.spectra, file = "./mutac_ms/data/mutspectra_domains.RData")
load("./mutac_ms/data/mutspectra_domains.RData")

#triplets <- aineisto.spectra$triplet

table(triplets) #Show table of all triplets and mutations

##Load trinucleotide frequencies
trinucfreq.wg <- read.table("Nc_trinuc_freq_whole_genome.txt", header = F)
trinucfreq.h3k9 <- read.table("Nc_trinuc_freq_H3K9.txt", header = F)
trinucfreq.eu <- read.table("Nc_trinuc_freq_euchromatin.txt", header = F)
trinucfreq.cent <- read.table("Nc_trinuc_freq_centromeres.txt", header = F)
trinucfreq.h3k27 <- read.table("Nc_trinuc_freq_H3K27_exK9.txt", header = F)

### **** Analysis of triplets mutations rate in different domains taking triplets into account

trip.eu <- trinucleotides(filter(aineisto.spectra, EUCHR == 1)$triplet)
trip.eu <- sort.trinucleotides(trinucfreq.eu, trip.eu)
trip.eu$EUCHR <- rep(1, 32)
trip.eu$cent <- rep(0, 32)
trip.eu$H3K9 <- rep(0,32)
trip.eu$H3K27 <- rep(0,32)

##Building the H3K9 dataset (need to subtract centromere number of mutations and trinuc counts from H3K9 counts)
trip.cent <- trinucleotides(filter(aineisto.spectra, centromere == 1)$triplet)
trip.cent <- sort.trinucleotides(trinucfreq.cent, trip.cent)
trip.cent$EUCHR <- rep(0, 32)
trip.cent$cent <- rep(1, 32)
trip.cent$H3K9 <- rep(1,32)
trip.cent$H3K27 <- rep(0,32)

trip.h3k9 <- trinucleotides(filter(aineisto.spectra, H3K9 == 1)$triplet)
trip.h3k9 <- sort.trinucleotides(trinucfreq.h3k9, trip.h3k9)
trip.h3k9$nmut <- trip.h3k9$nmut - trip.cent$nmut #Subtracting centromere counts and mutations
trip.h3k9$count <- trip.h3k9$count - trip.cent$count
trip.h3k9$EUCHR <- rep(0, 32)
trip.h3k9$cent <- rep(0, 32)
trip.h3k9$H3K9 <- rep(1,32)
trip.h3k9$H3K27 <- rep(0,32)

trip.h3k27 <- trinucleotides(filter(aineisto.spectra, H3K27exK9 == 1)$triplet)
trip.h3k27 <- sort.trinucleotides(trinucfreq.h3k27, trip.h3k27)
trip.h3k27$EUCHR <- rep(0, 32)
trip.h3k27$cent <- rep(0, 32)
trip.h3k27$H3K9 <- rep(0,32)
trip.h3k27$H3K27 <- rep(1,32)

#Combine all
trip.all <- rbind(trip.eu, trip.cent, trip.h3k9, trip.h3k27)
#Calculate the expected number of mutations provided that everything has equal rate
countfreq <- trip.all$count / sum(trip.all$count)
trip.all$expected <- countfreq*sum(trip.all$nmut)

## Maybe I should fit the trinucleotide effect as a random effect instead?

##Doing some model comparison
model.all.trip1 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip2 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + H3K9,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip3 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + H3K9 + class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip4 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + H3K9 + class + H3K9:class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip5 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip6 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip7 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip8 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + class + H3K9:class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip9 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + class + H3K9:class + cent:class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip10 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + H3K27 + class + H3K9:class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip11 <- brm(data = trip.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + H3K27 + class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip0 <- brm(data = trip.all, family = poisson,
                        nmut ~ 1 + offset(log(expected)),
                        prior = c(prior(normal(0, 10), class = Intercept)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

##Note: we did not attempt to fit H3K27:class interactions, as H3K27 has 0 mutations in so many triplets, there is not enough to estimate interactions.

#countfreq <- trip.all$count / sum(trip.all$count)
#trip.all$expected <- countfreq*sum(trip.all$nmut)

#sum(trip.all$nmut)

#model.all.trip12 <- brm(data = trip.all, family = poisson,
#                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + H3K27 + class,
#                        prior = c(prior(normal(0, 10), class = b)),
#                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

### Model comparisons  ###############################
#using WAIC instead
model.all.trip0 <- add_criterion(model.all.trip0, "waic")
model.all.trip1 <- add_criterion(model.all.trip1, "waic")
model.all.trip2 <- add_criterion(model.all.trip2, "waic")
model.all.trip3 <- add_criterion(model.all.trip3, "waic")
model.all.trip4 <- add_criterion(model.all.trip4, "waic")
model.all.trip5 <- add_criterion(model.all.trip5, "waic")
model.all.trip6 <- add_criterion(model.all.trip6, "waic")
model.all.trip7 <- add_criterion(model.all.trip7, "waic")
model.all.trip8 <- add_criterion(model.all.trip8, "waic")
model.all.trip9 <- add_criterion(model.all.trip9, "waic")
model.all.trip10 <- add_criterion(model.all.trip10, "waic")
model.all.trip11 <- add_criterion(model.all.trip11, "waic")


mcomp.trip <- loo_compare(model.all.trip0, model.all.trip1, model.all.trip2, model.all.trip3, model.all.trip4, model.all.trip5, model.all.trip6, model.all.trip7, model.all.trip8, model.all.trip9, model.all.trip10, model.all.trip11, criterion = "waic")

print(mcomp.trip, simplify = F)

mweights.trip <- model_weights(model.all.trip0, model.all.trip1, model.all.trip2, model.all.trip3, model.all.trip4, model.all.trip5, model.all.trip6, model.all.trip7, model.all.trip8, model.all.trip9, model.all.trip10, model.all.trip11, weights = "waic")

tripmodelcomp <- cbind(mcomp.trip, sort(mweights.trip, decreasing = T))
colnames(tripmodelcomp)[9] <- "weight"

###

######################################################
tripmodelres <- fixef(model.all.trip11)

mutation.rsquare.res <- bayes_R2(model.all.trip11)

#save(tripmodelcomp, tripmodelres, mutation.rsquare.res, file = "./mutac_ms/data/tripmodelcomp.RData")
#load("./mutac_ms/data/tripmodelcomp.RData")

##Getting posterior predictions from the model
post.trip <- posterior_samples(model.all.trip11)

#save(post.trip, file = "./mutac_ms/data/trip_posterior.RData")
#Can load the posterior distribution from above...
#load("./mutac_ms/data/trip_posterior.RData")

##Results for domains
cent.trip <- post.trip[,1] + post.trip[,2] #Mutation rate in centromeric regions is the sum of centromere effect and H3K9 effect
trip.domain.res <- apply(post.trip[,2:3], 2, exp)
trip.domain.res <- cbind(exp(cent.trip), trip.domain.res)
colnames(trip.domain.res) <- c("centromere", "H3K9", "H3K27")

trip.results.domains <- data.frame(estimate = apply(trip.domain.res,2,median), HPDinterval(as.mcmc(trip.domain.res)))
trip.results.domains$domain <- rownames(trip.results.domains)

trip.res <- apply(post.trip[,4:35], 2, exp) #Convert back to normal scale
colnames(trip.res) <- gsub("b_class", "", colnames(trip.res))
trip.res <- trip.res / median(trip.res) #Standardize effects of mutation rate (because there was no intercept) to relative mutation rates
trip.results <- data.frame(estimate = apply(trip.res,2,median), HPDinterval(as.mcmc(trip.res)))
trip.results$class <- factor(rownames(trip.results))

#This is the order of trinucleotides I want them in plots
#fct_relevel(class, "TCT:AGA", "TCC:GGA", "TCG:CGA", "TCA:TGA", "CCT:AGG", "CCC:GGG", "CCG:CGG", "CCA:TGG", "GCT:AGC", "GCC:GGC", "GCG:CGC", "GCA:TGC",  "ACT:AGT", "ACC:GGT", "ACG:CGT", "ACA:TGT", "TAT:ATA", "TAC:GTA", "TAG:CTA", "TAA:TTA", "CAT:ATG", "CAC:GTG", "CAG:CTG", "CAA:TTG", "GAT:ATC", "GAC:GTC", "GAG:CTC", "GAA:TTC", "AAT:ATT", "AAC:GTT", "AAG:CTT", "AAA:TTT")

triplot <- ggplot(trip.results, aes(x = fct_relevel(class, "TCT:AGA", "TCC:GGA", "TCG:CGA", "TCA:TGA", "CCT:AGG", "CCC:GGG", "CCG:CGG", "CCA:TGG", "GCT:AGC", "GCC:GGC", "GCG:CGC", "GCA:TGC",  "ACT:AGT", "ACC:GGT", "ACG:CGT", "ACA:TGT", "TAT:ATA", "TAC:GTA", "TAG:CTA", "TAA:TTA", "CAT:ATG", "CAC:GTG", "CAG:CTG", "CAA:TTG", "GAT:ATC", "GAC:GTC", "GAG:CTC", "GAA:TTC", "AAT:ATT", "AAC:GTT", "AAG:CTT", "AAA:TTT"), y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +
    ylab("Relative mutation rate") +
    coord_flip()

save_plot(filename = "./mutac_ms/fig/trinucleotides.pdf", plot = triplot, base_width = 3.71, base_height = 2.5*3.71)

epitri.effect <- ggplot(trip.results.domains, aes(y = domain, x = estimate, xmin = lower, xmax = upper)) +
    geom_pointrange() +
    geom_vline(xintercept = 1, lty = "dashed") +
    ylab("") +
    xlab("Mutation rate relative to euchromatin") +
    scale_x_continuous(limits = c(0, 18), breaks = seq(0,16,2))


#Analysis of transversions and transitions
trip.eu.tv <- trinucleotides(filter(aineisto.spectra, EUCHR == 1, point.type == "transversion")$triplet)
trip.eu.tv <- sort.trinucleotides(trinucfreq.eu, trip.eu.tv)
trip.eu.tv$EUCHR <- rep(1, 32)
trip.eu.tv$cent <- rep(0, 32)
trip.eu.tv$H3K9 <- rep(0,32)
trip.eu.tv$H3K27 <- rep(0,32)

trip.eu.ts <- trinucleotides(filter(aineisto.spectra, EUCHR == 1, point.type == "transition")$triplet)
trip.eu.ts <- sort.trinucleotides(trinucfreq.eu, trip.eu.ts)
trip.eu.ts$EUCHR <- rep(1, 32)
trip.eu.ts$cent <- rep(0, 32)
trip.eu.ts$H3K9 <- rep(0,32)
trip.eu.ts$H3K27 <- rep(0,32)

##Building the H3K9 dataset (need to subtract centromere number of mutations and trinuc counts from H3K9 counts)
trip.cent.tv <- trinucleotides(filter(aineisto.spectra, centromere == 1, point.type == "transversion")$triplet)
trip.cent.tv <- sort.trinucleotides(trinucfreq.cent, trip.cent.tv)
trip.cent.tv$EUCHR <- rep(0, 32)
trip.cent.tv$cent <- rep(1, 32)
trip.cent.tv$H3K9 <- rep(1,32)
trip.cent.tv$H3K27 <- rep(0,32)

trip.cent.ts <- trinucleotides(filter(aineisto.spectra, centromere == 1, point.type == "transition")$triplet)
trip.cent.ts <- sort.trinucleotides(trinucfreq.cent, trip.cent.ts)
trip.cent.ts$EUCHR <- rep(0, 32)
trip.cent.ts$cent <- rep(1, 32)
trip.cent.ts$H3K9 <- rep(1,32)
trip.cent.ts$H3K27 <- rep(0,32)

trip.h3k9.tv <- trinucleotides(filter(aineisto.spectra, H3K9 == 1, point.type == "transversion")$triplet)
trip.h3k9.tv <- sort.trinucleotides(trinucfreq.h3k9, trip.h3k9.tv)
trip.h3k9.tv$nmut <- trip.h3k9.tv$nmut - trip.cent.tv$nmut #Subtracting centromere counts and mutations
trip.h3k9.tv$count <- trip.h3k9.tv$count - trip.cent.tv$count
trip.h3k9.tv$EUCHR <- rep(0, 32)
trip.h3k9.tv$cent <- rep(0, 32)
trip.h3k9.tv$H3K9 <- rep(1,32)
trip.h3k9.tv$H3K27 <- rep(0,32)

trip.h3k9.ts <- trinucleotides(filter(aineisto.spectra, H3K9 == 1, point.type == "transition")$triplet)
trip.h3k9.ts <- sort.trinucleotides(trinucfreq.h3k9, trip.h3k9.ts)
trip.h3k9.ts$nmut <- trip.h3k9.ts$nmut - trip.cent.ts$nmut #Subtracting centromere counts and mutations
trip.h3k9.ts$count <- trip.h3k9.ts$count - trip.cent.ts$count
trip.h3k9.ts$EUCHR <- rep(0, 32)
trip.h3k9.ts$cent <- rep(0, 32)
trip.h3k9.ts$H3K9 <- rep(1,32)
trip.h3k9.ts$H3K27 <- rep(0,32)

trip.h3k27.tv <- trinucleotides(filter(aineisto.spectra, H3K27exK9 == 1, point.type == "transversion")$triplet)
trip.h3k27.tv <- sort.trinucleotides(trinucfreq.h3k27, trip.h3k27.tv)
trip.h3k27.tv$EUCHR <- rep(0, 32)
trip.h3k27.tv$cent <- rep(0, 32)
trip.h3k27.tv$H3K9 <- rep(0,32)
trip.h3k27.tv$H3K27 <- rep(1,32)

trip.h3k27.ts <- trinucleotides(filter(aineisto.spectra, H3K27exK9 == 1, point.type == "transition")$triplet)
trip.h3k27.ts <- sort.trinucleotides(trinucfreq.h3k27, trip.h3k27.ts)
trip.h3k27.ts$EUCHR <- rep(0, 32)
trip.h3k27.ts$cent <- rep(0, 32)
trip.h3k27.ts$H3K9 <- rep(0,32)
trip.h3k27.ts$H3K27 <- rep(1,32)

#Combine all transversions
trip.all.tv <- rbind(trip.eu.tv, trip.cent.tv, trip.h3k9.tv, trip.h3k27.tv)
#Calculate the expected number of mutations provided that everything has equal rate
countfreq.tv <- trip.all.tv$count / sum(trip.all.tv$count)
trip.all.tv$expected <- countfreq.tv*sum(trip.all.tv$nmut)

#Combine all transitions
trip.all.ts <- rbind(trip.eu.ts, trip.cent.ts, trip.h3k9.ts, trip.h3k27.ts)
#Calculate the expected number of mutations provided that everything has equal rate
countfreq.ts <- trip.all.ts$count / sum(trip.all.ts$count)
trip.all.ts$expected <- countfreq.ts*sum(trip.all.ts$nmut)

model.all.trip.tv <- brm(data = trip.all.tv, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + H3K27 + class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.trip.ts <- brm(data = trip.all.ts, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + cent + H3K9 + H3K27 + class,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

##Getting posterior predictions from the model
post.trip.tv <- posterior_samples(model.all.trip.tv)
post.trip.ts <- posterior_samples(model.all.trip.ts)

##Results for domains
cent.trip.tv <- post.trip.tv[,1] + post.trip.tv[,2] #Mutation rate in centromeric regions is the sum of centromere effect and H3K9 effect
trip.domain.res.tv <- apply(post.trip.tv[,2:3], 2, exp)
trip.domain.res.tv <- cbind(exp(cent.trip.tv), trip.domain.res.tv)
colnames(trip.domain.res.tv) <- c("centromere", "H3K9", "H3K27")

trip.results.domains.tv <- data.frame(estimate = apply(trip.domain.res.tv,2,median), HPDinterval(as.mcmc(trip.domain.res.tv)))
trip.results.domains.tv$domain <- rownames(trip.results.domains.tv)

cent.trip.ts <- post.trip.ts[,1] + post.trip.ts[,2] #Mutation rate in centromeric regions is the sum of centromere effect and H3K9 effect
trip.domain.res.ts <- apply(post.trip.ts[,2:3], 2, exp)
trip.domain.res.ts <- cbind(exp(cent.trip.ts), trip.domain.res.ts)
colnames(trip.domain.res.ts) <- c("centromere", "H3K9", "H3K27")

trip.results.domains.ts <- data.frame(estimate = apply(trip.domain.res.ts,2,median), HPDinterval(as.mcmc(trip.domain.res.ts)))
trip.results.domains.ts$domain <- rownames(trip.results.domains.ts)

trip.results.domains.all <- rbind(trip.results.domains, trip.results.domains.tv, trip.results.domains.ts)
trip.results.domains.all$type <- c(rep("All mutations", 3), rep("Transversions", 3), rep("Transitions", 3))

### Results for trinucleotides
trip.res.tv <- apply(post.trip.tv[,4:35], 2, exp) #Convert back to normal scale
colnames(trip.res.tv) <- gsub("b_class", "", colnames(trip.res.tv))
trip.res.tv <- trip.res.tv / median(trip.res.tv[,-20], na.rm = T) #Standardize effects of mutation rate (because there was no intercept) to relative mutation rates (drop GAT:ATC from calculation
trip.results.tv <- data.frame(estimate = apply(trip.res.tv,2,median), HPDinterval(as.mcmc(trip.res.tv)))
trip.results.tv$class <- factor(rownames(trip.results.tv))
trip.results.tv[20,1:3] <- NA #No mutations in GAT:ATC, set as missing data


trip.res.ts <- apply(post.trip.ts[,4:35], 2, exp) #Convert back to normal scale
colnames(trip.res.ts) <- gsub("b_class", "", colnames(trip.res.ts))
trip.res.ts <- trip.res.ts / median(trip.res.ts) #Standardize effects of mutation rate (because there was no intercept) to relative mutation rates
trip.results.ts <- data.frame(estimate = apply(trip.res.ts,2,median), HPDinterval(as.mcmc(trip.res.ts)))
trip.results.ts$class <- factor(rownames(trip.results.ts))

###Combine all
trip.results.all <- data.frame(rbind(trip.results, trip.results.tv, trip.results.ts), type = c(rep("All mutations", 32), rep("Transversions", 32), rep("Transitions", 32)))

triplot.all <-  ggplot(trip.results.all, aes(x = fct_relevel(class, "TCT:AGA", "TCC:GGA", "TCG:CGA", "TCA:TGA", "CCT:AGG", "CCC:GGG", "CCG:CGG", "CCA:TGG", "GCT:AGC", "GCC:GGC", "GCG:CGC", "GCA:TGC",  "ACT:AGT", "ACC:GGT", "ACG:CGT", "ACA:TGT", "TAT:ATA", "TAC:GTA", "TAG:CTA", "TAA:TTA", "CAT:ATG", "CAC:GTG", "CAG:CTG", "CAA:TTG", "GAT:ATC", "GAC:GTC", "GAG:CTC", "GAA:TTC", "AAT:ATT", "AAC:GTT", "AAG:CTT", "AAA:TTT"), y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange(fatten = 0.99) +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +
    ylab("Relative mutation rate") +
    coord_flip() +
    facet_grid(~ type)


epitri.effect.all <- ggplot(trip.results.domains.all, aes(y = domain, x = estimate, xmin = lower, xmax = upper)) +
    geom_pointrange(fatten = 0.99) +
    geom_vline(xintercept = 1, lty = "dashed") +
    ylab("") +
    xlab("Mutation rate relative to euchromatin") +
    #scale_x_continuous(limits = c(0, 18), breaks = seq(0,16,2)) +
    facet_grid( ~ type)

trisupplot <- plot_grid(triplot.all, epitri.effect.all, ncol = 1, labels = c("A", "B"), rel_heights = c(1, 0.2))
save_plot(filename = "./mutac_ms/fig/trinuc_sup.pdf", plot = trisupplot, base_width = 3.71*2.5, base_height = 3*3.71)

### Analysing the effects of 5' and 3' bases
#load("./mutac_ms/data/trip_posterior.RData")

#Setting up the data
trip.res <- apply(post.trip[,4:35], 2, exp) #Convert back to normal scale
colnames(trip.res) <- gsub("b_class", "", colnames(trip.res))
trip.res <- trip.res / median(trip.res) #Standardize effects of mutation rate (because there was no intercept) to relative mutation rates
trip.results <- data.frame(estimate = apply(trip.res,2,median), HPDinterval(as.mcmc(trip.res)), stdev = apply(trip.res,2,sd))
trip.results$class <- factor(rownames(trip.results))
trip.results$base <- ifelse(substr(rownames(trip.results),2,2) == "A" , "A:T", "C:G")
trip.results$fiveprime <- ifelse(substr(rownames(trip.results),1,1) == "A" | substr(rownames(trip.results),1,1) == "T" , "A:T", "C:G")
trip.results$threeprime <- ifelse(substr(rownames(trip.results),3,3) == "A" | substr(rownames(trip.results),3,3) == "T" , "A:T", "C:G")

#Running the model

#Initial values
ilist <- list(estimate = trip.results$estimate)
inits_list <- list(ilist, ilist, ilist, ilist) #Needs to be a list of lists of the size of nchains

#Calculate slopes for relationships
relmut.fit <- brm(data = trip.results, family = gaussian,
                  estimate | se(stdev, sigma = TRUE) ~ 1 + base + fiveprime + threeprime + fiveprime:base + threeprime:base,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4, inits = inits_list)


post.bases <- posterior_samples(relmut.fit)

pred.data <- data.frame(base = c(rep("A:T", 4), rep("C:G", 4)), fiveprime = rep(c(rep("A:T",2), rep("C:G",2)),2), threeprime = rep(c("A:T", "C:G"),4), stdev = rep(mean(trip.results$stdev),8))

bases.pred <- data.frame(fitted(relmut.fit, newdata = pred.data))
bases.pred <- cbind(pred.data, bases.pred)
bases.pred$flanks <- rep(c("A:T/A:T", "A:T/C:G", "C:G/A:T", "C:G/C:G"),2)

predictions.plot <- ggplot(bases.pred, aes(y = flanks, x = Estimate, xmin = Q2.5, xmax = Q97.5, color = base)) +
    geom_pointrange(position = position_dodge(width = 0.3)) +
    ylab("Flanking bases (5'/3')") +
    xlab("Relative mutation rate") +
    scale_color_discrete(name = "Focal base") +
    theme(legend.position = "top")    

bases.effects <- data.frame(fixef(relmut.fit))
bases.effects$coef <- c("alpha", "betab", "beta5", "beta3", "betaI5", "betaI3")

effects.plot <- ggplot(bases.effects[-1,], aes(y = fct_relevel(coef, "betaI3", "betaI5", "beta3", "beta5", "betab"), x = Estimate, xmin = Q2.5, xmax = Q97.5)) +
    geom_pointrange() +
    geom_vline(xintercept = 0, lty = "dashed") +
    ylab("") +
    scale_y_discrete(labels = c("betab" = parse(text = TeX("$\\beta_b$")), "beta5" = parse(text = TeX("$\\beta_5$")), "beta3" = parse(text = TeX("$\\beta_3$")), "betaI5" = parse(text = TeX("$\\beta_{I5}$")), "betaI3" = parse(text = TeX("$\\beta_{I3}$"))))
        
rcol <- plot_grid(predictions.plot, effects.plot, align = "v" , axis = "l", ncol = 1, rel_heights = c(1, 0.5), labels = c("B", "C"))

trifinal <- plot_grid(triplot, rcol, ncol = 2, labels = c("A", ""))

save_plot(filename = "./mutac_ms/fig/trinucleotides2.pdf", plot = trifinal, base_width = 2*3.71, base_height = 2.5*3.71)

##Effects of pyrimidine dimers
trip.results$npyr <- c(2,1,2,1,0,1,0,1,1,0,1,0,1,2,1,1,2,1,2,1,0,1,0,1,0,0,1,0,0,2,1,2)
trip.results$nTT <-  c(2,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
trip.results$nCT <-  c(0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,2,1,0,0,0,1,0,0,1,0,1,1,1,2)
trip.results$nCC <-  c(0,0,0,0,0,1,0,0,0,0,0,0,1,2,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0)

ggplot(trip.results, aes(x = npyr, y = estimate, ymin = lower, ymax = upper, color = base)) +
    geom_pointrange(position = position_dodge2(.2)) +
    geom_smooth(method = "lm")

ggplot(trip.results, aes(x = nTT, y = estimate, ymin = lower, ymax = upper, color = base)) +
    geom_pointrange(position = position_dodge2(.2)) +
    geom_smooth(method = "lm")

ggplot(trip.results, aes(x = nCT, y = estimate, ymin = lower, ymax = upper, color = base)) +
    geom_pointrange(position = position_dodge2(.2)) +
    geom_smooth(method = "lm")

ggplot(trip.results, aes(x = nCC, y = estimate, ymin = lower, ymax = upper, color = base)) +
    geom_pointrange(position = position_dodge2(.2)) +
    geom_smooth(method = "lm")

npymodel <- brm(data = trip.results, family = gaussian,
                  estimate | se(stdev, sigma = TRUE) ~ -1 + base + base:npyr,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4, inits = inits_list)

### **** Analysis of triplets (separately)

##Complementary trinucleotides need to be combined
### All mutations (whole genome)
tripletdata <- trinucleotides(aineisto.spectra$triplet)

tripletdata <- sort.trinucleotides(trinucfreq.wg, tripletdata)

#ATtrip <- tripletdata[1:16,] #Mutations that happened at A:T
#ATtrip$frequency <- ATtrip$count / sum(ATtrip$count)
#ATtrip$expected <- sum(ATtrip$nmut)*ATtrip$frequency

tripletdata$expected <- sum(tripletdata$nmut)*tripletdata$frequency #Calculate expected frequencies


model.trip <- brm(data = tripletdata, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res), split = "b_class")
trinuc.res$class <- unname(sapply(classes, '[[', 2))

### Transversions (whole genome)
tripletdata.tv <- trinucleotides(filter(aineisto.spectra, point.type == "transversion")$triplet)
tripletdata.tv <- sort.trinucleotides(trinucfreq.wg, tripletdata.tv)
tripletdata.tv$expected <- sum(tripletdata.tv$nmut)*tripletdata.tv$frequency

model.trip.tv <- brm(data = tripletdata.tv, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.tv)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.tv <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.tv), split = "b_class")
trinuc.res.tv$class <- unname(sapply(classes, '[[', 2))
trinuc.res.tv[20,1:3] <- NA #set those classes with too mutations to missing

### Transitions (whole genome)
tripletdata.ts <- trinucleotides(filter(aineisto.spectra, point.type == "transition")$triplet)
tripletdata.ts <- sort.trinucleotides(trinucfreq.wg, tripletdata.ts)
tripletdata.ts$expected <- sum(tripletdata.ts$nmut)*tripletdata.ts$frequency

model.trip.ts <- brm(data = tripletdata.ts, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.ts)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.ts <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.ts), split = "b_class")
trinuc.res.ts$class <- unname(sapply(classes, '[[', 2))

### All mutations (H3K9)
tripletdata.h3k9 <- trinucleotides(filter(aineisto.spectra, H3K9 == 1)$triplet)
tripletdata.h3k9 <- sort.trinucleotides(trinucfreq.h3k9, tripletdata.h3k9)
tripletdata.h3k9$expected <- sum(tripletdata.h3k9$nmut)*tripletdata.h3k9$frequency

#CGtrip.h3k9 <- tripletdata.h3k9[17:32,]
#CGtrip.h3k9$frequency <- CGtrip.h3k9$count / sum(CGtrip.h3k9$count)
#CGtrip.h3k9$expected <- sum(CGtrip.h3k9$nmut)*CGtrip.h3k9$frequency

model.trip.h3k9 <- brm(data = tripletdata.h3k9, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.h3k9)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.h3k9 <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.h3k9), split = "b_class")
trinuc.res.h3k9$class <- unname(sapply(classes, '[[', 2))
trinuc.res.h3k9[10,1:3] <- NA #set those classes with no mutations to missing

### Transversions (H3K9)
tripletdata.h3k9.tv <- trinucleotides(filter(aineisto.spectra, H3K9 == 1 & point.type == "transversion")$triplet)
tripletdata.h3k9.tv <- sort.trinucleotides(trinucfreq.h3k9, tripletdata.h3k9.tv)
tripletdata.h3k9.tv$expected <- sum(tripletdata.h3k9.tv$nmut)*tripletdata.h3k9.tv$frequency

model.trip.h3k9.tv <- brm(data = tripletdata.h3k9.tv, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.h3k9.tv)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.h3k9.tv <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.h3k9.tv), split = "b_class")
trinuc.res.h3k9.tv$class <- unname(sapply(classes, '[[', 2))
trinuc.res.h3k9.tv[c(10,20),1:3] <- NA

### Transitions (H3K9)
tripletdata.h3k9.ts <- trinucleotides(filter(aineisto.spectra, H3K9 == 1 & point.type == "transition")$triplet)
tripletdata.h3k9.ts <- sort.trinucleotides(trinucfreq.h3k9, tripletdata.h3k9.ts)
tripletdata.h3k9.ts$expected <- sum(tripletdata.h3k9.ts$nmut)*tripletdata.h3k9.ts$frequency

model.trip.h3k9.ts <- brm(data = tripletdata.h3k9.ts, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.h3k9.ts)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.h3k9.ts <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.h3k9.ts), split = "b_class")
trinuc.res.h3k9.ts$class <- unname(sapply(classes, '[[', 2))
trinuc.res.h3k9.ts[c(10),1:3] <- NA


### All mutations (Euchromatin)
tripletdata.eu <- trinucleotides(filter(aineisto.spectra, EUCHR == 1)$triplet)
tripletdata.eu <- sort.trinucleotides(trinucfreq.eu, tripletdata.eu)
tripletdata.eu$expected <- sum(tripletdata.eu$nmut)*tripletdata.eu$frequency

model.trip.eu <- brm(data = tripletdata.eu, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.eu)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.eu <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.eu), split = "b_class")
trinuc.res.eu$class <- unname(sapply(classes, '[[', 2))

### Transversions (Euchromatin)
tripletdata.eu.tv <- trinucleotides(filter(aineisto.spectra, EUCHR == 1 & point.type == "transversion")$triplet)
tripletdata.eu.tv <- sort.trinucleotides(trinucfreq.eu, tripletdata.eu.tv)
tripletdata.eu.tv$expected <- sum(tripletdata.eu.tv$nmut)*tripletdata.eu.tv$frequency

model.trip.eu.tv <- brm(data = tripletdata.eu.tv, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.eu.tv)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.eu.tv <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.eu.tv), split = "b_class")
trinuc.res.eu.tv$class <- unname(sapply(classes, '[[', 2))
trinuc.res.eu.tv[c(11,20,24,26),1:3] <- NA

### Transitions (Euchromatin)
tripletdata.eu.ts <- trinucleotides(filter(aineisto.spectra, EUCHR == 1 & point.type == "transition")$triplet)
tripletdata.eu.ts <- sort.trinucleotides(trinucfreq.eu, tripletdata.eu.ts)
tripletdata.eu.ts$expected <- sum(tripletdata.eu.ts$nmut)*tripletdata.eu.ts$frequency

model.trip.eu.ts <- brm(data = tripletdata.eu.ts, family = poisson,
                  nmut ~ -1 + offset(log(expected)) + class,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.trip.eu.ts)

post <- exp(post[,-33]) #Convert posterior samples back to normal scale

trinuc.res.eu.ts <- data.frame(estimate = apply(post, 2, median),HPDinterval(as.mcmc(post)), stdev = apply(post, 2, sd)) #Calculate HPD intervals
#Need to adjust class names
classes <- strsplit(rownames(trinuc.res.eu.ts), split = "b_class")
trinuc.res.eu.ts$class <- unname(sapply(classes, '[[', 2))
#trinuc.res.h3k9.ts[c(10),1:3] <- NA


##Plotting

trinuc.res.all <- rbind(trinuc.res, trinuc.res.tv, trinuc.res.ts, trinuc.res.h3k9, trinuc.res.h3k9.tv, trinuc.res.h3k9.ts, trinuc.res.eu, trinuc.res.eu.tv, trinuc.res.ts)
trinuc.res.all$type <- rep(c(rep("All mutations", 32), rep("Transversions", 32), rep("Transitions", 32)),3)
trinuc.res.all$domain <- c(rep("Whole genome", 96), rep("H3K9", 96), rep("Euchromatin", 96))

#Save the data
#save(trinuc.res.all, file = "./mutac_ms/data/trinuc_data.RData")

trinucplot <- ggplot(trinuc.res.all, aes(x = fct_relevel(class, "AAA:TTT", "AAC:GTT", "AAG:CTT", "AAT:ATT", "CAA:TTG", "GAA:TTC", "TAA:TTA", "CAC:GTG", "CAG:CTG", "CAT:ATG", "GAG:CTC", "GAC:GTC", "GAT:ATC", "TAC:GTA", "TAG:CTA", "TAT:ATA", "CCC:GGG", "CCA:TGG", "CCG:CGG", "CCT:AGG", "ACC:GGT", "GCC:GGC", "TCC:GGA", "ACA:TGT", "ACG:CGT", "ACT:AGT", "GCG:CGC", "GCA:TGC", "GCT:AGC", "TCT:AGA", "TCA:TGA", "TCG:CGA" ), y = estimate, ymin = lower, ymax = upper, colour = domain)) +
    geom_hline(yintercept = 1, lty = "dashed") +
    geom_pointrange(position = position_dodge2(width = 0.5)) +
    xlab("") +
    ylab("Relative mutation rate") +
    scale_y_continuous(expand = c(0,0), breaks = seq(0,6,1)) +
    scale_colour_manual(values = c("purple", "red", "black")) +
    coord_flip() +
    facet_wrap( ~ type) +
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.4))

save_plot("./mutac_ms/fig/trinuc_domains_new.pdf", trinucplot, base_height=3*3.71, base_width = 3.71*1.618*2)


### **** Trinucleotides and GC-content

#Tripletdata from earlier
#save(tripletdata, file = "./mutac_ms/data/triplets.RData")
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/triplets.RData")

#Load GC.data
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data.RData")

##Load trinucleotide frequencies
trinucfreq.wg <- read.table("Nc_trinuc_freq_whole_genome.txt", header = F)
trinucfreq.h3k9 <- read.table("Nc_trinuc_freq_H3K9.txt", header = F)
trinucfreq.eu <- read.table("Nc_trinuc_freq_euchromatin.txt", header = F)

#Calculate frequencies
colnames(trinucfreq.wg) <- c("trinucleotide", "count")
colnames(trinucfreq.h3k9) <- c("trinucleotide", "count")
colnames(trinucfreq.eu) <- c("trinucleotide", "count")
trinucfreq.wg$freq <- trinucfreq.wg$count / sum(trinucfreq.wg$count)
trinucfreq.h3k9$freq <- trinucfreq.h3k9$count / sum(trinucfreq.h3k9$count)
trinucfreq.eu$freq <- trinucfreq.eu$count / sum(trinucfreq.eu$count)

#Load called sites
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/calledsites.RData")

#tripletdata

#Calculate expected trinucleotide frequencies based on GC-content
gc.res.wg <- calc.gc.from.called(called.chr.matA)
gc.res.h3k9 <- calc.gc.from.called(called.H3K9.matA)
gc.res.eu <- calc.gc.from.called(called.euchr.matA)

trinucfreq.wg$expected <- calc.exp.trinuc.freq(trinucfreq.wg[,1], gc.res.wg[-1])
trinucfreq.h3k9$expected <- calc.exp.trinuc.freq(trinucfreq.h3k9[,1], gc.res.h3k9[-1])
trinucfreq.eu$expected <- calc.exp.trinuc.freq(trinucfreq.eu[,1], gc.res.eu[-1])

tripletfreqs.wg <- sort.trinucleotides.freq(tripletdata, trinucfreq.wg)[,-2]
tripletfreqs.h3k9 <- sort.trinucleotides.freq(tripletdata, trinucfreq.h3k9)[,-2]
tripletfreqs.eu <- sort.trinucleotides.freq(tripletdata, trinucfreq.eu)[,-2]

tripletfreqs <- rbind(tripletfreqs.wg, tripletfreqs.h3k9, tripletfreqs.eu)
tripletfreqs$domain <- c(rep("Whole genome", 32), rep("H3K9", 32), rep("Euchromatin", 32))

trifreq.plot <- ggplot(tripletfreqs, aes(x = fct_relevel(class, "TCT:AGA", "TCC:GGA", "TCG:CGA", "TCA:TGA", "CCT:AGG", "CCC:GGG", "CCG:CGG", "CCA:TGG", "GCT:AGC", "GCC:GGC", "GCG:CGC", "GCA:TGC",  "ACT:AGT", "ACC:GGT", "ACG:CGT", "ACA:TGT", "TAT:ATA", "TAC:GTA", "TAG:CTA", "TAA:TTA", "CAT:ATG", "CAC:GTG", "CAG:CTG", "CAA:TTG", "GAT:ATC", "GAC:GTC", "GAG:CTC", "GAA:TTC", "AAT:ATT", "AAC:GTT", "AAG:CTT", "AAA:TTT"), y = obsfreq / expfreq)) +
    geom_bar(stat = "identity", fill = "grey", colour = "black") +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +
    ylab("Observed / Expected") +
    scale_y_continuous(expand = c(0,0), limits = c(0,3)) +
    coord_flip() +
    facet_wrap( ~ domain) +
    theme(panel.spacing = unit(1, "lines"))

save_plot(filename = "./mutac_ms/fig/trifreqplot.pdf", trifreq.plot, base_height=3.71*1.5*1.618, base_width = 3.71*2) 


### **** Analysis of triplet and 6 different mutations togeter (96 mutation types)

aineisto.spectra$mutation <- classify.mutations(aineisto.spectra$anc.base, aineisto.spectra$sample.base)[,1]
aineisto.spectra$tripclass <- classify.trinucleotides(aineisto.spectra$triplet)

aineisto.spectra.mut.wg <- group_by(aineisto.spectra, mutation, tripclass, .drop = F)
aineisto.spectra.mut.wg <- summarise(aineisto.spectra.mut.wg, nmut = n())
#Add missing combinations...


### *** Insertions, deletions, and translocations

##Small indels
aineisto.indel <- filter(aineisto, type == "deletion" | type == "insertion")
aineisto.indel$type <- factor(aineisto.indel$type)

aineisto.indelsn <- group_by(aineisto.indel, Line, type, .drop = F)
aineisto.indelsn <- summarise(aineisto.indelsn, nmut = n())


##Mutation rate, including all deletions and insertions
model.indels <- brm(data = aineisto.indelsn, family = poisson,
                  nmut ~ -1 + type,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.indels)
deletion.rate <- exp(post[,1])/(mitoses*transfers)
insertion.rate <- exp(post[,2])/(mitoses*transfers)

del.rate <- c(median(deletion.rate),HPDinterval(as.mcmc(deletion.rate)))
ins.rate <- c(median(insertion.rate),HPDinterval(as.mcmc(insertion.rate)))

#Ratio of insertions to deletions
indelrat <- insertion.rate/deletion.rate
indelratio <- c(median(indelrat), HPDinterval(as.mcmc(indelrat)))

#Analysis of deletions and insertions excluding microsatellites and homopolymers
#Filtering microsatellites, homopolymers, and other repeats away
aineisto.indel.ex <- filter(aineisto.indel, !grepl("microsatellite|homopolymer|repeat", aineisto.indel$notes))

#Some summary statistics related to indel lenght
indel.length.ex <- group_by(aineisto.indel.ex, type)
indel.length.ex <- summarise(indel.length.ex, indel.mean = mean(length), indel.sd = sd(length), indel.n = n())

aineisto.indelsn.ex <- group_by(aineisto.indel.ex, Line, type, .drop = F)
aineisto.indelsn.ex <- summarise(aineisto.indelsn.ex, nmut = n())

model.indels.ex <- brm(data = aineisto.indelsn.ex, family = poisson,
                  nmut ~ -1 + type,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.indels.ex)
deletion.rate.ex <- exp(post[,1])/(mitoses*transfers)
insertion.rate.ex <- exp(post[,2])/(mitoses*transfers)

del.rate.ex <- c(median(deletion.rate.ex),HPDinterval(as.mcmc(deletion.rate.ex)))
ins.rate.ex <- c(median(insertion.rate.ex),HPDinterval(as.mcmc(insertion.rate.ex)))

#Ratio of insertions to deletions
indelrat.ex <- insertion.rate.ex/deletion.rate.ex
indelratio.ex <- c(median(indelrat.ex), HPDinterval(as.mcmc(indelrat.ex)))

#Only microsatellites, homopolymer and other repeat associated changes
aineisto.indel.repeats <- filter(aineisto.indel, grepl("microsatellite|homopolymer|repeat", aineisto.indel$notes))

#Some summary statistics related to indel lenght
indel.length.rep <- group_by(aineisto.indel.repeats, type)
indel.length.rep <- summarise(indel.length.rep, indel.mean = mean(length), indel.sd = sd(length), indel.n = n())

aineisto.indelsn.rep <- group_by(aineisto.indel.repeats, Line, type, .drop = F)
aineisto.indelsn.rep <- summarise(aineisto.indelsn.rep, nmut = n())

model.indels.rep <- brm(data = aineisto.indelsn.rep, family = poisson,
                  nmut ~ -1 + type,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.indels.rep)
deletion.rate.rep <- exp(post[,1])/(mitoses*transfers)
insertion.rate.rep <- exp(post[,2])/(mitoses*transfers)

del.rate.rep <- c(median(deletion.rate.rep),HPDinterval(as.mcmc(deletion.rate.rep)))
ins.rate.rep <- c(median(insertion.rate.rep),HPDinterval(as.mcmc(insertion.rate.rep)))

#Ratio of insertions to deletions
indelrat.rep <- insertion.rate.rep/deletion.rate.rep
indelratio.rep <- c(median(indelrat.rep), HPDinterval(as.mcmc(indelrat.rep)))

#Make the final results table
indel.table <- matrix(rep(0, 3*9), ncol = 9)
rownames(indel.table) <- c("all", "exclude repeats", "repeats only")
colnames(indel.table) <- c("del.rate.med", "del.rate.lower", "del.rate.upper", "ins.rate.med", "ins.rate.lower", "ins.rate.upper", "ratio.med", "ratio.low", "ratio.up")
indel.table[1,1:3] <- del.rate; indel.table[1,4:6] <- ins.rate; indel.table[1,7:9] <- indelratio
indel.table[2,1:3] <- del.rate.ex; indel.table[2,4:6] <- ins.rate.ex; indel.table[2,7:9] <- indelratio.ex
indel.table[3,1:3] <- del.rate.rep; indel.table[3,4:6] <- ins.rate.rep; indel.table[3,7:9] <- indelratio.rep

#Calculating summaries of 1 bp indels
indel.1bp <- filter(aineisto.indel, length == 1 | length == -1)
notinhp <-  filter(indel.1bp, !grepl("homopolymer", indel.1bp$notes))
indel.1bp.res <- c(nrow(indel.1bp), nrow(notinhp), nrow(indel.1bp) - nrow(notinhp))

###### Analysis of mutations that happened in homopolymers ###############################
#Check if there is any difference between repeats of A:T and C:G?
indel.hp <- filter(aineisto.indel, grepl("homopolymer", aineisto.indel$notes))
#Get repeated base for each hp
indel.hp$hpbase <- repeated.base(indel.hp$anc.base, indel.hp$sample.base, indel.hp$type)
indel.hp$hpbase <- factor(indel.hp$hpbase)
#Get length of homopolymers where mutations happened
#Get coordinates
#write.csv(indel.hp[,c(1,2,3)], file = "hp_coords.csv")
#Load homopolymer length data
hplength <- read.csv("hp_coords.csv", header = TRUE)
hplength <- arrange(hplength, ID)
indel.hp$hplength <- hplength$hplength

#Load homopolymer occurrances in the reference genome
hpref.counts <- read.csv("hp_genome.csv", header = T)

##Check numbers of mutations in different lenghts of hp, expected vs. observed
AT.hp <- sum(hpref.counts$A.T)
CG.hp <- sum(hpref.counts$C.G)

##Total number of mutations in hp
AT.mut <- nrow(filter(indel.hp, hpbase == "A:T"))
CG.mut <- nrow(filter(indel.hp, hpbase == "C:G"))

##Effect of homopolymer base
aineisto.indel.hp <- group_by(indel.hp, Line, hpbase, .drop = F)
aineisto.indel.hp <- summarise(aineisto.indel.hp, nmut = n())

model.indel.hpbase <- brm(data = aineisto.indel.hp, family = poisson,
                  nmut ~ -1 + hpbase,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.indel.hpbase)
AT.hp.rate <- exp(post[,1])/(mitoses*transfers*AT.hp)
CG.hp.rate <- exp(post[,2])/(mitoses*transfers*CG.hp)

AT.hp.mutrate <- c(median(AT.hp.rate),HPDinterval(as.mcmc(AT.hp.rate)))
CG.hp.mutrate <- c(median(CG.hp.rate),HPDinterval(as.mcmc(CG.hp.rate)))

#Calculates the ratio of AT homopolymer mutations to CG homopolymer mutations
hp.rate.ratio <- AT.hp.rate / CG.hp.rate
hp.mutrate.ratio <- c(median(hp.rate.ratio),HPDinterval(as.mcmc(hp.rate.ratio)))

##Effect of homopolymer length

##Is mutation rate similar across all hp lengths?
#Calculate number mutations for each repeat size
hp.muts <- group_by(indel.hp, hpbase, hplength, .drop = F)
hp.muts <- summarise(hp.muts, nmut = n())

obs.mut <- matrix(rep(0, nrow(hpref.counts)*2), ncol = 2)
colnames(obs.mut) <- c("obs.A:T", "obs.C:G")
hpref.counts <- cbind(hpref.counts, obs.mut)
#Check observed mutations and store them in correct place
for(i in 1:nrow(hp.muts)) {
    row.ind <- which(hpref.counts$Repeats == as.integer(hp.muts[i,2])) #Get row index for hpref.counts
    if(hp.muts[i,1] == "A:T") { hpref.counts[row.ind,4] <- hp.muts[i,3] }
    if(hp.muts[i,1] == "C:G") { hpref.counts[row.ind,5] <- hp.muts[i,3] }
}

hpref.counts$'exp.A:T' <- (AT.mut/AT.hp)*hpref.counts[,2]
hpref.counts$'exp.C:G' <- (CG.mut/CG.hp)*hpref.counts[,3]

hpref.counts$'rel.obs.A:T' <- hpref.counts[,4]/hpref.counts[,6]
hpref.counts$'rel.obs.C:G' <- hpref.counts[,5]/hpref.counts[,7]

#Probably makes sense to combine all mutations for the analysis of length
hpref.counts$genomecount <- hpref.counts[,2] + hpref.counts[,3]
hpref.counts$obs.all <- hpref.counts[,4] + hpref.counts[,5]
hpref.counts$exp.all <- (sum(hpref.counts$obs.all)/(sum(hpref.counts$genomecount)))*hpref.counts$genomecount
hpref.counts$rel.obs.all <- hpref.counts$obs.all/hpref.counts$exp.all

##Save results
hp.results <- c(sum(hpref.counts$obs.all), AT.hp.mutrate, CG.hp.mutrate, hp.mutrate.ratio)
names(hp.results) <- c("observed", "med.AT.hp", "low.AT.hp", "high.AT.hp", "med.CG.hp", "low.CG.hp", "high.CG.hp", "med.ratio.hp", "low.ratio.hp", "high.ratio.hp")

aineisto.indel.hp2 <- group_by(indel.hp, Line, hpbase, hplength, .drop = F)
aineisto.indel.hp2 <- summarise(aineisto.indel.hp2, nmut = n())

   

### 
### Include number of occurrances for each hp class and their frequencies as exposure coefs

###Making a input file for hp for poisson regression
hpinput <- filter(hpref.counts, A.T != 0 | C.G != 0)
hpfinal <- data.frame(Repeats = rep(hpinput$Repeats, 2), Counts = c(hpinput$A.T, hpinput$C.G), nmut = c(hpinput$'obs.A:T', hpinput$'obs.C:G'), hpbase = factor(c(rep("A:T", nrow(hpinput)), rep("C:G", nrow(hpinput)))))

hpfinal <- filter(hpfinal, Counts != 0)


model.indel.hp2 <- brm(data = hpfinal, family = poisson,
                  nmut ~ 1 + offset(log(Counts)) + hpbase + Repeats + hpbase:Repeats,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)
## Interactions effects not different from zero

model.indel.hp <- brm(data = hpfinal, family = poisson,
                  nmut ~ 1 + offset(log(Counts)) + hpbase + Repeats,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

##Getting posterior predictions from the model
post.hp <- posterior_samples(model.indel.hp)

hp.model.res <- fixef(model.indel.hp)

l.summary <- fitted(model.indel.hp)
l.summary <- l.summary/hpfinal$Counts
l.summary <- data.frame(l.summary)
l.summary$hpbase <- hpfinal$hpbase
l.summary$Repeats <- hpfinal$Repeats

l.summary2 <- fitted(model.indel.hp)
l.summary2 <- data.frame(l.summary2)
l.summary2$hpbase <- hpfinal$hpbase
l.summary2$Repeats <- hpfinal$Repeats

#Data used for prediction
#Note that need to predict for 1 count (to get the base estimates)
preddata.hp <- data.frame(Repeats = c(5:40, 5:40), Counts = rep(1, 36*2), hpbase = c(rep("A:T", 36), rep("C:G", 36)))

l.summary3 <- fitted(model.indel.hp, newdata = preddata.hp)
l.summary3 <- cbind(l.summary3, preddata.hp)

hp.model.plot <- ggplot(l.summary3, aes( x = Repeats, y = log(Estimate), ymin = log(Q2.5), ymax = log(Q97.5), fill = hpbase)) +
       geom_ribbon(alpha = 0.5) +
       geom_smooth(method = "lm", aes(colour = hpbase)) +
       xlab("Length of homopolymer (bp)") +
       ylab("log(model estimate)") +
       theme(legend.position = c(0.1, 0.8), legend.title = element_blank())    

ggplot(l.summary, aes(x = Repeats, y = Estimate, ymin = Q2.5, ymax = Q97.5, colour = hpbase)) +
    geom_smooth() +
    scale_x_continuous(limits = c(5, 34), breaks = seq(5,34,2))

ggplot(hpfinal, aes(x = Repeats, y = nmut/hpfinal$Counts, colour = hpbase)) +
    geom_point() +
    geom_smooth(data = l.summary2, aes(x = Repeats, y = Estimate/hpfinal$Counts, ymin = Q2.5/hpfinal$Counts, ymax = Q97.5/hpfinal$Counts, colour = hpbase)) +
    scale_x_continuous(limits = c(5, 34), breaks = seq(5,34,2))    
    #scale_y_continuous(limits = c(0, 0.025))

##
ggplot(hpfinal, aes(x = log(Repeats), y = log10(nmut/Counts), colour = hpbase)) +
    geom_point() +
    geom_smooth(data = l.summary, aes(x = log(Repeats), y = log10(Estimate), ymin = log10(Q2.5), ymax = log10(Q97.5), colour = hpbase)) 
    #scale_x_continuous(limits = c(5, 34), breaks = seq(5,34,2)) +
    #scale_y_continuous(limits = c(0, 0.1))

ggplot(hpfinal, aes(x = log(Repeats), y = log10(nmut/Counts), colour = hpbase)) +
    geom_point()

###Model with Repeats as a factor to estimate mutation rate
model.indel.hp2 <- brm(data = filter(hpfinal, Repeats < 23), family = poisson,
                  nmut ~ 1 + offset(log(Counts)) + hpbase + factor(Repeats) + hpbase:factor(Repeats),
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)
                       
l.summary3 <- fitted(model.indel.hp2)
l.summary3 <- data.frame(l.summary3)
l.summary3$hpbase <- filter(hpfinal, Repeats < 23)$hpbase
l.summary3$Repeats <- filter(hpfinal, Repeats < 23)$Repeats


##Effect of homopolymer base
koe <- filter(indel.hp, hplength != 2)
koe$hplength <- factor(koe$hplength)

koe <- group_by(koe, Line, hpbase, hplength, .drop = F)
koe <- summarise(koe, nmut = n())

koe <- data.frame(koe)
koe[,3] <- as.numeric(as.character(koe[,3]))

count <- rep(0, nrow(koe))
for(i in 1:nrow(koe)) { if(koe$hpbase[i] == "A:T") {
      ind <- which(hpref.counts$Repeats == koe$hplength[i])
      count[i] <- hpref.counts$A.T[ind] }
                       if(koe$hpbase[i] == "C:G") {
      ind <- which(hpref.counts$Repeats == koe$hplength[i])
      count[i] <- hpref.counts$C.G[ind] }
                    }
koe$count <- count
koe <- filter(koe, count > 0)

model.koe <- brm(data = koe, family = poisson,
                  nmut ~ 1 + offset(log(count)) + hpbase + hplength,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post.koe <- posterior_samples(model.koe)

#Making predictions using the model
pred.data <- data.frame(hplength = rep(seq(5,22,1),2), hpbase = factor(c(rep("A:T", 18), rep("C:G", 18))), estimate = rep(0,36), estmin = rep(0,36), estmax = rep(0,36), pred.mut = rep(0,36), mutmin = rep(0,36), mutmax = rep(0,36))
pred.data$obs.mut <- c(hpfinal[1:18,]$nmut, hpfinal[71:88,]$nmut)
pred.data$counts <- c(hpfinal[1:18,]$Counts, hpfinal[71:88,]$Counts)

for(i in 1:nrow(pred.data)) {
lambda <- exp(post.koe[,1] + ifelse(pred.data[i,2] == "C:G", 1, 0)*(post.koe[,2]) + pred.data[i,1]*(post.koe[,3]))
pred.data[i,3:5] <- c(median(lambda),HPDinterval(as.mcmc(lambda)))
#Predicted mutations
predmut <- lambda*39*pred.data$counts[i] #Predicted mutations are rate(per locus)*Lines*loci (transfer not in estimate, so not needed here
pred.data[i,6:8] <- c(median(predmut),HPDinterval(as.mcmc(predmut)))
}
#Still the same problem

#############

### Translocations
aineisto.transloc <- filter(aineisto, type == "translocation")
aineisto.transloc$type <- factor(aineisto.transloc$type)

aineisto.translocn <- group_by(aineisto.transloc, Line, type, .drop = F)
aineisto.translocn <- summarise(aineisto.translocn, nmut = n())

model.transloc <- brm(data = aineisto.translocn, family = poisson,
                  nmut ~ 1,
                  prior = c(prior(normal(0, 10), class = Intercept)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.transloc)
translocation.rate <- exp(post[,1])/(mitoses*transfers)

translocation.rate <- c(median(translocation.rate),HPDinterval(as.mcmc(translocation.rate)))

results.transloc <- c(nrow(aineisto.transloc), translocation.rate)
names(results.transloc) <- c("observed.count", "med.transloc.rate", "low.transloc.rate", "high.transloc.rate")

### Drawing a figure about distribution of indel lenghts ###
delinset <- ggplot(filter(aineisto.indel, type == "deletion"), aes(x = length)) +
    geom_histogram(fill = "red", colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks = c(-15000, -10000, -5000, 0), labels = c(-15, -10, -5, 0)) +    
    xlab("Length of deletion (kb)") +
    ylab("")    

#Plot deletions and insertions
p <- ggplot(filter(aineisto.indel, type == "deletion"), aes(x = length)) +
    geom_histogram(fill = "red", colour = "black", binwidth = 1) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = -c(75, 0), breaks = seq(from = -70, to = 0, by = 10)) +
    xlab("Length of deletion (bp)") +
    ylab("Number of deletions")    

p1 <- ggdraw(p) +
    draw_plot(delinset, .15, .45, .5, .5)
#Inset showing the whole range        

p2 <- ggplot(filter(aineisto.indel, type == "insertion"), aes(x = length)) +
    geom_histogram(fill = "blue", colour = "black", binwidth = 1) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(.01,0), breaks = seq(0, 125, 10)) +
    xlab("Length of insertion (bp)") +
    ylab("Number of insertions")

#Plot lengths of hp streches where mutations occurred
p3 <- ggplot(indel.hp, aes(x = hplength)) +
    geom_histogram(colour = "black", fill = "white") +
    xlab("Length of homopolymer (bp)") +
    ylab("Number of mutations") +
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(~ hpbase)

p4pre <- ggplot(hpref.counts, aes(x = Repeats, y = rel.obs.all)) +
    geom_bar(stat = "identity", colour = "black", fill = "white") +
    geom_text(aes(label = hpref.counts$obs.all), vjust = -0.75) +
    scale_y_continuous(limits = c(0,30), expand = c(0,0)) +
    scale_x_continuous(limits = c(5,22), breaks = seq(5,30,1)) +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Observed / Expected") +
    xlab("Length of homopolymer (bp)") 

indelplot <- plot_grid(p1, p2, p3, p4pre, labels = c("A", "B", "C", "D"), rel_heights = c(2,1), label_y = c(1,1,1.2,1.2))
save_plot("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/fig/indelplot.pdf", indelplot, nrow = 1, base_height=6, base_width=10)

indelplot2 <- plot_grid(p1, p2, p3, hp.model.plot, labels = c("A", "B", "C", "D"), rel_heights = c(2, 1.5), label_y = c(1,1,1.1, 1.1))
save_plot("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/fig/indelplot2.pdf", indelplot2, nrow = 1, base_height=6, base_width = 10)

pdf(file = "./mutac_ms/fig/hp_model.pdf")
ggplot(filter(hpfinal, Repeats < 23), aes(x = Repeats, y = nmut, colour = hpbase)) +
    geom_point() +
    geom_pointrange(data = l.summary3, aes(x = Repeats, y = Estimate, ymin = Q2.5, ymax = Q97.5, colour = hpbase), position = position_dodge(0.5), shape = 21, alpha = 0.5 ) +
    scale_x_continuous(limits = c(4, 23), breaks = seq(5,22,1)) +
    scale_colour_discrete(name = "Type") +
    ylab("Number of mutations") +
    xlab("Homopolymer length (bp)") +
    theme(legend.position=c(0.9,0.9))
      dev.off()

########### END PLOTTING indel stuff######################################################
##########################################################################################


#### Save the indel results  ###########################################################
save(indel.table, indel.length.rep, indel.length.ex, hp.results, indel.1bp.res, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/indels.RData")

save(results.transloc, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/translocations.RData")
### Done saving indel results ##########################################################

### **** New analysis of indels but now with taking chromosome domains into account

# 1. Load the data about genomic features
# 2. Load the data about called bases in different domains

##Small indels
aineisto.indel.domains <- filter(aineisto, type == "deletion" | type == "insertion")
aineisto.indel.domains$type <- factor(aineisto.indel.domains$type)
indel.domain <- rep("", nrow(aineisto.indel.domains))
for(i in 1:length(indel.domain)) {
    if(aineisto.indel.domains$centromere[i] == 1) {indel.domain[i] <- "Centromere"}
    if(aineisto.indel.domains$EUCHR[i] == 1) { indel.domain[i] <- "Euchromatin" }
    if(aineisto.indel.domains$H3K9[i] == 1 & aineisto.indel.domains$centromere[i] == 0) { indel.domain[i] <- "H3K9" }
    if(aineisto.indel.domains$H3K27exK9[i] == 1) {indel.domain[i] <- "H3K27" }
}
aineisto.indel.domains$domain <- factor(indel.domain)

#All insertion and deletions
aineisto.indelsn.domains <- group_by(aineisto.indel.domains, type, domain, .drop = F)
aineisto.indelsn.domains <- summarise(aineisto.indelsn.domains, nmut = n())

#Only those indels that do not occur in repeats
aineisto.indel.domains.exr <- filter(aineisto.indel.domains, !grepl("microsatellite|homopolymer|repeat", aineisto.indel.domains$notes))
aineisto.indeldom.exr.counts <- group_by(aineisto.indel.domains.exr, type, domain, .drop = F)
aineisto.indeldom.exr.counts <- summarise(aineisto.indeldom.exr.counts, nmut = n())

#Only those indels that occur in repeats
aineisto.indel.domains.rep <- filter(aineisto.indel.domains, grepl("microsatellite|homopolymer|repeat", aineisto.indel.domains$notes))
aineisto.indeldom.rep.counts <- group_by(aineisto.indel.domains.rep, type, domain, .drop = F)
aineisto.indeldom.rep.counts <- summarise(aineisto.indeldom.rep.counts, nmut = n())

#Load the number of sites called
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/calledsites.RData")

aineisto.indelsn.domains$sites <- rep( c(  mean(c(sum(unlist(called.centromeric.matA)), sum(unlist(called.centromeric.mata)))), mean(c(sum(unlist(called.euchr.matA)), sum(unlist(called.euchr.mata)))), mean(c(sum(unlist(called.H3K27.exK9.matA)), sum(unlist(called.H3K27.exK9.mata)))), mean(c(sum(unlist(called.H3K9.exc.matA)), sum(unlist(called.H3K9.exc.mata))))), 2)
aineisto.indeldom.exr.counts$sites <- aineisto.indelsn.domains$sites
aineisto.indeldom.rep.counts$sites <- aineisto.indelsn.domains$sites

#Saving the data, so can load it later
#save(aineisto, aineisto.indel.domains, aineisto.indelsn.domains, aineisto.indeldom.exr.counts, aineisto.indeldom.rep.counts, file = "./mutac_ms/data/indels.domains.RData")
#load the above
#load("./mutac_ms/data/indels.domains.RData")

##Mutation rate, including all deletions and insertions
model.indels.domains <- brm(data = aineisto.indelsn.domains, family = poisson,
                  nmut ~ -1 + offset(log(sites)) + type:domain,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

#Making posterior estimate for each type of mutation
post <- posterior_samples(model.indels.domains)[,-9]

#Calculating rates, relative to euchromatin
post.rates <- apply(post, 2, exp)
del.rates <- post.rates[,c(1,5,7)] #Euchromatin not included
del.rel.rates <- del.rates / post.rates[,3]

ins.rates <- post.rates[,c(2,6,8)] #Euchromating not included
ins.rel.rates <- ins.rates / post.rates[,4]

#Making a table of results
indel.domains.res <- data.frame(estimate = c(apply(del.rel.rates, 2, median), apply(ins.rel.rates, 2, median)), rbind(HPDinterval(as.mcmc(del.rel.rates)), HPDinterval(as.mcmc(ins.rel.rates))) )
indel.domains.res$type <- factor(c(rep("Deletion", 3), rep("Insertion", 3)))
indel.domains.res$domain <- factor(rep(c("Centromere", "H3K27", "H3K9"), 2))


##Mutation rate, including only indels that do not occur in repeats
model.indels.exr.domains <- brm(data = aineisto.indeldom.exr.counts, family = poisson,
                  nmut ~ -1 + offset(log(sites)) + type:domain,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

#Making posterior estimate for each type of mutation
post <- posterior_samples(model.indels.exr.domains)[,-9]

#Calculating rates, relative to euchromatin
post.rates <- apply(post, 2, exp)
del.rates <- post.rates[,c(1,5,7)] #Euchromatin not included
del.rel.rates <- del.rates / post.rates[,3]

ins.rates <- post.rates[,c(2,6,8)] #Euchromating not included
ins.rel.rates <- ins.rates / post.rates[,4]

#Making a table of results
indel.domains.exr.res <- data.frame(estimate = c(apply(del.rel.rates, 2, median), apply(ins.rel.rates, 2, median)), rbind(HPDinterval(as.mcmc(del.rel.rates)), HPDinterval(as.mcmc(ins.rel.rates))) )
indel.domains.exr.res$type <- factor(c(rep("Deletion", 3), rep("Insertion", 3)))
indel.domains.exr.res$domain <- factor(rep(c("Centromere", "H3K27", "H3K9"), 2))
indel.domains.exr.res[4,1:3] <- NA #No insertions occurred in Centromeres

### Mutation rate, including only indels that occur in repeats ###
model.indels.rep.domains <- brm(data = aineisto.indeldom.rep.counts, family = poisson,
                  nmut ~ -1 + offset(log(sites)) + type:domain,
                  prior = c(prior(normal(0, 10), class = b)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

#Making posterior estimate for each type of mutation
post <- posterior_samples(model.indels.rep.domains)[,-9]

#Calculating rates, relative to euchromatin
post.rates <- apply(post, 2, exp)
del.rates <- post.rates[,c(1,5,7)] #Euchromatin not included
del.rel.rates <- del.rates / post.rates[,3]

ins.rates <- post.rates[,c(2,6,8)] #Euchromating not included
ins.rel.rates <- ins.rates / post.rates[,4]

#Making a table of results
indel.domains.rep.res <- data.frame(estimate = c(apply(del.rel.rates, 2, median), apply(ins.rel.rates, 2, median)), rbind(HPDinterval(as.mcmc(del.rel.rates)), HPDinterval(as.mcmc(ins.rel.rates))) )
indel.domains.rep.res$type <- factor(c(rep("Deletion", 3), rep("Insertion", 3)))
indel.domains.rep.res$domain <- factor(rep(c("Centromere", "H3K27", "H3K9"), 2))

######################################################################################3
### Making the final results.table
indel.domains.res.all <- rbind(indel.domains.res, indel.domains.exr.res, indel.domains.rep.res)
indel.domains.res.all$class <- c(rep("All indels", 6), rep("Repeats excluded", 6), rep("Repeats only", 6))

##Save results for plot
save(indel.domains.res.all, file = "./mutac_ms/data/indel.domains.res.RData")

###Making the plot

#Load results table for plotting
#load("./mutac_ms/data/indel.domains.res.RData") #indel.domains.res.all contains results for plotting
indel.domains.res.all$signif <- ifelse(indel.domains.res.all$lower > 1, "red", "black")


indel.domains.plot <- ggplot(indel.domains.res.all, aes(x = domain, y = estimate, ymin = lower, ymax = upper, colour = signif)) +
    geom_hline(yintercept = 1, lty = "dashed") +
    geom_pointrange() +
    xlab("") +
    ylab("Mutation rate relative to euchromatin") +
    scale_colour_manual(values = c("black", "red")) +
    facet_grid(class ~ type) +
    theme(legend.position="none")

save_plot(filename = "./mutac_ms/fig/indel_domains.pdf", indel.domains.plot, base_height = 3.71*2, base_width = 6)


### *** Complex mutations

aineisto.complex <- filter(aineisto, type == "complex")

aineisto.complexmut <- group_by(aineisto.complex, Line, .drop = F)
aineisto.complexmut <- summarise(aineisto.complexmut, nmut = n())

model.complex <- brm(data = aineisto.complexmut, family = poisson,
                  nmut ~ 1,
                  prior = c(prior(normal(0, 10), class = Intercept)),
                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.complex)
#Calculate mutation rate
mutation.rate.complex <- exp(post[,1])/(mitoses*transfers*called.bases)
#get estimates
rate.complex <- c(median(mutation.rate.complex), HPDinterval(as.mcmc(mutation.rate.complex)))

##What the ratio of point mutations to complex mutations?
#This requires calculating point mutation rate from above
ratio.complex <- mutation.rate.point / mutation.rate.complex
ratio.complex.estimate <- c(median(ratio.complex), HPDinterval(as.mcmc(ratio.complex)))

ncomplex <- nrow(aineisto.complex)

### Save complex mutation results #################################################################
save(ncomplex, rate.complex, ratio.complex.estimate, file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/complex.RData")
###################################################################################################

### *** Genomic features and mutation rate

### **** Checking mutation qualities in euchromatin and H3K9 regions

#Mutations in the seven chromosomes only
point.mut <- filter(aineisto, type == "point")

#Overwhelming number of point mutations have genotype quality of 99
q1 <- ggplot(point.mut, aes(x = sample.GQ)) +
    geom_histogram(fill = "grey", colour = "black") +
    xlab("Mutation genotype quality") +
    scale_y_continuous(expand = c(0,0))

for(i in 1:nrow(point.mut)) {
    if(point.mut$H3K9[i] == 1 & point.mut$centromere[i] == 0) { point.mut$domain[i] <- "H3K9 exc. cent." }
    if(point.mut$centromere[i] == 1) { point.mut$domain[i] <- "Centromeric" }
    if(point.mut$EUCHR[i] == 1) { point.mut$domain[i] <- "Euchromatin" }
}

point.mut$domain <- forcats::fct_relevel(point.mut$domain, "H3K9 exc. cent.", "Centromeric", "Euchromatin")

q2 <- ggplot(point.mut, aes(x = sample.GQ)) +
    geom_histogram(fill = "grey", colour = "black") +
    xlab("Mutation genotype quality") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 400)) +
    facet_wrap( ~ domain, ncol = 2)

#Make final plot
qfinal <- plot_grid(q1, q2, nrow = 1, labels = c("A", "B"), align = "h", axis = "b", rel_widths = c(1,3))
save_plot(filename = "./mutac_ms/fig/mutquality.pdf", plot = qfinal, base_height = 3.71, base_width = 3.71*4)

#New version of the plot
qp1 <- ggplot(point.mut, aes(x = sample.GQ)) +
    geom_histogram(fill = "grey", colour = "black") +
    xlab("Mutation genotype quality") +
    ggtitle("Whole genome") +    
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(29, 101), breaks = seq(30,100,10)) +
    theme(text = element_text(size = 20))    

qp2 <- ggplot(filter(point.mut, domain == "Euchromatin"), aes(x = sample.GQ)) +
    geom_histogram(fill = "grey", colour = "black") +
    xlab("Mutation genotype quality") +
    ggtitle("Euchromatin") +    
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(29, 101), breaks = seq(30,100,10)) +
    theme(text = element_text(size = 20))    

qp3 <- ggplot(filter(point.mut, domain == "H3K9 exc. cent."), aes(x = sample.GQ)) +
    geom_histogram(fill = "grey", colour = "black") +
    xlab("Mutation genotype quality") +
    ggtitle("H3K9 exc. cent.") +    
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(29, 101), breaks = seq(30,100,10)) +
    theme(text = element_text(size = 20))    

qp4 <- ggplot(filter(point.mut, domain == "Centromeric"), aes(x = sample.GQ)) +
    geom_histogram(fill = "grey", colour = "black") +
    xlab("Mutation genotype quality") +
    ggtitle("Centromeric") +    
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(29, 101), breaks = seq(30,100,10)) +
    theme(text = element_text(size = 20))    

qpfinal <- plot_grid(qp1, qp2, qp3, qp4, nrow = 2, align = "h", axis = "b")
save_plot(filename = "./mutac_ms/fig/mutquality.pdf", plot = qpfinal, base_height = 3.71*2, base_width = 3.71*2*1.68)


#Save point mutations to check which have been verified by Sanger
write.table(point.mut, file = "point_mutations_domains.csv", sep = ",", row.names = F, quote = F)

#How many (%) mutations have genotype quality less than 80?
sum(point.mut$sample.GQ < 80) / length(point.mut$sample.GQ) * 100

sum(point.mut$sample.GQ < 50) / length(point.mut$sample.GQ) * 100

### *** Comparing the different chromosomes
bychr <- filter(aineisto, Chromosome != "Supercontig_12.8" & Chromosome != "Supercontig_12.9" & Chromosome != "Supercontig_12.10" & Chromosome != "Supercontig_12.11" & Chromosome != "Supercontig_12.12" & Chromosome != "Supercontig_12.13" & Chromosome != "Supercontig_12.14" & Chromosome != "Supercontig_12.15" & Chromosome != "Supercontig_12.17")
bychr$Chromosome <- factor(bychr$Chromosome)

#Get line numbers and mating type for each line
line <- strsplit(as.character(bychr[,3]), "L")
line <- sapply(line, '[[', 2)
line <- strsplit(line, "G")
line <- as.numeric(sapply(line, '[[', 1))
bychr$mat <- ifelse(line <= 20, "matA", "mata") #Line numbers <= 20 are mat A, >= 21 are mat a

#bychrgroups <- group_by(bychr, Chromosome, mat, duplicatedreg, centromere, H3K27, H3K9)
##Since number of called sites is slightly different for mat A and mat a lines, I need to adjust for that
##Mutation by chromosome
bychrgroups <- group_by(bychr, Chromosome, mat)
bychr.all <- summarise(bychrgroups, nmut = n())

#unlist(lapply(called.chr.matA, sum)) #Called sites for each chromosome for mat A lines
#unlist(lapply(called.chr.mata, sum)) #Called sites for each chromosome for mat a lines
bychr.all$sites <- c(rbind(unlist(lapply(called.chr.matA, sum)), unlist(lapply(called.chr.mata, sum))))
bychr.all$lines <- c(rbind(rep(20,7), rep(19,7))) #Number of lines, since line 32 was dropped there are only 19 mat a lines.

#Number of mutations by mating type
matAtotal <- as.numeric(summarise(group_by(bychr.all, mat), total = sum(nmut))[1,2])
matatotal <- as.numeric(summarise(group_by(bychr.all, mat), total = sum(nmut))[2,2])

#Expected number of mutations for each chromosome
#(Total number of mutations / Total number of called sites, mat A) * called sites per chromosome
exp.matA <- (matAtotal / sum(unlist(lapply(called.chr.matA, sum)))) * unlist(lapply(called.chr.matA, sum))
exp.mata <- (matatotal / sum(unlist(lapply(called.chr.mata, sum)))) * unlist(lapply(called.chr.mata, sum))

bychr.all$expected <- c(rbind(exp.matA, exp.mata))

#Model to calculate mutation rate for each chromosome
model.chromosome <- brm(data = bychr.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + Chromosome,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.chromosome)

chromosome.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-8,]

rate.bychr.all <- data.frame(Chromosome = c("Chr I", "Chr II", "Chr III", "Chr IV", "Chr V", "Chr VI", "Chr VII"), chromosome.rel.rates)
colnames(rate.bychr.all) <- c("Chromosome", "estimate", "lower", "upper")

##Plot the results
ggplot(rate.bychr.all, aes(x = Chromosome, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Relative mutation rate")    


### Same analysis but for point mutations only ###
pointmutations <- filter(bychr, type == "point")

bychrgroups.point <- group_by(pointmutations, Chromosome, mat)
bychr.point <- summarise(bychrgroups.point, nmut = n())

bychr.point$sites <- c(rbind(unlist(lapply(called.chr.matA, sum)), unlist(lapply(called.chr.mata, sum))))
bychr.point$lines <- c(rbind(rep(20,7), rep(19,7)))

#Number of mutations by mating type
matAtotal.point <- as.numeric(summarise(group_by(bychr.point, mat), total = sum(nmut))[1,2])
matatotal.point <- as.numeric(summarise(group_by(bychr.point, mat), total = sum(nmut))[2,2])

#Expected number of mutations for each chromosome
#(Total number of mutations / Total number of called sites, mat A) * called sites per chromosome
exp.matA.point <- (matAtotal.point / sum(unlist(lapply(called.chr.matA, sum)))) * unlist(lapply(called.chr.matA, sum))
exp.mata.point <- (matatotal.point / sum(unlist(lapply(called.chr.mata, sum)))) * unlist(lapply(called.chr.mata, sum))

bychr.point$expected <- c(rbind(exp.matA.point, exp.mata.point))

#bychr <- filter(pointmutations, Chromosome != "Supercontig_12.8" & Chromosome != "Supercontig_12.9" & Chromosome != "Supercontig_12.10" & Chromosome != "Supercontig_12.11" & Chromosome != "Supercontig_12.12" & Chromosome != "Supercontig_12.13" & Chromosome != "Supercontig_12.14" & Chromosome != "Supercontig_12.15" & Chromosome != "Supercontig_12.17")
#bychr$Chromosome <- factor(bychr$Chromosome)
#bychrgroups <- group_by(bychr, Chromosome, Line)
#bychr.point <- summarise(bychrgroups, nmut = n())

model.chromosome.point <- brm(data = bychr.point, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + Chromosome, #No intercept 
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.chromosome.point)

chromosome.rel.rates.point <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-8,]

rate.bychr.point <- data.frame(Chromosome = c("Chr I", "Chr II", "Chr III", "Chr IV", "Chr V", "Chr VI", "Chr VII"), chromosome.rel.rates.point)
colnames(rate.bychr.point) <- c("Chromosome", "estimate", "lower", "upper")

##Plot the results
ggplot(rate.bychr.point, aes(x = Chromosome, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Relative mutation rate") 

### *** Effects of epigenetic domains etc.

#Need to use the data in bychr
groups.h3k9 <- group_by(bychr, H3K9, mat)
h3k9.all <- summarise(groups.h3k9, nmut = n())

#How many sites called in total, and how many in H3K9?
#called.chr.matA
#called.H3K9.matA

h3k9.all$sites <- c(sum(unlist(called.chr.matA)) - sum(unlist(called.H3K9.matA)), sum(unlist(called.chr.mata)) - sum(unlist(called.H3K9.mata)), sum(unlist(called.H3K9.matA)), sum(unlist(called.H3K9.mata)) )

#Expected number of mutations for each domain
#(Total number of mutations / Total number of called sites, mat A) * called sites per domain

h3k9.all$expected <- c(matAtotal / sum(unlist(called.chr.matA)) * (sum(unlist(called.chr.matA)) - sum(unlist(called.H3K9.matA))), matatotal / sum(unlist(called.chr.mata)) * (sum(unlist(called.chr.mata)) - sum(unlist(called.H3K9.mata))), matAtotal / sum(unlist(called.chr.matA)) * sum(unlist(called.H3K9.matA)), matatotal / sum(unlist(called.chr.mata)) * sum(unlist(called.H3K9.mata))  )

h3k9.all$domain <- c("euchromatin", "euchromatin", "H3K9", "H3K9")

#Model
model.H3K9.all <- brm(data = h3k9.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + domain,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.H3K9.all)

#h3k9.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-3,]
h3k9.rel.rates <- quantile(exp(post[,1])/exp(post[,2]), probs = c(0.5, 0.025, 0.975))

rate.h3k9.all <- data.frame(Domain = "H3K9", t(h3k9.rel.rates))
colnames(rate.h3k9.all) <- c("Domain", "estimate", "lower", "upper")

##Plot the results
ggplot(rate.h3k9.all, aes(x = Domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Relative mutation rate")

groups.h3k27 <- group_by(bychr, H3K27exK9, mat) #Mutations in H3K27 ex K9 regions
h3k27.all <- summarise(groups.h3k27, nmut = n())

groups.euchr <- group_by(bychr, EUCHR, mat)
euchr.all <- summarise(groups.euchr, nmut = n())

h3k27.all <- rbind(euchr.all[3:4,], h3k27.all[3:4,])

h3k27.all$sites <- c(sum(unlist(called.euchr.matA)), sum(unlist(called.euchr.mata)), sum(unlist(called.H3K27.exK9.matA)), sum(unlist(called.H3K27.exK9.mata)))

#Expected number of mutations for each domain
#(Total number of mutations / Total number of called sites, mat A) * called sites per domain

h3k27.all$expected <- unlist(c(sum(h3k27.all[1,3] + h3k27.all[3,3]) / sum(h3k27.all[1,5] + h3k27.all[3,5]) * h3k27.all[1,5], sum(h3k27.all[2,3] + h3k27.all[4,3]) / sum(h3k27.all[2,5] + h3k27.all[4,5]) * h3k27.all[2,5], sum(h3k27.all[1,3] + h3k27.all[3,3]) / sum(h3k27.all[1,5] + h3k27.all[3,5]) * h3k27.all[3,5], sum(h3k27.all[2,3] + h3k27.all[4,3]) / sum(h3k27.all[2,5] + h3k27.all[4,5]) * h3k27.all[4,5]))

h3k27.all$domain <- c("euchromatin", "euchromatin", "H3K27", "H3K27")

model.H3K27.all <- brm(data = h3k27.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + domain,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.H3K27.all)

#h3k27.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-3,]
h3k27.rel.rates <- quantile(exp(post[,1])/exp(post[,2]), probs = c(0.5, 0.025, 0.975))

rate.h3k27.all <- data.frame(Domain = c("H3K27"), t(h3k27.rel.rates))
colnames(rate.h3k27.all) <- c("Domain", "estimate", "lower", "upper")

##Plot the results
ggplot(rate.h3k27.all, aes(x = Domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Relative mutation rate")

#Centrome regions
groups.cent <- group_by(bychr, centromere, mat)
cent.all <- summarise(groups.cent, nmut = n())

cent.all <- rbind(h3k9.all[1:2,], cent.all[3:4,])

cent.all$sites <- c(sum(unlist(called.chr.matA)) - sum(unlist(called.centromeric.matA)), sum(unlist(called.chr.mata)) - sum(unlist(called.centromeric.mata)), sum(unlist(called.centromeric.matA)), sum(unlist(called.centromeric.mata)) )

cent.all$expected <- c(matAtotal / sum(unlist(called.chr.matA)) * (sum(unlist(called.chr.matA)) - sum(unlist(called.centromeric.matA))), matatotal / sum(unlist(called.chr.mata)) * (sum(unlist(called.chr.mata)) - sum(unlist(called.centromeric.mata))), matAtotal / sum(unlist(called.chr.matA)) * sum(unlist(called.centromeric.matA)), matatotal / sum(unlist(called.chr.mata)) * sum(unlist(called.centromeric.mata))  )

cent.all$domain <- c("euchromatin", "euchromatin", "centromeric", "centromeric")

#Note that need to take into account that when comparing to euchromatin, other H3K9 mutations need to be excluded

cent.all <- rbind(h3k9.all[1:2,], cent.all[3:4,])

cent.all$expected <- c(sum(cent.all$nmut[1],cent.all$nmut[3])/sum(cent.all$sites)*cent.all$sites[1], sum(cent.all$nmut[2],cent.all$nmut[4])/sum(cent.all$sites)*cent.all$sites[2], sum(cent.all$nmut[1],cent.all$nmut[3])/sum(cent.all$sites)*cent.all$sites[3], sum(cent.all$nmut[2],cent.all$nmut[4])/sum(cent.all$sites)*cent.all$sites[4])

model.cent.all <- brm(data = cent.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + domain,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.cent.all)

#cent.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-3,]
cent.rel.rates <- quantile(exp(post[,1])/exp(post[,2]), probs = c(0.5, 0.025, 0.975))

rate.cent.all <- data.frame(Domain = c("Centromeric"), t(cent.rel.rates))
colnames(rate.cent.all) <- c("Domain", "estimate", "lower", "upper")

##Plot the results
ggplot(rate.cent.all, aes(x = Domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Relative mutation rate")

#Duplicated regions
groups.dup <- group_by(bychr, duplicatedreg, mat)
dup.all <- summarise(groups.dup, nmut = n())

dup.all$sites <- c(sum(unlist(called.chr.matA)) - sum(unlist(called.duplicated.matA)), sum(unlist(called.chr.mata)) - sum(unlist(called.duplicated.mata)), sum(unlist(called.duplicated.matA)), sum(unlist(called.duplicated.mata)) )

dup.all$expected <- c(matAtotal / sum(unlist(called.chr.matA)) * (sum(unlist(called.chr.matA)) - sum(unlist(called.duplicated.matA))), matatotal / sum(unlist(called.chr.mata)) * (sum(unlist(called.chr.mata)) - sum(unlist(called.duplicated.mata))), matAtotal / sum(unlist(called.chr.matA)) * sum(unlist(called.duplicated.matA)), matatotal / sum(unlist(called.chr.mata)) * sum(unlist(called.duplicated.mata))  )

dup.all$domain <- c("euchromatin", "euchromatin", "duplicated", "duplicated")

model.dup.all <- brm(data = dup.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + domain,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.dup.all)

#dup.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-3,]
dup.rel.rates <- quantile(exp(post[,1])/exp(post[,2]), probs = c(0.5, 0.025, 0.975))

rate.dup.all <- data.frame(Domain = c("Duplicated"), t(dup.rel.rates))
colnames(rate.dup.all) <- c("Domain", "estimate", "lower", "upper")

##Plot the results
ggplot(rate.dup.all, aes(x = Domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Relative mutation rate")


#Mutations that were in excluded regions (H3K9 domains that overlapped with centromeric regions)
#Number of mutations by mating type
bychr.ex <- filter(bychr, ex.cent == 0) #Mutations in excluded regions excluded ; )

bychrgroups.ex <- group_by(bychr.ex, mat)
bychr.all.ex <- summarise(bychrgroups.ex, nmut = n())

matAtotal.ex.rm <- as.numeric(bychr.all.ex[1,2])
matatotal.ex.rm <- as.numeric(bychr.all.ex[2,2])

#Need to use the data in bychr
groups.h3k9.ex.cent <- group_by(bychr.ex, H3K9.ex.cent, mat)
h3k9.ex.cent.all <- summarise(groups.h3k9.ex.cent, nmut = n())

#called.H3K9.exc.matA #Contains sites for H3K9 domains with centrometic regions excluded
#called.exc.matA #Contains sites that have to be removed from total

h3k9.ex.cent.all$sites <- c(sum(unlist(called.chr.matA)) - sum(unlist(called.exc.matA)) - sum(unlist(called.H3K9.exc.matA)), sum(unlist(called.chr.mata)) - sum(unlist(called.exc.mata)) - sum(unlist(called.H3K9.exc.mata)), sum(unlist(called.H3K9.exc.matA)), sum(unlist(called.H3K9.exc.mata)) )

#Expected number of mutations for each domain
#(Total number of mutations / Total number of called sites, mat A) * called sites per domain
#- sum(unlist(called.exc.matA))

h3k9.ex.cent.all$expected <- c(matAtotal.ex.rm / (sum(unlist(called.chr.matA)) - sum(unlist(called.exc.matA)))  * (sum(unlist(called.chr.matA)) - sum(unlist(called.exc.matA)) - sum(unlist(called.H3K9.exc.matA))), matatotal.ex.rm / (sum(unlist(called.chr.mata)) - sum(unlist(called.exc.mata))) * (sum(unlist(called.chr.mata))  - sum(unlist(called.exc.mata)) - sum(unlist(called.H3K9.exc.mata))), matAtotal.ex.rm / (sum(unlist(called.chr.matA)) - sum(unlist(called.exc.matA))) * sum(unlist(called.H3K9.exc.matA)), matatotal.ex.rm / (sum(unlist(called.chr.mata)) - sum(unlist(called.exc.mata))) * sum(unlist(called.H3K9.exc.mata))  )

h3k9.ex.cent.all$domain <- c("euchromatin", "euchromatin", "H3K9 ex. centromeric", "H3K9 ex. centromeric")

model.H3K9.ex.cent.all <- brm(data = h3k9.ex.cent.all, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + domain,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.H3K9.ex.cent.all)

#h3k9.ex.cent.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-3,]
h3k9.ex.cent.rel.rates <- quantile(exp(post[,1])/exp(post[,2]), probs = c(0.5, 0.025, 0.975))

rate.h3k9.ex.cent.all <- data.frame(Domain = c("H3K9 ex. centromeric"), t(h3k9.ex.cent.rel.rates))
colnames(rate.h3k9.ex.cent.all) <- c("Domain", "estimate", "lower", "upper")

##Plot the results
ggplot(rate.h3k9.ex.cent.all, aes(x = Domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    ylab("Relative mutation rate")

##All domain results together
all.domain.results <- rbind(rate.h3k9.all, rate.h3k27.all, rate.cent.all, rate.dup.all, rate.h3k9.ex.cent.all)
#all.domain.results <- filter(all.domain.results, Domain != "Euchromatin")
#all.domain.results$Domain <- factor(all.domain.results$Domain)

#Results of relative mutation rate to euchromatin for different domains
save(all.domain.results, file = "./mutac_ms/data/alldomain_relrates.RData")

##Load the data to make the plot
##Plot the results
my.labels <- c("H3K9", "H3K27", "Centromeric", "Duplicated", "H3K9 ex.\ncentromeric")
domain.res.plot <- ggplot(all.domain.results, aes(x = Domain, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("") +    
    scale_x_discrete(labels = my.labels) +
    scale_y_continuous(breaks = seq(1,11,2)) +
    ylab("Mutation rate relative \n to euchromatin")

###For making the GC agains mutation rate plot

##For each mutation get the GC-content of local region (nearest 200 bp window midpoint)
gc.data$Mid <- gc.data$Start + 100
#gc.data <- filter(gc.data, Chromosome %in% chr.sizes[,1])
bychr$GCcont <- find.gc.window(bychr, gc.data)*100

genomegc <- nucstats.plot[,4]*100
#Can get bin counts by table(genomegc), then compare that to rest of the genome
#If GC content alone determines mutation rate, then relationship with GC-content is expected within all domains (Euchromatin, H3K9me, H3K27me etc.)

euchrom.gc <- filter(gc.data, Domain == "Euchromatic") #GC data for euchromatic regions
aineisto.euchrom <- filter(bychr, centromere == 0 & H3K27 == 0 & H3K9 == 0) #Mutations in euchrom

mutbygc.eu <- table(cut(aineisto.euchrom$GCcont, breaks = 19)) #Mutations in different bins of GC

#Prepare data for model fitting, calculate expected numbers of mutation for different bins
windows.euchrom.gc <- data.frame(GCcont = names(table(euchrom.gc$GCcont*100)), count = as.vector(table(euchrom.gc$GCcont)))

data.euchrom.gc <- data.frame(GCcont = names(mutbygc.eu), nmut = as.vector(mutbygc.eu))
data.euchrom.gc$count <- rep(0, nrow(data.euchrom.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.euchrom.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.euchrom.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.euchrom.gc$GCcont)) >= inte[1] & as.numeric(as.character(windows.euchrom.gc$GCcont)) <= inte[2]
    data.euchrom.gc$count[i] <- sum(windows.euchrom.gc$count[index])
}

#Calculate expected numbers of mutations
data.euchrom.gc$expected <- (sum(data.euchrom.gc$nmut) / sum(data.euchrom.gc$count)) * data.euchrom.gc$count

data.euchrom.gc$GC <- unname(sapply(strsplit(as.character(data.euchrom.gc$GCcont), ",|]"), '[[', 2))

model.euchrom.gc <- brm(data = data.euchrom.gc, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.euchrom.gc)

euchrom.gc.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-20,]
euchrom.gc.rel.rates[3,] <- c(NA, NA, NA) #No mutations observed here

rate.euchrom.gc.all <- data.frame(GCcont = data.euchrom.gc$GC, euchrom.gc.rel.rates)
colnames(rate.euchrom.gc.all) <- c("GCcont", "estimate", "lower", "upper")

##Plot the results
myylab <- c(TeX("Log_{10}(relative mutation rate)"))
rate.gc.euchrom.plot <- ggplot(rate.euchrom.gc.all, aes(x = GCcont, y = log10(estimate), ymin = log10(lower), ymax = log10(upper))) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_pointrange(colour = "deepskyblue") +
    ylab(myylab) +
    xlab("GC content (%)")

##Mutation rate and GC content for H3K9 regions
h3k9.gc <- filter(gc.data, Domain == "H3K9") #GC data for H3K9 regions
aineisto.h3k9 <- filter(bychr, H3K9 == 1) #Mutations in H3K9 regions

mutbygc.h3k9 <- table(cut(aineisto.h3k9$GCcont, breaks = 19))  #Mutations in different bins of GC

#Prepare data for model fitting, calculate expected numbers of mutation for different bins
windows.h3k9.gc <- data.frame(GCcont = names(table(h3k9.gc$GCcont*100)), count = as.vector(table(h3k9.gc$GCcont)))

data.h3k9.gc <- data.frame(GCcont = names(mutbygc.h3k9), nmut = as.vector(mutbygc.h3k9))
data.h3k9.gc$count <- rep(0, nrow(data.h3k9.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.h3k9.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.h3k9.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.h3k9.gc$GCcont)) >= inte[1] & as.numeric(as.character(windows.h3k9.gc$GCcont)) <= inte[2]
    data.h3k9.gc$count[i] <- sum(windows.h3k9.gc$count[index])
}

#Calculate expected numbers of mutations
data.h3k9.gc$expected <- (sum(data.h3k9.gc$nmut) / sum(data.h3k9.gc$count)) * data.h3k9.gc$count

data.h3k9.gc$GC <- unname(sapply(strsplit(as.character(data.h3k9.gc$GCcont), ",|]"), '[[', 2))

model.h3k9.gc <- brm(data = data.h3k9.gc, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.h3k9.gc)

h3k9.gc.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-20,]
h3k9.gc.rel.rates[18,] <- c(NA, NA, NA) #No mutations observed here

rate.h3k9.gc.all <- data.frame(GCcont = data.h3k9.gc$GC, h3k9.gc.rel.rates)
colnames(rate.h3k9.gc.all) <- c("GCcont", "estimate", "lower", "upper")

##Plot the results
myylab <- c(TeX("Log_{10}(relative mutation rate)"))
rate.gc.h3k9.plot <- ggplot(rate.h3k9.gc.all, aes(x = GCcont, y = log10(estimate), ymin = log10(lower), ymax = log10(upper))) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_pointrange(colour = "red") +
    ylab(myylab) +
    xlab("GC content (%)")


##Mutation rate and GC content for H3K27 regions
h3k27.gc <- filter(gc.data, Domain == "H3K27exK9") #GC data for H3K27 regions
aineisto.h3k27 <- filter(bychr, H3K27exK9 == 1) #Mutations in H3K27 regions

mutbygc.h3k27 <- table(cut(aineisto.h3k27$GCcont, breaks = 19))  #Mutations in different bins of GC

#Prepare data for model fitting, calculate expected numbers of mutation for different bins
windows.h3k27.gc <- data.frame(GCcont = names(table(h3k27.gc$GCcont*100)), count = as.vector(table(h3k27.gc$GCcont)))

data.h3k27.gc <- data.frame(GCcont = names(mutbygc.h3k27), nmut = as.vector(mutbygc.h3k27))
data.h3k27.gc$count <- rep(0, nrow(data.h3k27.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.h3k27.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.h3k27.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.h3k27.gc$GCcont)) >= inte[1] & as.numeric(as.character(windows.h3k27.gc$GCcont)) <= inte[2]
    data.h3k27.gc$count[i] <- sum(windows.h3k27.gc$count[index])
}

#Calculate expected numbers of mutations
data.h3k27.gc$expected <- (sum(data.h3k27.gc$nmut) / sum(data.h3k27.gc$count)) * data.h3k27.gc$count

data.h3k27.gc$GC <- unname(sapply(strsplit(as.character(data.h3k27.gc$GCcont), ",|]"), '[[', 2))

model.h3k27.gc <- brm(data = data.h3k27.gc, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.h3k27.gc)

h3k27.gc.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-20,]
h3k27.gc.rel.rates[c(2,3,4,5,6,8),] <- c(NA, NA, NA) #No mutations observed here

rate.h3k27.gc.all <- data.frame(GCcont = data.h3k27.gc$GC, h3k27.gc.rel.rates)
colnames(rate.h3k27.gc.all) <- c("GCcont", "estimate", "lower", "upper")

##Plot the results
myylab <- c(TeX("Log_{10}(relative mutation rate)"))
rate.gc.h3k27.plot <- ggplot(rate.h3k27.gc.all, aes(x = GCcont, y = log10(estimate), ymin = log10(lower), ymax = log10(upper))) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_pointrange(colour = "blue") +
    ylab(myylab) +
    xlab("GC content (%)")

##Mutation rate and GC content for centromeric regions
cent.gc <- filter(gc.data, Domain == "Centromeric") #GC data for centromeric regions
aineisto.cent <- filter(bychr, centromere == 1) #Mutations in centromeric regions

mutbygc.cent <- table(cut(aineisto.cent$GCcont, breaks = 19))  #Mutations in different bins of GC

#Prepare data for model fitting, calculate expected numbers of mutation for different bins
windows.cent.gc <- data.frame(GCcont = names(table(cent.gc$GCcont*100)), count = as.vector(table(cent.gc$GCcont)))

data.cent.gc <- data.frame(GCcont = names(mutbygc.cent), nmut = as.vector(mutbygc.cent))
data.cent.gc$count <- rep(0, nrow(data.cent.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.cent.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.cent.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.cent.gc$GCcont)) >= inte[1] & as.numeric(as.character(windows.cent.gc$GCcont)) <= inte[2]
    data.cent.gc$count[i] <- sum(windows.cent.gc$count[index])
}

#Calculate expected numbers of mutations
data.cent.gc$expected <- (sum(data.cent.gc$nmut) / sum(data.cent.gc$count)) * data.cent.gc$count

data.cent.gc$GC <- unname(sapply(strsplit(as.character(data.cent.gc$GCcont), ",|]"), '[[', 2))

model.cent.gc <- brm(data = data.cent.gc, family = poisson,
                        nmut ~ -1 + offset(log(expected)) + GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.cent.gc)

cent.gc.rel.rates <- t(apply(exp(post), 2, quantile, probs = c(0.5, 0.025, 0.975)))[-20,]
cent.gc.rel.rates[c(15,18),] <- c(NA, NA, NA) #No mutations observed here

rate.cent.gc.all <- data.frame(GCcont = data.cent.gc$GC, cent.gc.rel.rates)
colnames(rate.cent.gc.all) <- c("GCcont", "estimate", "lower", "upper")

##Plot the results
myylab <- c(TeX("Log_{10}(relative mutation rate)"))
rate.gc.cent.plot <- ggplot(rate.cent.gc.all, aes(x = GCcont, y = log10(estimate), ymin = log10(lower), ymax = log10(upper))) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_pointrange(colour = "grey") +
    ylab(myylab) +
    xlab("GC content (%)")

##Plot the results
gc.plot.final <- plot_grid(rate.gc.euchrom.plot, rate.gc.h3k9.plot, rate.gc.h3k27.plot, rate.gc.cent.plot, nrow = 4, labels = c("A", "B", "C", "D"), label_y = 1.02, align = "v")

save_plot("./mutac_ms/fig/relrate_gc.pdf", gc.plot.final, base_height=3.71*4, base_width=3.71*2*1.618)

##Need an analysis to check whether there is an independent effect of H3K9, centromeres and GC-content

#both.gc <- filter(gc.data, Domain == "Euchromatic" | Domain == "H3K27") #GC data for euchromatic and H3K27 regions
#aineisto.both <- filter(bychr, centromere == 0 & H3K27 == 0 | H3K27 == 1 & H3K9 == 0) #Mutations in euchrom and H3K27

euchr.gc <- filter(gc.data, Domain == "Euchromatic")
aineisto.euchr <- filter(bychr, centromere == 0 & H3K27 == 0 & H3K9 == 0)
aineisto.h3k27 <- filter(bychr, H3K27exK9 == 1)

mutbygc.euchr <- table(cut(aineisto.euchr$GCcont, breaks = c(9.5, 12, 14.5, 17, 19.5, 22, 24.5, 27, 29.5, 32, 34.5, 37, 39.5, 42, 44.5, 47, 49.5, 52, 54.5, 57, 59.5, 62, 64.5, 67, 69.5)))

mutbygc.h3k27 <- table(cut(aineisto.h3k27$GCcont, breaks = c(9.5, 12, 14.5, 17, 19.5, 22, 24.5, 27, 29.5, 32, 34.5, 37, 39.5, 42, 44.5, 47, 49.5, 52, 54.5, 57, 59.5, 62, 64.5, 67, 69.5)))

#mutbygc.both <- table(cut(aineisto.both$GCcont, breaks = c(9.5, 12, 14.5, 17, 19.5, 22, 24.5, 27, 29.5, 32, 34.5, 37, 39.5, 42, 44.5, 47, 49.5, 52, 54.5, 57, 59.5, 62, 64.5, 67, 69.5)))

mutbygc.h3k9 <- table(cut(aineisto.h3k9$GCcont, breaks = c(9.5, 12, 14.5, 17, 19.5, 22, 24.5, 27, 29.5, 32, 34.5, 37, 39.5, 42, 44.5, 47, 49.5, 52, 54.5, 57, 59.5, 62, 64.5, 67, 69.5)))

mutbygc.cent <- table(cut(aineisto.cent$GCcont, breaks = c(9.5, 12, 14.5, 17, 19.5, 22, 24.5, 27, 29.5, 32, 34.5, 37, 39.5, 42, 44.5, 47, 49.5, 52, 54.5, 57, 59.5, 62, 64.5, 67, 69.5)))

#Prepare data for model fitting, calculate expected numbers of mutation for different bins
#windows.both.gc <- data.frame(GCcont = names(table(both.gc$GCcont*100)), count = as.vector(table(both.gc$GCcont)))

#data.both.gc <- data.frame(GCcont = names(mutbygc.both), nmut = as.vector(mutbygc.both))
#data.both.gc$count <- rep(0, nrow(data.both.gc))

#Sum counts of intervals 
#for(i in 1:nrow(data.both.gc)) {
#    inte <- as.numeric(unlist(strsplit(as.character(data.both.gc[i,1]), "\\(|,|]"))[-1])
#    index <- as.numeric(as.character(windows.both.gc$GCcont)) > inte[1] & as.numeric(as.character(windows.both.gc$GCcont)) <= inte[2]
#    data.both.gc$count[i] <- sum(windows.both.gc$count[index])
#}

#data.both.gc$domain <- "Euchrom" #Note capitalize name here so that this is default level

#Euchromatin 
data.euchr.gc <- data.frame(GCcont = names(mutbygc.euchr), nmut = as.vector(mutbygc.euchr))
data.euchr.gc$count <- rep(0, nrow(data.euchr.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.euchr.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.euchr.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.euchrom.gc$GCcont)) > inte[1] & as.numeric(as.character(windows.euchrom.gc$GCcont)) <= inte[2]
    data.euchr.gc$count[i] <- sum(windows.euchrom.gc$count[index])
}

data.euchr.gc$domain <- "Euchrom"

data.h3k9.gc <- data.frame(GCcont = names(mutbygc.h3k9), nmut = as.vector(mutbygc.h3k9))
data.h3k9.gc$count <- rep(0, nrow(data.h3k9.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.h3k9.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.h3k9.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.h3k9.gc$GCcont)) > inte[1] & as.numeric(as.character(windows.h3k9.gc$GCcont)) <= inte[2]
    data.h3k9.gc$count[i] <- sum(windows.h3k9.gc$count[index])
}

data.h3k9.gc$domain <- "H3K9"

#H3K27
data.h3k27.gc <- data.frame(GCcont = names(mutbygc.h3k27), nmut = as.vector(mutbygc.h3k27))
data.h3k27.gc$count <- rep(0, nrow(data.h3k27.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.h3k27.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.h3k27.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.h3k27.gc$GCcont)) > inte[1] & as.numeric(as.character(windows.h3k27.gc$GCcont)) <= inte[2]
    data.h3k27.gc$count[i] <- sum(windows.h3k27.gc$count[index])
}

data.h3k27.gc$domain <- "H3K27"

data.cent.gc <- data.frame(GCcont = names(mutbygc.cent), nmut = as.vector(mutbygc.cent))
data.cent.gc$count <- rep(0, nrow(data.cent.gc))

#Sum counts of intervals 
for(i in 1:nrow(data.cent.gc)) {
    inte <- as.numeric(unlist(strsplit(as.character(data.cent.gc[i,1]), "\\(|,|]"))[-1])
    index <- as.numeric(as.character(windows.cent.gc$GCcont)) > inte[1] & as.numeric(as.character(windows.cent.gc$GCcont)) <= inte[2]
    data.cent.gc$count[i] <- sum(windows.cent.gc$count[index])
}

data.cent.gc$domain <- "Centromeric"

#Since centromeric and H3K9 overlap completely, need to adjust the numbers in H3K9 to avoid double counting mutations
data.h3k9.gc$nmut <- data.h3k9.gc$nmut - data.cent.gc$nmut
data.h3k9.gc$count <- data.h3k9.gc$count - data.cent.gc$count

data.all.gc <- rbind(data.euchr.gc, data.h3k27.gc, data.h3k9.gc, data.cent.gc)
data.all.gc$H3K9 <- c(rep(0, 48), rep(1, 48)) #H3K9 present (also in centromeric regions)
data.all.gc$Centromere <- c(rep(0,72), rep(1, 24)) #Centromeric regions
data.all.gc$H3K27 <- c(rep(0,24), rep(1,24), rep(0,48))

#Calculate expected numbers of mutations
data.all.gc$expected <- (sum(data.all.gc$nmut) / sum(data.all.gc$count)) * data.all.gc$count

data.all.gc$GC <- unname(sapply(strsplit(as.character(data.all.gc$GCcont), ",|]"), '[[', 2))
data.all.gc$GC <- as.numeric(data.all.gc$GC)

###Save the data file, so that it can be easily loaded later
#save(data.all.gc, file = "./mutac_ms/data/data.all.gc.RData")

##Doing some model comparison
model.all.gc1 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc2 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + H3K9,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc3 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + H3K9 + GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc4 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + H3K9 + GC + H3K9:GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc5 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + Centromere,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc6 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + Centromere + H3K9,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc7 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + Centromere + H3K9 + GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc8 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + Centromere + H3K9 + GC + H3K9:GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc9 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + Centromere + H3K9 + GC + H3K9:GC + Centromere:GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc10 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc11 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC + H3K27:GC,
                        prior = c(prior(normal(0, 10), class = b)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.all.gc0 <- brm(data = data.all.gc, family = poisson,
                        nmut ~ 1 + offset(log(count)),
                        prior = c(prior(normal(0, 10), class = Intercept)),
                        iter = 3000, warmup = 1000, chains = 4, cores = 4)


### Model comparisons  ###############################
##For some reason loo does not work with these models(?)
#model.all.gc0 <- add_criterion(model.all.gc0, "loo", reloo = T)
#model.all.gc1 <- add_criterion(model.all.gc1, "loo", reloo = T)
#model.all.gc2 <- add_criterion(model.all.gc2, "loo", reloo = T)
#model.all.gc3 <- add_criterion(model.all.gc3, "loo", reloo = T)
#model.all.gc4 <- add_criterion(model.all.gc4, "loo", reloo = T)
#model.all.gc5 <- add_criterion(model.all.gc5, "loo", reloo = T)
#model.all.gc6 <- add_criterion(model.all.gc6, "loo", reloo = T)
#model.all.gc7 <- add_criterion(model.all.gc7, "loo", reloo = T)
#model.all.gc8 <- add_criterion(model.all.gc8, "loo", reloo = T)
#model.all.gc9 <- add_criterion(model.all.gc9, "loo", reloo = T)

#using WAIC instead
model.all.gc0 <- add_criterion(model.all.gc0, "waic")
model.all.gc1 <- add_criterion(model.all.gc1, "waic")
model.all.gc2 <- add_criterion(model.all.gc2, "waic")
model.all.gc3 <- add_criterion(model.all.gc3, "waic")
model.all.gc4 <- add_criterion(model.all.gc4, "waic")
model.all.gc5 <- add_criterion(model.all.gc5, "waic")
model.all.gc6 <- add_criterion(model.all.gc6, "waic")
model.all.gc7 <- add_criterion(model.all.gc7, "waic")
model.all.gc8 <- add_criterion(model.all.gc8, "waic")
model.all.gc9 <- add_criterion(model.all.gc9, "waic")
model.all.gc10 <- add_criterion(model.all.gc10, "waic")
model.all.gc11 <- add_criterion(model.all.gc11, "waic")


mcomp.gc <- loo_compare(model.all.gc0, model.all.gc1, model.all.gc2, model.all.gc3, model.all.gc4, model.all.gc5, model.all.gc6, model.all.gc7, model.all.gc8, model.all.gc9, model.all.gc10, model.all.gc11, criterion = "waic")

print(mcomp.gc, simplify = F)

mweights.gc <- model_weights(model.all.gc0, model.all.gc1, model.all.gc2, model.all.gc3, model.all.gc4, model.all.gc5, model.all.gc6, model.all.gc7, model.all.gc8, model.all.gc9, model.all.gc10, model.all.gc11, weights = "waic")

GCmodelcomp <- cbind(mcomp.gc, sort(mweights.gc, decreasing = T))
colnames(GCmodelcomp)[9] <- "weight"

######################################################
GCmodelres <- fixef(model.all.gc10)

mutation.rsquare.res <- bayes_R2(model.all.gc10)


##Getting posterior predictions from the model
post.gc <- posterior_samples(model.all.gc10)

gc.summary <- data.frame(fitted(model.all.gc10))
gc.summary$GC <- data.all.gc$GC
gc.summary$domain <- data.all.gc$domain
gc.summary$nmut <- data.all.gc$nmut

#Comparing model estimates and data
ggplot(gc.summary, aes(x = GC, y = Estimate, ymin = Q2.5, ymax = Q97.5, colour = domain)) +
    geom_pointrange(shape = 21) +
    geom_point(data = gc.summary, aes(x = GC, y = nmut, colour = domain)) +
    ylab("Number of mutations")

###Rate estimates for Rest of genome and H3K9 related to GC-content
### Need to predict with the same count of intervals
pred.data <- data.frame(GC = rep(seq(15, 70, length.out = 100), 4), H3K27 = c(rep(0,100), rep(1,100), rep(0,200)), H3K9 = c(rep(0,200), rep(1,200)), Centromere = c(rep(0,300), rep(1,100)),  count = rep(1, 400))
gc.pred <- data.frame(fitted(model.all.gc10, newdata = pred.data))
gc.pred$GC <- pred.data$GC
gc.pred$H3K27 <- pred.data$H3K27
gc.pred$H3K9 <- pred.data$H3K9
gc.pred$Centromere <- pred.data$Centromere
gc.pred$domain <- c(rep("Euchrom", 100), rep("H3K27", 100), rep("H3K9", 100), rep("Centromere", 100))

gc.pred.plot <- ggplot(gc.pred, aes(x = GC, y = log(Estimate), ymin = log(Q2.5), ymax = log(Q97.5), colour = domain)) +
    geom_ribbon(aes(fill = domain), colour = NA, alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE) +    
    scale_colour_manual(breaks = c("Euchrom", "H3K27", "H3K9", "Centromere"), values = c(Euchrom = "deepskyblue", H3K27 = "blue", H3K9 = "red", Centromere = "grey"), labels = c("Euchromatin", "H3K27 domains", "H3K9 domains", "Centromeric"),  guide = guide_legend(label.position = "left")) +
    scale_fill_manual(breaks = c("Euchrom", "H3K27", "H3K9", "Centromere"), values = c(Euchrom = "deepskyblue", H3K27 = "blue", H3K9 = "red", Centromere = "grey"), labels = c("Euchromatin", "H3K27 domains", "H3K9 domains", "Centromeric"), guide = guide_legend(label.position = "left") )   +
    scale_x_continuous(expand = c(0,0)) +
    xlab("GC-content (%)") +
    ylab("log(model estimate)") +
    theme(legend.title = element_blank(), legend.position = c(0.45, 0.90))
        

### Saving the GC-model and model comparison results ####################
save(GCmodelcomp, GCmodelres, mutation.rsquare.res, file = "./mutac_ms/data/GCmodel.RData")
#save(

##Making a new plot in Finnish for a presentation
predFI <- ggplot(gc.pred, aes(x = GC, y = log(Estimate), ymin = log(Q2.5), ymax = log(Q97.5), colour = domain)) +
    geom_ribbon(aes(fill = domain), colour = NA, alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE) +    
    scale_colour_manual(breaks = c("Euchrom", "H3K27", "H3K9", "Centromere"), values = c(Euchrom = "deepskyblue", H3K27 = "blue", H3K9 = "red", Centromere = "grey"), labels = c("Eukromatiini", "H3K27 alueet", "H3K9 alueet", "Sentromeeri"),  guide = guide_legend(label.position = "left")) +
    scale_fill_manual(breaks = c("Euchrom", "H3K27", "H3K9", "Centromere"), values = c(Euchrom = "deepskyblue", H3K27 = "blue", H3K9 = "red", Centromere = "grey"), labels = c("Eukromatiini", "H3K27 alueet", "H3K9 alueet", "Sentromeeri"), guide = guide_legend(label.position = "left") )   +
    scale_x_continuous(expand = c(0,0)) +
    xlab("GC-pitoisuus (%)") +
    ylab("log(vaikutus mutaatiovauhtiin)") +
    theme(legend.title = element_blank(), legend.position = c(0.45, 0.90))

save_plot("./mutac_ms/fig/predFI.pdf", predFI)

#########################################################################

##Making a compound plot for mutations, domains, and GC-content

bottomrow <- plot_grid(domain.res.plot, gc.domains, gc.pred.plot, nrow = 1, rel_widths = c(1.3, 1.2, 0.8), align = "h", labels = c("B", "C", "D"), label_y = 1.02)

finalmutplot <- plot_grid(mutchr.plot, bottomrow, nrow = 2, labels = c("A", ""), rel_heights = c(6,3))

save_plot("./mutac_ms/fig/mutchrfinal.pdf", finalmutplot, base_height=9.5, base_width=9*1.618)

### *** Mutation spectrum in different domains
### **** Basefreqs in different domains
As <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "A"])
Ts <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "T"])
Cs <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "C"])
Gs <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "G"])

As.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "A"])
Ts.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "T"])
Cs.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "C"])
Gs.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "G"])

total.h3k9 <- c(As + As.a + Ts + Ts.a + Cs + Cs.a + Gs + Gs.a)
basefreq.h3k9 <- c(As + As.a, Ts + Ts.a, Cs + Cs.a, Gs + Gs.a)/total.h3k9

#Base frequencies in euchromatin
As.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "A"])
Ts.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "T"])
Cs.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "C"])
Gs.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "G"])

As.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "A"])
Ts.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "T"])
Cs.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "C"])
Gs.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "G"])

total.eu <- As.eu + As.eu.a + Ts.eu + Ts.eu.a + Cs.eu + Cs.eu.a + Gs.eu + Gs.eu.a
basefreq.eu <- c(As.eu + As.eu.a,  Ts.eu + Ts.eu.a,  Cs.eu + Cs.eu.a,  Gs.eu + Gs.eu.a)/total.eu

#Base freqs in total
As.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "A"])
Ts.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "T"])
Cs.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "C"])
Gs.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "G"])

As.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "A"])
Ts.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "T"])
Cs.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "C"])
Gs.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "G"])

total.chrs <- As.t + As.t2 + Ts.t + Ts.t2 + Cs.t + Cs.t2 + Gs.t + Gs.t2
basefreqs <- c(As.t + As.t2, Ts.t + Ts.t2, Cs.t + Cs.t2, Gs.t + Gs.t2) / total.chrs

### **** Mutation spectra in different domains
### Load the data from
#load(file = "./mutac_ms/data/triplet_domains.RData") #Load triplet data separated by domains
#contains aineisto.spectra and base frequencies for different domains
mutspec.temp <- rep(0, nrow(aineisto.spectra))

#Loop over all mutations and assign them to correct category
for(i in 1:nrow(aineisto.spectra)) {
    anc.base <- as.character(aineisto.spectra$anc.base[i])
    sample.base <- as.character(aineisto.spectra$sample.base[i])
    if((anc.base == "A" & sample.base == "T") | (anc.base == "T" & sample.base == "A")) {mutspec.temp[i] <- "A:T → T:A"}
    if((anc.base == "C" & sample.base == "G") | (anc.base == "G" & sample.base == "C")) {mutspec.temp[i] <- "C:G → G:C"}
    if((anc.base == "C" & sample.base == "A") | (anc.base == "G" & sample.base == "T")) {mutspec.temp[i] <- "C:G → A:T"}
    if((anc.base == "A" & sample.base == "C") | (anc.base == "T" & sample.base == "G")) {mutspec.temp[i] <- "A:T → C:G"}
    if((anc.base == "A" & sample.base == "G") | (anc.base == "T" & sample.base == "C")) {mutspec.temp[i] <- "A:T → G:C"}
    if((anc.base == "C" & sample.base == "T") | (anc.base == "G" & sample.base == "A")) {mutspec.temp[i] <- "C:G → T:A"}
}
                                                                                         
aineisto.spectra$mutation <- mutspec.temp #Store mutations

#Mutations by line
aineisto.spectra.indmut <- group_by(aineisto.spectra, Line, mutation, .drop = F)
spectra.res.lines <- summarise(aineisto.spectra.indmut, nmut = n())

#Mutations as aggregated
#For the whole genome
aineisto.spectra.agg <- group_by(aineisto.spectra, mutation, .drop = F)
aineisto.spectra.agg <- summarise(aineisto.spectra.agg, nmut = n())

#For H3K9 only
aineisto.spectra.agg.h3k9 <- group_by(filter(aineisto.spectra, H3K9 == 1), mutation, .drop = F)
aineisto.spectra.agg.h3k9 <- summarise(aineisto.spectra.agg.h3k9, nmut = n())

#Euchromatin
aineisto.spectra.agg.eu <- group_by(filter(aineisto.spectra, EUCHR == 1), mutation, .drop = F)
aineisto.spectra.agg.eu <- summarise(aineisto.spectra.agg.eu, nmut = n())

#Calculating expected frequencies
aineisto.spectra.agg$expected.freq <- c(rep((basefreqs[1]+basefreqs[2])/3,3), rep((basefreqs[3]+basefreqs[4])/3,3))
aineisto.spectra.agg.h3k9$expected.freq <- c(rep((basefreq.h3k9[1]+basefreq.h3k9[2])/3,3), rep((basefreq.h3k9[3]+basefreq.h3k9[4])/3,3))
aineisto.spectra.agg.eu$expected.freq <- c(rep((basefreq.eu[1]+basefreq.eu[2])/3,3), rep((basefreq.eu[3]+basefreq.eu[4])/3,3))

#expected numbers
aineisto.spectra.agg$expected <- aineisto.spectra.agg$expected.freq*sum(aineisto.spectra.agg$nmut)
aineisto.spectra.agg.h3k9$expected <- aineisto.spectra.agg.h3k9$expected.freq*sum(aineisto.spectra.agg.h3k9$nmut)
aineisto.spectra.agg.eu$expected <- aineisto.spectra.agg.eu$expected.freq*sum(aineisto.spectra.agg.eu$nmut)

#Fit models
model.spectra.all <- brm(data = aineisto.spectra.agg, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.spectra.all)
#Convert posterior samples back to normal scale and calculate intervals
spectra <- apply(exp(post[,-7]), 2, quantile, probs = c(0.5, 0.025, 0.975))

aineisto.spectra.agg$estimate <- t(spectra)[,1]
aineisto.spectra.agg$lower <- t(spectra)[,2]
aineisto.spectra.agg$upper <- t(spectra)[,3]
aineisto.spectra.agg$type <- c("transversion", "transition", "transversion", "transversion", "transversion", "transition")

#H3K9
model.spectra.h3k9 <- brm(data = aineisto.spectra.agg.h3k9, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.spectra.h3k9)
#Convert posterior samples back to normal scale and calculate intervals
spectra <- apply(exp(post[,-7]), 2, quantile, probs = c(0.5, 0.025, 0.975))
post.h3k9 <- exp(post[,-7])

aineisto.spectra.agg.h3k9$estimate <- t(spectra)[,1]
aineisto.spectra.agg.h3k9$lower <- t(spectra)[,2]
aineisto.spectra.agg.h3k9$upper <- t(spectra)[,3]
aineisto.spectra.agg.h3k9$type <- c("transversion", "transition", "transversion", "transversion", "transversion", "transition")

#Euchromatin
model.spectra.eu <- brm(data = aineisto.spectra.agg.eu, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.spectra.eu)
#Convert posterior samples back to normal scale and calculate intervals
spectra <- apply(exp(post[,-7]), 2, quantile, probs = c(0.5, 0.025, 0.975))
post.eu <- exp(post[,-7])
#save(post.h3k9, post.eu, file = "./mutac_ms/data/spectra_asex_posterior.RData")

aineisto.spectra.agg.eu$estimate <- t(spectra)[,1]
aineisto.spectra.agg.eu$lower <- t(spectra)[,2]
aineisto.spectra.agg.eu$upper <- t(spectra)[,3]
aineisto.spectra.agg.eu$type <- c("transversion", "transition", "transversion", "transversion", "transversion", "transition")

#Making a dataframe that combines all of the data
aineisto.spectra.all <- rbind(aineisto.spectra.agg, aineisto.spectra.agg.h3k9, aineisto.spectra.agg.eu)
aineisto.spectra.all$domain <- c(rep("Whole genome", 6), rep("H3K9", 6), rep("Euchromatin",6))

##Need to calculate differences in relative mutation rates between H3K9 and H3K9 excluded
post.h3k9 <- posterior_samples(model.spectra.h3k9)
post.eu <- posterior_samples(model.spectra.eu)

ratios <- exp(post.h3k9[,-7]) / exp(post.eu[,-7])
rat.res <- t(apply(ratios, 2, quantile, probs = c(0.5, 0.025, 0.975)))

ratio.2plot <- data.frame(mutation = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"), estimate = rat.res[c(3,1,5,4,2,6),1], lower = rat.res[c(3,1,5,4,2,6),2], upper = rat.res[c(3,1,5,4,2,6),3])
ratio.2plot$mutation <- factor(ratio.2plot$mutation, levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"))
ratio.2plot$type <- c(rep("Transversion", 4), rep("Transition", 2))

save(aineisto.spectra.all, ratio.2plot, file = "./mutac_ms/data/spectra_results.RData")

###Calculating Ts/Tv ratios for euchromatin and H3K9

#called bases
called.H3K9 <- (sum(unlist(lapply(called.H3K9.matA, sum))) + sum(unlist(lapply(called.H3K9.mata, sum))))/2

euchrom.matA <- sum(unlist(lapply(called.euchr.matA, sum)))
euchrom.mata <- sum(unlist(lapply(called.euchr.mata, sum)))
called.eu <- (euchrom.matA + euchrom.mata)/2

#Ts/tv ratios
#For euchromatin
euchrom.tstv <- summarise(group_by(aineisto.spectra.agg.eu, type), nmut = sum(nmut))
euchrom.tstv$bases <- called.eu

model.tstv.eu <- brm(data = euchrom.tstv, family = poisson,
                     nmut ~ -1 + offset(log(bases)) + type,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post.eu <- posterior_samples(model.tstv.eu)

tstv.eu <- exp(post.eu[,1])/exp(post.eu[,2])
tstv.eu.res <- quantile(tstv.eu, probs = c(0.5, 0.025, 0.975))

#For H3K9
h3k9.tstv <- summarise(group_by(aineisto.spectra.agg.h3k9, type), nmut = sum(nmut))
h3k9.tstv$bases <- called.H3K9

model.tstv.h3k9 <- brm(data = h3k9.tstv, family = poisson,
                     nmut ~ -1 + offset(log(bases)) + type,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post.h3k9 <- posterior_samples(model.tstv.h3k9)

tstv.h3k9 <- exp(post.h3k9[,1])/exp(post.h3k9[,2])
tstv.h3k9.res <- quantile(tstv.h3k9, probs = c(0.5, 0.025, 0.975))

save(tstv.eu.res, tstv.h3k9.res, file = "./mutac_ms/data/tstv_domains.RData")


##Plotting
load(file = "./mutac_ms/data/spectra_results.RData")
load(file = "./mutac_ms/data/spectra_domains_ratio_trinuc.RData") #In variable ratio.res

ratio.res$type <- c(rep("Transversion", 4), rep("Transition", 2))
ratio.res <- ratio.res[,c(4,1,2,3,5)] #Rearrange columns

ratio.2plot <- rbind(ratio.2plot, ratio.res)
ratio.2plot$correction <- c(rep("nucleotide", 6), rep("trinucleotide", 6)) 

library(latex2exp)

mylabels <- c(TeX("A:T $\\rightarrow$ T:A"), TeX("A:T $\\rightarrow$ C:G"), TeX("C:G $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ A:T"), TeX("A:T $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ T:A") )
aineisto.spectra.all$mutation <- factor(aineisto.spectra.all$mutation, levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"))

#Plot of the mutation spectra in different regions of the genome
p1 <- ggplot(aineisto.spectra.all, aes(y = estimate, ymin = lower, ymax = upper, x = mutation, fill = type)) +
    geom_bar(stat = "identity", colour = "black") +
    geom_errorbar(width = 0.1) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.5), breaks = c(seq(0,2.25, by = 0.25))) +
    scale_x_discrete(labels = mylabels) +
    ylab("Relative mutation rate") +
    xlab("Mutation") +
    geom_hline(yintercept = 1, lty = "dashed") +
    facet_wrap( ~ domain) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5))

#Plot of ratios between H3K9 and H3K9 excluded
p2 <- ggplot(data = ratio.2plot, aes(x = mutation, y = estimate, ymin = lower, ymax = upper, colour = type, shape = correction)) +
    geom_hline(yintercept = 1, lty = "dashed") +
    geom_pointrange(position = position_dodge2(width =0.5)) +
    scale_x_discrete(labels = mylabels) +
    scale_y_continuous(expand = c(0,0), limits = c(0,4)) +
    scale_shape_manual(values = c("nucleotide" = 20, "trinucleotide" = 23)) +
    xlab("Mutation") +
    ylab("Ratio of H3K9 / Euchromatin") +
    guides(color = guide_legend(override.aes = list(linetype = 0, size=1))) +   
    theme(legend.title = element_blank())

grDevices::cairo_pdf("./mutac_ms/fig/spectra_domains.pdf", width = 10, height = 8)
plot_grid(p1, p2, labels = c("A", "B"), align = "v", axis = "l", ncol = 1, rel_heights = c(1, 0.7))
dev.off()

### Mutation rate in in euchromatic and H3K9 domains
bp.eu <- total.eu/2
bp.h3k9 <- total.h3k9/2

mutperline.domains <- group_by(aineisto.spectra, Line, H3K9, EUCHR, .drop = F)
mutperline.domains <- summarise(mutperline.domains, nmut = n())
mutperline.domains <- filter(mutperline.domains, H3K9 == 1 | EUCHR == 1)
mutperline.domains$Domain <- factor(mutperline.domains$H3K9, levels = c(0,1), labels = c("euchromatin", "H3K9"))

model.rate.domains <- brm(data = mutperline.domains, family = poisson,
                     nmut ~ -1 + Domain,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post.rate.domains <- posterior_samples(model.rate.domains)

rate.euchromatin <- exp(post.rate.domains[,1]) / (bp.eu*transfers*mitoses)
rate.h3k9 <- exp(post.rate.domains[,2]) / (bp.h3k9*transfers*mitoses)

rate.euchromatin.res <- quantile(rate.euchromatin, probs = c(0.5, 0.025, 0.975))
rate.h3k9.res <- quantile(rate.h3k9, probs = c(0.5, 0.025, 0.975))

save(rate.euchromatin, rate.h3k9, file = "./mutac_ms/data/rates_domains_post.RData")
save(rate.euchromatin.res, rate.h3k9.res, file = "./mutac_ms/data/rates_domains.RData")

### **** Mutation spectra by taking trinucleotide frequencies into account

load("./mutac_ms/data/mutspectra.RData") #Load the data about point mutations, triplets included
load("./mutac_ms/data/mutspectra_domains.RData")

##Load trinucleotide frequencies
trinucfreq.wg <- read.table("Nc_trinuc_freq_whole_genome.txt", header = F)
trinucfreq.h3k9 <- read.table("Nc_trinuc_freq_H3K9.txt", header = F)
trinucfreq.eu <- read.table("Nc_trinuc_freq_euchromatin.txt", header = F)
trinucfreq.cent <- read.table("Nc_trinuc_freq_centromeres.txt", header = F)
trinucfreq.h3k27 <- read.table("Nc_trinuc_freq_H3K27_exK9.txt", header = F)

aineisto.spectra$mutation <- classify.mutations(aineisto.spectra$anc.base, aineisto.spectra$sample.base)[,1]
aineisto.spectra$tripclass <- classify.trinucleotides(aineisto.spectra$triplet)

#Over the whle genome
mutations96.wg <- classify.96mutations(aineisto.spectra, trinucfreq.wg)

#In euchromatin
mutations96.eu <- classify.96mutations(filter(aineisto.spectra, EUCHR == 1), trinucfreq.eu)

#In H3K9 domains
mutations96.h3k9 <- classify.96mutations(filter(aineisto.spectra, H3K9 == 1), trinucfreq.h3k9)

#Combine euchromatin and H3K9
mutations96.domains <- rbind(mutations96.eu, mutations96.h3k9)
mutations96.domains$domain <- c( rep("Euchromatin", 96), rep("H3K9me3", 96))

#Trying to estimate trinucleotide effects for all mutations separately is difficult
#Estimating them separately
mutations.domains.AT2TA <- filter(mutations96.domains, mutation == "A:T -> T:A")
mutations.domains.AT2CG <- filter(mutations96.domains, mutation == "A:T -> C:G")
mutations.domains.CG2GC <- filter(mutations96.domains, mutation == "C:G -> G:C")
mutations.domains.CG2AT <- filter(mutations96.domains, mutation == "C:G -> A:T")
mutations.domains.AT2GC <- filter(mutations96.domains, mutation == "A:T -> G:C")
mutations.domains.CG2TA <- filter(mutations96.domains, mutation == "C:G -> T:A")

model.domains.AT2TA <- brm(data = mutations.domains.AT2TA, family = poisson,
                     nmut ~ -1 + offset(log(count)) + domain + class,
                     prior = c(prior(normal(0, 4), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.domains.AT2CG <- brm(data = mutations.domains.AT2CG, family = poisson,
                     nmut ~ -1 + offset(log(count)) + domain + class,
                     prior = c(prior(normal(0, 4), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.domains.CG2GC <- brm(data = mutations.domains.CG2GC, family = poisson,
                     nmut ~ -1 + offset(log(count)) + domain + class,
                     prior = c(prior(normal(0, 4), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.domains.CG2AT <- brm(data = mutations.domains.CG2AT, family = poisson,
                     nmut ~ -1 + offset(log(count)) + domain + class,
                     prior = c(prior(normal(0, 4), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.domains.AT2GC <- brm(data = mutations.domains.AT2GC, family = poisson,
                     nmut ~ -1 + offset(log(count)) + domain + class,
                     prior = c(prior(normal(0, 4), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.domains.CG2TA <- brm(data = mutations.domains.CG2TA, family = poisson,
                     nmut ~ -1 + offset(log(count)) + domain + class,
                     prior = c(prior(normal(0, 4), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

##Posteriors
post.eu <- cbind(posterior_samples(model.domains.AT2TA)[,1], posterior_samples(model.domains.AT2CG)[,1], posterior_samples(model.domains.CG2GC)[,1], posterior_samples(model.domains.CG2AT)[,1], posterior_samples(model.domains.AT2GC)[,1], posterior_samples(model.domains.CG2TA)[,1])
colnames(post.eu) <- c("AT2TA", "AT2CG", "CG2GC", "CG2AT", "AT2GC", "CG2TA")

post.h3k9 <- cbind(posterior_samples(model.domains.AT2TA)[,2], posterior_samples(model.domains.AT2CG)[,2], posterior_samples(model.domains.CG2GC)[,2], posterior_samples(model.domains.CG2AT)[,2], posterior_samples(model.domains.AT2GC)[,2], posterior_samples(model.domains.CG2TA)[,2])
colnames(post.h3k9) <- c("AT2TA", "AT2CG", "CG2GC", "CG2AT", "AT2GC", "CG2TA")

mut.eu <- apply(post.eu, 2,  exp)
mut.eu <- mut.eu / mean(mut.eu) #Calculate relative mutation rates

mutrate.eu <- data.frame(estimate = apply(mut.eu, 2,  median), HPDinterval(as.mcmc(mut.eu)), mutation = c("AT2TA", "AT2CG", "CG2GC", "CG2AT", "AT2GC", "CG2TA"), type = c(rep("transversion", 4), rep("transition", 2)), domain = rep("Euchromatin",6))

mut.h3k9 <- apply(post.h3k9, 2,  exp)
mut.h3k9 <- mut.h3k9 / mean(mut.h3k9) #Calculate relative mutation rates

mutrate.h3k9 <- data.frame(estimate = apply(mut.h3k9, 2,  median), HPDinterval(as.mcmc(mut.h3k9)), mutation = c("AT2TA", "AT2CG", "CG2GC", "CG2AT", "AT2GC", "CG2TA"), type = c(rep("transversion", 4), rep("transition", 2)), domain = rep("H3K9me3",6))

mutrate.all <- rbind(mutrate.eu, mutrate.h3k9)

#Can look at the data...
ggplot(mutrate.all, aes(x = mutation, y = estimate, ymin = lower, ymax = upper, fill = type)) +
    geom_bar(stat = "identity") +
    facet_wrap( ~ domain)


#Calculate ratios
ratios <- mut.h3k9 / mut.eu

#mylabels <- c(TeX("A:T $\\rightarrow$ T:A"), TeX("A:T $\\rightarrow$ C:G"), TeX("C:G $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ A:T"), TeX("A:T $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ T:A") )

aineisto.spectra.all$mutation <- factor(aineisto.spectra.all$mutation, levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"))

grDevices::cairo_pdf("./mutac_ms/fig/spectra_domains.pdf", width = 10, height = 8)
plot_grid(p1, p2, labels = c("A", "B"), align = "v", axis = "l", ncol = 1, rel_heights = c(1, 0.7))
dev.off()

ratio.res <- data.frame(estimate = apply(ratios, 2,  median), HPDinterval(as.mcmc(ratios)), mutation = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"), type = c(rep("transversion", 4), rep("transition", 2)))

#Save results file
save(ratio.res, file = "./mutac_ms/data/spectra_domains_ratio_trinuc.RData")

#Making the plot
grDevices::cairo_pdf("./mutac_ms/fig/spectra_domains_trinuc.pdf", width = 10, height = 4)
ggplot(ratio.res, aes(x = fct_relevel(mutation, "A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A" ), y = estimate, ymin = lower, ymax = upper, colour = type)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = "dashed") +
    xlab("Mutation") +
    ylab("Ratio of H3K9 / Euchromatin") +
    theme(legend.title = element_blank())
dev.off()


### **** Mutation spectrum in line 21 versus rest of the MA lines

#Line 21 had an excess number of mutations

load("./mutac_ms/data/mutspectra.RData") #Load the data about point mutations, triplets included
load("./mutac_ms/data/mutspectra_domains.RData")
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/calledsites.RData")
#Nucleotide frequencies

#Base freqs in total
As.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "A"])
Ts.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "T"])
Cs.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "C"])
Gs.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "G"])

As.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "A"])
Ts.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "T"])
Cs.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "C"])
Gs.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "G"])

total.chrs <- As.t + As.t2 + Ts.t + Ts.t2 + Cs.t + Cs.t2 + Gs.t + Gs.t2
basefreqs <- c(As.t + As.t2, Ts.t + Ts.t2, Cs.t + Cs.t2, Gs.t + Gs.t2) / total.chrs


aineisto.line21 <- filter(aineisto.spectra, Line == "L21G40")
aineisto.rest <- filter(aineisto.spectra, Line != "L21G40")

aineisto.line21 <- cbind(aineisto.line21, classify.mutations(aineisto.line21$anc.base, aineisto.line21$sample.base))
aineisto.rest <- cbind(aineisto.rest, classify.mutations(aineisto.rest$anc.base, aineisto.rest$sample.base))

mutations.line21 <- group_by(aineisto.line21, mutation, .drop = F)
mutations.line21 <- summarise(mutations.line21, nmut = n())

mutations.rest <- group_by(aineisto.rest, mutation, .drop = F)
mutations.rest <- summarise(mutations.rest, nmut = n())

#Base frequencies
mutations.line21$basefreq <- c(rep((basefreqs[1]+basefreqs[2])/3,3), rep((basefreqs[3]+basefreqs[4])/3,3))
mutations.rest$basefreq <- c(rep((basefreqs[1]+basefreqs[2])/3,3), rep((basefreqs[3]+basefreqs[4])/3,3))

#Expected numbers of mutations
mutations.line21$expected <- mutations.line21$basefreq*sum(mutations.line21$nmut)
mutations.rest$expected <- mutations.rest$basefreq*sum(mutations.rest$nmut)

##Fit models to calculate relative mutation spectra
model.rest <- brm(data = mutations.rest, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 5), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.line21 <- brm(data = mutations.line21, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 5), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post.rest <- posterior_samples(model.rest)[,-7]
post.line21 <- posterior_samples(model.line21)[,-7]

rate.rest <- apply(post.rest, 2, exp)
rate.line21 <- apply(post.line21, 2, exp)

c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A")

#Spectrum for rest of the lines
res.rest <- data.frame(mutation = c("A:T → C:G", "A:T → G:C", "A:T → T:A", "C:G → A:T", "C:G → G:C", "C:G → T:A"), estimate = apply(rate.rest, 2, median), HPDinterval(as.mcmc(rate.rest)), type = c("Transversion", "Transition", "Transversion", "Transversion", "Transversion", "Transition"))

#Spectrum for line 21
res.line21 <- data.frame(mutation =c("A:T → C:G", "A:T → G:C", "A:T → T:A", "C:G → A:T", "C:G → G:C", "C:G → T:A"), estimate = apply(rate.line21, 2, median), HPDinterval(as.mcmc(rate.line21)), type = c("Transversion", "Transition", "Transversion", "Transversion", "Transversion", "Transition"))

#Ratios of mutations rates Line 21 / Rest of the MA lines
rate.ratio <- rate.line21 / rate.rest

res.ratio <- data.frame(mutation = c("A:T → C:G", "A:T → G:C", "A:T → T:A", "C:G → A:T", "C:G → G:C", "C:G → T:A"), estimate = apply(rate.ratio, 2, median), HPDinterval(as.mcmc(rate.ratio)), type = c("Transversion", "Transition", "Transversion", "Transversion", "Transversion", "Transition"))

#Making a plot, with spectra of line 21 and comparison to rest
mylabels <- c(TeX("A:T $\\rightarrow$ T:A"), TeX("A:T $\\rightarrow$ C:G"), TeX("C:G $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ A:T"), TeX("A:T $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ T:A") )
aineisto.spectra.all$mutation <- factor(aineisto.spectra.all$mutation, levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"))

p1 <- ggplot(res.line21, aes(y = estimate, ymin = lower, ymax = upper, x = fct_relevel(mutation, "A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"), fill = type)) +
    geom_bar(stat = "identity", colour = "black") +
    geom_errorbar(width = 0.1) +
    scale_y_continuous(expand = c(0,0), limits = c(0,3)) +
    scale_x_discrete(labels = mylabels) +    
    ylab("Relative mutation rate") +
    xlab("Mutation") +
    geom_hline(yintercept = 1, lty = "dashed") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5))    

p2 <- ggplot(res.ratio, aes(y = estimate, ymin = lower, ymax = upper, x = fct_relevel(mutation, "A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"), colour = type)) +
           geom_pointrange() +
           ylab("Ratio of Line 21 / Other MA lines")  +
           xlab("Mutation") +
           geom_hline(yintercept = 1, lty = "dashed") +
           theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

grDevices::cairo_pdf("./mutac_ms/fig/spectra_domains_L21.pdf", width = 9, height = 5)
plot_grid(p1, p2, labels = c("A", "B"), align = "h", axis = "b", nrow = 1, rel_widths = c(0.5,1))
dev.off()



### * Genetic diversity in natural strains

### ** Load the data

library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

#Check instructions for using the data.table from:
#https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html

#Read data (using the fread function)

#This load the first chromosome
natvar.chr1 <- fread("/mnt/HDD2/NGS/natpop/natpop_allsites_chr1.txt", header = T, sep = "\t", na.strings = "NA")
#2nd
natvar.chr2 <- fread("/mnt/HDD2/NGS/natpop/natpop_allsites_chr2.txt", header = T, sep = "\t", na.strings = "NA")
#3rd
natvar.chr3 <- fread("/mnt/HDD2/NGS/natpop/natpop_allsites_chr3.txt", header = T, sep = "\t", na.strings = "NA")
#4th
natvar.chr4 <- fread("/mnt/HDD2/NGS/natpop/natpop_allsites_chr4.txt", header = T, sep = "\t", na.strings = "NA")
#5th
natvar.chr5 <- fread("/mnt/HDD2/NGS/natpop/natpop_allsites_chr5.txt", header = T, sep = "\t", na.strings = "NA")
#6th
natvar.chr6 <- fread("/mnt/HDD2/NGS/natpop/natpop_allsites_chr6.txt", header = T, sep = "\t", na.strings = "NA")
#7th
natvar.chr7 <- fread("/mnt/HDD2/NGS/natpop/natpop_allsites_chr7.txt", header = T, sep = "\t", na.strings = "NA")



#Calculate theta over the chromosome
theta.chr1 <- summary.over.intervals(natvar.chr1, interval = 200)
theta.chr2 <- summary.over.intervals(natvar.chr2, interval = 200)
theta.chr3 <- summary.over.intervals(natvar.chr3, interval = 200)
theta.chr4 <- summary.over.intervals(natvar.chr4, interval = 200)
theta.chr5 <- summary.over.intervals(natvar.chr5, interval = 200)
theta.chr6 <- summary.over.intervals(natvar.chr6, interval = 200)
theta.chr7 <- summary.over.intervals(natvar.chr7, interval = 200)

#Calculating theta over a different window size (400 bp)
theta.chr1 <- summary.over.intervals(natvar.chr1, interval = 400)
theta.chr2 <- summary.over.intervals(natvar.chr2, interval = 400)
theta.chr3 <- summary.over.intervals(natvar.chr3, interval = 400)
theta.chr4 <- summary.over.intervals(natvar.chr4, interval = 400)
theta.chr5 <- summary.over.intervals(natvar.chr5, interval = 400)
theta.chr6 <- summary.over.intervals(natvar.chr6, interval = 400)
theta.chr7 <- summary.over.intervals(natvar.chr7, interval = 400)

#Calculating theta over a different window size (600 bp)
theta.chr1 <- summary.over.intervals(natvar.chr1, interval = 600)
theta.chr2 <- summary.over.intervals(natvar.chr2, interval = 600)
theta.chr3 <- summary.over.intervals(natvar.chr3, interval = 600)
theta.chr4 <- summary.over.intervals(natvar.chr4, interval = 600)
theta.chr5 <- summary.over.intervals(natvar.chr5, interval = 600)
theta.chr6 <- summary.over.intervals(natvar.chr6, interval = 600)
theta.chr7 <- summary.over.intervals(natvar.chr7, interval = 600)

#Calculating theta over a different window size (800 bp)
theta.chr1 <- summary.over.intervals(natvar.chr1, interval = 800)
theta.chr2 <- summary.over.intervals(natvar.chr2, interval = 800)
theta.chr3 <- summary.over.intervals(natvar.chr3, interval = 800)
theta.chr4 <- summary.over.intervals(natvar.chr4, interval = 800)
theta.chr5 <- summary.over.intervals(natvar.chr5, interval = 800)
theta.chr6 <- summary.over.intervals(natvar.chr6, interval = 800)
theta.chr7 <- summary.over.intervals(natvar.chr7, interval = 800)

#Calculating theta over a different window size (1000 bp)
theta.chr1 <- summary.over.intervals(natvar.chr1, interval = 1000)
theta.chr2 <- summary.over.intervals(natvar.chr2, interval = 1000)
theta.chr3 <- summary.over.intervals(natvar.chr3, interval = 1000)
theta.chr4 <- summary.over.intervals(natvar.chr4, interval = 1000)
theta.chr5 <- summary.over.intervals(natvar.chr5, interval = 1000)
theta.chr6 <- summary.over.intervals(natvar.chr6, interval = 1000)
theta.chr7 <- summary.over.intervals(natvar.chr7, interval = 1000)


theta.chr.all <- rbind(theta.chr1, theta.chr2, theta.chr3, theta.chr4, theta.chr5, theta.chr6, theta.chr7)
theta.chr.all <- theta.chr.all[,-3] #Drop the anx column which was not used
interval <- 1000 #CHANGE VALUE HERE!
theta.chr.all$mid <- theta.chr.all$start + interval/2 #Coordinate of window mid point
theta.chr.all$Chromosome <- c(rep("Supercontig_12.1", nrow(theta.chr1)), rep("Supercontig_12.2", nrow(theta.chr2)), rep("Supercontig_12.3", nrow(theta.chr3)), rep("Supercontig_12.4", nrow(theta.chr4)), rep("Supercontig_12.5", nrow(theta.chr5)), rep("Supercontig_12.6", nrow(theta.chr6)), rep("Supercontig_12.7", nrow(theta.chr7)))


#Load the GC-data for correct interval
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data.RData") #200 bp interval
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data400.RData")
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data600.RData")
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data800.RData")
load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/gc.data1000.RData")

#Check that CG-content, presence of different domains are in the data correct
#The variable gc.data already contains GC content and domains
gc.data$Mid <- gc.data$Start + interval/2 #midpoint of window
gc.data$thetaW <- combine.gc.theta(gc.data, theta.chr.all) #Get theta values for each window

#Get information about GC content and chromatin domains
theta.stats <- combine.theta.all(theta.chr.all, nucstats.plot, cent, h3k9.plot, h3k27.plot)

theta.chr.all <- cbind(theta.chr.all, theta.stats)

#Need to fix H3K27, for regions where H3K9 and H3K27 overlap
theta.chr.all$H3K27exK9 <- rep(0, nrow(theta.chr.all))
for(i in 1:nrow(theta.chr.all)) {
    theta.chr.all$H3K27exK9[i] <- ifelse(theta.chr.all$H3K27[i] == 1 & theta.chr.all$H3K9[i] == 0, 1, 0) }


#Save the theta.W data
#save(theta.chr.all, gc.data, file = "/mnt/HDD2/NGS/natpop/thetaW.RData") #This has also been saved to /.../MA_WGS/
#save(theta.chr.all, gc.data, file = "/mnt/HDD2/NGS/natpop/thetaW_400.RData") #This has also been saved to /.../MA_WGS/
#save(theta.chr.all, gc.data, file = "/mnt/HDD2/NGS/natpop/thetaW_600.RData") #This has also been saved to /.../MA_WGS/
#save(theta.chr.all, gc.data, file = "/mnt/HDD2/NGS/natpop/thetaW_800.RData") #This has also been saved to /.../MA_WGS/
save(theta.chr.all, gc.data, file = "/mnt/HDD2/NGS/natpop/thetaW_1000.RData") #This has also been saved to /.../MA_WGS/

### ** Analysis of natural genetic variation and mutation rate

#load(file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/thetaW.RData")

#Get only variable windows
#gc.data.var <- filter(gc.data, thetaW > 0) #Removing windows where thetaW = 0
gc.data.var <- gc.data

myylab <- c(TeX("$\\theta_W$"))
#Plotting theta.W in different domains
theta.domains <- ggplot(filter(gc.data.var, Domain != "H3K27"), aes(x = Domain, y = thetaW)) +
    geom_violin() +
    geom_boxplot(width = 0.15, alpha = 0.1, outlier.shape = NA) +
    ylab(myylab) +
    xlab("") +
    scale_x_discrete(labels = c("Centromeric", "Euchromatic", "H3K27", "H3K9", "H3K9 ex. centromeric")) +
    coord_flip()

#Calculate theta.W for different domains:
gc.data.groups <- group_by(filter(gc.data.var, Domain != "H3K27"), Domain)
#theta.summary <- summarise(gc.data.groups, nobs = n(), theta = quantile(thetaW, probs = 0.5), theta.min = quantile(thetaW, probs = 0.025), theta.max = quantile(thetaW, probs = 0.975))

theta.summary <- summarise(gc.data.groups, nobs = n(), theta = quantile(thetaW, probs = 0.5, na.rm = T), theta.min = quantile(thetaW, probs = 0.025, na.rm = T), theta.max = quantile(thetaW, probs = 0.975, na.rm = T))

#save(theta.summary, file = "./mutac_ms/data/thetaSummary.RData")

##Relationship between thetaW and GC-content
theta.gc <- ggplot(filter(theta.chr.all, theta.W > 0), aes(x = GC*100, y = theta.W)) +
   geom_hex() +
   scale_fill_viridis_c() +
   geom_smooth(method = "lm") +
   #geom_density_2d(colour = "black", bins = 100) +
   ylab(myylab) +
   xlab("GC-content (%)")

save_plot("./mutac_ms/fig/thetagc.pdf", theta.gc)

##Predicted mutation rate and theta
load("./mutac_ms/data/GCmodel.RData")

#Calculate predicted mutation rates for each window
#Should maybe consider the error associated with model predictions (but so many points, does it make a difference?)
#nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC,
#GC content was from 0 to 100 in the model fit
predmutrate <- GCmodelres[1,1] + GCmodelres[2,1]*theta.chr.all[,7] + GCmodelres[3,1]*theta.chr.all[,8] + GCmodelres[4,1]*theta.chr.all[,10] + GCmodelres[5,1]*theta.chr.all[,6]*100 + GCmodelres[6,1]*theta.chr.all[,6]*100*theta.chr.all[,8]

theta.chr.all$predmut <- predmutrate
#theta.plot <- filter(theta.chr.all, theta.W > 0)
theta.plot <- theta.chr.all

colnames(theta.plot)[3] <- "thetaW"

#This takes very long because the dataset is so large
model.theta <- brm(data = theta.plot, family = gaussian,
                    thetaW ~ 1 + predmut,
                    prior = c(prior(normal(0, 10), class = b)),
                    iter = 3000, warmup = 1000, chains = 4, cores = 4)

model.theta.res <- fixef(model.theta) #Model results

##This requires lot of memory!######
#theta.rsquare.res <- bayes_R2(model.theta)
####################################



#Model prediction for drawing slope
pred.data <- data.frame(predmut = seq(-10, 0, length.out = 100))
theta.pred <- data.frame(fitted(model.theta, newdata = pred.data))
theta.pred$predmut <- pred.data$predmut

##Load data calculated with cluster
#load("./mutac_ms/data/thetaModel.RData")
load("./mutac_ms/data/thetaModel_new200.RData")
    
theta.pred.plot <- ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    #geom_point(alpha = 0.05)
    geom_hex() +
    scale_fill_viridis_c() +
    geom_abline(intercept = model.theta.res[1,1], slope = model.theta.res[2,1], lwd = 1, color = "blue") +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    #geom_smooth(method = "lm") +
    ylab(myylab) +
    xlab("log (predicted mutation rate)")    
#Note that error in the slope is not visible, because it is so small...

ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    geom_point(alpha = 0.05)
    #geom_hex() +
    #scale_fill_viridis_c() +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    geom_smooth(method = "lm") +
    ylab(myylab) +
    xlab("log (predicted mutation rate)")

ggplot(theta.plot, aes(x = predmut)) +
    geom_histogram()

theta.plot$Domain <- rep(0, nrow(theta.plot))
for(i in 1:nrow(theta.plot)) {
    if(theta.plot$H3K9[i] == 1) {theta.plot$Domain[i] <- "H3K9"}
    if(theta.plot$centromere[i] == 1) {theta.plot$Domain[i] <- "Centromeric"}
    if(theta.plot$H3K27exK9[i] == 1) {theta.plot$Domain[i] <- "H3K27"}
    if(theta.plot$H3K9[i] == 0 & theta.plot$H3K27exK9[i] == 0) {theta.plot$Domain[i] <- "Euchromatic"}
}

ggplot(theta.plot, aes(x = Domain, y = predmut)) +
    geom_boxplot()

ggplot(theta.plot, aes(x = Domain, y = thetaW)) +
    geom_boxplot()

ggplot(theta.plot, aes(x = GC*100, y = thetaW)) +
    geom_point() +
    facet_grid(~ Domain)

ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    geom_point() +
    facet_grid(~ Domain)

ggplot(theta.plot, aes(x = predmut, y = GC*100)) +
    geom_point() +
    facet_grid(~ Domain)

### Ordinary linear model, to please the reviewers ################
test.theta <- lm(thetaW ~ predmut, data = theta.plot)
summary(test.theta) #R^2 = 0.22, seems quite high even!
cor.test(theta.plot$predmut, theta.plot$thetaW, use = "complete")
###################################################################


### Does the model amount of natural variation within different domains?
myylab <- c(TeX("$\\theta_W$"))
ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    #geom_point(alpha = 0.05)
    geom_hex() +
    scale_fill_viridis_c() +
    #geom_abline(intercept = model.theta.res[1,1], slope = model.theta.res[2,1], lwd = 1, color = "blue") +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    geom_smooth(method = "lm") +
    ylab(myylab) +
    xlab("log (predicted mutation rate)") +
    facet_wrap(~ Domain)

##Yes! Except that within centromeric domains relationship is weak.

##Are differences in thetaW in different domains significant?
##Ordinary linear model for thetaW in different domains###
###Need to set Euchromatin as the first
theta.plot$Domain <- fct_relevel(theta.plot$Domain, "Euchromatic", "Centromeric", "H3K9", "H3K27")
test.thetadomain <- lm(thetaW ~ Domain, data = theta.plot)
summary(test.thetadomain)
summary(lm(thetaW ~ -1 + Domain, data = theta.plot))
##########################################################

##Run brms model
### This takes very long because the dataset is so large ###
model.domaintheta <- brm(data = theta.plot, family = gaussian,
                    thetaW ~ 1 + Domain,
                    prior = c(prior(normal(0, 10), class = b)),
                    iter = 3000, warmup = 1000, chains = 4, cores = 4)

domainthetares <- fixef(model.domaintheta)
save(domainthetares, file = "./mutac_ms/data/domainthetares.RData")
###################################################################



final.theta <- plot_grid(theta.domains, theta.pred.plot, ncol = 2, labels = c("A", "B"), align = "h")

#save_plot(filename = "./mutac_ms/fig/thetares.pdf", final.theta, base_height=3.71, base_width=3.71*1.8*1.618)
save_plot(filename = "./mutac_ms/fig/thetares_new.pdf", final.theta, base_height=3.71, base_width=3.71*1.8*1.618)

#Plot theta.W across chromosomes
#need to load cent, h3k9.plot, h3k27.plot from earlier

labels <- c(Supercontig_12.1 = "Chr I", Supercontig_12.2 = "Chr II", Supercontig_12.3 = "Chr III", Supercontig_12.4 = "Chr IV", Supercontig_12.5 = "Chr V", Supercontig_12.6 = "Chr VI", Supercontig_12.7 = "Chr VII")
myylab <- c(TeX("$\\theta_W$"))

theta_chr.plot <- ggplot(theta.chr.all, aes()) +
    geom_line(aes(x = mid, y = theta.W)) +
     geom_rect(data = cent, aes(xmin = start, xmax = end, ymin = -0.01, ymax = -0.05, fill = "Centromeric")) +    
    geom_rect(data = h3k9.plot, aes(xmin = Start, xmax = End, ymax = -0.06, ymin = -0.1, fill = "H3K9me")) +
    geom_rect(data = h3k27.plot, aes(xmin = Start, xmax = End, ymax = -0.11, ymin = -0.15, fill = "H3K27me3")) +
    scale_y_continuous(breaks = seq(0,0.3, 0.1)) +    
    scale_fill_manual(breaks = c("Centromeric", "H3K27me3", "H3K9me"), values = c(Centromeric = "grey", H3K27me3 = "blue", H3K9me = "red")) +    
    xlab("Position (bp)") +
    ylab(myylab) +    
    facet_grid(Chromosome ~ ., labeller = labeller(Chromosome = labels)) +
    theme(legend.position = "top", legend.justification = c(0.5, 0.5), legend.title = element_blank())


save_plot("./mutac_ms/fig/theta_chr.pdf", theta_chr.plot, base_height = 7, base_width = 9*1.618)

### *** Genetic diversity and mutation rate correlation for different window sizes

load("./mutac_ms/data/GCmodel.RData")


##Load data one window size at a time

### 200 bp #####################################################################
load(file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/thetaW.RData")
load("./mutac_ms/data/thetaModel_new200.RData") #Load data for 200 bp window

#Calculate predicted mutation rates for each window
#nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC,
#GC content was from 0 to 100 in the model fit
predmutrate <- GCmodelres[1,1] + GCmodelres[2,1]*theta.chr.all[,7] + GCmodelres[3,1]*theta.chr.all[,8] + GCmodelres[4,1]*theta.chr.all[,10] + GCmodelres[5,1]*theta.chr.all[,6]*100 + GCmodelres[6,1]*theta.chr.all[,6]*100*theta.chr.all[,8]

theta.chr.all$predmut <- predmutrate
theta.plot <- theta.chr.all
colnames(theta.plot)[3] <- "thetaW"

myylab <- c(TeX("$\\theta_W$"))
blabel <- paste0("$\\beta = ", round(model.theta.res[2,1],4), "\\, \\[", round(model.theta.res[2,3],4), ", \\, ", round(model.theta.res[2,4],4), "\\]$")
r2label <- paste0("$R^2 = ", round(theta.rsquare.res[1,1],2), "\\, \\[", round(theta.rsquare.res[1,3],2), ", \\, ", round(theta.rsquare.res[1,4],2), "\\]$")
theta.plot200 <- ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    #geom_point(alpha = 0.05)
    geom_hex() +
    scale_fill_viridis_c() +
    geom_abline(intercept = model.theta.res[1,1], slope = model.theta.res[2,1], lwd = 1, color = "blue") +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    #geom_smooth(method = "lm") +
    ggtitle("200 bp windows") +
    annotate("text", x = -4, y = 0.3 , label = paste("n = ", nobs.theta, sep =""), cex = 4.5) +
    annotate("text", x = -4, y = 0.28, label = TeX(blabel), cex = 4.5) +
    annotate("text", x = -4, y = 0.26, label = TeX(r2label), cex = 4.5) +     
    ylab(myylab) +
    xlab("log (predicted mutation rate)")

########################################################################################

### 400 bp #############################################################################
load(file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/thetaW_400.RData")
load("./mutac_ms/data/thetaModel_new400.RData") #Load data for 200 bp window

#Calculate predicted mutation rates for each window
#nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC,
#GC content was from 0 to 100 in the model fit
predmutrate <- GCmodelres[1,1] + GCmodelres[2,1]*theta.chr.all[,7] + GCmodelres[3,1]*theta.chr.all[,8] + GCmodelres[4,1]*theta.chr.all[,10] + GCmodelres[5,1]*theta.chr.all[,6]*100 + GCmodelres[6,1]*theta.chr.all[,6]*100*theta.chr.all[,8]

theta.chr.all$predmut <- predmutrate
theta.plot <- theta.chr.all
colnames(theta.plot)[3] <- "thetaW"

myylab <- c(TeX("$\\theta_W$"))
blabel <- paste0("$\\beta = ", round(model.theta.res[2,1],4), "\\, \\[", round(model.theta.res[2,3],4), ", \\, ", round(model.theta.res[2,4],4), "\\]$")
r2label <- paste0("$R^2 = ", round(theta.rsquare.res[1,1],2), "\\, \\[", round(theta.rsquare.res[1,3],2), ", \\, ", round(theta.rsquare.res[1,4],2), "\\]$")
theta.plot400 <- ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    #geom_point(alpha = 0.05)
    geom_hex() +
    scale_fill_viridis_c() +
    geom_abline(intercept = model.theta.res[1,1], slope = model.theta.res[2,1], lwd = 1, color = "blue") +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    #geom_smooth(method = "lm") +
    ggtitle("400 bp windows") +
    annotate("text", x = -4, y = 0.3 , label = paste("n = ", nobs.theta, sep =""), cex = 4.5) +
    annotate("text", x = -4, y = 0.28, label = TeX(blabel), cex = 4.5) +
    annotate("text", x = -4, y = 0.26, label = TeX(r2label), cex = 4.5) +     
    ylab(myylab) +
    xlab("log (predicted mutation rate)")

############################################################################################

### 600 bp ##################################################################################
load(file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/thetaW_600.RData")
load("./mutac_ms/data/thetaModel_new600.RData") #Load data for 200 bp window

#Calculate predicted mutation rates for each window
#nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC,
#GC content was from 0 to 100 in the model fit
predmutrate <- GCmodelres[1,1] + GCmodelres[2,1]*theta.chr.all[,7] + GCmodelres[3,1]*theta.chr.all[,8] + GCmodelres[4,1]*theta.chr.all[,10] + GCmodelres[5,1]*theta.chr.all[,6]*100 + GCmodelres[6,1]*theta.chr.all[,6]*100*theta.chr.all[,8]

theta.chr.all$predmut <- predmutrate
theta.plot <- theta.chr.all
colnames(theta.plot)[3] <- "thetaW"

myylab <- c(TeX("$\\theta_W$"))
blabel <- paste0("$\\beta = ", round(model.theta.res[2,1],4), "\\, \\[", round(model.theta.res[2,3],4), ", \\, ", round(model.theta.res[2,4],4), "\\]$")
r2label <- paste0("$R^2 = ", round(theta.rsquare.res[1,1],2), "\\, \\[", round(theta.rsquare.res[1,3],2), ", \\, ", round(theta.rsquare.res[1,4],2), "\\]$")
theta.plot600 <- ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    #geom_point(alpha = 0.05)
    geom_hex() +
    scale_fill_viridis_c() +
    geom_abline(intercept = model.theta.res[1,1], slope = model.theta.res[2,1], lwd = 1, color = "blue") +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    #geom_smooth(method = "lm") +
    ggtitle("600 bp windows") +
    annotate("text", x = -4, y = 0.3 , label = paste("n = ", nobs.theta, sep =""), cex = 4.5) +
    annotate("text", x = -4, y = 0.28, label = TeX(blabel), cex = 4.5) +
    annotate("text", x = -4, y = 0.26, label = TeX(r2label), cex = 4.5) +     
    ylab(myylab) +
    xlab("log (predicted mutation rate)")

###############################################################################################

### 800 bp ##################################################################################
load(file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/thetaW_800.RData")
load("./mutac_ms/data/thetaModel_new800.RData") #Load data for 800 bp window

#Calculate predicted mutation rates for each window
#nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC,
#GC content was from 0 to 100 in the model fit
predmutrate <- GCmodelres[1,1] + GCmodelres[2,1]*theta.chr.all[,7] + GCmodelres[3,1]*theta.chr.all[,8] + GCmodelres[4,1]*theta.chr.all[,10] + GCmodelres[5,1]*theta.chr.all[,6]*100 + GCmodelres[6,1]*theta.chr.all[,6]*100*theta.chr.all[,8]

theta.chr.all$predmut <- predmutrate
theta.plot <- theta.chr.all
colnames(theta.plot)[3] <- "thetaW"

myylab <- c(TeX("$\\theta_W$"))
blabel <- paste0("$\\beta = ", round(model.theta.res[2,1],4), "\\, \\[", round(model.theta.res[2,3],4), ", \\, ", round(model.theta.res[2,4],4), "\\]$")
r2label <- paste0("$R^2 = ", round(theta.rsquare.res[1,1],2), "\\, \\[", round(theta.rsquare.res[1,3],2), ", \\, ", round(theta.rsquare.res[1,4],2), "\\]$")
theta.plot800 <- ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    #geom_point(alpha = 0.05)
    geom_hex() +
    scale_fill_viridis_c() +
    geom_abline(intercept = model.theta.res[1,1], slope = model.theta.res[2,1], lwd = 1, color = "blue") +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    #geom_smooth(method = "lm") +
    ggtitle("800 bp windows") +
    annotate("text", x = -4, y = 0.3 , label = paste("n = ", nobs.theta, sep =""), cex = 4.5) +
    annotate("text", x = -4, y = 0.28, label = TeX(blabel), cex = 4.5) +
    annotate("text", x = -4, y = 0.26, label = TeX(r2label), cex = 4.5) +     
    ylab(myylab) +
    xlab("log (predicted mutation rate)")

###############################################################################################

### 1000 bp ##################################################################################
load(file = "~/Documents/tutkijatohtori/epimutation/MA_WGS/thetaW_1000.RData")
load("./mutac_ms/data/thetaModel_new1000.RData") #Load data for 1000 bp window

#Calculate predicted mutation rates for each window
#nmut ~ offset(log(count)) + Centromere + H3K9 + H3K27 + GC + H3K9:GC,
#GC content was from 0 to 100 in the model fit
predmutrate <- GCmodelres[1,1] + GCmodelres[2,1]*theta.chr.all[,7] + GCmodelres[3,1]*theta.chr.all[,8] + GCmodelres[4,1]*theta.chr.all[,10] + GCmodelres[5,1]*theta.chr.all[,6]*100 + GCmodelres[6,1]*theta.chr.all[,6]*100*theta.chr.all[,8]

theta.chr.all$predmut <- predmutrate
theta.plot <- theta.chr.all
colnames(theta.plot)[3] <- "thetaW"

myylab <- c(TeX("$\\theta_W$"))
blabel <- paste0("$\\beta = ", round(model.theta.res[2,1],4), "\\, \\[", round(model.theta.res[2,3],4), ", \\, ", round(model.theta.res[2,4],4), "\\]$")
r2label <- paste0("$R^2 = ", round(theta.rsquare.res[1,1],2), "\\, \\[", round(theta.rsquare.res[1,3],2), ", \\, ", round(theta.rsquare.res[1,4],2), "\\]$")
theta.plot1000 <- ggplot(theta.plot, aes(x = predmut, y = thetaW)) +
    #geom_point(alpha = 0.05)
    geom_hex() +
    scale_fill_viridis_c() +
    geom_abline(intercept = model.theta.res[1,1], slope = model.theta.res[2,1], lwd = 1, color = "blue") +
    #geom_ribbon(data = theta.pred, aes(x = predmut, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    #geom_smooth(method = "lm") +
    ggtitle("1000 bp windows") +
    annotate("text", x = -4, y = 0.3 , label = paste("n = ", nobs.theta, sep =""), cex = 4.5) +
    annotate("text", x = -4, y = 0.28, label = TeX(blabel), cex = 4.5) +
    annotate("text", x = -4, y = 0.26, label = TeX(r2label), cex = 4.5) +     
    ylab(myylab) +
    xlab("log (predicted mutation rate)")

###############################################################################################

finalthetaplot <- plot_grid(theta.plot200, theta.plot400, theta.plot600, theta.plot800, theta.plot1000)
save_plot(filename = "./mutac_ms/fig/thetares_windows.pdf", finalthetaplot, base_height=3.71*2.5, base_width=3.71*1.618*3)


save_plot(filename = "./mutac_ms/fig/thetares_new.pdf", final.theta, base_height=3.71, base_width=3.71*1.8*1.618)
          
### * Reanalysis of data from Wang et al. 2020

### ** Load the data

##Load the mutation data from Wang et al.

wangdata <- read.table("wang_mutations.csv", header = T, sep = ",")
#Mutations only in the seven chromosomes...
#Change the order of columns, so that chromosome is the 1st and position the 2nd column
wangdata <- wangdata[,c(2:7,1)]
#Get anc.base and sample.base
wangdata$anc.base <- substr(as.character(wangdata$Change),1,1)
wangdata$sample.base <- substr(as.character(wangdata$Change),3,3)
##Get the cross where mutations happened
mylist <- strsplit(as.character(wangdata$Ascospores),"-")
wangdata$tetrad <- unname(sapply(mylist, '[[', 1))
wangdata$tetrad <- factor(wangdata$tetrad) #This is the tetrad (represents one meiosis)
#There were 5 different crosses and a number of replicates within them
#In total 67 meioses were sampled

mutspec.temp <- rep(0, nrow(wangdata))

#Loop over all mutations and assign them to correct category
for(i in 1:nrow(wangdata)) {
    anc.base <- as.character(wangdata$anc.base[i])
    sample.base <- as.character(wangdata$sample.base[i])
    if((anc.base == "A" & sample.base == "T") | (anc.base == "T" & sample.base == "A")) {mutspec.temp[i] <- "A:T → T:A"}
    if((anc.base == "C" & sample.base == "G") | (anc.base == "G" & sample.base == "C")) {mutspec.temp[i] <- "C:G → G:C"}
    if((anc.base == "C" & sample.base == "A") | (anc.base == "G" & sample.base == "T")) {mutspec.temp[i] <- "C:G → A:T"}
    if((anc.base == "A" & sample.base == "C") | (anc.base == "T" & sample.base == "G")) {mutspec.temp[i] <- "A:T → C:G"}
    if((anc.base == "A" & sample.base == "G") | (anc.base == "T" & sample.base == "C")) {mutspec.temp[i] <- "A:T → G:C"}
    if((anc.base == "C" & sample.base == "T") | (anc.base == "G" & sample.base == "A")) {mutspec.temp[i] <- "C:G → T:A"}
}
                                                                                         
wangdata$mutation <- mutspec.temp #Store mutations

#Then load the Chip-seq domain data and amount of bases
### *** Load data about genomic features

### Duplicated regions

##Data about duplicated regions from Wang et al. 2021
dupregions <- read.csv("duplicated.csv", header = T)

#Consolidate duplicated regions, select only those detected by Blast   
dupreg <- consolidate.dupregions(dupregions, method = "Blast")        
sum(dupreg$End - (dupreg$Start - 1)) #Calculate length of all duplicated regions: 6 585 945 bp

dupcheck <- in.dup(wangdata, dupreg) #Check which mutations happened in duplicated regions
##Little over half of the mutations occur in duplicated regions (seems quite a lot)
## What is the number of called bases in duplicated regions?
wangdata$duplicatedreg <- dupcheck

### Centromeric regions

cent <- read.csv("centromeres.csv", header = T)

##For each mutation look whether it occurred in a centromeric region
wangdata$centromere <- in.centromeric(wangdata, cent)

sum(wangdata$centromere) #1073 mutations happened in centromeric regions. This seems quite a lot


### Regions in H3K27me domains (duplicates removed from chip-seq)
h3k27 <- read.table("2489.H3K27.domains.duprm.bed", header = FALSE, sep = "\t")
colnames(h3k27) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k27[,2] <- h3k27[,2] + 1 #BED files have zero based start coordinate, need to fix this
h3k27 <- filter(h3k27, Chromosome != "mtDNA") #There are reads from mtDNA, but it should not have histones?, but drop it from the dataset
sum(h3k27$End - (h3k27$Start - 1)) #Lenght of all H3K27 domains: 4 538 600 bp

wangdata$H3K27 <- in.dup(wangdata, h3k27) #in.dup function should work for this
sum(wangdata$H3K27) # 735 mutations in H3K27

#H3K27 regions with overlapping H3K9 regions excluded
h3k27exk9 <- read.table("2489.H3K27_exK9.bed", header = FALSE, sep = "\t")
colnames(h3k27exk9) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k27exk9[,2] <- h3k27exk9[,2] + 1 #BED files have zero based start coordinate, need to fix this
sum(h3k27exk9$End - (h3k27exk9$Start - 1)) #Length of H3K27 domains with no K9 3 839 000 bp

wangdata$H3K27exK9 <- in.dup(wangdata, h3k27exk9)
sum(wangdata$H3K27exK9) # 349 mutations in H3K27 without K9


### Regions in H3K9me domtains (duplicates removed from chip-seq)
h3k9 <- read.table("2489.H3K9.domains.duprm.bed", header = FALSE, sep = "\t")
colnames(h3k9) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k9[,2] <- h3k9[,2] + 1 #BED files have zero based start coordinate, need to fix this
h3k9 <- filter(h3k9, Chromosome != "mtDNA") #There are reads from mtDNA, but it should not have histones?, but drop it from the dataset
sum(h3k9$End - (h3k9$Start - 1)) #Length of all H3K9 domains: 7 513 000 bp

wangdata$H3K9 <- in.dup(wangdata, h3k9) #in.dup function should work for this
sum(wangdata$H3K9) # 8691 mutations in H3K9 (this is almost all of the mutations, seems a lot)


### Regions in H3K36me domains (duplicates removed from chip-seq)
h3k36 <- read.table("2489.H3K36.domains.duprm.bed", header = FALSE, sep = "\t")
colnames(h3k36) <- c("Chromosome", "Start", "End", "change", "avg.count", "domain.score", "strand")
h3k36[,2] <- h3k36[,2] + 1 #BED files have zero based start coordinate, need to fix this
h3k36 <- filter(h3k36, Chromosome != "mtDNA") #There are reads from mtDNA, but it should not have histones?, but drop it from the dataset
sum(h3k36$End - (h3k36$Start - 1)) #Length of all H3K36 domains: 32 307 600 bp

wangdata$H3K36 <- in.dup(wangdata, h3k36)
sum(wangdata$H3K36) #1775 mutations in H3K36


### Euchromatic regions
euchromatin <- read.table("2489.euchromatin.bed", header = FALSE, sep = "\t")
colnames(euchromatin) <-  c("Chromosome", "Start", "End")
euchromatin[,2] <- euchromatin[,2] + 1 #BED files have zero based start coordinate, need to fix this
euchromatin <- filter(euchromatin, Chromosome != "mtDNA") #filter out mtDNA
sum(euchromatin$End - (euchromatin$Start - 1)) #Length 29 692 084 bp
euchromatin <- euchromatin[1:260,] #Only the seven chromosomes

wangdata$EUCHR <- in.dup(wangdata, euchromatin)
sum(wangdata$EUCHR) #1453 mutations in euchromatin

#Mutations in euchromatin out of total mutations is:
sum(wangdata$EUCHR)/nrow(wangdata)*100  #13.8% of mutations

### *** Called bases and base frequencies

load("~/Documents/tutkijatohtori/epimutation/MA_WGS/mutac_ms/data/calledsites.RData")

As <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "A"])
Ts <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "T"])
Cs <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "C"])
Gs <- sum(unlist(called.H3K9.matA)[names(unlist(called.H3K9.matA)) == "G"])

As.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "A"])
Ts.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "T"])
Cs.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "C"])
Gs.a <- sum(unlist(called.H3K9.mata)[names(unlist(called.H3K9.mata)) == "G"])

total.h3k9 <- c(As + As.a + Ts + Ts.a + Cs + Cs.a + Gs + Gs.a)
basefreq.h3k9 <- c(As + As.a, Ts + Ts.a, Cs + Cs.a, Gs + Gs.a)/total.h3k9

#Base frequencies in euchromatin
As.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "A"])
Ts.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "T"])
Cs.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "C"])
Gs.eu <- sum(unlist(called.euchr.matA)[names(unlist(called.euchr.matA)) == "G"])

As.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "A"])
Ts.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "T"])
Cs.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "C"])
Gs.eu.a <- sum(unlist(called.euchr.mata)[names(unlist(called.euchr.mata)) == "G"])

total.eu <- As.eu + As.eu.a + Ts.eu + Ts.eu.a + Cs.eu + Cs.eu.a + Gs.eu + Gs.eu.a
basefreq.eu <- c(As.eu + As.eu.a,  Ts.eu + Ts.eu.a,  Cs.eu + Cs.eu.a,  Gs.eu + Gs.eu.a)/total.eu

#Base freqs in total
As.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "A"])
Ts.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "T"])
Cs.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "C"])
Gs.t <- sum(unlist(called.chr.matA)[names(unlist(called.chr.matA)) == "G"])

As.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "A"])
Ts.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "T"])
Cs.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "C"])
Gs.t2 <- sum(unlist(called.chr.mata)[names(unlist(called.chr.mata)) == "G"])

total.chrs <- As.t + As.t2 + Ts.t + Ts.t2 + Cs.t + Cs.t2 + Gs.t + Gs.t2
basefreqs <- c(As.t + As.t2, Ts.t + Ts.t2, Cs.t + Cs.t2, Gs.t + Gs.t2) / total.chrs

bp.eu <- total.eu/2
bp.h3k9 <- total.h3k9/2
### ** Calculate mutation rates and spectra in different domains

mutpercross.domains <- group_by(wangdata, tetrad, H3K9, EUCHR, .drop = F)
mutpercross.domains <- summarise(mutpercross.domains, nmut = n())
mutpercross.domains <- filter(mutpercross.domains, H3K9 == 1 | EUCHR == 1)
mutpercross.domains$Domain <- factor(mutpercross.domains$H3K9, levels = c(0,1), labels = c("Euchromatin", "H3K9"))
mutpercross.domains$cross <- factor(substr(as.character(mutpercross.domains$tetrad),1,1)) #cross

mutpercross.euchromatin <- filter(mutpercross.domains, Domain == "Euchromatin")
mutpercross.h3k9 <- filter(mutpercross.domains, Domain == "H3K9")

median(mutpercross.euchromatin$nmut) #22
median(mutpercross.h3k9$nmut) #38

p1 <- ggplot(mutpercross.domains, aes(x = nmut)) +
    geom_histogram(fill = "grey", colour = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    ylab("Number of tetrads") +
    xlab("Number of mutations") +    
    facet_wrap( ~ Domain, scales = "free")
##There is some heterogeinity in the data that needs to be modelled

p2 <- ggplot(mutpercross.domains, aes(y = nmut, x = cross)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ylab("Number of mutations per tetrad") +
    xlab("Cross") +    
    facet_wrap( ~Domain, scales = "free")

#Estimate mutation rate per meiosis for euchromatin
model.meiosisrate.eu <- brm(data = mutpercross.euchromatin, family = poisson,
                     nmut ~ 1 + (1 | cross),
                     prior = c(prior(normal(0, 10), class = Intercept)),
                     iter = 6000, warmup = 2000, chains = 4, cores = 4,
                     control = list(adapt_delta = 0.9) )

pp_check(model.meiosisrate.eu) #Okay, this is much better fit to the data

post.meiosis.eu <- posterior_samples(model.meiosisrate.eu)
rate.meiosis.eu <- exp(post.meiosis.eu[,1]) / (bp.eu*67) #Calculate mutation rate
#Note that the mutation rate here is the mean mutation rate of all the crosses
res.meiosis.eu <- c(median(rate.meiosis.eu), HPDinterval(as.mcmc(rate.meiosis.eu)))

#Estimate mutation rate per meiosis for H3K9 domains
model.meiosisrate.h3k9 <- brm(data = mutpercross.h3k9, family = poisson,
                     nmut ~ 1 + (1 | cross),
                     prior = c(prior(normal(0, 10), class = Intercept)),
                     iter = 6000, warmup = 2000, chains = 4, cores = 4,
                     control = list(adapt_delta = 0.9) )

pp_check(model.meiosisrate.h3k9) #Still bad fit, even if little better...

##Try gamma-poisson model...
model.meiosisrate.h3k9.gamma <- brm(data = mutpercross.h3k9, family = negbinomial,
                                    nmut ~ 1 + (1 | cross),
                                    prior =c(prior(normal(0, 10), class = Intercept)),
                                    iter = 6000, warmup = 2000, chains = 4, cores = 4,
                                    control = list(adapt_delta = 0.9) )

pp_check(model.meiosisrate.h3k9.gamma)

post.meiosis.h3k9 <- posterior_samples(model.meiosisrate.h3k9.gamma)
rate.meiosis.h3k9 <- exp(post.meiosis.h3k9[,1]) / (bp.h3k9*67) #Calculate mutation rate
#Note that the mutation rate here is the mean mutation rate of all the crosses
res.meiosis.h3k9 <- c(median(rate.meiosis.h3k9), HPDinterval(as.mcmc(rate.meiosis.h3k9)))

#Ratio of h3k9 to euchromatin in meiosis
ratio <- rate.meiosis.h3k9 / rate.meiosis.eu

res.ratio.meiosis <- c(median(ratio), HPDinterval(as.mcmc(ratio)))

#Save the results
save(res.meiosis.eu, res.meiosis.h3k9, res.ratio.meiosis, file = "./mutac_ms/data/rate_meiosis.RData")

### Calculate spectra of mutation in different domains
#Mutations as aggregated
#For the whole genome
wangdata.agg <- group_by(wangdata, mutation, .drop = F)
wangdata.agg <- summarise(wangdata.agg, nmut = n())

#For H3K9 only
wangdata.agg.h3k9 <- group_by(filter(wangdata, H3K9 == 1), mutation, .drop = F)
wangdata.agg.h3k9 <- summarise(wangdata.agg.h3k9, nmut = n())

#Euchromatin
wangdata.agg.eu <- group_by(filter(wangdata, EUCHR == 1), mutation, .drop = F)
wangdata.agg.eu <- summarise(wangdata.agg.eu, nmut = n())

#Calculating expected frequencies
wangdata.agg$expected.freq <- c(rep((basefreqs[1]+basefreqs[2])/3,3), rep((basefreqs[3]+basefreqs[4])/3,3))
wangdata.agg.h3k9$expected.freq <- c(rep((basefreq.h3k9[1]+basefreq.h3k9[2])/3,3), rep((basefreq.h3k9[3]+basefreq.h3k9[4])/3,3))
wangdata.agg.eu$expected.freq <- c(rep((basefreq.eu[1]+basefreq.eu[2])/3,3), rep((basefreq.eu[3]+basefreq.eu[4])/3,3))

#expected numbers
wangdata.agg$expected <- wangdata.agg$expected.freq*sum(wangdata.agg$nmut)
wangdata.agg.h3k9$expected <- wangdata.agg.h3k9$expected.freq*sum(wangdata.agg.h3k9$nmut)
wangdata.agg.eu$expected <- wangdata.agg.eu$expected.freq*sum(wangdata.agg.eu$nmut)

#Fit models
model.wang.all <- brm(data = wangdata.agg, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.wang.all)

#Convert posterior samples back to normal scale and calculate intervals
wangdata.spectra.agg <- data.frame(estimate = apply(exp(post[,-7]), 2, median), HPDinterval(as.mcmc(exp(post[,-7]))))
wangdata.spectra.agg$type <- c("transversion", "transition", "transversion", "transversion", "transversion", "transition")


#H3K9
model.wang.h3k9 <- brm(data = wangdata.agg.h3k9, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.wang.h3k9)
#Convert posterior samples back to normal scale and calculate intervals
wangdata.spectra.agg.h3k9 <- data.frame(estimate = apply(exp(post[,-7]), 2, median), HPDinterval(as.mcmc(exp(post[,-7]))))
wangdata.spectra.agg.h3k9$type <- c("transversion", "transition", "transversion", "transversion", "transversion", "transition")
post.h3k9.wang <- exp(post[,-7])

#Euchromatin
model.wang.eu <- brm(data = wangdata.agg.eu, family = poisson,
                     nmut ~ -1 + offset(log(expected)) + mutation,
                     prior = c(prior(normal(0, 10), class = b)),
                     iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(model.wang.eu)
#Convert posterior samples back to normal scale and calculate intervals
wangdata.spectra.agg.eu <- data.frame(estimate = apply(exp(post[,-7]), 2, median), HPDinterval(as.mcmc(exp(post[,-7]))))
wangdata.spectra.agg.eu$type <- c("transversion", "transition", "transversion", "transversion", "transversion", "transition")
post.eu.wang <- exp(post[,-7])
#save(post.h3k9.wang, post.eu.wang, file = "./mutac_ms/data/spectra_sex_posterior.RData")

#Making a dataframe that combines all of the data
wangdata.spectra.all <- rbind(wangdata.spectra.agg, wangdata.spectra.agg.h3k9, wangdata.spectra.agg.eu)
wangdata.spectra.all$domain <- c(rep("Whole genome", 6), rep("H3K9", 6), rep("Euchromatin",6))

wangdata.spectra.all$mutation <- factor(rep(c("A:T → C:G", "A:T → G:C", "A:T → T:A", "C:G → A:T", "C:G → G:C", "C:G → T:A"),3), levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"))

mylabels <- c(TeX("A:T $\\rightarrow$ T:A"), TeX("A:T $\\rightarrow$ C:G"), TeX("C:G $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ A:T"), TeX("A:T $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ T:A") )


#Plot of the mutation spectra in different regions of the genome
 p3 <- ggplot(wangdata.spectra.all, aes(y = estimate, ymin = lower, ymax = upper, x = mutation, fill = type)) +
    geom_bar(stat = "identity", colour = "black") +
    geom_errorbar(width = 0.1) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    #scale_y_continuous(expand = c(0,0), limits = c(0,2.5), breaks = c(seq(0,2.25, by = 0.25))) +
    scale_x_discrete(labels = mylabels) +
    ylab("Relative mutation rate") +
    xlab("Mutation") +
    geom_hline(yintercept = 1, lty = "dashed") +
    facet_wrap( ~ domain, scales="free_y") +
    theme(legend.position = "top", legend.justification = 0.5, legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

meiosis.plot <- plot_grid(p1, p2, p3, nrow = 3, labels = c("A", "B", "C"), rel_heights = c(1/5, 1/5, 1/4))

grDevices::cairo_pdf("./mutac_ms/fig/meiosis.pdf", width = 10, height = 12)
plot_grid(p1, p2, p3, nrow = 3, labels = c("A", "B", "C"), rel_heights = c(1/5, 1/5, 1/4))
dev.off()

### ** Mutations in genes and genetic load

mutpercross.domains.gen <- group_by(wangdata, H3K9, EUCHR, genic_regions, .drop = F)
mutpercross.domains.gen <- summarise(mutpercross.domains.gen, nmut = n())
mutpercross.domains.gen <- filter(mutpercross.domains.gen, H3K9 == 1 | EUCHR == 1)

euchrom.gen <- filter(mutpercross.domains.gen, EUCHR == 1)
h3k9.gen <- filter(mutpercross.domains.gen, H3K9 == 1)

#Let's compare the amount of non-synonymous mutations
eu.nsyn <- 308+2+9+3
eu.tot <- sum(euchrom.gen$nmut)
h3k9.nsyn <- 8+1+6
h3k9.tot <- sum(h3k9.gen$nmut)

load.results <- c((eu.nsyn/eu.tot)*100, (h3k9.nsyn/h3k9.tot)*100)

save(load.results, file = "./mutac_ms/data/load.res.RData")
### ** Comparing spectra of asexual and sexual mutations

#load("./mutac_ms/data/spectra_asex_posterior.RData")
#load("./mutac_ms/data/spectra_sex_posterior.RData")

#Calculate ratio of sexual mutation spectra to asexual mutation spectra for euchromatin and H3K9

post.rat.eu <- post.eu.wang / post.eu
ratio.spectra.eu <- data.frame(estimate = apply(post.rat.eu, 2, median), HPDinterval(as.mcmc(post.rat.eu)))

post.rat.h3k9 <- post.h3k9.wang / post.h3k9
ratio.spectra.h3k9 <- data.frame(estimate = apply(post.rat.h3k9, 2, median), HPDinterval(as.mcmc(post.rat.h3k9)))

ratio.spectra.all <- rbind(ratio.spectra.eu, ratio.spectra.h3k9)
ratio.spectra.all$mutation <- factor(rep(c("A:T → C:G", "A:T → G:C", "A:T → T:A", "C:G → A:T", "C:G → G:C", "C:G → T:A"),2), levels = c("A:T → T:A", "A:T → C:G", "C:G → G:C", "C:G → A:T", "A:T → G:C", "C:G → T:A"))
ratio.spectra.all$type <- rep(c("transversion", "transition", "transversion", "transversion", "transversion", "transition"),2)
ratio.spectra.all$domain <- c(rep("Euchromatin",6), rep("H3K9", 6))

mylabels <- c(TeX("A:T $\\rightarrow$ T:A"), TeX("A:T $\\rightarrow$ C:G"), TeX("C:G $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ A:T"), TeX("A:T $\\rightarrow$ G:C"), TeX("C:G $\\rightarrow$ T:A") )

sex.asex.ratio <- ggplot(ratio.spectra.all, aes(y = estimate, ymin = lower, ymax = upper, x = mutation, color = type)) +
    geom_pointrange() +
    scale_x_discrete(labels = mylabels) +
    scale_y_continuous(limits = c(0,6), breaks = seq(0,6,1)) +
    ylab("Ratio of meiosis / mitosis relative rates") +
    xlab("Mutation") +
    geom_hline(yintercept = 1, lty = "dashed") +
    facet_wrap( ~ domain) +
    theme(legend.position = "top", legend.justification = 0.5, legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

grDevices::cairo_pdf("./mutac_ms/fig/reproduction_ratio.pdf", width = 8, height = 5)
sex.asex.ratio
dev.off()
### * Plot genome coverage stats


#aineisto <- read.table(paste(getwd(),"/sequencing/2489anc.gc.report", sep = ""), sep = "\t", dec = ",", skip = 5, header = T)

aineisto <- read.table("ANC1A.gc.report", sep = "\t", dec = ",", skip = 5, header = T)
aineisto <- mutate(aineisto, GCcont = WINDOWS/sum(WINDOWS)*100)

#To get two different y-axes, we plot the data on the same scale (and alter the other dataset by multiplying by 2.5
#pdf("GCbias.pdf")
gcplot <- ggplot(aineisto) +
    geom_bar(mapping = aes(x = GC, y = GCcont), stat = "identity", fill = "grey") +
    geom_pointrange(mapping = aes(x = GC, y = NORMALIZED_COVERAGE*2.5, ymin = (NORMALIZED_COVERAGE*2.5-ERROR_BAR_WIDTH*2.5), ymax = (NORMALIZED_COVERAGE*2.5+ERROR_BAR_WIDTH*2.5) ), size = 0.5, colour = "blue") +
    scale_y_continuous(name = "% of genome", limits = c(0,5), expand = c(0,0), sec.axis = sec_axis(~./2.5, name = "Normalized coverage")) +
    scale_x_continuous(name = "%GC content in 100bp windows", limits = c(0, 90), expand = c(0,0)) +
    theme(axis.title.y = element_text(color = "grey"), axis.title.y.right = element_text(color = "blue"))
#dev.off()

save_plot(filename = "./mutac_ms/fig/gcbiasplot.pdf", gcplot)
