### This is R scripts for estimating the number of mitoses that happened during experiment

#Load libraries

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(brms)

### * Load data

### ** Data about test tube diameter

ttdiam <- read.table("tt_diameters.csv", header = T)

#Diameter of test tube
mean(ttdiam$diameter) #9.5 mm

### ** Data about conidial numbers

conidia <- read.csv("colonies_on_plates2.csv", header = T)

conidia$genotype <- paste(conidia$Line, conidia$mat, conidia$Generation, sep = " ")

### *** Analysis of conidial numbers

ggplot(conidia, aes(x = conidia)) +
    geom_histogram(colour = "black", fill = "white") +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(genotype ~ operator)


ggplot(conidia, aes(x = conidia/10^6)) +
    geom_histogram(colour = "black", fill = "white") +
    scale_y_continuous(expand = c(0,0))

conidia2 <- conidia
conidia2$conidia <- conidia2$conidia/10^6

conidiamodel2 <- brm(conidia ~ 1 + (1|genotype), data = conidia2, family = gaussian,
                                  #prior = c(prior(normal(0,240^6), class = Intercept)),
                                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

#post <- posterior_samples(conidiamodel)
post2 <- posterior_samples(conidiamodel2)


#conidia.post <- exp(post$b_Intercept)*10^6
conidia.post <- post2$b_Intercept*10^6 #Data was earlier divided by 10^6, so it needs to be multiplied by 10^6 because of convergence is easier with the smaller numbers

### ** Data about colony size on sorbose plates

colonysize <- read.csv("sorbose_colonies.csv", header = T, dec = ",", sep = ";")

### *** Analysis of colony sizes on sorbose plates

ggplot(colonysize, aes(x = Area)) +
    geom_histogram(colour = "black", fill = "white") +
    scale_y_continuous(expand = c(0,0), limits = c(0,50)) +
    facet_grid(Line ~ mat)

##Looks like size has a distribution that is multimodal
##Need to take line into account

#Make a new factor that combines, line, mat and generation
colonysize$genotype <- paste(colonysize$Line, colonysize$mat, colonysize$generation, sep = " ")

ggplot(colonysize, aes(x = Area)) +
    geom_histogram(colour = "black", fill = "white") +
    scale_y_continuous(expand = c(0,0), limits = c(0,50)) +
    facet_wrap( ~ genotype)

##Model to estimate size, while taking genotype into account

colsizemodel <- brm(Area ~ 1 + (1|genotype), data = colonysize, family = gaussian,
                                  prior = c(prior(normal(0, 3), class = Intercept)),
                                  iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(colsizemodel)
colsize.post <- post$b_Intercept #Posterior distribution of colony size on sorbose plates

### ** Data about nuclei counts

nuclei <- read.csv("nuclei_counts.csv", header = T)

### *** Analysis of nuclei counts

#VM is at least bimodel distribution of hyphal width
ggplot(nuclei, aes(x = hyphalwidth)) +
    geom_histogram(colour = "black", fill = "white") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 85)) +
    facet_wrap( ~ media)

#Count the number of nuclei per area
nuclei$nucleiperarea <- nuclei$nucleicount / nuclei$area

##Need to check if there is a string effect of hyphal width and nuclei per area
ggplot(nuclei, aes(x = hyphalwidth, y = nucleiperarea)) +
    geom_point() +
    facet_wrap( ~ media)

#relationship between hyphal width and nuclei per area
summary(lm(data = nuclei, nucleiperarea ~ hyphalwidth + media + hyphalwidth:media))
#yes there is relationship between hyphal width and area in both
#negative in VM, and positive in sorbose

ggplot(nuclei, aes(x = nucleiperarea)) +
    geom_histogram(colour = "black", fill = "white") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 65)) +    
    facet_wrap( ~ media)

#Nuclei per area
mean(filter(nuclei, media == "sorbose")$nucleiperarea) #0.0518 nuclei per square micro meter

mean(filter(nuclei, media == "VM")$nucleiperarea) #0.0648 nuclei per square micro meter

nucleiareamodel <- brm(bf(nucleiperarea ~ 1 + media, sigma ~ media), data = nuclei, family = gaussian,
                       prior = c(prior(normal(0, 0.1), class = Intercept)),
                       iter = 3000, warmup = 1000, chains = 4, cores = 4)

post <- posterior_samples(nucleiareamodel)
nuclei.area.VM <- post$b_Intercept #Posterior distribution of nuclei per area in VM
nuclei.area.sorbose <- post$b_Intercept + post$b_mediasorbose #Post. of nuclei per area in sorbose

### Calculating the number of nuclei

##Number of nuclei in sorbose colonies
nuclei.sorbose.col <- colsize.post*(10^6)*nuclei.area.sorbose

##Number of nuclei in slant

#Slant area
slant.area <- pi*(mean(ttdiam$diameter)/2)^2

#Nuclei in VM slant
nuclei.VM.slant <- slant.area*(10^6)*nuclei.area.VM

#Conidia in VM slant (2 nuclei per conidia)
nuclei.conidia <- conidia.post*2

###Number of mitotic divisions

mitoses <- log2(nuclei.sorbose.col/2) + log2((nuclei.VM.slant+nuclei.sorbose.col)/nuclei.sorbose.col) + log2((nuclei.conidia+nuclei.VM.slant+nuclei.sorbose.col)/(nuclei.VM.slant+nuclei.sorbose.col))

### *** Save the data

save(file = "mitoses.RData", mitoses)

