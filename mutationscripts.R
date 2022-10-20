##Some functions that are required for analysing mutation data

#Classifies point mutations
classify.mutations <- function(anc.base, new.base) {
    anc.base <- as.character(anc.base)
    new.base <- as.character(new.base)
    nmutations <- length(anc.base)
    res.mat <- matrix(rep(0, nmutations*2), ncol = 2)
    colnames(res.mat) <- c("mutation", "mut.type")
    for(i in 1:nmutations) {
        if( (anc.base[i] == "A" & new.base[i] == "T") | (anc.base[i] == "T" & new.base[i] == "A") ) { res.mat[i,] <- c("A:T -> T:A", "transversion") }
        if( (anc.base[i] == "A" & new.base[i] == "C") | (anc.base[i] == "T" & new.base[i] == "G") ) { res.mat[i,] <- c("A:T -> C:G", "transversion") }
        if( (anc.base[i] == "C" & new.base[i] == "G") | (anc.base[i] == "G" & new.base[i] == "C") ) { res.mat[i,] <- c("C:G -> G:C", "transversion") }
        if( (anc.base[i] == "C" & new.base[i] == "A") | (anc.base[i] == "G" & new.base[i] == "T") ) { res.mat[i,] <- c("C:G -> A:T", "transversion") }
        if( (anc.base[i] == "A" & new.base[i] == "G") | (anc.base[i] == "T" & new.base[i] == "C") ) { res.mat[i,] <- c("A:T -> G:C", "transition") }
        if( (anc.base[i] == "C" & new.base[i] == "T") | (anc.base[i] == "G" & new.base[i] == "A") ) { res.mat[i,] <- c("C:G -> T:A", "transition") }
    }
    return(res.mat)
}

#Classifies triplets
classify.trinucleotides <- function(triplet) {
    classes <- c("AAA:TTT", "AAC:GTT", "AAG:CTT", "AAT:ATT", "CAA:TTG", "GAA:TTC", "TAA:TTA", "CAC:GTG", "CAG:CTG", "CAT:ATG", "GAG:CTC", "GAC:GTC", "GAT:ATC", "TAC:GTA", "TAG:CTA", "TAT:ATA", "CCC:GGG", "CCA:TGG", "CCG:CGG", "CCT:AGG", "ACC:GGT", "GCC:GGC", "TCC:GGA", "ACA:TGT", "ACG:CGT", "ACT:AGT", "GCG:CGC", "GCA:TGC", "GCT:AGC", "TCT:AGA", "TCA:TGA", "TCG:CGA")
    ntrip <- length(triplet)
    tripclass <- rep("", ntrip)
    for(i in 1:ntrip) {
        tripclass[i] <- classes[grep(triplet[i], classes)]
    }
    return(tripclass)
}

#This function classifies the 96-mutations
classify.96mutations <- function(aineisto, counts) {
    #Build the table
    classes <- c("AAA:TTT", "AAC:GTT", "AAG:CTT", "AAT:ATT", "CAA:TTG", "GAA:TTC", "TAA:TTA", "CAC:GTG", "CAG:CTG", "CAT:ATG", "GAG:CTC", "GAC:GTC", "GAT:ATC", "TAC:GTA", "TAG:CTA", "TAT:ATA", "CCC:GGG", "CCA:TGG", "CCG:CGG", "CCT:AGG", "ACC:GGT", "GCC:GGC", "TCC:GGA", "ACA:TGT", "ACG:CGT", "ACT:AGT", "GCG:CGC", "GCA:TGC", "GCT:AGC", "TCT:AGA", "TCA:TGA", "TCG:CGA")
    classes <- c(sort(classes[1:16]), sort(classes[17:32])) #Fixing some order
    c1 <- c(rep(classes[1:16],3), rep(classes[17:32],3))
    m1 <- c(rep("A:T -> T:A", 16), rep("A:T -> C:G", 16), rep("A:T -> G:C", 16), rep("C:G -> G:C", 16), rep("C:G -> A:T", 16), rep("C:G -> T:A", 16))
    res.mat <- data.frame(mutation = m1, class = c1, count = rep(0,96), nmut = rep(0, 96) )
    #Loop over the data and fill the table with mutations
    for(i in 1:nrow(aineisto)) {
        index <- which(res.mat$mutation == aineisto$mutation[i] & res.mat$class == aineisto$tripclass[i])
        res.mat$nmut[index] <- res.mat$nmut[index] + 1 #Increment mutation count
    }

    #Add counts for trinucleotide classes
    for(i in 1:nrow(res.mat)) {
        seqs <- unlist(strsplit(as.character(res.mat$class[i]), ":")) #Get both sequences
        index1 <- which(seqs[1] == counts[,1]) #Get index of sequence in counts
        index2 <- which(seqs[2] == counts[,1])
        res.mat$count[i] <- counts[index1,2] + counts[index2,2]
    }
    return(res.mat)
}

#Classify mutations by type
mutation.types <- function(anc.base, new.base) {
    anc.base <- as.character(anc.base)
    new.base <- as.character(new.base)
    nchar.new <- nchar(new.base)
    nchar.anc <- nchar(anc.base)
    nmutations <- length(anc.base)
    res.mat <- matrix(rep(0, nmutations*2), ncol = 2)
    colnames(res.mat) <- c("mutation", "type")
    res.mat <- as.data.frame(res.mat)
    for(i in 1:length(anc.base)) {
        if(is.na(anc.base[i]) == F & is.na(new.base[i]) == F) {
            if(nchar.anc[i] == 1 & nchar.new[i] == 1) {
                res.mat[i,1] <- "point"
                if( (anc.base[i] == "G" & new.base[i] == "A") | (anc.base[i] == "A" & new.base[i] == "G") | (anc.base[i] == "C" & new.base[i] == "T") | (anc.base[i] == "T" & new.base[i] == "C")) {
                    res.mat[i,2] <- "transition" } else {res.mat[i,2] <- "transversion"} }
            if(nchar.anc[i] > nchar.new[i]) {res.mat[i,1] <- "deletion"; res.mat[i,2] <- NA}
            if(nchar.anc[i] < nchar.new[i]) {res.mat[i,1] <- "insertion"; res.mat[i,2] <- NA}
        } else { res.mat[i,1] <- NA; res.mat[i,2] <- NA }
    }
    #    
    return(res.mat)
}

#Get repeated base for each homopolymer
repeated.base <- function(anc.base, new.base, type) {
    anc.base <- as.character(anc.base)
    new.base <- as.character(new.base)
    nbases <- length(anc.base) #This is the number of bases
    results <- rep(0, nbases) #Initialize results
    for(i in 1:nbases) { #Loop over bases to be checked
        if(type[i] == "deletion") { base <- stringr::str_sub(anc.base[i],-1) }
        if(type[i] == "insertion") { base <- stringr::str_sub(new.base[i],-1) }
        if(base == "A" | base == "T") { results[i] <- "A:T" }
        if(base == "C" | base == "G") { results[i] <- "C:G" }
    }

    return(results)
}

##This function retrieves the adjacent bases based on a coordinate
##Return a triplet, where the middle base is the coordinate base (based on genome ref)
##This function calls a python scripts that uses wormtable to retrieve the coordinates
retrieve.triplet <- function(chromosome, position) {
    #First make the command string
    komento <- paste("python", "~/Genomics/Neurospora/mutacc/retrieve.py", chromosome, position, sep = " ")
    system(komento, intern = TRUE)
}


##This function check that ref base and corresponding base in the fasta reference are identical
##Uses samtools (assumes that index file exists in the same folder as reference)
##"fasta" needs to be a path to the fasta file 
check.ref.fasta <- function(fasta, chromosome, position, refbase) {
    #First make the command string
    komento <- paste("samtools faidx", fasta, paste(chromosome,":",position,"-",position, sep = ""), sep = " ")
    fastabase <- system(komento, intern = TRUE)[2] #retrieve the ref base in fasta reference
    return(fastabase == refbase)
}

##This function retrieves the base context of a mutation based on a given chromosome and coordinate
##Uses samtools as above
### kmer can be 3 bases or 5 bases
retrieve.kmer.fasta <- function(fasta, chromosome, position, kmer) {
    if(kmer == 3) {
        #First make the command string
        komento <- paste("samtools faidx", fasta, paste(chromosome,":",position-1,"-",position+1, sep = ""), sep = " ")
    }
    if(kmer == 5) {
        komento <- paste("samtools faidx", fasta, paste(chromosome,":",position-2,"-",position+2, sep = ""), sep = " ")
    }
    basecontext <- system(komento, intern = TRUE)[2] #retrieve the base context of the mutation
    return(basecontext)
}

##Trinucleotides
#This function combines trinucleotides 
trinucleotides <- function(trinucleotide) {
#make a results table
    res.mat <- data.frame(class = c("AAA:TTT", "AAC:GTT", "AAG:CTT", "AAT:ATT", "CAA:TTG", "GAA:TTC", "TAA:TTA", "CAC:GTG", "CAG:CTG", "CAT:ATG", "GAG:CTC", "GAC:GTC", "GAT:ATC", "TAC:GTA", "TAG:CTA", "TAT:ATA", "CCC:GGG", "CCA:TGG", "CCG:CGG", "CCT:AGG", "ACC:GGT", "GCC:GGC", "TCC:GGA", "ACA:TGT", "ACG:CGT", "ACT:AGT", "GCG:CGC", "GCA:TGC", "GCT:AGC", "TCT:AGA", "TCA:TGA", "TCG:CGA"), nmut = rep(0, 32))
    res.mat[1,2] <- sum(trinucleotide == "AAA" | trinucleotide == "TTT")
    res.mat[2,2] <- sum(trinucleotide == "AAC" | trinucleotide == "GTT") #reverse complement
    res.mat[3,2] <- sum(trinucleotide == "AAG" | trinucleotide == "CTT")
    res.mat[4,2] <- sum(trinucleotide == "AAT" | trinucleotide == "ATT")
    res.mat[5,2] <- sum(trinucleotide == "CAA" | trinucleotide == "TTG")
    res.mat[6,2] <- sum(trinucleotide == "GAA" | trinucleotide == "TTC")
    res.mat[7,2] <- sum(trinucleotide == "TAA" | trinucleotide == "TTA")
    res.mat[8,2] <- sum(trinucleotide == "CAC" | trinucleotide == "GTG")
    res.mat[9,2] <- sum(trinucleotide == "CAG" | trinucleotide == "CTG")
    res.mat[10,2] <- sum(trinucleotide == "CAT" | trinucleotide == "ATG")
    res.mat[11,2] <- sum(trinucleotide == "GAG" | trinucleotide == "CTC")
    res.mat[12,2] <- sum(trinucleotide == "GAC" | trinucleotide == "GTC")
    res.mat[13,2] <- sum(trinucleotide == "GAT" | trinucleotide == "ATC")
    res.mat[14,2] <- sum(trinucleotide == "TAC" | trinucleotide == "GTA")
    res.mat[15,2] <- sum(trinucleotide == "TAG" | trinucleotide == "CTA")
    res.mat[16,2] <- sum(trinucleotide == "TAT" | trinucleotide == "ATA")
    res.mat[17,2] <- sum(trinucleotide == "CCC" | trinucleotide == "GGG")
    res.mat[18,2] <- sum(trinucleotide == "CCA" | trinucleotide == "TGG")
    res.mat[19,2] <- sum(trinucleotide == "CCG" | trinucleotide == "CGG")
    res.mat[20,2] <- sum(trinucleotide == "CCT" | trinucleotide == "AGG")
    res.mat[21,2] <- sum(trinucleotide == "ACC" | trinucleotide == "GGT")
    res.mat[22,2] <- sum(trinucleotide == "GCC" | trinucleotide == "GGC")
    res.mat[23,2] <- sum(trinucleotide == "TCC" | trinucleotide == "GGA")
    res.mat[24,2] <- sum(trinucleotide == "ACA" | trinucleotide == "TGT")
    res.mat[25,2] <- sum(trinucleotide == "ACG" | trinucleotide == "CGT")
    res.mat[26,2] <- sum(trinucleotide == "ACT" | trinucleotide == "AGT")
    res.mat[27,2] <- sum(trinucleotide == "GCG" | trinucleotide == "CGC")
    res.mat[28,2] <- sum(trinucleotide == "GCA" | trinucleotide == "TGC")
    res.mat[29,2] <- sum(trinucleotide == "GCT" | trinucleotide == "AGC")
    res.mat[30,2] <- sum(trinucleotide == "TCT" | trinucleotide == "AGA")
    res.mat[31,2] <- sum(trinucleotide == "TCA" | trinucleotide == "TGA")
    res.mat[32,2] <- sum(trinucleotide == "TCG" | trinucleotide == "CGA")
    #return results table
    return(res.mat)
}

#Function checks if number is even
is.even <- function(x) x %% 2 == 0

#This function sorts the trinucleotides outputted by emboss into the same classes as mutations
sort.trinucleotides <- function(trinucleotides, mutations) {
    n.tri <- nrow(trinucleotides) #Should be 64
    mutations$count <- rep(0, nrow(mutations)) #Initialise results
    for(i in 1:n.tri) { #Loop over all trinucleotides
        index <- grep(as.character(trinucleotides[i,1]), as.character(mutations$class)) #get index
        mutations[index,3] <- mutations[index,3] + trinucleotides[i,2] #Add counts to corresponding calss
    }
    mutations$frequency <- mutations$count/sum(mutations$count) #Calculate frequency
    return(mutations) #return results
}

sort.trinucleotides.freq <- function(trinucleotides, frequencies) {
    n.tri <- nrow(frequencies) #Should be 64
    trinucleotides$obsfreq <- rep(0, nrow(trinucleotides))
    trinucleotides$expfreq <- rep(0, nrow(trinucleotides))
    for(i in 1:n.tri) { #Loop over all trinucleotides
        index <- grep(as.character(frequencies[i,1]), as.character(trinucleotides$class)) #get index
        trinucleotides$obsfreq[index] <- trinucleotides$obsfreq[index] + frequencies[i,3] #Add counts to corresponding calss
        trinucleotides$expfreq[index] <- trinucleotides$expfreq[index] + frequencies[i,4]
    }
    return(trinucleotides)
}
    
    

#Function to put triplet data into correct boxes
#Input is table(my triplets), type (either "AT" or "CG"), basefrequencies (f(A), f(T), f(C), f(G)) 
triplet.sorter <- function(taulukko, type, basefreqs) {
    table.len <- length(taulukko) #Should be 32
    #First check which kind of triplets are being sorted: AT or CG (in the middle)
    if(type == "AT") {
        #Make a results table
        datamat <- data.frame(matrix(rep(0, 10*5), ncol = 5))
        colnames(datamat) <- c("Class", "expected.freq", "observed", "expected", "rel.obs")
        datamat$Class <- factor(c("AAA:TTT", "TAA:ATT,AAT:TTA", "TAT:ATA", "CAA:GTT,AAC:TTG", "CAC:GTG", "CAG:GTC,GAC:CTG", "GAA:CTT,AAG:TTC", "GAG:CTC", "CAT:GTA,TAC:ATG", "GAT:CTA,TAG:ATC"))
        trips <- names(taulukko)
        #Store numbers of mutations for different triplet classes
        datamat[1,3] <- sum(taulukko[trips == "AAA" | trips == "TTT"])
        datamat[2,3] <- sum(taulukko[trips == "TAA" | trips == "ATT" | trips == "AAT" | trips == "TTA"])
        datamat[3,3] <- sum(taulukko[trips == "TAT" | trips == "ATA"])
        datamat[4,3] <- sum(taulukko[trips == "CAA" | trips == "GTT" | trips == "AAC" | trips == "TTG"])
        datamat[5,3] <- sum(taulukko[trips == "CAC" | trips == "GTG"])
        datamat[6,3] <- sum(taulukko[trips == "CAG" | trips == "GTC" | trips == "GAC" | trips == "CTG"])
        datamat[7,3] <- sum(taulukko[trips == "GAA" | trips == "CTT" | trips == "AAG" | trips == "TTC"])
        datamat[8,3] <- sum(taulukko[trips == "GAG" | trips == "CTC"])
        datamat[9,3] <- sum(taulukko[trips == "CAT" | trips == "GTA" | trips == "TAC" | trips == "ATG"])
        datamat[10,3] <- sum(taulukko[trips == "GAT" | trips == "CTA" | trips == "TAG" | trips == "ATC"])
        total <- sum(taulukko) #Total number of mutations
        #Calculate expected and relative observed values etc.
        #Expected frequencies of triplets occurring calculated from observed base frequencies
        datamat[1,2] <- basefreqs[1]^3 + basefreqs[2]^3
        datamat[2,2] <- 2*(basefreqs[2]*basefreqs[1]^2) + 2*(basefreqs[1]*basefreqs[2]^2)
        datamat[3,2] <- basefreqs[1]*basefreqs[2]^2 + basefreqs[2]*basefreqs[1]^2
        datamat[4,2] <- 2*(basefreqs[3]*basefreqs[1]^2) + 2*(basefreqs[4]*basefreqs[2]^2)
        datamat[5,2] <- basefreqs[1]*basefreqs[3]^2 + basefreqs[2]*basefreqs[4]^2
        datamat[6,2] <- 2*(basefreqs[1]*basefreqs[3]*basefreqs[4]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[3])
        datamat[7,2] <- 2*(basefreqs[4]*basefreqs[1]^2) + 2*(basefreqs[3]*basefreqs[2]^2)
        datamat[8,2] <- basefreqs[1]*basefreqs[4]^2 + basefreqs[2]*basefreqs[3]^2
        datamat[9,2] <- 2*(basefreqs[3]*basefreqs[1]*basefreqs[2]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[1])
        datamat[10,2] <- 2*(basefreqs[1]*basefreqs[2]*basefreqs[4]) + 2*(basefreqs[3]*basefreqs[2]*basefreqs[1])
        #
        #Normalize frequencies so that they sum to 1
        datamat[,2] <- datamat[,2]/sum(datamat[,2])
        #Expected number of mutations
        datamat[,4] <- datamat[,2]*total
        #Relative observed / expected
        datamat[,5] <- datamat[,3]/datamat[,4]
    }
    if(type == "CG") {
        #Make a results table
        datamat <- data.frame(matrix(rep(0, 10*5), ncol = 5))
        colnames(datamat) <- c("Class", "expected.freq", "observed", "expected", "rel.obs")
        datamat$Class <- factor(c("CCC:GGG", "GCC:CGG,CCG:GGC", "GCG:CGC", "ACC:TGG,CCA:GGT", "ACA:TGT", "ACG:TGC,GCA:CGT", "TCC:AGG,CCT:GGA", "TCT:AGA", "ACT:TGA,TCA:AGT", "TCG:AGC,GCT:CGA"))
        trips <- names(taulukko)
        #Store numbers of mutations for different triplet classes
        datamat[1,3] <- sum(taulukko[trips == "CCC" | trips == "GGG"])
        datamat[2,3] <- sum(taulukko[trips == "GCC" | trips == "CGG" | trips == "CCG" | trips == "GGC"])
        datamat[3,3] <- sum(taulukko[trips == "GCG" | trips == "CGC"])
        datamat[4,3] <- sum(taulukko[trips == "ACC" | trips == "TGG" | trips == "CCA" | trips == "GGT"])
        datamat[5,3] <- sum(taulukko[trips == "ACA" | trips == "TGT"])
        datamat[6,3] <- sum(taulukko[trips == "ACG" | trips == "TGC" | trips == "GCA" | trips == "CGT"])
        datamat[7,3] <- sum(taulukko[trips == "TCC" | trips == "AGG" | trips == "CCT" | trips == "GGA"])
        datamat[8,3] <- sum(taulukko[trips == "TCT" | trips == "AGA"])
        datamat[9,3] <- sum(taulukko[trips == "ACT" | trips == "TGA" | trips == "TCA" | trips == "AGT"])
        datamat[10,3] <- sum(taulukko[trips == "TCG" | trips == "AGC" | trips == "GCT" | trips == "CGA"])
        total <- sum(taulukko) #Total number of mutations
        #Calculate expected and relative observed values etc.
        #Expected frequencies of triplets occurring calculated from observed base frequencies
        datamat[1,2] <- basefreqs[3]^3 + basefreqs[4]^3
        datamat[2,2] <- 2*(basefreqs[4]*basefreqs[3]^2) + 2*(basefreqs[3]*basefreqs[4]^2)
        datamat[3,2] <- basefreqs[3]*basefreqs[4]^2 + basefreqs[4]*basefreqs[3]^2
        datamat[4,2] <- 2*(basefreqs[1]*basefreqs[3]^2) + 2*(basefreqs[2]*basefreqs[4]^2)
        datamat[5,2] <- basefreqs[3]*basefreqs[1]^2 + basefreqs[4]*basefreqs[2]^2
        datamat[6,2] <- 2*(basefreqs[1]*basefreqs[3]*basefreqs[4]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[3])
        datamat[7,2] <- 2*(basefreqs[2]*basefreqs[3]^2) + 2*(basefreqs[1]*basefreqs[4]^2)
        datamat[8,2] <- basefreqs[3]*basefreqs[2]^2 + basefreqs[4]*basefreqs[1]^2
        datamat[9,2] <- 2*(basefreqs[3]*basefreqs[1]*basefreqs[2]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[1])
        datamat[10,2] <- 2*(basefreqs[2]*basefreqs[3]*basefreqs[4]) + 2*(basefreqs[3]*basefreqs[2]*basefreqs[4])
        #
        #Normalize frequencies so that they sum to 1
        datamat[,2] <- datamat[,2]/sum(datamat[,2])
        #Expected number of mutations
        datamat[,4] <- datamat[,2]*total
        #Relative observed / expected
        datamat[,5] <- datamat[,3]/datamat[,4]
    }
    return(list(resultstable = datamat, totalmutations = total))
}


#Function to calculate triplet frequencies (order is the same as in triplet sorter function
calc.tripfreq <- function(basefreqs) {
    tripfreq <- c(basefreqs[1]^3 + basefreqs[2]^3, 2*(basefreqs[2]*basefreqs[1]^2) + 2*(basefreqs[1]*basefreqs[2]^2), basefreqs[1]*basefreqs[2]^2 + basefreqs[2]*basefreqs[1]^2, 2*(basefreqs[3]*basefreqs[1]^2) + 2*(basefreqs[4]*basefreqs[2]^2), basefreqs[1]*basefreqs[3]^2 + basefreqs[2]*basefreqs[4]^2, 2*(basefreqs[1]*basefreqs[3]*basefreqs[4]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[3]), 2*(basefreqs[4]*basefreqs[1]^2) + 2*(basefreqs[3]*basefreqs[2]^2), basefreqs[1]*basefreqs[4]^2 + basefreqs[2]*basefreqs[3]^2, 2*(basefreqs[3]*basefreqs[1]*basefreqs[2]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[1]), 2*(basefreqs[1]*basefreqs[2]*basefreqs[4]) + 2*(basefreqs[3]*basefreqs[2]*basefreqs[1]), basefreqs[3]^3 + basefreqs[4]^3, 2*(basefreqs[4]*basefreqs[3]^2) + 2*(basefreqs[3]*basefreqs[4]^2), basefreqs[3]*basefreqs[4]^2 + basefreqs[4]*basefreqs[3]^2, 2*(basefreqs[1]*basefreqs[3]^2) + 2*(basefreqs[2]*basefreqs[4]^2), basefreqs[3]*basefreqs[1]^2 + basefreqs[4]*basefreqs[2]^2, 2*(basefreqs[1]*basefreqs[3]*basefreqs[4]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[3]), 2*(basefreqs[2]*basefreqs[3]^2) + 2*(basefreqs[1]*basefreqs[4]^2), basefreqs[3]*basefreqs[2]^2 + basefreqs[4]*basefreqs[1]^2, 2*(basefreqs[3]*basefreqs[1]*basefreqs[2]) + 2*(basefreqs[4]*basefreqs[2]*basefreqs[1]), 2*(basefreqs[2]*basefreqs[3]*basefreqs[4]) + 2*(basefreqs[3]*basefreqs[2]*basefreqs[4]))
    #I guess there is rounding error or something... make sure that frequencies sum to 1
    tripfreq <- tripfreq/sum(tripfreq)
    return(tripfreq)
}


#Function that consolidates the duplicated regions
#There are many regions that end and begin with adjacent base pairs but have different methods only
consolidate.dupregions <- function(dupregions, method = "Blast") {
    #First take only those regions that have been detected with the chosen method
    ind <- grepl(method, dupregions$Methods)
    dupregions <- dupregions[ind,] #Filter regions
    #
    ### Consolidate regions ###
    ##Make results matrix
    res.mat <- data.frame(matrix(rep(NA, nrow(dupregions)*3), ncol = 3))
    colnames(res.mat) <- c("Chromosome", "Start", "End")
    #
    #Initialize current chromosome and start
    current.chr <- dupregions[1,1]
    current.start <- dupregions[1,2]
    current.end <- dupregions[1,3]
    #Loop over all duplicated regions
    for(i in 2:nrow(dupregions)) {
        #Check is the next start right after current end
        if(current.chr == dupregions[i,1] & current.end + 1 == dupregions[i,2]) { current.end <- dupregions[i,3] }
        #If next start is further away or in a different chromosome, or end has been reached
        #Update the current positions
        if(current.end + 1 < dupregions[i,2] | current.chr != dupregions[i,1] | i == nrow(dupregions)) {
            res.mat[i,1] <- as.character(current.chr)
            res.mat[i,2:3] <- c(current.start, current.end)
            current.chr <- dupregions[i,1]
            current.start <- dupregions[i,2]
            current.end <- dupregions[i,3]
        }
    }
    #Remove useless stuff from res.mat
    res.mat <- res.mat[!is.na(res.mat[,1]),] #Remove all rows that are NA
    return(res.mat)
}


##Function to check whether a mutation is in duplicated regions
in.dup <- function(aineisto, dupreg) {
    res <- rep(0, nrow(aineisto)) #Store results here
    for( i in 1:nrow(aineisto)) { #Loop over all mutations
        #Check in which chromosome we are in
        chr <- as.character(aineisto[i,1]) #For compatibility
        current.dupreg <- filter(dupreg, Chromosome == chr) #Get the duplicated regions of current chr
        #Is there any region where position of mutation is larger than region start
        #but smaller than region end?
        check <- aineisto[i,2] > current.dupreg$Start & aineisto[i,2] < current.dupreg$End
        if(any(check) == T) {res[i] <- TRUE} else {res[i] <- FALSE}
    } #Done looping over mutations
    return(res)
}



#This function checks whether a mutation is in a centromeric region
##Note aineisto has to contain only the seven chromosomes
in.centromeric <- function(aineisto, cent) {
    res <- rep(0, nrow(aineisto)) #Store results here
    for(i in 1:nrow(aineisto)) { #Loop over all mutations
        #Check in which chromosome we are in
        chr <- aineisto$Chromosome[i]
        ind <- which(chr == cent[,1]) #Row index
        if(aineisto$Position[i] >= cent[ind,2] & aineisto$Position[i] <= cent[ind,3]) { res[i] <- TRUE } else { res[i] <- FALSE }
    }
    return(res)
}

##This function gets the number of called sites according the certain coordinates
called.sites <- function(calledsites, criteria, type = "cent") {
    #Make results list
    results <- rep(list(0),7)
    #
    if(type == "cent") {
        ##Loop over chromosomes
        for(i in 1:7) {
            chrind <- calledsites$Chromosome == criteria[i,1] & calledsites$Position >= criteria[i,2] & calledsites$Position <= criteria[i,3]
            results[[i]] <- table(calledsites[chrind,3])
        }
        return(results)
    } #End processing centromeric regions
    if(type == "dup") {
        calledsites[,3] <- factor(calledsites[,3]) #As factor so that table does not drop levels
        ##Loop over chromosomes
        for(i in 1:7) {
            chr <- paste("Supercontig_12.",i, sep = "")
            called.chr <- calledsites[calledsites$Chromosome == chr,] #Filter current chromosome
            current.coord <- criteria[criteria$Chromosome == chr,] #Filter current dupregions
            ##Loop over coordinates
            for(j in 1:nrow(current.coord)) {
                chrind <- called.chr$Position >= current.coord[j,2] & called.chr$Position <= current.coord[j,3]
                results[[i]] <- results[[i]] + table(called.chr[chrind,3])
            }
        }
        return(results)
    } #End processing duplicate regions       
}

#This functions excludes centromeric regions from chip-seq domains
ex.centro <- function(domains, cent) {
    domains$Chromosome <- factor(domains$Chromosome)
    check <- rep(0,nrow(domains)) #Vector for checking
    for(i in 1:7) { #Loop over all chromosomes
        sc <- cent[i,2]
        ec <- cent[i,3]
        check <- check + (domains$Chromosome == cent[i,1] & domains[,2] <= sc & domains[,3] >= sc |
                  domains$Chromosome == cent[i,1] & domains[,2] <= ec & domains[,3] >= ec | domains$Chromosome == cent[i,1] & domains[,2] >= sc & domains[,3] <= ec)
    }
    results <- domains[!check,] #Drop domains that overlap centromeric regions
    return(list(results, domains[as.logical(check),])) #Need as logical to convert to boolean
}

#This function calculates GC-content of given bases (returns %)
calculate.gc <- function(bases) { (sum(bases == "C") + sum(bases == "G") ) / (sum(bases == "A") + sum(bases == "G") + sum(bases == "C") + sum(bases == "T")) * 100 }

find.gc.window <- function(aineisto, gcstats) {
    #Set up results frame
    results <- rep(0, nrow(aineisto))
    #Assumes that aineisto and gcstats only contain the same chromosomes
    #Loop over all mutations
    for(i in 1:nrow(aineisto)) {
        current.chr <- as.character(aineisto[i,1]) #Take current chromosome
        current.gc <- filter(gcstats, Chromosome == current.chr)
        ind <- which.min(abs(current.gc$Mid - aineisto[i,2])) #Find index of nearest window
        results[i] <- current.gc$GCcont[ind] #Store GC content of window
    }
    return(results)
}


##Okay write a function that implements Watterson's theta with missing data and calculate it over 100 bp (or some interval) windows

#Function to calculate the a_nx value
#a_n = sigma^(n-1)_i=1 1/i
anx.calc <- function(n) { sum( 1 / (1:(n-1))) }

#Then for a given interval
theta.W <- function(aineisto) {
    #
    S <- sum(aineisto$type == "P") #Count number of polymorphic sites
    #Count how many bases inds have been sampled for a given site
    nx <- rowSums(!is.na(aineisto[,-c(1:4)]))
    anx <- sapply(nx, anx.calc) #Calculate a_nx for each site
    theta <- S / (sum(anx)) #Theta
    return(theta)
}


###Function to loop over using intervals of n bp and calculate theta
summary.over.intervals <- function(DT, interval = 100) {
    #Total number of bases
    nbases <- nrow(DT)
    #Number of whole intervals
    nint <- floor(nbases/interval)
    #Note, this leaves the final interval at chromosome end missing. If interval is small (< 2 kb) I don't think this is a problem since there are no good quality genotypes at chromosome ends anyway. So, data is missing no matter what. #Alternative would be to make the last interval shorter
    #
    #Prepare the results file
    res.mat <- data.frame(matrix(rep(0, nint*4), ncol = 4))
    colnames(res.mat) <- c("start", "end", "anx", "theta.W")
    res.mat[,1] <- (1:nint)*interval - (interval-1) #Interval start coordinates
    res.mat[,2] <- (1:nint)*interval #Interval end coordinates
    #
    #Loop over all intervals
    for(i in 1:nint) {
        aineisto <- DT[res.mat[i,1]:res.mat[i,2],] #Take the data from the interval
        aineisto <- aineisto[type != "N"] #Remove all alignment gaps
        #Removing all alignment gaps can results in an empty data table
        if(nrow(aineisto) > 5) { #Need to have more than 5 high quality bases in an interval
            #
            #totalsites <- nrow(aineisto) #Total number of genotyped sites in interval
            theta <- theta.W(aineisto) #Calculate theta.W
        }
        else { theta <- NA } #Assign theta to NA if there were not enough sites to calculate it
        res.mat[i,4] <- theta #Store results
    }
    #
    return(res.mat)
}


#Combine gc content and theta
combine.gc.theta <- function(gc.data, theta.chr.all) {
    #Make results vector
    theta <- rep(0, nrow(gc.data))
    #loop over all windows in gc.data
    for(i in 1:nrow(gc.data)) {
        cur.chr <- gc.data[i,2]
        cur.mid <- gc.data[i,6]
        #Find current window from theta.chr.all
        cur.ind <- which(cur.chr == theta.chr.all$Chromosome & cur.mid == theta.chr.all$mid)
        if(identical(cur.ind, integer(0)) == FALSE ) { theta[i] <- theta.chr.all[cur.ind,3] }    else    { theta[i] <- NA }
    } #Done looping over all windows
    return(theta)
}

##Write a function that gets GC content, and H3K9, H3K27, centromere domains
combine.theta.all <- function(theta.chr.all, nucstats.plot, cent, h3k9.plot, h3k27.plot ) {
    res.mat <- data.frame(matrix(rep(0, nrow(theta.chr.all)*4), ncol = 4))
    colnames(res.mat) <- c("GC", "centromere", "H3K9", "H3K27")
    #
    for(i in 1:nrow(theta.chr.all)) { #Loop over all windows
        cur.chr <- theta.chr.all[i,5] #Current chromosome
        cur.mid <- theta.chr.all[i,4] #Current window mid-point
        ws <- theta.chr.all[i,1] #window start
        we <- theta.chr.all[i,2] #window end
        #Find current window from nuctstats.plot for GC
        cur.ind <- which(cur.chr == nucstats.plot$Chromosome & cur.mid == nucstats.plot$Mid)
        if(identical(cur.ind, integer(0)) == FALSE ) { res.mat[i,1] <- nucstats.plot[cur.ind,4] }    else    { res.mat[i,1] <- NA }
        #Find the domains
        #centromere
        res.mat[i,2] <- ifelse(any(cur.chr == cent$Chromosome & ws < cent$start & we > cent$start | cur.chr == cent$Chromosome & ws >= cent$start & we <= cent$end | cur.chr == cent$Chromosome & ws < cent$end & we > cent$end), 1, 0)
        #h3k9
        res.mat[i,3] <- ifelse(any(cur.chr == h3k9.plot$Chromosome & ws < h3k9.plot$Start & we > h3k9.plot$Start | cur.chr == h3k9.plot$Chromosome & ws >= h3k9.plot$Start & we <= h3k9.plot$End | cur.chr == h3k9.plot$Chromosome & ws < h3k9.plot$End & we > h3k9.plot$End), 1, 0)
        #h3k27
        res.mat[i,4] <- ifelse(any(cur.chr == h3k27.plot$Chromosome & ws < h3k27.plot$Start & we > h3k27.plot$Start | cur.chr == h3k27.plot$Chromosome & ws >= h3k27.plot$Start & we <= h3k27.plot$End | cur.chr == h3k27.plot$Chromosome & ws < h3k27.plot$End & we > h3k27.plot$End), 1, 0)
    }
    return(res.mat) #return results
}


#This function takes a list of called sites and calculates GC-content and nucleotide frequencies
calc.gc.from.called <- function(sites.data) {
    A <- sum(unname(sapply(sites.data, '[[', 1)))
    C <- sum(unname(sapply(sites.data, '[[', 2)))
    G <- sum(unname(sapply(sites.data, '[[', 3)))
    T <- sum(unname(sapply(sites.data, '[[', 4)))
    #
    total <- A + C + G + T
    #
    GC.cont <- (G + C) / total
    fA <- A / total
    fT <- T / total
    fC <- C / total
    fG <- G / total
    tulokset <- c(GC.cont, fA, fT, fC, fG)
    names(tulokset) <- c("GC-cont", "f(A)", "f(T)", "f(C)", "f(G)")
    return(tulokset)
}

#This function calculates expected trinucleotide frequencies for trinucleotides
calc.exp.trinuc.freq <- function(trinucleotides, basefreq) {
    ntri <- length(trinucleotides) #Should be 64
    res.vec <- rep(0, ntri)
    #
    for(i in 1:ntri) {
        trinuc <- as.character(trinucleotides[i])
        bases <- unlist(strsplit(trinuc, ""))
        bases <- replace(bases, bases == "A", basefreq[1])
        bases <- replace(bases, bases == "T", basefreq[2])
        bases <- replace(bases, bases == "C", basefreq[3])
        bases <- replace(bases, bases == "G", basefreq[4])
        res.vec[i] <- prod(as.numeric(bases))
    }
    return(res.vec)
}
