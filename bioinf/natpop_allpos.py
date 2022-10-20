#Use wormtable to call high quality SNPs from the vcf file and output every position
#Make separate files for each chromosome
from __future__ import division
import wormtable as wt
import sys
import numpy
import pdb

#Define a function to check if record is an indel or not
def indel(record):
    for i in record.ALT:
        if i:
            if len(i) > 1:
                return True
        else:
            return False
    if len(record.REF) > 1:
        return True
    else:
        return False

#Define a function to check if 
def is_het(genotype, separator):
    return (genotype.split(separator)[0] != genotype.split(separator)[1])

#Define a function to check if samples share the same genotype
def is_shared(genotypes):
    return(all(genotypes[i] == '1/1' for i in range(len(genotypes))) or all(genotypes[i] == '0/0' for i in range(len(genotypes))))

#Define a function to check if genotypes are equal
def GT_diff(genotype1, separator1, genotype2, separator2):
    return (genotype1.split(separator1) != genotype2.split(separator2))

#Define a function to check if genotype is alt
#def GT_alt(genotype, separator):

#This function checks that an array of bases is polymorphic
def is_polymorphic(bases):
    if 'NA' in bases:
        un = numpy.unique(bases)
        myind = numpy.where(un == 'NA') #are there NA's in the array?
        result = numpy.delete(un, myind)
        if len(result) >= 2: #two or more alleles
            return (True)
        else:
            return (False)
    else:
        un = numpy.unique(bases)
        if len(un) >= 2:
            return (True)
        else:
            return (False)

def all_NAs(bases):
    numna = sum(bases == 'NA')
    total = len(bases)
    if numna/total >= 0.9: #If 90% or more samples are NA's
        return (True)
    else:
        return (False)

aineisto = wt.open_table('./allsamples.wt') # open the wormtable
chrompos_index = aineisto.open_index('CHROM+POS')
#aineisto = wt.open_table('./small.wt')
f = open('natpop_allsites.txt', 'w') #Open file mutations.txt for writing
#f2 = open('natpop_info.txt', 'w')

samples = ['10948', '10886', '10932', '1165', '8816', '3223', '8845', '10908', '10904', '851', '1131', '8850', '8819', '4708', '4712', '6203', '4824', '8783', '8790', '3975', '10928', '10912', '3210', '10923', '10950', '10951', '10946', '3211', '10906', 'P4452', 'P4463', 'P4468', 'P4471', 'P4476', 'P4479', '10882', '10883', '10884', '10892', '10907', '10914', '10915', '10918', '10925', '10926', '10927', '10935', '10937', '10943', '10983', '3943', 'P4489', '5910', '4730', '1133', '4716'] #samples for natural populations (remove the ones I don't need)

##Need to write a header for the hapmap file
f.write('alleles' + '\t' + 'type' +'\t' + 'chrom' + '\t' + 'pos' + '\t' + "\t".join(samples) + '\n')

#samples.sort() #The ancestor should always be the first sample [0]
if len(samples) >2:
    print samples
    print 'Analysing', len(samples), 'samples'
else:
    print "Problem with sample list"
    sys.exit()

#Columns to be retrieved from wormtable, note that genotype are 'samplename.GT' so all these are generated
values = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO.MQ'] + [i+'.GT' for i in samples] + [i+'.DP' for i in samples] + [i+'.GQ' for i in samples] + [i+'.RGQ' for i in samples]
empty = ['', "./.", None, ".", 'None', '*'] #missing_value_character
chrs = ['1', '2', '3', '4', '5', '6', '7'] #Chromosomes to genotype

#For each chromosome
# Chr1 #9798893 + 1 since otherwise last base is not evaluated otherwise
for row in chrompos_index.cursor(values, start=("Supercontig_12.7",1), stop=("Supercontig_12.7",4255303)):
    #pdb.set_trace()
    chromosome, position, ref, alt, qual, filter, MQ = row[:7]
    
    #print chromosome, position #For debugging
    #Reference genotype qualities
    RGQs = row[7+3*len(samples):7+4*len(samples)]
    DPs = row[7+len(samples):7+2*len(samples)]
    DPs = [0L if i is None else i for i in DPs] #Converts None to 0
    
    #Initialize the array where to write bases
    bases = numpy.empty(len(samples), dtype = object)
    
    if '.' in chromosome:
        chromosome = chromosome.split(".")[1] #Take chromosome_number

    #pdb.set_trace()
    
    #1. Check that position is polymorphic in the sample
    if alt not in empty and ref not in empty and chromosome in chrs:
        alt = alt.split(",")
        #if type(alt) != list:
        #    print "Error!"
        #    alt = [alt]
        #2. Check that polymorphic site is a SNP, and not something else    
        if len(alt) == 1 and len(alt[0]) == 1 and len(ref) == 1:
            alt=alt[0]
            #position, qual, MQ = int(position), float(qual), float(MQ) #This line ???
            GTs = row[7:7+len(samples)]
            GQs = row[7+2*len(samples):7+3*len(samples)]
            GQs = [0L if i is None else i for i in GQs]
            #3. Check that overall quality of the site is OK
            if numpy.mean(DPs) >= 5 and numpy.mean(GQs) >= 30 and MQ >= 40:
                #4. Check each genotype and write bases
                for i in range(len(samples)):
                    if DPs[i] >= 5 and GQs[i] >= 30 and (is_het(GTs[i], GTs[i][1]) == False) and GTs[i] not in empty:
                        if GTs[i][0] == '1':
                            bases[i] = alt
                        if GTs[i][0] == '0':
                            bases[i] = ref
                    else:
                        bases[i] = 'NA' #If genotype fails quality of heterozygosity check, it is NA
                #5. Check again that there are not too many NA's and site is still polymorphic, if not count it is a monomorphic site
                if is_polymorphic(bases) == True and all_NAs(bases) == False:
                    #Write genotypes to output file
                    alleles = ''.join([ref, alt]) #Write the two alleles as ref alt
                    #rs = ''.join(['snp', snp_counter])
                    #write bases to output file
                    f.write(alleles + '\t' + 'P' + '\t' + str(chromosome) + '\t' + str(position) + '\t' + "\t".join(bases) + '\n')
                if is_polymorphic(bases) == False: #If only polymorphic base turned out to be bad...
                    #5.2 Check each genotype and write bases
                    f.write(ref + '\t' + 'M' + '\t' + str(chromosome) + '\t' + str(position) + '\t' + "\t".join(bases) + '\n')
            else:
                bases = numpy.repeat('NA', len(samples))
                f.write(ref + '\t' + 'N' + '\t' + str(chromosome) + '\t' + str(position) + '\t' + "\t".join(bases) + '\n')
        else:
            bases = numpy.repeat('NA', len(samples))
            f.write(ref + '\t' + 'N' + '\t' + str(chromosome) + '\t' + str(position) + '\t' + "\t".join(bases) + '\n')
    #1. If site is monomorphic
    if alt in empty:
        #For monomorphic sites check RGQ of ancestor
        RGQs = [0L if i is None else i for i in RGQs] #Converts None to 0
        #pdb.set_trace()
        #Alternative is to drop 0's [x for x in RGQs if z != None]
        if numpy.mean(RGQs) >= 30 and numpy.mean(DPs) >= 5 and len(ref) == 1:
            #bases = numpy.repeat(ref, len(samples))
            for i in range(len(samples)):
                if DPs[i] >= 5 and RGQs[i] >= 30:
                    bases[i] = ref
                else:
                    bases[i] = 'NA'
            f.write(ref + '\t' + 'M' + '\t' + str(chromosome) + '\t' + str(position) + '\t' + "\t".join(bases) + '\n')
        else:
            bases = numpy.repeat('NA', len(samples))
            if len(ref) > 1:
                ref = 'N'
            f.write(ref + '\t' + 'N' + '\t' + str(chromosome) + '\t' + str(position) + '\t' + "\t".join(bases) + '\n')
#Store information about called sites
