#Implementing the mutation finder that I did previously with PyVCF in wormtable
import wormtable as wt
import sys

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

aineisto = wt.open_table('./genotypes3/allsamples.wt') # open the wormtable
f = open('sitess_MA.txt', 'w') #Open file mutations.txt for writing
f2 = open('mut_info_MA.txt', 'w')
f3 = open('called_sites_MA.txt', 'w')

#Need to define all samples that are investigated
samples = ['ANC1A', 'L1G40', 'L2G40', 'L3G40', 'L4G40', 'L5G40', 'L6G40', 'L7G40', 'L8G40', 'L9G40', 'L10G40', 'L11G40', 'L12G40', 'L13G40', 'L14G40', 'L15G40', 'L16G40', 'L17G40', 'L18G40', 'L19G40', 'L20G40'] #mat A samples

#samples = ['ANC2a', 'L21G40', 'L22G40', 'L23G40',  'L24G40', 'L25G40', 'L26G40', 'L27G40', 'L28G40', 'L29G40', 'L30G40', 'L31G40', 'L32G40', 'L33G40', 'L34G40', 'L35G40', 'L36G40', 'L37G40', 'L38G40', 'L39G40', 'L40G40'] #mat a samples

##Need to write a header for f2 file

#samples.sort() #The ancestor should always be the first sample [0]
if len(samples) >2:
    print samples
else:
    print "Problem with sample list"
    sys.exit()

#Columns to be retrieved from wormtable, note that genotype are 'samplename.GT' so all these are generated
values = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO.MQ'] + [i+'.GT' for i in samples] + [i+'.DP' for i in samples] + [i+'.GQ' for i in samples] + [i+'.RGQ' for i in samples]
empty = ['', "./.", None, ".", 'None'] #missing_value_character

called_sites = 0
called_A = 0
called_T = 0
called_C = 0
called_G = 0

#For each row of the wormtable
for row in aineisto.cursor(values):
    chromosome, position, ref, alt, qual, filter, MQ = row[:7]
    
    #Reference genotype qualities
    RGQs = row[7+3*len(samples):7+4*len(samples)]

    #Check that position is polymorphic in the sample
    if alt not in empty:
        alt = alt.split(",")
        if type(alt) != list:
            print "Error!"
            alt = [alt]
        #position, qual, MQ = int(position), float(qual), float(MQ) #This line ???
        GTs = row[7:7+len(samples)]
        DPs = row[7+len(samples):7+2*len(samples)]
        GQs = row[7+2*len(samples):7+3*len(samples)]
        
        #Filter positions where all samples share the same genetic background SNP
        #For some samples the ancestor has not been sequenced and thus genetic background SNPs have to be distinguished if all samples share the same SNP
        #Note, May need to do some manual curation since in theory if the same mutation has happened twice and third sample has low GQ etc.
        if is_shared(GTs[1:]) == False:
            #Loop over all samples and test for filters
            for i in range(len(samples)):
                #Check that ancestor genotype is different from sample genotype
                if (GTs[0] and GTs[i]) not in empty \
                and GT_diff(GTs[0], GTs[0][1], GTs[i], GTs[i][1]): #and GTs[0] != GTs[i]:
                    #Then do some quality filtering, on depth, heterozygosity and genotype quality
                    if DPs[0] >= 5 and DPs[i] >= 5 \
                    and (is_het(GTs[0], GTs[0][1]) == False) and (is_het(GTs[i], GTs[i][1]) == False) \
                    and GQs[0] >= 30 and GQs[i] >= 30:
                        #Get mutation, (get alt bases as pure strings not list)
                        if int(GTs[0][0]) == 0: anc_base = ref
                        else: anc_base = alt[int(GTs[0][0])-1]
                        if int(GTs[i][0]) == 0: new_base = ref
                        else: new_base = alt[int(GTs[i][0])-1]
                        if ref == "A":
                            called_A += 1 #Increment number of called sites
                        if ref == "T":
                            called_T += 1
                        if ref == "G":
                            called_G += 1
                        if ref == "C":
                            called_C += 1
                        #Print some information to the screen
                        print chromosome, position, "in individual", samples[i], \
                        "Ancestor GT=", GTs[0], "GQ =", GQs[0], \
                        "Sample GT=", GTs[i], "GQ =", GQs[i], \
                        "Mutation is from", anc_base, "to", new_base
                        #f.write(chromosome + '\t' + str(position) + '\n')
                        f2.write(chromosome + '\t' + str(position) + '\t' + samples[i] + \
                        '\t' + GTs[0] + '\t' + str(GQs[0]) + '\t' + GTs[i] + '\t' + str(GQs[i]) + \
                        '\t' + anc_base + '\t' + new_base + '\n')
                        f3.write(chromosome + '\t' + str(position) + '\t' + ref + '\n')
    else:
        #For monomorphic sites check RGQ of ancestor
        if RGQs[0] >= 30:
            if ref == "A":
                called_A += 1 #Increment number of called sites
            if ref == "T":
                called_T += 1
            if ref == "G":
                called_G += 1
            if ref == "C":
                called_C += 1
            f3.write(chromosome + '\t' + str(position) + '\t' + ref + '\n')
called_sites = called_A + called_T + called_G + called_C #Total number of called sites
#Store information about called sites
f.write("Number of called sites: " + str(called_sites) + '\n' + "Number of called A's: " + str(called_A) + '\n' + "Number of called T's: " + str(called_T) + '\n' + \
"Number of called C's: " + str(called_C) + '\n' + "Number of called G's: " + str(called_G) )
#Note need to store also number of A's, T's , C's and G's called

f.close()
f2.close()
f3.close()
### Done        
