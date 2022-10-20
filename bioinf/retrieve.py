#This script retrieves the base and adjacent bases from the wormtable
import wormtable as wt
import sys
#import os
chromosome = sys.argv[1]
position = int(sys.argv[2])

aineisto = wt.open_table('/home/ililkron/Genomics/Neurospora/mutacc/genotypes3/allsamples.wt') # open the wormtable
chrompos_index = aineisto.open_index('CHROM+POS')

#samples = ['ANC1A', 'ANC2a']

values = ['CHROM', 'POS', 'REF']
empty = ['', "./.", None, ".", 'None'] #missing_value_character

#debug
#print(chromosome, position)
seq = ""
#For given chromosome and bases of the wormtable
for c, p, r in chrompos_index.cursor(values, start=(chromosome, position-1), stop=(chromosome, position+2)):
    if(len(r) > 1): r = r[0] #Get the first character of a string 
    seq+=r
    #print(c, p, r) #Debug purposes

print(seq)


