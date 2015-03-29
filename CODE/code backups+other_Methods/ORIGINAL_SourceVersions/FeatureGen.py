#!/sw/bin/python3.3
#! E:\Python33\python
#Read FASTA files from current directory, generate output txt file with values for features.
#30.6.2013 . Edited order of features generated, and bigrams from absolute to relative freq. + added length
'4.4.2014 - New version begun to use in aa molec.info proj. '
#import Bio

import re
from itertools import product
from math import log

import numpy as np
from Bio.SeqUtils import ProtParam as pp
from numba import autojit


unwanted_residues = ['U','X','Z','B']


'Reduced AA alph - for later conversion , https://github.com/ddofer/epitopes/blob/master/epitopes/reduced_alphabet.py'
AA_20= ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
		"M", "F", "P", "S", "T", "W", "Y", "V"]

' Assumes an "alpha_transformer class from the github..'
def dict_from_list(groups):
    result = {}
    for i, group in enumerate(groups):
        for c in group:
            result[c.upper()] = i
            result[c.lower()] = i
    return result

gbmr4 = dict_from_list(["ADKERNTSQ", "YFLIVMCWH", "G", "P"])

sdm12 = dict_from_list(
  ["A", "D", "KER", "N",  "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"]
)

hsdm17 = dict_from_list(
  ["A", "D", "KE", "R", "N", "T", "S", "Q", "Y", "F", "LIV", "M", "C", "W", "H", "G", "P"])


#hydrophilic vs. hydrophobic
hp2 = dict_from_list(["AGTSNQDEHRKP", "CMFILVWY"])

murphy10 = dict_from_list(
  ["LVIM", "C", "A", "G", "ST", "P", "FYW", "EDNQ", "KR", "H"])

alex6 = dict_from_list(["C", "G", "P", "FYW", "AVILM", "STNQRHKDE"])

aromatic2 = dict_from_list(["FHWY", "ADKERNTSQLIVMCGP"])

hp_vs_aromatic = dict_from_list(["H", "CMILV", "FWY", "ADKERNTSQGP"])



#Used tofilter out ,U,B,Z,X non standrd AAs from a given sequence/string. Returns true if illegals present
#@autojit
def contain_illegals(str, illegals):
  for c in illegals:
    if c in str:
      return True
    else:  return False


# Data is Imported from a FASTA sequence file:
# list - each entry is appended strings/the sequence
def parse_fasta(filename) :
   #f = open(sys.argv[1]+'/'+filename, 'r')
   f = open(filename, 'r')
   sequences =[]
   i = 0
   for line in f:
         if line.startswith('>'):
            i = (i+1)
         else:
            if not contain_illegals(line,unwanted_residues):
               if (len(sequences) - 1 == i):
                  sequences[i] += line.rstrip('\n')
               else:
                  sequences.append(line.rstrip('\n'))
   return sequences

#Writes out results to (same file each time) file with name "outname".txt, param is the key/values dict to print
#Modify: - 'for Key in dictionary, write out k[V], /n .... (after iterating over all key's values, close file)
## This loop syntax accesses the whole dict by looping over the .items() tuple list
def outWrite(param, outName) :
   out = open('./'+outName + '.txt', "w")
   for k, v in param.items() :
      out.write('\t'.join(map(str, v)))
      out.write('\n')
      #print(' Values: \n' + ('\t'.join(map(str, v))))
   out.close()
   return




@autojit
def KROne(seq) :
   seq = re.sub("\D", '0', seq)
   return seq
# GKR = 1 (Glycine)
@autojit
def GKROne(seq) :
   seq = seq.replace('G', '1')
   seq = re.sub("\D", '0', seq)
   return seq
#Hydrophibicity - "Charged AA" (DERHK) = '1'
@autojit
def chargeOne(seq) :
   seq = seq.replace('D', '1').replace('E', '1').replace('H', '1')
   seq = re.sub("\D", '0', seq)
   return seq

#@autojit
def MerCount(s) :
#Combo [list] holds all the 2^5 binary combinations
#Transpose list of permutations (combo) into a new defaultdict (key,0)
   d = dict.fromkeys(combo, 0)
   for i in xrange(len(s) - 4) :
      d[s[i :i + 5]] += 1
   return d.values()

# a and b are  relative volume of valine and Leu/Ile side chains to side chain of alanine.
# http://stackoverflow.com/questions/991350/counting-repeated-characters-in-a-string-in-python
@autojit
def aliphaticness(seq) :
   a = 2.9
   b = 3.9
   length = float(len(seq))
   alanine_per = (seq.count('A') / length )
   valine_per = (seq.count('V') / length )
   isoleucine_per = (seq.count('I') / length )
   leucine_per = (seq.count('L') / length )
   # Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
   aliphatic_index = (100 * (alanine_per + a * valine_per + b * (isoleucine_per + leucine_per )))
   return aliphatic_index

@autojit
def Autocorrellation(seq,loc) :
   seq = seq.replace('RR', '1').replace('KK', '1').replace('RK', '1').replace('KR', '1')
   seq = seq.replace('R', '1').replace('K', '1')
   seq = re.sub("\D", '0', seq)
   seq = map(int, seq)
   selfCor = np.correlate(seq, seq, 'full')
   #Avoid divide by zero error:
   if sum(seq)==0:
      return 0
   #Normalization - By Sum ("1's") or seq.length?
   autoCor = sorted(selfCor)[loc] / float(len(seq))
   #Second highest (NOT "100%" Autocorrelation : loc=-2
   return autoCor

#@autojit
# Counts percentage of occurences of biGrams (from bigramDict) for  a given seq
'Would be MUCH more efficient to first get all the counts (one run of the sequence, to def.dict), then %..!'
def bigramsFreq(seq, bigramDict) :
   length=(len(seq))
   for Aa in bigramDict.keys() :
      bigramDict[Aa] = ((seq.count(str(Aa)))/length) #TODO: Change here for speed
      #print bigramDict
   return bigramDict.values()

#@autojit
def seq_Entropy(seq) :
   length = float(len(seq))
   letters = list('ACDEFGHIKLMNPQRSTVWY')
   amino_acids = dict.fromkeys(letters, 0)
   for Aa in amino_acids :
      hits = []
      hits = [a.start() for a in list(re. finditer(Aa, seq))]
      p_prev = 0
      p_next = 1
      sum = 0
      while p_next < len(hits) :
         distance = (hits[p_next] - hits[p_prev]) / length
         sum += distance * log(distance, 2)
         p_prev = p_next
         p_next += 1
      amino_acids[Aa] = -sum
   return amino_acids.values()


def NSites (seq,length):
	''' N Glycosylation sites'''
	NSites = len(re.findall(r'N[^P][ST][^P]', seq))
	return (NSites/length)


   #counts # of suspected cleavage sites according to known motif model
   # Xxx-Xxx-Lys-Lys# , Xxx-Xxx-Lys-Arg# , Xxx-Xxx-Arg-Arg# ,  Arg-Xxx-Xxx-Lys# , # # Arg-Xxx-Xxx-Arg#
   # lysine: K. arginine: R.   #"Not Proline" = [^P]
@autojit
def cleavageCounts(seq) :
   count1 = len(re.findall('R.[^P][RK]', seq))
   #Arg-Xxx-Xxx-Arg|Lys
   count2 = (len(re.findall('.[^P][RK][RK]', seq)))
   return (count1 + count2)

#@autojit
# Bigramsdict
'USe COMBO code! bi_comb = (''.join(x) for x in product("ACDEFGHIKLMNPQRSTVWY", repeat=2)) '
"Could just do bigramDict=dict.fromkeys(bi_comb, 0)) "
def gen_BigramDict():
    bigramsAll = []
    for i in (product('ACDEFGHIKLMNPQRSTVWY', 2)) :
       bigramsAll.append(i[0] + i[1])
    bigramDict = dict.fromkeys(bigramsAll, 0)
    return bigramDict
# bigramDict = dict containing all the valid 2 letter bigrams as keys


#==============================================================================
# #main CODE:
#==============================================================================
combo = [''.join(x) for x in product('01', repeat=5)] #combo = list of all  possible [0/1] ,length 5 combinations
sampled_proteins = {}
sampled_seq = []
sequences = {} #will contain sequences as values (rather than JUST as "keys" as is case with 'sampled_proteins')
aa_groups = ('FYW', 'P', 'C', 'RHK', 'DE' , 'CSTMNQ','RK',
'ST','LASGVTIPMC','EKRDNQH')
bigramDict=gen_BigramDict()

#Read FASTA files from current directory
#for f in os.listdir(sys.argv[1]) :
# files = [f for f in os.listdir(os.curdir) if (os.path.isfile(f) and f.endswith(".fasta"))]
# for f in files:
# #for f in os.listdir(os.curdir) :
#    if (negative_set): #If negative sets marker set to true due to user input
#       if f.endswith(".fasta") and not f.startswith("_") and f.startswith('NEG') or f.startswith('Neg') :
#          Fasta_seq = parse_fasta(f)
# #num_samples = How many samples to sample at random from each file (in the given directory).
#          sampled_seq += random.sample(Fasta_seq, num_samples)
#          sampled_proteins = dict.fromkeys(sampled_seq,0)
#          #NP+ Positive  case
#    elif f.endswith(".fasta") and not f.startswith("_"):
Fasta_seq = parse_fasta(f)
sampled_seq +=Fasta_seq
sampled_proteins = dict.fromkeys(sampled_seq,0)
sequences = dict(zip(sampled_seq, sampled_seq))

#http://biopython.org/wiki/ProtParam
for seq in sampled_proteins :
   length = float(len(seq))
   Z = pp.ProteinAnalysis(sequences[seq].replace('X', '').replace('Z', ''))
   sampled_proteins[seq] = []
   window_mer = sequences[seq].replace('R', '1').replace('K', '1')

   sampled_proteins[seq].append(length)
   sampled_proteins[seq].append(Z.isoelectric_point())
   sampled_proteins[seq].append(Z.molecular_weight())
   sampled_proteins[seq].append(Z.gravy())
   sampled_proteins[seq].append(Z.aromaticity())
   sampled_proteins[seq].append(Z.instability_index())
   # (Z.flexibility())
   #protparam AA% returns a dict. of K(Aa):V(%) pairs
   sampled_proteins[seq].append(Autocorrellation(sequences[seq],-2))
   #sampled_proteins[seq].append(Autocorrellation(sequences[seq],-3))
   sampled_proteins[seq].append(aliphaticness(sequences[seq]))

   # N Glycosylation sites
   #NSites = len(re.findall(r'N[^P][ST][^P]', sequences[seq]))
   sampled_proteins[seq].append(NSites(sequences[seq],length))
   sampled_proteins[seq].append(hydroxSites(seq,length))

   #Optional - counts of suspected cleavage sites
   sampled_proteins[seq].append(cleavageCounts(sequences[seq]) / length)

   sampled_proteins[seq].extend(Z.get_amino_acids_percent().values())
   sampled_proteins[seq].extend(MerCount(KROne(window_mer)))
   sampled_proteins[seq].extend(MerCount(GKROne(window_mer)))
   sampled_proteins[seq].extend(MerCount(chargeOne(window_mer)))
   sampled_proteins[seq].extend(seq_Entropy(sequences[seq]))
      #AA Bigrams (400) frequencies:
   #sampled_proteins[seq].extend(bigramsFreq(sequences[seq], bigramDict))
   for Aa in aa_groups :
      sampled_proteins[seq].append(PrefixCount(Aa, sequences[seq]))
      sampled_proteins[seq].append(SuffixCount(Aa, sequences[seq]))
#Finally Write out results (with seperate lines per key/sequence) to a new file:
outWrite(sampled_proteins, outPut_filename)

#print 'length'
#print len(sampled_proteins)
print ('Done')


