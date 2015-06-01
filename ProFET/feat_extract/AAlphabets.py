"""
Check to make  alphabets, dicts, strings
 are persistant  and not recalculated each time this method called!!

Amino acid groupings from
'Reduced amino acid alphabets improve the sensitivity...' by
Peterson, Kondev, et al.
http://www.rpgroup.caltech.edu/publications/Peterson2008.pdf

Other alphabets from
http://bio.math-inf.uni-greifswald.de/viscose/html/alphabets.html

"""

'TODO:'
'Add AA Propensities (From BioPython, articles, comp.profiler, etc) - eg, AAindex, http://bioinf.icm.uu.se/kbib/project13/convertAAstoProperties/'


from collections import defaultdict

# ambiguous amina acids: [ 'aspartic acid or asparagine', 'leucine or isoleucine',
#  'glutamic acid[E] or glutamine[Q]'] :
ambiguous_aa = 'BJZX'
# special amino acids - 'selenocysteine', 'pyrralysine'
aa_special_alph = 'UO'
UNKNOWN_AA = "Z" #'unknown amino acid',

'''
ILLEGALS = [c for c in ambiguous_aa+aa_special_alph+'Z']
'''
ILLEGALS = ['B', 'J', 'Z', 'X', 'U', 'O', 'Z']
# print(ILLEGALS)

def TransDict_from_list(groups):
    '''
    Given a list of letter groups, returns a dict mapping each group to a
    single letter from the group - for use in translation.
    >>> alex6=["C", "G", "P", "FYW", "AVILM", "STNQRHKDE"]
    >>> trans_a6 = TransDict_from_list(alex6)
    >>> print(trans_a6)
    {'V': 'A', 'W': 'F', 'T': 'D', 'R': 'D', 'S': 'D', 'P': 'P',
     'Q': 'D', 'Y': 'F', 'F': 'F',
     'G': 'G', 'D': 'D', 'E': 'D', 'C': 'C', 'A': 'A',
      'N': 'D', 'L': 'A', 'M': 'A', 'K': 'D', 'H': 'D', 'I': 'A'}
    '''
    transDict = dict()

    result = {}
    for group in groups:
        g_members = sorted(group) #Alphabetically sorted list
        for c in g_members:
            # print('c' + str(c))
            # print('g_members[0]' + str(g_members[0]))
            result[c] = str(g_members[0]) #K:V map, use group's first letter as represent.
    # print(result)
    return result

def translate_sequence (seq, TranslationDict):
    '''
    Given (seq) - a string/sequence to translate,
    Translates into a reduced alphabet, using a translation dict provided
    by the TransDict_from_list() method.
    Returns the string/sequence in the new, reduced alphabet.
    Remember - in Python string are immutable..

    '''
    from_list = []
    to_list = []
    for k,v in TranslationDict.items():
        from_list.append(k)
        to_list.append(v)
    # TRANS_seq = seq.translate(str.maketrans(zip(from_list,to_list)))
    TRANS_seq = seq.translate(str.maketrans(TranslationDict))
    return TRANS_seq

def Get_Letters (TranslationDict):
    '''
    Given a TranslationDict,
    return, as string,  the letters retained after translation
    by that dict.
    '''
    e = set(TranslationDict.values())
    res = sorted (e)
    return ("".join(res))


AA20 = 'ACDEFGHIKLMNPQRSTVWY'  #"Standard alphabet"

'Invented, based roughly on Ofer8, for Dibasic cleavage prediction'
OferKR = TransDict_from_list(["C", "G", "P", "FYW", "AVILM", "R","K","H", "DE", "STNQ"])


ofer14=TransDict_from_list(["A", "D", "KR","E", "N", "TS","Q",
 "YF", "LIVM", "C", "W", "H", "G", "P"])
ofer13=TransDict_from_list(["A", "DE", "KR", "N", "TS","Q",
 "YF", "LIVM", "C", "W", "H", "G", "P"])

"modifed from wang-wang, Clustering of the Protein Design Alphabets by Using Hierarchical SOM "
ofer_w8 = TransDict_from_list(["FIL", "CY", "MVW", "HAT", "GP", "RK", "QSN", "DE"])
# ofer14=TransDict_from_list(['LIVM', 'D', 'G', 'A', 'C', 'N', 'H', 'KE','R', 'W', 'P', 'TSQ', 'YF'])

# Ofer7=TransDict_from_list(["C", "G", "P", "FYW", "AVILM","KR", "STNQHDE"])

'Look at: http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=1594927&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D1594927'
ofer_tail =TransDict_from_list(["FAILV","TS","C","G","P","KR","DE","MWY","NQH"])

ofer8=TransDict_from_list(["C", "G", "P", "FYW", "AVILM", "RKH", "DE", "STNQ"])

ofer_gbm5 = TransDict_from_list(["ANTSQ", "YFLIVMCWH","DKER" "G", "P"])
gbm4 = TransDict_from_list(["ADKERNTSQ", "YFLIVMCWH", "G", "P"])
sdm12 =TransDict_from_list(
    ["A", "D", "KER", "N",  "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"] )

hsdm17 =TransDict_from_list(
  ["A", "D", "KE", "R", "N", "T", "S", "Q", "Y", "F", "LIV",
  "M", "C", "W", "H", "G", "P"])

alex6=TransDict_from_list(["C", "G", "P", "FYW", "AVILM", "STNQRHKDE"])

shen7 =TransDict_from_list(["AGV","ILFP","YMTS","HNQW","RK","DE","C"])
'''
"Shen 7" From: "Predicting protein-protein interactions based only on sequences information.",
Shen J,Jiang H. et al. PNAS.  2007.
    Suggested ese as trimers, and/or with RNA 4-mers for predicting Protein-interaction,
(Protein-RNA idea, from: "Predicting RNA-Protein Interactions Using Only Sequence Information",
    BMC Bioinformatics. 2011; Dobbs et al)
'''


#hydrophilic vs. hydrophobic
hp2 =TransDict_from_list(["AGTSNQDEHRKP", "CMFILVWY"])
#Hydrophilic, Hydrophobic, Charged. (Custom Ofer)
hp3 = TransDict_from_list(["AGTSNQP", "CMFILVWY", "RKHED"])
#Hydrophilic, Hydrophobic, Positively Charged. (Custom Ofer)
hp3_Plus = TransDict_from_list(["AGTSNQPHED", "CMFILVWY", "RK"])

murphy10 =TransDict_from_list(  ["LVIM", "C", "A", "G", "ST",
 "P", "FYW", "EDNQ", "KR", "H"])

aromatic2 =TransDict_from_list(["FHWY", "ADKERNTSQLIVMCGP"])

hp_aroma_4 =TransDict_from_list(["H", "CMILV", "FWY", "ADKERNTSQGP"])

# https://github.com/biopython/biopython/blob/master/Bio/Alphabet/Reduced.py
murphy15 = {"L": "L",             "V": "L",             "I": "L",
             "M": "L",             "C": "C",             "A": "A",
             "G": "G",             "S": "S",             "T": "T",
             "P": "P",             "F": "F",             "Y": "F",
             "W": "W",             "E": "E",
             "D": "D",             "N": "N",             "Q": "Q",
             "K": "K",             "R": "K",             "H": "H"}

murphy_8 = {"L": "L",             "V": "L",             "I": "L",
            "M": "L",
             "C": "L",             "A": "A",             "G": "A",
             "S": "S",             "T": "S",
             "P": "P",             "F": "F",             "Y": "F",
             "W": "F",             "E": "E",
             "D": "E",             "N": "E",             "Q": "E",
             "K": "K",             "R": "K",             "H": "H"}
pc5 = {"I": "A", # Aliphatic
         "V": "A",         "L": "A",
         "F": "R", # Aromatic
         "Y": "R",         "W": "R",         "H": "R",
         "K": "C", # Charged
         "R": "C",         "D": "C",         "E": "C",
         "G": "T", # Tiny
         "A": "T",         "C": "T",         "S": "T",
         "T": "D", # Diverse
         "M": "D",         "Q": "D",         "N": "D",
         "P": "D"}



### ProFEAT propensity based scales: ####
# modified from ProFEAT + CTD.  (Intended for letter: number use there)
Disorder_3=TransDict_from_list(['ARSQEGKP','ILNCFYVW', 'DHMT'])
Hydrophobicity_3 = TransDict_from_list(['RKEDQN','GASTPHY','CLVIMFW'])
# #'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity
Polarity_3 = TransDict_from_list(['LIFWCMVY','PATGS','HQRKNED']) #ProFeat based
# #'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)
Polarizability_3 = TransDict_from_list(['GASDT','CPNVEQIL','KMHFRYW'])
# #'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)
Charge_3 = TransDict_from_list(['KR','ANCQGHILMFPSTWYV','DE'])
# #'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative
SecondaryStr_3 = TransDict_from_list(['EALMQKRH','VIYCWFT','GNPSD']) #Orig
# #1'stand for Helix; '2'stand for Strand, '3' stand for coil
NormVDWV_3 = TransDict_from_list(['GASTPDC','NVEQIL','MHKFRYW'])
# #1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)
SolventA_3 = TransDict_from_list(['ALFCGIVW','RKQEND','MPSTHY'])
# #1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate
SurfaceTension_3 = TransDict_from_list(['GQDNAHR','KTSEC','ILMFPWYV'])
# Hierarchical Classification of Protein Folds Using a Novel Ensemble Classifier. PLoS ONE

THREE_LETTER_ALPH_NAMES = ['Disorder_3','Hydrophobicity_3',
'Polarity_3','Polarizability_3','Charge_3','SecondaryStr_3',
'NormVDWV_3','SolventA_3','SurfaceTension_3']

'Call alphabet by name from this dict, then feed value into translator func:'
REDUCED_ALPHABETS_TRANSDICTS = {
'ofer14':(ofer14),
'ofer_w8':ofer_w8,
'ofer13':(ofer13),
'ofer8':ofer8,
'ofer_tail':ofer_tail,
'gbm4':(gbm4),
'murphy10':(murphy10),
'hp_aroma_4':(hp_aroma_4),
'hp2':(hp2),
'hp3':(hp3),
'alex6':(alex6),
'sdm12':(sdm12),
'hsdm17':(hsdm17),
'murphy15':murphy15,
'pc5':pc5,
'Disorder_3':Disorder_3,
'Hydrophobicity_3':Hydrophobicity_3,
'Polarity_3':Polarity_3,
'Polarizability_3':Polarizability_3,
'Charge_3':Charge_3,
'SecondaryStr_3':SecondaryStr_3,
'NormVDWV_3':NormVDWV_3,
'SolventA_3':SolventA_3,
'hp3_Plus':hp3_Plus,
'ofer_gbm5':ofer_gbm5,
'shen7':shen7
}


def Get_Alph_Letters(REDUCED_ALPHABETS_TRANSDICTS):
    REDUCED_ALPHABETS_LETTERS = defaultdict(str)
    for k,v in REDUCED_ALPHABETS_TRANSDICTS.items():
        REDUCED_ALPHABETS_LETTERS[k]=Get_Letters(v)
    REDUCED_ALPHABETS_LETTERS['AA20'] = 'ACDEFGHIKLMNPQRSTVWY' #Include full, nonreduced alphabet.
    return REDUCED_ALPHABETS_LETTERS

'Make this run once! Not every time method is called! (Potentially)'
REDUCED_ALPHABETS_LETTERS = Get_Alph_Letters(REDUCED_ALPHABETS_TRANSDICTS)

##############################################################################

if __name__=="__main__":

    '''Check this all works..'''
    # print(Reduced_Alphabets)
    protein="MQNEEDACLEAGYCLGTTLSSWRLHFMEEQSQSTMLMGIGIGALLTLAFVGIFFFVYRRVRRLRRAEDQQGTDDESDYQTEYEEELPAIPKETYADFQSTGIELDSDSEYEPSMLQGPPSLTSPEQSQDSFPWLPNQDDQGPRLEHPS"
    print(REDUCED_ALPHABETS_TRANSDICTS['gbm4'])
    print(translate_sequence(protein,REDUCED_ALPHABETS_TRANSDICTS['gbm4']))
    print(REDUCED_ALPHABETS_LETTERS)
    print(REDUCED_ALPHABETS_LETTERS['ofer14'])
    # for k,v in REDUCED_ALPHABETS_LETTERS.items():
    #     print (str(k), str(len(set(v))))
    print(translate_sequence(protein,REDUCED_ALPHABETS_TRANSDICTS['Charge_3']))



    '''
    #Internet:
    import string
    s='abracadabra'
    from_list='abcdr'
    to_list='?*!@|'
    print s.translate(string.maketrans(from_list,to_list)),
    # ?*|?!?@?*|?
    '''
