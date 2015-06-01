'''
AA Propensity Scales.

TODO: (Add/Note "combined" metrics: Georgiev scales. Kidera factors.. )
Some from BioPython.
(BioPython stored them as dictionaries, e.g: Bio.SeqUtils.ProtParamData.kd).

May need to be SCALED to 0-1 range ?!?

Data initially acquired from BioPython:
https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParamData.py
    Bio.SeqUtils.ProtParamData

    Some more descriptors:
    https://github.com/ddofer/Protein-Descriptors/blob/master/src/csdsML/Descriptors.py
'''

# Kyte & Doolittle {kd} index of hydrophobicity
hp = {'A': 1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C': 2.5,
      'Q':-3.5, 'E':-3.5, 'G':-0.4, 'H':-3.2, 'I': 4.5,
      'L': 3.8, 'K':-3.9, 'M': 1.9, 'F': 2.8, 'P':-1.6,
      'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V': 4.2 }

# Flexibility
# Normalized flexibility parameters (B-values), average (Vihinen et al., 1994)
Flex= {'A': 0.984, 'C': 0.906, 'E': 1.094, 'D': 1.068,
       'G': 1.031, 'F': 0.915, 'I': 0.927, 'H': 0.950,
       'K': 1.102, 'M': 0.952, 'L': 0.935, 'N': 1.048,
       'Q': 1.037, 'P': 1.049, 'S': 1.046, 'R': 1.008,
       'T': 0.997, 'W': 0.904, 'V': 0.931, 'Y': 0.929}

# Hydrophilicity
# 1 Hopp & Wood
# Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).
hw = {'A':-0.5, 'R': 3.0, 'N': 0.2, 'D': 3.0, 'C':-1.0,
      'Q': 0.2, 'E': 3.0, 'G': 0.0, 'H':-0.5, 'I':-1.8,
      'L':-1.8, 'K': 3.0, 'M':-1.3, 'F':-2.5, 'P': 0.0,
      'S': 0.3, 'T':-0.4, 'W':-3.4, 'Y':-2.3, 'V':-1.5 }

# Surface accessibility {"em"}
# 1 Emini Surface fractional probability
sa = {'A': 0.815, 'R': 1.475, 'N': 1.296, 'D': 1.283, 'C': 0.394,
      'Q': 1.348, 'E': 1.445, 'G': 0.714, 'H': 1.180, 'I': 0.603,
      'L': 0.603, 'K': 1.545, 'M': 0.714, 'F': 0.695, 'P': 1.236,
      'S': 1.115, 'T': 1.184, 'W': 0.808, 'Y': 1.089, 'V': 0.606 }

# 2 Janin Interior to surface transfer energy scale
ja = {'A': 0.28, 'R':-1.14, 'N':-0.55, 'D':-0.52, 'C': 0.97,
      'Q':-0.69, 'E':-1.01, 'G': 0.43, 'H':-0.31, 'I': 0.60,
      'L': 0.60, 'K':-1.62, 'M': 0.43, 'F': 0.46, 'P':-0.42,
      'S':-0.19, 'T':-0.32, 'W': 0.29, 'Y':-0.15, 'V': 0.60 }

# Disorder Propensity scale
#"TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder"
#Campen, Uversky, Dunker et al. Protein Pept Lett. 2008

# Positive  values indicate protein (or windows) are likely to be ordered, etc' /
TOP_IDP ={'A':0.06, 'R':0.180, 'N':0.007, 'D':0.192, 'C': 0.02,
      'Q':0.318, 'E':0.736, 'G': 0.166, 'H':0.303, 'I': -0.486,
      'L': -0.326, 'K':0.586, 'M': -0.397, 'F': -0.697, 'P':0.987,
      'S':0.341, 'T':0.059, 'W': -0.884, 'Y':-0.510, 'V': -0.121 }

# https://github.com/ddofer/Protein-Descriptors/blob/master/src/csdsML/Descriptors.py
polarizability= {'A':0.046,'R':0.291,'N':0.134,'D':0.105,'C': 0.128,'Q':0.180,
                                         'E':0.151,'G':0.000,'H':0.230,'I':0.186,'L':0.186,'K':0.219,'M':0.221,
                                          'F':0.290,'P':0.131,'S':0.062,'T':0.108,'W':0.409,'Y':0.298,'V':0.140}

ASAInTripeptide = {'A':115,'R':225,'N':160,'D':150,'C':135,'Q':180,
                                         'E':190,'G':75,'H':195,'I':175,'L':170,'K':200,'M':185,
                                          'F':210,'P':145,'S':115,'T':140,'W':255,'Y':230,'V':155}
Volume = {'A':52.6,'R':109.1,'N':75.7,'D':68.4,'C':68.3,'Q':89.7,
                                         'E':84.7,'G':36.3,'H':91.9,'I':102.0,'L':102.0,'K':105.1,'M':97.7,
                                          'F':113.9,'P':73.6,'S':54.9,'T':71.2,'W':135.4,'Y':116.2,'V':85.1}

StericParam = {'A':0.52,'R':0.68,'N':0.76,'D':0.76,'C':0.62,'Q':0.68,
                                         'E':0.68,'G':0.00,'H':0.70,'I':1.02,'L':0.98,'K':0.68,'M':0.78,
                                          'F':0.70,'P':0.36,'S':0.53,'T':0.50,'W':0.70,'Y':0.70,'V':0.76}
Mutability = {'A':100,'R':65,'N':134,'D':106,'C':20,'Q':93,
                                         'E':102,'G':49,'H':66,'I':96,'L':40,'K':56,'M':94,
                                          'F':41,'P':56,'S':120,'T':97,'W':18,'Y':41,'V':74}

# Hydrophobicity_kd_TMD = Bio.SeqUtils.ProtParamData.kd
# Hydrophilicity = Bio.SeqUtils.ProtParamData.hw
# Surface_access  = Bio.SeqUtils.ProtParamData.em
# Ja_transfer_energy = Bio.SeqUtils.ProtParamData.ja
# flexibility = Bio.SeqUtils.ProtParamData.Flex


'GeorgievScales:'
#Acquired from georgiev's paper of AAscales using helper script "GetTextData.py". + RegEx cleaning
gg_1 = {'Q': -2.54, 'L': 2.72, 'T': -0.65, 'C': 2.66, 'I': 3.1, 'G': 0.15, 'V': 2.64, 'K': -3.89, 'M': 1.89, 'F': 3.12, 'N': -2.02, 'R': -2.8, 'H': -0.39, 'E': -3.08, 'W': 1.89, 'A': 0.57, 'D': -2.46, 'Y': 0.79, 'S': -1.1, 'P': -0.58}
gg_2 = {'Q': 1.82, 'L': 1.88, 'T': -1.6, 'C': -1.52, 'I': 0.37, 'G': -3.49, 'V': 0.03, 'K': 1.47, 'M': 3.88, 'F': 0.68, 'N': -1.92, 'R': 0.31, 'H': 1, 'E': 3.45, 'W': -0.09, 'A': 3.37, 'D': -0.66, 'Y': -2.62, 'S': -2.05, 'P': -4.33}
gg_3 = {'Q': -0.82, 'L': 1.92, 'T': -1.39, 'C': -3.29, 'I': 0.26, 'G': -2.97, 'V': -0.67, 'K': 1.95, 'M': -1.57, 'F': 2.4, 'N': 0.04, 'R': 2.84, 'H': -0.63, 'E': 0.05, 'W': 4.21, 'A': -3.66, 'D': -0.57, 'Y': 4.11, 'S': -2.19, 'P': -0.02}
gg_4 = {'Q': -1.85, 'L': 5.33, 'T': 0.63, 'C': -3.77, 'I': 1.04, 'G': 2.06, 'V': 2.34, 'K': 1.17, 'M': -3.58, 'F': -0.35, 'N': -0.65, 'R': 0.25, 'H': -3.49, 'E': 0.62, 'W': -2.77, 'A': 2.34, 'D': 0.14, 'Y': -0.63, 'S': 1.36, 'P': -0.21}
gg_5 = {'Q': 0.09, 'L': 0.08, 'T': 1.35, 'C': 2.96, 'I': -0.05, 'G': 0.7, 'V': 0.64, 'K': 0.53, 'M': -2.55, 'F': -0.88, 'N': 1.61, 'R': 0.2, 'H': 0.05, 'E': -0.49, 'W': 0.72, 'A': -1.07, 'D': 0.75, 'Y': 1.89, 'S': 1.78, 'P': -8.31}
gg_6 = {'Q': 0.6, 'L': 0.09, 'T': -2.45, 'C': -2.23, 'I': -1.18, 'G': 7.47, 'V': -2.01, 'K': 0.1, 'M': 2.07, 'F': 1.62, 'N': 2.08, 'R': -0.37, 'H': 0.41, 'E': 0, 'W': 0.86, 'A': -0.4, 'D': 0.24, 'Y': -0.53, 'S': -3.36, 'P': -1.82}
gg_7 = {'Q': 0.25, 'L': 0.27, 'T': -0.65, 'C': 0.44, 'I': -0.21, 'G': 0.41, 'V': -0.33, 'K': 4.01, 'M': 0.84, 'F': -0.15, 'N': 0.4, 'R': 3.81, 'H': 1.61, 'E': -5.66, 'W': -1.07, 'A': 1.23, 'D': -5.15, 'Y': -1.3, 'S': 1.39, 'P': -0.12}
gg_8 = {'Q': 2.11, 'L': -4.06, 'T': 3.43, 'C': -3.49, 'I': 3.45, 'G': 1.62, 'V': 3.93, 'K': -0.01, 'M': 1.85, 'F': -0.41, 'N': -2.47, 'R': 0.98, 'H': -0.6, 'E': -0.11, 'W': -1.66, 'A': -2.32, 'D': -1.17, 'Y': 1.31, 'S': -1.21, 'P': -1.18}
gg_9 = {'Q': -1.92, 'L': 0.43, 'T': 0.34, 'C': 2.22, 'I': 0.86, 'G': -0.47, 'V': -0.21, 'K': -0.26, 'M': -2.05, 'F': 4.2, 'N': -0.07, 'R': 2.43, 'H': 3.55, 'E': 1.49, 'W': -5.87, 'A': -2.01, 'D': 0.73, 'Y': -0.56, 'S': -2.83, 'P': 0}
gg_10 = {'Q': -1.67, 'L': -1.2, 'T': 0.24, 'C': -3.78, 'I': 1.98, 'G': -2.9, 'V': 1.27, 'K': -1.66, 'M': 0.78, 'F': 0.73, 'N': 7.02, 'R': -0.99, 'H': 1.52, 'E': -2.26, 'W': -0.66, 'A': 1.31, 'D': 1.5, 'Y': -0.95, 'S': 0.39, 'P': -0.66}
gg_11 = {'Q': 0.7, 'L': 0.67, 'T': -0.53, 'C': 1.98, 'I': 0.89, 'G': -0.98, 'V': 0.43, 'K': 5.86, 'M': 1.53, 'F': -0.56, 'N': 1.32, 'R': -4.9, 'H': -2.28, 'E': -1.62, 'W': -2.49, 'A': -1.14, 'D': 1.51, 'Y': 1.91, 'S': -2.92, 'P': 0.64}
gg_12 = {'Q': -0.27, 'L': -0.29, 'T': 1.91, 'C': -0.43, 'I': -1.67, 'G': -0.62, 'V': -1.71, 'K': -0.06, 'M': 2.44, 'F': 3.54, 'N': -2.44, 'R': 2.09, 'H': -3.12, 'E': -3.97, 'W': -0.3, 'A': 0.19, 'D': 5.61, 'Y': -1.26, 'S': 1.27, 'P': -0.92}
gg_13 = {'Q': -0.99, 'L': -2.47, 'T': 2.66, 'C': -1.03, 'I': -1.02, 'G': -0.11, 'V': -2.93, 'K': 1.38, 'M': -0.26, 'F': 5.25, 'N': 0.37, 'R': -3.08, 'H': -1.45, 'E': 2.3, 'W': -0.5, 'A': 1.66, 'D': -3.85, 'Y': 1.57, 'S': 2.86, 'P': -0.37}
gg_14 = {'Q': -1.56, 'L': -4.79, 'T': -3.07, 'C': 0.93, 'I': -1.21, 'G': 0.15, 'V': 4.22, 'K': 1.78, 'M': -3.09, 'F': 1.73, 'N': -0.89, 'R': 0.82, 'H': -0.77, 'E': -0.06, 'W': 1.64, 'A': 4.39, 'D': 1.28, 'Y': 0.2, 'S': -1.88, 'P': 0.17}
gg_15 = {'Q': 6.22, 'L': 0.8, 'T': 0.2, 'C': 1.43, 'I': -1.78, 'G': -0.53, 'V': 1.06, 'K': -2.71, 'M': -1.39, 'F': 2.14, 'N': 3.13, 'R': 1.32, 'H': -4.18, 'E': -0.35, 'W': -0.72, 'A': 0.18, 'D': -1.98, 'Y': -0.76, 'S': -2.42, 'P': 0.36}
gg_16 = {'Q': -0.18, 'L': -1.43, 'T': -2.2, 'C': 1.45, 'I': 5.71, 'G': 0.35, 'V': -1.31, 'K': 1.62, 'M': -1.02, 'F': 1.1, 'N': 0.79, 'R': 0.69, 'H': -2.91, 'E': 1.51, 'W': 1.75, 'A': -2.6, 'D': 0.05, 'Y': -5.19, 'S': 1.75, 'P': 0.08}
gg_17 = {'Q': 2.72, 'L': 0.63, 'T': 3.73, 'C': -1.15, 'I': 1.54, 'G': 0.3, 'V': -1.97, 'K': 0.96, 'M': -4.32, 'F': 0.68, 'N': -1.54, 'R': -2.62, 'H': 3.37, 'E': -2.29, 'W': 2.73, 'A': 1.49, 'D': 0.9, 'Y': -2.56, 'S': -2.77, 'P': 0.16}
gg_18 = {'Q': 4.35, 'L': -0.24, 'T': -5.46, 'C': -1.64, 'I': 2.11, 'G': 0.32, 'V': -1.21, 'K': -1.09, 'M': -1.34, 'F': 1.46, 'N': -1.71, 'R': -1.49, 'H': 1.87, 'E': -1.47, 'W': -2.2, 'A': 0.46, 'D': 1.38, 'Y': 2.87, 'S': 3.36, 'P': -0.34}
gg_19 = {'Q': 0.92, 'L': 1.01, 'T': -0.73, 'C': -1.05, 'I': -4.18, 'G': 0.05, 'V': 4.77, 'K': 1.36, 'M': 0.09, 'F': 2.33, 'N': -0.25, 'R': -2.57, 'H': 2.17, 'E': 0.15, 'W': 0.9, 'A': -4.22, 'D': -0.03, 'Y': -3.43, 'S': 2.67, 'P': 0.04}

Atch_1 = {'A': 0.591, 'C': 1.343, 'E': 1.357, 'D': 1.05, 'G': 0.384, 'F': 1.006, 'I': 1.239, 'H': 0.336, 'K': 1.831, 'M': 0.663, 'L': 1.019, 'N': 0.945, 'Q': 0.931, 'P': 0.189, 'S': 0.228, 'R': 1.538, 'T': 0.032, 'W': 0.595, 'V': 1.337, 'Y': 0.26}
Atch_2 = {'A': 1.302, 'C': 0.465, 'E': 1.453, 'D': 0.302, 'G': 1.652, 'F': 0.59, 'I': 0.547, 'H': 0.417, 'K': 0.561, 'M': 1.524, 'L': 0.987, 'N': 0.828, 'Q': 0.179, 'P': 2.081, 'S': 1.399, 'R': 0.055, 'T': 0.326, 'W': 0.009, 'V': 0.279, 'Y': 0.83}
Atch_3 = {'A': 0.733, 'C': 0.862, 'E': 1.477, 'D': 3.656, 'G': 1.33, 'F': 1.891, 'I': 2.131, 'H': 1.673, 'K': 0.533, 'M': 2.219, 'L': 1.505, 'N': 1.299, 'Q': 3.005, 'P': 1.628, 'S': 4.76, 'R': 1.502, 'T': 2.213, 'W': 0.672, 'V': 0.544, 'Y': 3.097}
Atch_4 = {'A': 1.57, 'C': 1.02, 'E': 0.113, 'D': 0.259, 'G': 1.045, 'F': 0.397, 'I': 0.393, 'H': 1.474, 'K': 0.277, 'M': 1.005, 'L': 1.266, 'N': 0.169, 'Q': 0.503, 'P': 0.421, 'S': 0.67, 'R': 0.44, 'T': 0.908, 'W': 2.128, 'V': 1.242, 'Y': 0.838}
Atch_5 = {'A': 0.146, 'C': 0.255, 'E': 0.837, 'D': 3.242, 'G': 2.064, 'F': 0.412, 'I': 0.816, 'H': 0.078, 'K': 1.648, 'M': 1.212, 'L': 0.912, 'N': 0.933, 'Q': 1.853, 'P': 1.392, 'S': 2.647, 'R': 2.897, 'T': 1.313, 'W': 0.184, 'V': 1.262, 'Y': 1.512}

MinScales_Dict = {'hp':hp, 'hw':hw,
        'sa':sa, 'TOP_IDP':TOP_IDP,
                  'Atch_1':Atch_1,'Atch_2':Atch_2,
                  'Atch_3':Atch_3,'Atch_4':Atch_4,
                  'Atch_5':Atch_5}
#Some scales removed from "full" scales dict, due ot redundnacy, partic if minScales dict is used on subsegments of sequence.
#If MinScales dict is NOT used, then it's HGIHLY recomended to re-add these features!
Scales_Dict = {'hp':hp, 'ja':ja,
'polarizability':polarizability,'Mutability':Mutability,'Volume':Volume,
'ASAInTripeptide':ASAInTripeptide,
        'gg_1' : gg_1,'gg_2' : gg_2,'gg_3' : gg_3,'gg_4' : gg_4,'gg_5' : gg_5,
'gg_6' : gg_6,'gg_7' : gg_7,'gg_8' : gg_8,'gg_9' : gg_9,'gg_10' : gg_10,'gg_11' : gg_11}
#,'gg_12' : gg_12
#,'gg_13' : gg_13,'Atch_1':Atch_1,'Atch_2':Atch_2,'Atch_3':Atch_3,'Atch_4':Atch_4,'Atch_5':Atch_5,
#,'gg_14' : gg_14,'gg_15' : gg_15,
#,'gg_16' : gg_16,'gg_17' : gg_17,'gg_18' : gg_18,'gg_19' : gg_19}





########################################################################################
'''
From PyPro,
Authors: Dongsheng Cao and Yizeng Liang.
:
'''
# def _mean(listvalue):
#     """
#     ########################################################################################
#     The mean value of the list data.

#     Usage:

#     result=_mean(listvalue)
#     ########################################################################################
#     """
#     return sum(listvalue)/len(listvalue)
# ##############################################################################################
# def _std(listvalue,ddof=1):
#     """
#     ########################################################################################
#     The standard deviation of the list data.

#     Usage:

#     result=_std(listvalue)
#     ########################################################################################
#     """
#     mean=_mean(listvalue)
#     temp=[math.pow(i-mean,2) for i in listvalue]
#     res=math.sqrt(sum(temp)/(len(listvalue)-ddof))
#     return res
# ##############################################################################################
"TODO: Fix to use proper way of normalizing, AND scaling. (Maybe sci-kit learn's preprocessor?"
def NormalizeAAP(AAP):
    """
    ########################################################################################
    Centralize and normalize amino acid indices (Scales) before calculations.

    Usage:

    result=NormalizeEachAAP(AAP)

    Input: AAP is a dict containing the properties of 20 amino acids.

    Output: result is the a dict form containing the normalized properties.
    ########################################################################################
    """
    if len(AAP.values())!=20:
        print ('Some Amino Acids are missing')
    else:
        Result={}
        for i,j in AAP.items():
            Result[i]=(j-_mean(AAP.values()))/_std(AAP.values(),ddof=0)

    return Result
########################################################################################
'''GetAAindex1 Requires the GetAAIndex.py from PyPro: '''
    # def GetAAindex1(self,name,path='.'):
    #     """
    #     Get the amino acid property values from aaindex1

    #     Usage:

    #     result=GetAAIndex1(name)

    #     Input: name is the name of amino acid property (e.g., KRIW790103)

    #     Output: result is a dict form containing the properties of 20 amino acids
    #     """

    #     return GetAAIndex1(name,path=path)

