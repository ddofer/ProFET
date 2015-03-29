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
# Disordered_State = —(< Top — IDP > −0.542) ; < Top — IDP > is the average TOP-IDP value for a protein (or window).
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

# Hydrophobicity_kd_TMD = Bio.SeqUtils.ProtParamData.kd
# Hydrophilicity = Bio.SeqUtils.ProtParamData.hw
# Surface_access  = Bio.SeqUtils.ProtParamData.em
# Ja_transfer_energy = Bio.SeqUtils.ProtParamData.ja
# flexibility = Bio.SeqUtils.ProtParamData.Flex

Scales_Dict = {'hp':hp, 'Flex':Flex, 'hw':hw,
        'sa':sa,'ja':ja, 'TOP_IDP':TOP_IDP        }






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

