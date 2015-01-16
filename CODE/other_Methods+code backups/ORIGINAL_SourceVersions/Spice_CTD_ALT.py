# -*- coding: utf-8 -*-
'''
#####################################################################################

This module is used for computing the composition, transition and distribution

descriptors based on the different properties of AADs. Based on SPICE - Sequtil.py, + protein_class. CTD
(Returns a list rather than dict). (Has much Cleaner and more modular code)

References:

[1]: Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim. Prediction

of protein folding class using global description of amino acid sequence. Proc.Natl.

Acad.Sci.USA, 1995, 92, 8700-8704.

[2]:Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou Kim.

Recognition of a Protein Fold in the Context of the SCOP classification. Proteins:

Structure, Function and Genetics,1999,35,401-407.

[3] Composition profiler
http://www.cprofiler.org/help.html
[4] PROFEAT  (Table 2)
? - ProBias  http://lcg.rit.albany.edu/ProBias/help.html#ref1

#####################################################################################
'''

import string, math, copy, numpy

# AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]


# _DisorderPropensity={'1':'ARSQEGKP','2':'ILNCFYVW'}

# 'Horribly kludgy - Would be better to convert to a dict?:'


# ### ProtPy propensity based scales: ####
# # #Orig ProtPy scales - seem buggy! Letters  were lacking.. Some based on Profeat. I fixed manually
# _Hydrophobicity={'1':'RKEDQN','2':'GASTPHY','3':'CLVIMFW'}
# # #'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

# _NormalizedVDWV={'1':'GASTPDC','2':'NVEQIL','3':'MHKFRYW'}
# # #'1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

# _Polarity={'1':'LIFWCMVY','2':'PATGS','3':'HQRKNED'} #ProFeat based
# # _Polarity={'1':'LIFWCMVY','2':'CPNVEQIL','3':'KMHFRYW'} # Orig - diff. (Extra V?)
# # #'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

# _Polarizability={'1':'GASDT','2':'CPNVEQIL','3':'KMHFRYW'}
# # #'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)

# _Charge={'1':'KR','2':'ANCQGHILMFPSTWYV','3':'DE'}
# # #'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

# _SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'} #Orig
# # #'1'stand for Helix; '2'stand for Strand, '3' stand for coil

# _SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND','3':'MPSTHY'}
# # #'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate




# ##You can  add other properties of AADs to compute descriptors of protein sequence.

# _AATProperty=(_Hydrophobicity,_NormalizedVDWV,_Polarity,_Charge,_SecondaryStr,_SolventAccessibility,_Polarizability)

# _AATPropertyName=('_Hydrophobicity','_NormalizedVDWV','_Polarity','_Charge','_SecondaryStr','_SolventAccessibility','_Polarizability')



##################################################################################################
##################################################################################################

aa_ambiguous_alph = 'BJZX'
aa_unambiguous_alph = 'ARNDCEQGHILKMFPSTWYV'
aa_special_alph = 'UO'
aa_ter_alph = '*'
aa_alph = aa_unambiguous_alph +\
    aa_ambiguous_alph + aa_special_alph + aa_ter_alph

aa_subset_dict = {
    'aliphatic_hydrophobic': 'AVLIMPFW',
    'polar_uncharged': 'GSYNQC',
    'acidic': 'ED',
    'basic': 'KRH',
    'aliphatic': 'ILV',
    'aromatic': 'FYWH',
    'charged': 'HKRED',
    'polar': 'YWHKRDETCSNQ',
    'small': 'VCAGTPSDN',
    'tiny': 'AGCST',
    'helix': 'MALEK',
    'sheet': 'YFWTVI'}
aa_subsets = sorted(aa_subset_dict.keys())


# amino acids subdivided into three clusters per 7 properties as obtained
# from PROFEAT paper
# 2 adjustments...
AA_PROPERTIES_DIVISIONS = {
    'hydrophobicity': ['RKEDQN', 'GASTPHY', 'CLVIMFW'],
    'normvdw': ['GACSTPD', 'NVEQIL', 'MHKFRYW'],
    'polarity': ['LIFWCMVY', 'PATGS', 'HQRKNED'],
    'polarizability': ['GASDT', 'CPNVEQIL', 'KMHFRYW'],
    'charge': ['KR', 'ANCQGHILMFPSTWYV', 'DE'],
    'ss': ['EALMQKRH', 'VIYCWFT', 'GNPSD'],
    'sa': ['ALFCGIVW', 'PKQEND', 'MRSTHY']
}

#########################################
def letter_composition(seq, alph):
    '''
    This function returns the letter composition of seq for the letters in
    alph.

    Args:
        seq (str):
        alph (str):

    Raises:
        ValueError: if the sequence is empty.

    If seq contains only letters that are in alph, than the returned
    list of floats adds to one. Otherwise the sum of the numbers is between
    0.0 and 1.0

    >>> letter_composition('AABBCBBACB', 'ABC')
    array([ 0.3,  0.5,  0.2])
    >>> sum(letter_composition('AABBCBBACB', 'ABC'))
    1.0
    >>> letter_composition('AABBCBBACB', 'AB')
    array([ 0.3,  0.5])
    >>> letter_composition('AABBCBBACB', 'AD')
    array([ 0.3,  0. ])
    >>> letter_composition('AAAAAAAAAA', 'A')
    array([ 1.])
    >>> letter_composition('AAAAAAAAAA', '')
    array([], dtype=float64)
    '''

    if(len(seq) == 0):
        raise ValueError('Cannot calculate composition of empty sequence.')

    return letter_count(seq, alph) / float(len(seq))
def letter_count(seq, alph):
    '''
    This function counts letter occurances in seq for each letter in alph.

    Args:
        seq (str): The sequence of which the letters will be counted.
        alph (str): The letters that will be counted in seq

    Returns:
        numpy.array List with letter counts in the order of alph.

    >>> letter_count('AABBCBBACB', 'ABC')
    array([3, 5, 2])
    >>> letter_count('', 'ABC')
    array([0, 0, 0])
    >>> letter_count('ABC', '')
    array([], dtype=int64)
    '''
    return numpy.array([seq.count(l) for l in alph], dtype=int)


def property_ctd(seq, property):
    '''
    Based on the given property, the amino acid alphabet is subdivided into
    three groups. The amino acid sequence (seq) is mapped to this three-letter
    alphabet.

    The composition, transition, and distribution (ctd) of the mapped sequence
    is calculated and returned.

    property must be one of: 'hydrophobicity', 'normvdw', 'polarity',
    'polarizability', 'charge', 'ss', 'sa'.

    hyd
    vdw
    plr
    plz
    chr
    ss
    sa

    >>> s = 'ACACACACAC'
    >>> pctd = property_ctd(s, 'hydrophobicity')
    >>> pctd[:3]
    (0.0, 0.5, 0.5)
    >>> pctd[3:6]
    (0.0, 0.0, 1.0)
    >>> pctd[6:11]
    (0.0, 0.0, 0.0, 0.0, 0.0)
    >>> pctd[11:16]
    (0.1, 0.1, 0.5, 0.7, 0.9)
    >>> pctd[16:21]
    (0.2, 0.2, 0.6, 0.8, 1.0)

    >>> pctd = property_ctd('A', 'hydrophobicity')
    >>> pctd[:3]
    (0.0, 1.0, 0.0)
    >>> pctd[3:6]
    (0.0, 0.0, 0.0)
    >>> pctd[6:11]
    (0.0, 0.0, 0.0, 0.0, 0.0)
    >>> pctd[11:16]
    (1.0, 1.0, 1.0, 1.0, 1.0)
    >>> pctd[16:21]
    (0.0, 0.0, 0.0, 0.0, 0.0)

    >>> property_ctd('', 'hydrophobicity')
    Traceback (most recent call last):
     ...
    ValueError: Cannot calculate composition of empty sequence.
    '''

    # get mapping from amino acids to the three property clusters
    letter_mapping = property_division_mapping(property)

    # map aa protein sequence to property sequence (3-letter alphabet)
    state_seq = ''.join([letter_mapping[l] for l in seq])

    # composition features (letter counts normalized by sequence length)
    c0, c1, c2 = letter_composition(state_seq, 'ABC')

    # transition features (transition counts normalized by total number of
    # transitions)
    # TODO add separate transition count function

    # check if there is at least one transition, to avoid division by zero
    if(len(state_seq) < 2):
        t0 = 0.0
        t1 = 0.0
        t2 = 0.0
    else:
        seq_length = float(len(state_seq))
        t0 = (state_seq.count('AB') + state_seq.count('BA')) / (seq_length - 1)
        t1 = (state_seq.count('AC') + state_seq.count('CA')) / (seq_length - 1)
        t2 = (state_seq.count('BC') + state_seq.count('CB')) / (seq_length - 1)

    # distribution
    fractions = [0.25, 0.5, 0.75, 1.0]

    d0 = distribution(state_seq, 'A', fractions)
    d1 = distribution(state_seq, 'B', fractions)
    d2 = distribution(state_seq, 'C', fractions)

    return (c0, c1, c2,
            t0, t1, t2,
            d0[0], d0[1], d0[2], d0[3], d0[4],
            d1[0], d1[1], d1[2], d1[3], d1[4],
            d2[0], d2[1], d2[2], d2[3], d2[4])


def distribution(seq, letter, fractions=[0.25, 0.5, 0.75, 1.0]):
    '''
    This function returns at what fractions of the sequence, the given
    fractions of the letter are reached.

    >>> s = 'AABBABABABAABBAAAAAB'

    >>> distribution(s, 'B')
    [0.15, 0.2, 0.4, 0.65, 1.0]

    >>> distribution(s, 'C')
    [0.0, 0.0, 0.0, 0.0, 0.0]

    >>> distribution('', 'A')
    [0.0, 0.0, 0.0, 0.0, 0.0]
    '''

    #'Dan: First, last occurence ? '

    # count how often letter occurs in seq
    num_letter = seq.count(letter)

    if(num_letter == 0):
        return [0.0] * (len(fractions) + 1)

    # get letter indices where fraction of the letters is reached
    letter_positions = [max(1, int(round(f * num_letter))) for f in fractions]

    seq_fractions = []
    letter_count = 0
    for index, l in enumerate(seq):
        if(l == letter):
            letter_count += 1
            # the first occurance
            if letter_count == 1:
                seq_fractions.append((index + 1.0) / len(seq))
            # the fraction occurences
            if letter_count in letter_positions:

                for i in range(letter_positions.count(letter_count)):
                    seq_fractions.append((index + 1.0) / len(seq))

    return seq_fractions


def property_division_mapping(property, extra_letters=True):
    '''
    This function returns a mapping from amino acid to property 'index': A, B,
    or C. Other than unambiguous amino acids are mapped to D if extra_letters
    is set to True.
    '''

    default_letters = 'ABC'
    extra_letter = 'D'

    clusters = AA_PROPERTIES_DIVISIONS[property]
    assert(len(default_letters) == len(clusters))

    d = {}
    for letter, cluster in zip(default_letters, clusters):
        for aa in cluster:
            d[aa] = letter

    if(extra_letters):
        for aa in aa_ambiguous_alph + aa_special_alph + aa_ter_alph:
            d[aa] = extra_letter

    if(extra_letters):
        assert(sorted(d.keys()) == sorted(aa_alph))
    else:
        assert(sorted(d.keys()) == sorted(aa_unambiguous_alph))

    return d

##################################################################################################
##################################################################################################


##################################################################################################

def StringtoNum(ProteinSequence,AAProperty):
    """
    ###############################################################################################
    Tranform the protein sequence into the string form such as 32123223132121123.

    Usage:

    result=StringtoNum(protein,AAProperty)

    Input: protein is a pure protein sequence.

    AAProperty is a dict form containing classifiation of amino acids such as _Polarizability.

    Output: result is a string such as 123321222132111123222
    ###############################################################################################
    """

    hardProteinSequence=copy.deepcopy(ProteinSequence)
    for k,m in AAProperty.items():
        for index in m:
            hardProteinSequence=str.replace(hardProteinSequence,index,k)
    TProteinSequence=hardProteinSequence

    return TProteinSequence


def CalculateComposition(ProteinSequence,AAProperty,AAPName):
    """
    ###############################################################################################
    A method used for computing composition descriptors.

    Usage:

    result=CalculateComposition(protein,AAProperty,AAPName)

    Input: protein is a pure protein sequence.

    AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

    AAPName is a string used for indicating a AAP name.

    Output: result is a dict form containing composition descriptors based on the given property.
    ###############################################################################################
    """
    TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
    Result={}
    Num=len(TProteinSequence)
    Result[AAPName+'C'+'1']=round(float(TProteinSequence.count('1'))/Num,3)
    Result[AAPName+'C'+'2']=round(float(TProteinSequence.count('2'))/Num,3)
    Result[AAPName+'C'+'3']=round(float(TProteinSequence.count('3'))/Num,3)
    return Result

def CalculateTransition(ProteinSequence,AAProperty,AAPName):
    """
    ###############################################################################################
    A method used for computing transition descriptors

    Usage:

    result=CalculateTransition(protein,AAProperty,AAPName)

    Input:protein is a pure protein sequence.

    AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

    AAPName is a string used for indicating a AAP name.

    Output:result is a dict form containing transition descriptors based on the given property.
    ###############################################################################################
    """

    TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
    Result={}
    Num=len(TProteinSequence)
    CTD=TProteinSequence
    Result[AAPName+'T'+'12']=round(float(CTD.count('12')+CTD.count('21'))/(Num-1),3)
    Result[AAPName+'T'+'13']=round(float(CTD.count('13')+CTD.count('31'))/(Num-1),3)
    Result[AAPName+'T'+'23']=round(float(CTD.count('23')+CTD.count('32'))/(Num-1),3)
    return Result



def CalculateDistribution(ProteinSequence,AAProperty,AAPName):

    """
    ###############################################################################################
    A method used for computing distribution descriptors.

    Usage:

    result=CalculateDistribution(protein,AAProperty,AAPName)

    Input:protein is a pure protein sequence.

    AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.

    AAPName is a string used for indicating a AAP name.

    Output:result is a dict form containing Distribution descriptors based on the given property.
    ###############################################################################################
    """
    TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
    Result={}
    Num=len(TProteinSequence)
    temp=('1','2','3')
    for i in temp:
        num=TProteinSequence.count(i)
        ink=1
        indexk=0
        cds=[]
        while ink<=num:
            #indexk=string.find(TProteinSequence,i,indexk)+1 #Orig
            indexk=str.find(TProteinSequence,i,indexk)+1
            cds.append(indexk)
            ink=ink+1

        if cds==[]:
            Result[AAPName+'D'+i+'001']=0
            Result[AAPName+'D'+i+'025']=0
            Result[AAPName+'D'+i+'050']=0
            Result[AAPName+'D'+i+'075']=0
            Result[AAPName+'D'+i+'100']=0
        else:

            Result[AAPName+'D'+i+'001']=round(float(cds[0])/Num*100,3)
            Result[AAPName+'D'+i+'025']=round(float(cds[int(math.floor(num*0.25))-1])/Num*100,3)
            Result[AAPName+'D'+i+'050']=round(float(cds[int(math.floor(num*0.5))-1])/Num*100,3)
            Result[AAPName+'D'+i+'075']=round(float(cds[int(math.floor(num*0.75))-1])/Num*100,3)
            Result[AAPName+'D'+i+'100']=round(float(cds[-1])/Num*100,3)

    return Result

##################################################################################################


def CalculateC(ProteinSequence):
    '''
    ###############################################################################################
    Calculate all composition descriptors based seven different properties of AADs.
    Usage:

    result=CalculateC(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing all composition descriptors.
    ###############################################################################################
    '''
    result={}
    result.update(CalculateCompositionPolarizability(ProteinSequence))
    result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    result.update(CalculateCompositionCharge(ProteinSequence))
    result.update(CalculateCompositionPolarity(ProteinSequence))
    result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    return result

def CalculateT(ProteinSequence):
    """
    ###############################################################################################
    Calculate all transition descriptors based seven different properties of AADs.

    Usage:

    result=CalculateT(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing all transition descriptors.
    ###############################################################################################
    """
    result={}
    result.update(CalculateTransitionPolarizability(ProteinSequence))
    result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    result.update(CalculateTransitionCharge(ProteinSequence))
    result.update(CalculateTransitionPolarity(ProteinSequence))
    result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    return result

def CalculateD(ProteinSequence):
    """
    ###############################################################################################
    Calculate all distribution descriptors based seven different properties of AADs.

    Usage:

    result=CalculateD(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing all distribution descriptors.
    ###############################################################################################
    """
    result={}
    result.update(CalculateDistributionPolarizability(ProteinSequence))
    result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    result.update(CalculateDistributionCharge(ProteinSequence))
    result.update(CalculateDistributionPolarity(ProteinSequence))
    result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    return result


def CalculateCTD(ProteinSequence):
    """
    ###############################################################################################
    Calculate all CTD descriptors based seven different properties of AADs.

    Usage:

    result=CalculateCTD(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing all CTD descriptors.
    ###############################################################################################
    """
    result={}
    result.update(CalculateCompositionPolarizability(ProteinSequence))
    result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    result.update(CalculateCompositionCharge(ProteinSequence))
    result.update(CalculateCompositionPolarity(ProteinSequence))
    result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    result.update(CalculateTransitionPolarizability(ProteinSequence))
    result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    result.update(CalculateTransitionCharge(ProteinSequence))
    result.update(CalculateTransitionPolarity(ProteinSequence))
    result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    result.update(CalculateDistributionPolarizability(ProteinSequence))
    result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    result.update(CalculateDistributionCharge(ProteinSequence))
    result.update(CalculateDistributionPolarity(ProteinSequence))
    result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    return result
##################################################################################################

if __name__=="__main__":

#   import scipy,string

#   result=scipy.zeros((268,147))
#   f=file('protein1.txt','r')
#   for i,j in enumerate(f:
#       temp=CalculateCTD(string.strip(j))
#       result[i,:]=temp.values()
#   scipy.savetxt('ResultNCTRER.txt', result, fmt='%15.5f',delimiter='')
#
    protein="ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    pctd = property_ctd('ACACACACAC', 'hydrophobicity')
    print(pctd)
    print(property_ctd(protein,'ss'))
#   print StringtoNum(protein,_Hydrophobicity)
#   print CalculateComposition(protein,_Hydrophobicity,'_Hydrophobicity')
#   print CalculateTransition(protein,_Hydrophobicity,'_Hydrophobicity')
#   print CalculateDistribution(protein,_Hydrophobicity,'_Hydrophobicity')
#   print CalculateDistributionSolventAccessibility(protein)
#   print len(CalculateCTD(protein))
#   print len(CalculateC(protein))
#   print len(CalculateT(protein))
#   print len(CalculateD(protein))
    # print (CalculateCTD(protein))


