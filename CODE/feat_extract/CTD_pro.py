# -*- coding: utf-8 -*-
'''
#####################################################################################
 - ProtPy version.
This module is used for computing the composition, transition and distribution

descriptors based on the different properties of AADs. The AADs with the same

properties are marked as the same number.
You can get 147 descriptors from the classic descriptors.

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

#####################################################################################
'''
from collections import defaultdict
import math
import copy

'TODO: Implement a better, partial getting of different C/T/Ds  from CTD Pro..'
'TODO: replace StringtoNum with something decent, then update rest of methods'

'TODO: Merge ROFEAT/ProtPy propensities with/into AALphabets part (contains 3 letter alphs)!'
AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]


### PROFEAT/ProtPy propensity based scales: ####
' TODO: take from AAlphabets translations! - For Consistency! (Requires backmapping keys to NUM)'
'Horribly kludgy - Would be better to convert to a dict?:'
_DisorderPropensity={'1':'ARSQEGKP','2':'ILNCFYVW', '3':'DHMT'}
# #Orig ProtPy scales - seem buggy! Letters  were lacking.. Some based on Profeat. I fixed manually
_Hydrophobicity={'1':'RKEDQN','2':'GASTPHY','3':'CLVIMFW'}
# #'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

_Polarity={'1':'LIFWCMVY','2':'PATGS','3':'HQRKNED'} #ProFeat based
# _Polarity={'1':'LIFWCMVY','2':'CPNVEQIL','3':'KMHFRYW'} # Orig - diff. (Extra V?)
# #'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

_Polarizability={'1':'GASDT','2':'CPNVEQIL','3':'KMHFRYW'}
# #'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)

_Charge={'1':'KR','2':'ANCQGHILMFPSTWYV','3':'DE'}
# #'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

_SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'} #Orig
# #'1'stand for Helix; '2'stand for Strand, '3' stand for coil

_NormalizedVDWV={'1':'GASTPDC','2':'NVEQIL','3':'MHKFRYW'}
# #'1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

_SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND','3':'MPSTHY'}
# #'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate


'Todo: Zip these two into a nice list or sets of tuples for pitys sake -'

##You can add other properties of AADs to compute descriptors of protein sequence.
#_AATProperties
AAG_Properties=[_Hydrophobicity,_NormalizedVDWV,
_Polarity,_Charge,_SecondaryStr,
_SolventAccessibility,_Polarizability,
_DisorderPropensity]

# AAG_Names=['_Hydrophobicity','_NormalizedVDWV',
# '_Polarity','_Charge',
# '_SecondaryStr','_SolventAccessibility',
# '_Polarizability','_DisorderPropensity']
AAG_Names=['Hydrophobicity','Normalized VDWV',
'Polarity','Charge',
'Secondary Str','Solvent Accessibility',
'Polarizability','Disorder Propensity']

# AAG_groups = zip(AAG_Names,AAG_Properties) #Make a list of tuples

# ag=[]
# for a in AAG_Properties:
#   ag.append((a,str(a)))
# print(ag)
##################################################################################################

'TODO: See to implementing something efficient - translator from Aalphabet.py, or Str.maketrans , or RegEx..'
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

    # # reversed_AA_dict = dict (zip(AAProperty.values(),AAProperty.keys()))
    # reversed_AA_dict= dict(((value, key) for key, value in AAProperty.items()))
    # print(reversed_AA_dict)

    # TProteinSequence = translate_sequence(ProteinSequence,reversed_AA_dict)

    hardProteinSequence=copy.deepcopy(ProteinSequence)
    for k,m in AAProperty.items():
      # print ("AAProperty: %s" %(AAProperty))
      # print ("AAProperty.items(): %s" %(AAProperty.items()))
      # print ("k,m = %s,%s, in AAProperty.items() " % (k,m))
      for index in m:
          # print ("index in m: %s" %(index))
          hardProteinSequence=str.replace(hardProteinSequence,index,k)
    return hardProteinSequence
    # TProteinSequence=hardProteinSequence
    # return TProteinSequence



def CalculateCTD(ProteinSequence,ctd_call='CTD'):
    """
    ###############################################################################################
    Calculate all CTD descriptors based on all saved AAG_Properties.


    Usage:

    result=CalculateCTD(protein,CTD)
    composition_results =CalculateCTD(protein,C)

    Input:ProteinSequence is a pure sequence.
    ctd = String of which properties (C,T,D) should be calculated and returned.

    Output:result is a dict containing all CTD descriptors.
    ###############################################################################################
    """
    # result={}
    result=defaultdict(float)
    #May want to change this later into a dict, or a more elegant format..

    #Check if only part of the properties were requested:
    get_C = True
    get_T = True
    get_D = True
    ctd_call = ctd_call.lower()
    if 'c' not in ctd_call:
        get_C=False
    if 't' not in ctd_call:
        get_T=False
    if 'd' not in ctd_call:
        get_D=False

    for i in range(len(AAG_Names)):
        AAProperty=AAG_Properties[i]
        AAPName = AAG_Names[i]
        if get_C:
            result.update(CalculateComposition(ProteinSequence, AAProperty, AAPName))
        if get_T:
            result.update(CalculateTransition(ProteinSequence, AAProperty, AAPName))
        if get_D:
            result.update(CalculateDistribution(ProteinSequence, AAProperty, AAPName))
    return result


def CalculateComposition(ProteinSequence,AAProperty,AAPName):
    """
    ###############################################################################################
    A method used for computing composition descriptors.

    Usage:

    result=CalculateComposition(protein,AAProperty,AAPName)

    Input: protein is a pure protein sequence.

    AAProperty is a dict form containing classification of amino acids such as _Polarizability.

    AAPName is a string used for indicating a AAP name.

    Output: result is a dict form containing composition descriptors based on the given property.
    ###############################################################################################
    """
    TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
    Result={}
    Num=len(TProteinSequence)
    Result[AAPName +' Composition:'+'1']=round(float(TProteinSequence.count('1'))/Num,3)
    Result[AAPName +' Composition:'+'2']=round(float(TProteinSequence.count('2'))/Num,3)
    Result[AAPName +' Composition:'+'3']=round(float(TProteinSequence.count('3'))/Num,3)
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
    # Result={}
    Result=defaultdict(float)
    Num=len(TProteinSequence)
    CTD=TProteinSequence
    Result[AAPName +' Transitions:'+'12']=round(float(CTD.count('12')+CTD.count('21'))/(Num-1),3)
    Result[AAPName +' Transitions:'+'13']=round(float(CTD.count('13')+CTD.count('31'))/(Num-1),3)
    Result[AAPName +' Transitions:'+'23']=round(float(CTD.count('23')+CTD.count('32'))/(Num-1),3)
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
    # Result={}
    Result=defaultdict(float)
    Num=len(TProteinSequence)
    temp=('1','2','3')
    for i in temp:
        num=TProteinSequence.count(i)
        ink=1
        indexk=0
        cds=[]
        while ink<=num:
            indexk=str.find(TProteinSequence,i,indexk)+1
            cds.append(indexk)
            ink=ink+1

        if cds==[]:
            Result[AAPName +' Distribution'+i+'001']=0
            Result[AAPName +' Distribution'+i+'025']=0
            Result[AAPName +' Distribution'+i+'050']=0
            Result[AAPName +' Distribution'+i+'075']=0
            Result[AAPName +' Distribution'+i+'100']=0
        else:

            Result[AAPName +' Distribution'+i+'001']=round(float(cds[0])/Num*100,3)
            Result[AAPName +' Distribution'+i+'025']=round(float(cds[int(math.floor(num*0.25))-1])/Num*100,3)
            Result[AAPName +' Distribution'+i+'050']=round(float(cds[int(math.floor(num*0.5))-1])/Num*100,3)
            Result[AAPName +' Distribution'+i+'075']=round(float(cds[int(math.floor(num*0.75))-1])/Num*100,3)
            Result[AAPName +' Distribution'+i+'100']=round(float(cds[-1])/Num*100,3)

    return Result


#####################################################################################



##################################################################################################
def CalculateCompositionHydrophobicity(ProteinSequence):

    """
    ###############################################################################################
    A method used for calculating composition descriptors based on Hydrophobicity of

    AADs.

    Usage:

    result=CalculateCompositionHydrophobicity(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Composition descriptors based on Hydrophobicity.
    ###############################################################################################
    """

    result=CalculateComposition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result

def CalculateCompositionNormalizedVDWV(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating composition descriptors based on NormalizedVDWV of

    AADs.

    Usage:

    result=CalculateCompositionNormalizedVDWV(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Composition descriptors based on NormalizedVDWV.
    ###############################################################################################
    """
    result=CalculateComposition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
    return result

def CalculateCompositionPolarity(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating composition descriptors based on Polarity of

    AADs.

    Usage:

    result=CalculateCompositionPolarity(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Composition descriptors based on Polarity.
    ###############################################################################################
    """

    result=CalculateComposition(ProteinSequence,_Polarity,'_Polarity')
    return result

def CalculateCompositionCharge(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating composition descriptors based on Charge of

    AADs.

    Usage:

    result=CalculateCompositionCharge(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Composition descriptors based on Charge.
    ###############################################################################################
    """

    result=CalculateComposition(ProteinSequence,_Charge,'_Charge')
    return result

def CalculateCompositionSecondaryStr(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating composition descriptors based on SecondaryStr of

    AADs.

    Usage:

    result=CalculateCompositionSecondaryStr(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Composition descriptors based on SecondaryStr.
    ###############################################################################################
    """

    result=CalculateComposition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
    return result

def CalculateCompositionSolventAccessibility(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating composition descriptors based on SolventAccessibility

    of  AADs.

    Usage:

    result=CalculateCompositionSolventAccessibility(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Composition descriptors based on SolventAccessibility.
    ###############################################################################################
    """

    result=CalculateComposition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
    return result

def CalculateCompositionPolarizability(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating composition descriptors based on Polarizability of

    AADs.

    Usage:

    result=CalculateCompositionPolarizability(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Composition descriptors based on Polarizability.
    ###############################################################################################
    """

    result=CalculateComposition(ProteinSequence,_Polarizability,'_Polarizability')
    return result

##################################################################################################


##################################################################################################
def CalculateTransitionHydrophobicity(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Transition descriptors based on Hydrophobicity of

    AADs.

    Usage:

    result=CalculateTransitionHydrophobicity(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Transition descriptors based on Hydrophobicity.
    ###############################################################################################
    """

    result=CalculateTransition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result

def CalculateTransitionNormalizedVDWV(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Transition descriptors based on NormalizedVDWV of

    AADs.

    Usage:

    result=CalculateTransitionNormalizedVDWV(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Transition descriptors based on NormalizedVDWV.
    ###############################################################################################
    """

    result=CalculateTransition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
    return result

def CalculateTransitionPolarity(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Transition descriptors based on Polarity of

    AADs.

    Usage:

    result=CalculateTransitionPolarity(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Transition descriptors based on Polarity.
    ###############################################################################################
    """

    result=CalculateTransition(ProteinSequence,_Polarity,'_Polarity')
    return result

def CalculateTransitionCharge(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Transition descriptors based on Charge of

    AADs.

    Usage:

    result=CalculateTransitionCharge(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Transition descriptors based on Charge.
    ###############################################################################################
    """

    result=CalculateTransition(ProteinSequence,_Charge,'_Charge')
    return result

def CalculateTransitionSecondaryStr(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Transition descriptors based on SecondaryStr of

    AADs.

    Usage:

    result=CalculateTransitionSecondaryStr(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Transition descriptors based on SecondaryStr.
    ###############################################################################################
    """

    result=CalculateTransition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
    return result

def CalculateTransitionSolventAccessibility(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Transition descriptors based on SolventAccessibility

    of  AADs.

    Usage:

    result=CalculateTransitionSolventAccessibility(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Transition descriptors based on SolventAccessibility.
    ###############################################################################################
    """

    result=CalculateTransition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
    return result

def CalculateTransitionPolarizability(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Transition descriptors based on Polarizability of

    AADs.

    Usage:

    result=CalculateTransitionPolarizability(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Transition descriptors based on Polarizability.
    ###############################################################################################
    """

    result=CalculateTransition(ProteinSequence,_Polarizability,'_Polarizability')
    return result

##################################################################################################
##################################################################################################
def CalculateDistributionHydrophobicity(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Distribution descriptors based on Hydrophobicity of

    AADs.

    Usage:

    result=CalculateDistributionHydrophobicity(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Distribution descriptors based on Hydrophobicity.
    ###############################################################################################
    """

    result=CalculateDistribution(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result

def CalculateDistributionNormalizedVDWV(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Distribution descriptors based on NormalizedVDWV of

    AADs.

    Usage:

    result=CalculateDistributionNormalizedVDWV(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Distribution descriptors based on NormalizedVDWV.
    ###############################################################################################
    """

    result=CalculateDistribution(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
    return result

def CalculateDistributionPolarity(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Distribution descriptors based on Polarity of

    AADs.

    Usage:

    result=CalculateDistributionPolarity(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Distribution descriptors based on Polarity.
    ###############################################################################################
    """

    result=CalculateDistribution(ProteinSequence,_Polarity,'_Polarity')
    return result

def CalculateDistributionCharge(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Distribution descriptors based on Charge of

    AADs.

    Usage:

    result=CalculateDistributionCharge(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Distribution descriptors based on Charge.
    ###############################################################################################
    """

    result=CalculateDistribution(ProteinSequence,_Charge,'_Charge')
    return result

def CalculateDistributionSecondaryStr(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Distribution descriptors based on SecondaryStr of

    AADs.

    Usage:

    result=CalculateDistributionSecondaryStr(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Distribution descriptors based on SecondaryStr.
    ###############################################################################################
    """

    result=CalculateDistribution(ProteinSequence,_SecondaryStr,'_SecondaryStr')
    return result

def CalculateDistributionSolventAccessibility(ProteinSequence):

    """
    ###############################################################################################
    A method used for calculating Distribution descriptors based on SolventAccessibility

    of  AADs.

    Usage:

    result=CalculateDistributionSolventAccessibility(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Distribution descriptors based on SolventAccessibility.
    ###############################################################################################
    """

    result=CalculateDistribution(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
    return result

def CalculateDistributionPolarizability(ProteinSequence):
    """
    ###############################################################################################
    A method used for calculating Distribution descriptors based on Polarizability of

    AADs.

    Usage:

    result=CalculateDistributionPolarizability(protein)

    Input:protein is a pure protein sequence.

    Output:result is a dict form containing Distribution descriptors based on Polarizability.
    ###############################################################################################
    """

    result=CalculateDistribution(ProteinSequence,_Polarizability,'_Polarizability')
    return result

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
    for i in range(len(AAG_Names)):
        AAProperty=AAG_Properties[i]
        AAPName = AAG_Names[i]
        result.update(CalculateComposition(ProteinSequence, AAProperty, AAPName))
    return result
    # result.update(CalculateCompositionPolarizability(ProteinSequence))
    # result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    # result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    # result.update(CalculateCompositionCharge(ProteinSequence))
    # result.update(CalculateCompositionPolarity(ProteinSequence))
    # result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    # result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    # return result

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
    for i in range(len(AAG_Names)):
        AAProperty=AAG_Properties[i]
        AAPName = AAG_Names[i]
        result.update(CalculateTransition(ProteinSequence, AAProperty, AAPName))

    return result
    # result.update(CalculateTransitionPolarizability(ProteinSequence))
    # result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    # result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    # result.update(CalculateTransitionCharge(ProteinSequence))
    # result.update(CalculateTransitionPolarity(ProteinSequence))
    # result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    # result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    # return result

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
    for i in range(len(AAG_Names)):
        AAProperty=AAG_Properties[i]
        AAPName = AAG_Names[i]
        result.update(CalculateDistribution(ProteinSequence, AAProperty, AAPName))
    return result
    # result.update(CalculateDistributionPolarizability(ProteinSequence))
    # result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    # result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    # result.update(CalculateDistributionCharge(ProteinSequence))
    # result.update(CalculateDistributionPolarity(ProteinSequence))
    # result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    # result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    # return result



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
    protein="ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    print(StringtoNum(protein,_Hydrophobicity))
    print(CalculateComposition(protein,_Hydrophobicity,'_Hydrophobicity'))
    # print(CalculateTransition(protein,_Hydrophobicity,'_Hydrophobicity'))
    # print CalculateDistribution(protein,_Hydrophobicity,'_Hydrophobicity')
    print (len(CalculateCTD(protein)))
    # print len(CalculateT(protein))
    # print len(CalculateD(protein))
    # print (CalculateCTD(protein))
    # print (CalculateCTD(protein,'c'))
    # print (CalculateCTD(protein,'Td'))


