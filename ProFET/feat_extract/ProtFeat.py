# -*- coding: utf-8 -*-
"""
A class used for computing different types of protein descriptors!

##TODO:
#Autocorrellation not quite working properly. (Also, add method to give mult. potent. letters)
# PI, PH..
# Entropy (both) - Sequence, letter - maybe wrong. (check that values returned is really log..)

# FUTURE IMPROVEMENT: Replace Dict datastructure storage of feature names:values
 with Numpy/PANDAS Series?:
    Allows Vectorization of mathematical operations. (Much FASTER! Add vectorization where possible).
    More elegant.
    Easy concatenation.
    Saves work when merging into a pd.dataframe later. (I.E , "GetSimpleFeatures")
# Future:
    Have called methods also update a self.value (store features of a protein?).
    (I.e, protfeat.f.getAProp() => f.self.AProp = {prop value}.
    - This makes more sense from a OO perspective, and is useful
     when performing multiple calls on the same sequence (vs just one call)

 Issues:
 PI, PH, PKA - using n/c terminal. (Needs fixing)
"""
# import AAlphabets
import time
from AAlphabets import *
from AAScales import *
from collections import Counter, defaultdict
from itertools import product
import CTD_pro
# from AAComposition import CalculateAAComposition, CalculateDipeptideComposition, GetSpectrumDict
import numpy as np
from Bio.SeqUtils import ProtParam as pp
import Bio.SeqUtils.ProtParamData
import re
from math import log
from scipy import stats
# from GetSubSeq import GetSubSequence

from Disorder import getDisordered

"""
AA20 = 'ACDEFGHIKLMNPQRSTVWY' #standard 20 AA Alphabet"
#'AA_14_List = Ofer_14 - custom list derived from SDM12. 14^2 = 196'
"""
# print(REDUCED_ALPHABETS_LETTERS)
# AA20=REDUCED_ALPHABETS_LETTERS['AA20']

"ofer14=REDUCED_ALPHABETS_LETTERS['ofer14'] #Be careful not to overload variables from alph!!"

# print(AA20)
# print(REDUCED_ALPHABETS_LETTERS['AA20'])

'Used to calc protein net charge'
'TODO: Modify PI, Netcharge calcs, to use mofidied deep copuies of these  when no N or C termini'
pKa     = {'D':3.9, 'E':4.3, 'H':6.1, 'C':8.3, 'Y':10.1, 'K':10.5, 'R':12, 'N-term':8, 'C-term':3.1}
charges = {'D':-1,  'E':-1,  'H':+1,  'C':-1,  'Y':-1,   'K':1,    'R':1,  'N-term':1, 'C-term':-1}


'We will want to add try&pickling for this! (OR for more classes):'
'TODO: OR - Memoize. (If carried across classes/objects?)'

def NGram_combinations(alphabet='AA20',k=2): # (alphabet=(str(AA20)),k=2):
    '''
    This function returns all possible ordered letter pairs for the provided
    alphabet alph, (I.E possible Bigrams), of length k.

    This gets "Permutations with Replacement" using Python's itertools.product.
    (It gets ALL the combs, and later this dict retains also "0" counts.. )

    TODO: Pickle or Memoize n-grams so as not to recalcualte for every call!

    alph should be a non-empty string without duplicate characters. (The letters
        of an alphabet)

    ******
    NOTE:
    This ENTIRE Step could be skipped IF features are PRESERVED from the
    training set, and TEST Set features are selected/filtered to include ONLY
    those features PRESENT in the TRAINING SET!
    (IF so - then one can simply count only those N-Grams found in the initial
    training set of inputs, without needing to get all possible combs here)
    ******
    Args:
        alph (str): The letter alphabet
        k (int) : number of repeats & output n-grams length
    Returns:
        [str] List with length-k strings.
    >>> NGram_combinations('AB')
    ['AA', 'AB', 'BA', 'BB']
    '''

    alph = REDUCED_ALPHABETS_LETTERS[str(alphabet)]
    # print('alph',str(alph))
    def get_combs(alph,k):
        bi_comb = [''.join(x) for x in product(str(alph), repeat=int(k))]
        return bi_comb

    bi_comb=get_combs(alph,k)
    CombinationsDict = dict.fromkeys(bi_comb,0)
    # print(CombinationsDict)
    return CombinationsDict
# comb_generator = (''.join(x) for x in product('AB', repeat=2))

'Store some common combinations of reduced alphabets, save on computation/Better to MEMOIZE!'
bigrams_20 = NGram_combinations('AA20',2)
bigrams_14 = NGram_combinations('ofer14',2)
bigrams_sdm_12 = NGram_combinations("sdm12",2)
grams_hp2_5 = NGram_combinations("hp2",5)

"NOTE: We could do a list here using the 3-letter alphs names, more neat"
# THREE_LETTER_ALPH_NAMES
trigrams_disorder_3 = NGram_combinations("Disorder_3",3)
trigrams_charge_3 = NGram_combinations("Charge_3",3)
trigrams_hp3_3 = NGram_combinations("Hydrophobicity_3",3)
trigrams_polarity_3 = NGram_combinations("Polarity_3",3)
trigrams_polarize_3 = NGram_combinations("Polarizability_3",3)
trigrams_SS_3 = NGram_combinations("SecondaryStr_3",3)
trigrams_NormVDWV_3 = NGram_combinations("NormVDWV_3",3)
trigrams_SolventA_3 = NGram_combinations("SolventA_3",3)

# trigrams_4 = NGram_combinations('gbm4',3)
# print(bigrams_20['MQ'])
# print('bigrams_14:',str(bigrams_14))



def autocorrelation_mb(sequence, scale, lag):
    '''
    Code copied from SPICE:
    "SPiCE: a web-based tool for sequence-based protein classification and exploration
    Bastiaan A van den Berg*, Marcel JT Reinders, Johannes A Roubos and Dick de Ridder"
    https://github.com/basvandenberg/spice
    ---------------------

    This function uses the provided scale to transform the sequence seq into a
    signal and returns the autocorrelation of this signal for the given lag.
    Normalized Moreau-Broto autocorrelation as given in Li (2006) PROFEAT.
    TODO formula
    This autocorrelation evaluates to 0.0 if lag is equal to larger than the
    sequence length.
    >>> autocorrelation_mb('ABACABAC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    0.5
    >>> autocorrelation_mb('ABACABAC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 2)
    -0.5
    >>> autocorrelation_mb('BBBBCCCC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    -1.0
    >>> autocorrelation_mb('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    1.0
    >>> autocorrelation_mb('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 8)
    0.0
    >>> autocorrelation_mb('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 12)
    0.0
    '''

    if(lag < 1):
        raise ValueError('The provided lag should be a positive integer.')

    # TODO check if this is the correct approach
    if(len(sequence) <= lag):
        return 0.0

    # transform sequence to signal using the provided scale
    signal = numpy.array(seq_signal_raw(sequence, scale))

    # calculate autocorrelation
    autocorr = sum(signal[:-lag] * signal[lag:])

    # return normalized autocorrelation
    return autocorr / float(len(sequence) - lag)




class ProtFeat:
    """
    The GetProDes/ProtFeat class collects most sequence features  calculation modules in a class.
    Different features  are returned in dictionary-(like) forms.

    """

    def __init__(self,ProteinSequence='',alph_letters='AA20',HAS_N=True,HAS_C=True):
        """
        Input a protein sequence.

        ProteinSequence = string, protein's sequence.
        (alph) - alph_letters = string of all the amino acids in the alphabet.
         Default is 20 (normal).Required when getting N-grams/K-mer counts,
         so as to count consistently.

         Feature keys/Values must be consistent in training and test! (Or MADE so post-hoc)

         HAS_N = "Does sequence have N_terminal end" (Or is it a 'subsequence').
         HAS_N = "Does sequence have C_terminal end".

         We may want to have called features update the class instance,
         when the feature is "unchanging" ,i.e  AAC, PI.., in addition to
         returning the feature and its values as a dict.


        """
        if len(ProteinSequence)==0:
            print ("You must input a protein sequence as a string when constructing!")
        else:
            self.seq=ProteinSequence
            'alph_letters = Total letters in biological alphabet used'
            self.alph_set = alph_letters #AA20, ofer14, etc'
            'self.alph on the other hand, holds the letters defined by the alphabet set'
            if alph_letters=='AA20':
                self.IsStandardAlph = True
            else:
                self.IsStandardAlph = False
            'Get Alphabet letters from AAlph'
            # self.alph = alph_letters #default is the 20 letter AA"
            self.alph=REDUCED_ALPHABETS_LETTERS[alph_letters]
            # print ('self.alph = alphabet',str(self.alph))
            # print ('alph_set ',str(self.alph_set))
            self.HAS_N = HAS_N #Does seq have "free" N-terminal tail end.
            self.HAS_C = HAS_C
            # self.length = float(len(ProteinSequence))
            self.length = len(ProteinSequence)

            """
            def Get_AA_Counts(self,AA_Types=AA20):
                '''
                Gets the total counts of letters/symbols in the sequence/self.

                Defaults to self.alph, but can be supplied with another alphabet,
                for example if one wishes to count specific groups of AA.
                '''
                letter_counts = defaultdict(int)
                for letter in AA_Types:
                    letter_counts[letter] = AA_Counter[letter]
                return letter_counts
            # self.AA_Counts=Get_AA_Counts(alph)
            """
            self.AA_Counts=Counter(ProteinSequence)

            self.Bio_ProtParam = pp.ProteinAnalysis(self.seq)

    'Remove one of the two (update calls).'
    ' Done: Updated aa counting for when alphabet is reduced..  '
    # @jit
    def GetAA_Freq(self):
        """
        amino acid composition (20).
        Usage:

        result = GetAA_Freq()
        """

        res = defaultdict(float)
        counts=self.AA_Counts
        for a in counts:
            _aa_ = str('AA: '+a+' Frequency')  #MOD Dan 12.1.2015
            res[_aa_]=round((float(self.AA_Counts[a]) / self.length) * 100, 3)

        return self.alphabet_prefix(res)

    def countToFreq(self,dict_counts):
        '''
        Gets percentage from a provided dict of AA counts.
        Could be optimized!.. (Use a list..)
        '''
        length = self.length
        # Result = defaultdict(float)
        Result = Counter()
        # print(dict_counts)
        for k,v in dict_counts.items():
            Result[k]=round(v / (length - 1) * 100, 2)
        return Result


    'Makes Seq of bigram chunks:'
    'Play with  this later to make discrete chunks of majority 2-3 mers, in base 3'
    def ordered_seq_pairs(self, distance=1):
        '''
        Returns sequence as pairs of bigrams. (Useful for counting).
        This function returns (sequentially ordered) letter pairs from seq that
        occur at the given distance from each other.

        The distance should be an integer larger then 0.

        An empty list will be returned if distance is equal to or larger than the
        sequence length.

        Args:
            seq (str): The sequence of which the pairs are returned
            distance (int): distance between the letter pairs
        Returns:
            [(str, str)] List of tuples, each tuple containing two letters of one
                         pair. The lengt of the list is max(0, len(seq) - 1).

        >>> ordered_seq_pairs('ABCDE', 1)
        ['AB', 'BC', 'CD', 'DE']
        >>> ordered_seq_pairs('ABCDE', 3)
        ['AD', 'BE']
        >>> ordered_seq_pairs('ABCDE', 5)
        []
        '''
        seq=self.seq
        if(distance <= 0):
            raise ValueError('Distance should be larger than 0.')

        return [''.join(t) for t in zip(seq[:-distance], seq[distance:])]
    'TODO - Implement this with a DECORATOR, for all the features returned!!'
    def alphabet_prefix (self,res):
        '''
        if non standard/reduced alphabet is used,
        then modify feature names (keys in given dictionary)
         of results to have appropiate prefix.

         res = dict of features. Keys = feature names (strings), values: features val.
        '''
        if self.IsStandardAlph == False:
            return (self.Dict_Keys_prefix(res,str(self.alph_set)))
        else:
            return res

    def MerCount(self,kmer_counts,k) :
        '''
        Counts all occurences of given k-mers (From the
            kmers dictionary's keys)
        in the sequence.
        Stores counts as values.
        '''
        s=self.seq
        length=int(self.length)
        # print(kmer_counts)
        for i in range(length - (k-1)) :
            kmer_counts[s[i :i + k]] += 1
        return kmer_counts

    """
    def GetBigramCount(self,k=2):

        'Currently using mercoutn code - also possible to order seq as bigrams and count over them..-'
        # 'BIGSpeed up possible!  - vectorizing;NumPy would give big speedup. But names/keys needed..'
        # # print("Bigram_dict-Bigram count")
        # # print(bigram_dict)
        #     # count pairs while walking over the sequence
        # for seq_pair in self.ordered_seq_pairs():
        #     # if(seq_pair in pair_counts.keys()):
        #     bigram_dict[seq_pair] += 1
        # 'BIG Possible optimization here - Numpy Array - see "diletter_count" sequtil- Spice..'

        # return bigram_dict
        return self.MerCount(bigram_dict,2)
    """

    'BIG potential speed up - use a 2d numpy recarray maybe, instead of dict for counting'
    'Also - storage vs calcing of different, common NGram_combinations'
    'TODO: If MEMOIZE implemented, then the "IF checks" for saved alph may not be needed'
    'TODO: IgnoreAlphabetType not implemented. Examine down the pipe "mercount" to see if needed.'
    'TODO: Implement UseVariableCombinations / varied counts. (At least minimally, and assuming pandas later). AND For KGram Freq, KMirrors..'
    def GetkgramCounts(self,k=2,UseVariableCombinations=False,IgnoreAlphabetType=False):
        """
        k-mer peptide composition counts (alphSize^k)

        If UseVariableCombinations=False, then calls on a dictionary
        of ALL possible k_mer combinations using def NGram_combinations,
        (In order to ensure feature count consistency).
        This could be skipped, with CARE; (If feature matrix is later kept 'fixed')

        Usage:

        result = GetkmerFreq()

        NOTE: Ideally - n-gram multiLetterComb output dictionary should be saved/pickled,
        rather than requiring that it be recalcualted.
        By default, length 2 bigram combinations are available for some alphabets.
        """
        # res=GetSpectrumDict(self.ProteinSequence)
        if k <= 0:
            print ("Please enter a positive k !")
            return
        # Don't count zero counts... Not implemented yet
        if UseVariableCombinations==True:
            # print ('Using VariableCombinations')
            kmers_dict = Counter()
        else:
            'Get stored dicts if available'
            if (self.alph_set == 'AA20') & (k==2):
                kmers_dict = bigrams_20.copy()
            elif  (self.alph_set == 'ofer14') & (k==2): #- re.sub alph first..
                kmers_dict = bigrams_14.copy() #Make a DEEP copy of the dict
            elif (self.alph_set == 'Disorder_3') & (k==3):
                kmers_dict = trigrams_disorder_3.copy()
            elif (self.alph_set == 'Charge_3') & (k==3):
                kmers_dict = trigrams_charge_3.copy()
            elif (self.alph_set == 'Hydrophobicity_3') & (k==3):
                kmers_dict = trigrams_hp3_3.copy()
            elif (self.alph_set == 'Polarity_3') & (k==3):
                kmers_dict = trigrams_polarity_3.copy()
            elif (self.alph_set == 'Polarizability_3') & (k==3):
                kmers_dict = trigrams_polarize_3.copy()
            elif (self.alph_set == 'SecondaryStr_3') & (k==3):
                kmers_dict = trigrams_SS_3.copy()
            elif (self.alph_set == 'NormVDWV_3') & (k==3):
                kmers_dict = trigrams_NormVDWV_3.copy()
            elif (self.alph_set == 'SolventA_3') & (k==3):
                kmers_dict = trigrams_SolventA_3.copy()
            elif (self.alph_set == 'hp2') & (k==5):
                kmers_dict = grams_hp2_5.copy()
            else: #Non stored AA.alph -> Gen. possible Bigrams using method NGram_combinations
                # kmers_dict = NGram_combinations (self.alph,k)
                kmers_dict = NGram_combinations (self.alph_set,k)

        # print(kmers_dict)
        'This step(NGram_combinations) Could potentially be skipped with care! (See note)'
        return self.MerCount(kmers_dict,k)

    def GetkgramFreq(self,k,IgnoreAlphabetType=False):
        '''
        Get K-mer frequency/composition.

        If IgnoreAlphabetType=True, (IS False by default) then don't add
        a prefix of the type of reduced aa alph to the feature keys.
        (Useful when further manipulating by key/string )
        '''
        res = self.countToFreq(self.GetkgramCounts(k))
        if IgnoreAlphabetType==True:
            return res
        else:
            return self.alphabet_prefix(res)

        # if self.IsStandardAlph == False:
        #     return (self.Dict_Keys_prefix(res,str(self.alph))
        # else:
        #     return res


    'TODO! - Ensure that Working currently'
    'BUG: Not working properly curretnly with reduced AA alphs due to prefixes..'
    'TODO: Add Get Kmirror COUNTS, not just freq'
    def GetKMirrorsFreq (self,k=2,AddMirrorsPrefix=True,getFreq=True):
        '''
        Get kmer mirror frequencies: if we assume mirror images of letters to be identical.
        Merges "mirror images" and sums their values into the remaining 'key'/image.
        (Palindromes are unaffected)
        I.E: we treat '01' == '10'. but, "11","00",  not affected. But 'KR' and 'RK' would be.

        This gives us a number of combinations acc. to Combination with Replacement, vs "Permutations with Replacement".
        https://docs.python.org/3.3/library/itertools.html#itertools.combinations_with_replacement
        http://www.calculatorsoup.com/calculators/discretemathematics/combinationsreplacement.php

        So: (Pseudocode):   #Function does not currently support passing it kgramFreq directly.
            original_kfreq=self.GetkgramFreq(k=2)
            print(original_kfreq)
            {'01':0.2,'10':0.7,'00':0.0,'11':0.1}
            k_mirrors = GetKgramMirrorFreq(original_kfreq)
            print(k_mirrors)
            {'10':6,'00':0.0,'11':0.1}

        For k=4, eg: 0111 would me merged with 1110..

        Can be configured to get counts rather than frequencies. (getFreq)

        '''
        # if k != 2: #Or %2 . .
        #     print ("Warning, mirror reduction is meant for 2-mers")

        #Be careful here not to mess up the "original" k_freq, if looking for efficiency..
        if getFreq==True:
            k_dict = self.GetkgramFreq(k=k,IgnoreAlphabetType=True)
        else:
            k_dict = self.GetkgramCounts(k=k,IgnoreAlphabetType=True)
        mirrored_k = self.KMirrors(k_dict)
        if AddMirrorsPrefix==True:
             mirrored_k=self.Dict_Keys_prefix(mirrored_k,' Mirror K-mer ')
        return self.alphabet_prefix(mirrored_k)

    def KMirrors(self,k_dict):
        '''
        Given a kmer (k_freq?) dict, returns it with mirror images merged and summed.
        (Palindromes unaffected).
         = Implementation part of GetKMirrorsFreq
        '''
        # d=Counter()
        d = defaultdict(float)
        # k_dict.items()
        old_keys = list(k_dict.keys())
        old_keys=sorted(old_keys)
        for k in old_keys:
            k=str(k)
            if k_dict.get(k) is not None:
                d[k] = k_dict[k]
                'Based on Oneliner palindrome checking trick'
                rev_k = k[::-1] # '01' -> '10'.
                if k != rev_k: #Avoid double counting palindromes; "11" + "11"..
                    d[k] += k_dict[rev_k] #Sum value.
                    # print ('summed (new)k: %s' %(d[k]))
                del(k_dict[rev_k]) #remove the "duplicate"/"mirror image" key+val.
        return d


    def GetCTD(self,ctd='CTD'):
        """
        Get Composition, Transition, Distribution descriptors.
        (Currently uses CTD_pro - external module).
        ctd = string of letters, according to which C,T,D will be called.

        Usage:
        result = GetCTD() #Returns all C,T,D
        composition = GetCTD('c') #Returns all C.
        Comp_And_Distribution =  = GetCTD('CD') #Returns all C and D.
        """
        res=CTD_pro.CalculateCTD(self.seq,str(ctd))
        return res

    # import protein_class
    'TODO -fix PI for cleaved or reduced aa sequences'
    def Get_SimpleProtParam (self):
        '''
        Get basic physical properties as in BioPython/ExPasy ProtParam
        module.
        Returns: PI, MW, GRAVY, aromaticity, Instability index.

        Note: These methods all assume a standard AA Alphabet.
        Warning! Returned PI is INNACCURATE For a parsed (Tail(s) removed) subseq.
        (BioPy-ProtParam.isoelectric_point assumes N,C terminii!)
        '''

        Bio_ProtParam = self.Bio_ProtParam #BioPython SequenceAnalysis object from str
        length=self.length
        PI = self.Get_PI()
        MW=Bio_ProtParam.molecular_weight()
        GRAVY=Bio_ProtParam.gravy()
        aromaticity=Bio_ProtParam.aromaticity()
        instability_index=Bio_ProtParam.instability_index()
        flex = self.GetFlex() #Returns 3 keys/values

        prot_pp = {'PI':PI,
        'Molecular Weight':round(MW,1), 'GRAVY':round(GRAVY,2),
        'Aromaticity':round(aromaticity,2),
         'Instability index':round(instability_index,2),
         'Sequence Length':length
         }
        prot_pp.update(flex)

        return prot_pp

    def GetLength (self):
        return {'Length':self.length}

    def GetFlex(self):
        '''
        Get parameters for B-values / flexibility. Built for a 9 window window;
        alters values from a BioPy ProtParam module.
        '''
        PP = self.Bio_ProtParam
        flex = PP.flexibility()
        Flex_mean =np.mean(flex)
        Flex_max=np.amax(flex)
        Flex_min=np.amin(flex)
        res = {'Flexibility (B-values):mean':Flex_mean,'Flexibility (B-values):Max':Flex_max,'Flexibility (B-values): Min':Flex_min}
        return res

    'TODO: Fix to work with subsequences (Remove N,C terminii charge in protparams calculation'
    def Get_PI (self):
        '''
        Get predicted PI. (Uses BioPython.ProtParam util)
        NOTE! For a "parsed"/cleaved subsequence,
        This is innaccurate unless fixed! (Doesn't take lack of N/C terminal into account)
        '''
        Bio_ProtParam = self.Bio_ProtParam
        'Maybe print a warning if self.HAS_N|C == True? '
        PI=Bio_ProtParam.isoelectric_point()
        return round(PI,2)

    def calculateProteinCharge(self, pH=7.2): # Has_N_Terminal=True,Has_C_Terminal=True):
        '''
        Get Net-Protein charge at given PH. Can alse be Used to predict PI of a protein.
        Be careful if used on a "partial" seq - +-N/C termini..
        If using only part of a sequence (Just the N-tail etc'), then modify call appropaitely.
        http://www.petercollingridge.co.uk/sites/files/peter/predictPI.txt
        We could also try using:
        http://pythonhosted.org/pyteomics/_modules/pyteomics/electrochem.html#charge
        '''
        sequence = self.seq
        AA_Counter = self.AA_Counts
        Has_N_Terminal = self.HAS_N
        Has_C_Terminal = self.HAS_N
        N_charge = 0
        C_charge = 0

        # @jit  #autojitCauses slow down for some reason
        def calculateAminoAcidCharge(amino_acid, pH):
            ratio = 1 / (1 + 10**(pH - pKa[amino_acid]))

            if charges[amino_acid] == 1:
                return ratio
            else:
                return ratio - 1

        if Has_N_Terminal:
            N_charge = calculateAminoAcidCharge('N-term', pH)
        if Has_C_Terminal:
            C_charge = calculateAminoAcidCharge('C-term', pH)
        protein_charge = N_charge+C_charge

        for amino_acid in pKa.keys():
            # protein_charge += sequence.count(amino_acid) * calculateAminoAcidCharge(amino_acid, pH)
            protein_charge += AA_Counter[amino_acid] * calculateAminoAcidCharge(amino_acid, pH)
        return protein_charge

    def get_netCharge(self,PH_ranges = [3.5,4.5,5.5,6.81,7.2,7.36,7.5,8.1]):
            '''
            SLOW!
            Return net charge of sequence for range of PH.
            Some (Wiki) PH ranges - Lysosomes - 4.5, Granules of chromaffin cells - 5.5,
            Cytosol 7.2, Blood 7.34–7.45, Mitochondrial matrix - 7.5,
            Pancreas secretions 8.1
            http://en.wikipedia.org/wiki/PH#Living_systems

            PH Of Stress response/oxidization (which can trigger unfolding in some proteins) ?
            '''
            st = "Net Charge-PH: "
            res = defaultdict(float)
            # 'add a map(round(x,3),x)'
            for PH in PH_ranges:
                res [st+str(PH)] = self.calculateProteinCharge(PH)
            return res

    def GetAliphaticness(self) :
        a = 2.9
        b = 3.9
        seq = self.seq
        length = float(self.length) #Otherwise, we 're dividing ints..
        AA_Counts = self.AA_Counts
        alanine_per = (AA_Counts['A'] / length )
        valine_per = (AA_Counts['V'] / length )
        isoleucine_per = (AA_Counts['I'] / length )
        leucine_per = (AA_Counts['L'] / length )
        # Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
        aliphatic_index = (100 * (alanine_per + a * valine_per + b * (isoleucine_per + leucine_per )))
        return {'Aliphaticness':aliphatic_index}

    '''
    TODO: Add Motifs.
    EG: Dibasic sites. Glycine amidation. Transit peptide motifs (from ends). Phosphorylation
    UTURE: Add option to select which motifs.
    AND-Stronger filtering -
    '''
    'TODO: option to filter protein for subsequence regions, and search for motifs ONLY in those regions (disordered, exposed etc'
    def GetPTMMotifs (self):
        '''
        Counts occurences of possible PTM motifs.
        This is a very weak feature, should be improved using
        local sequence properties and external subsequence learning
        (average regional flexibility, hydrophobicity, exposure, etc').

        We will likely want to replace or augment this with EXTERNAL Predictors.

        List of PTM motifs - ELM. (Can also do prediction, accounts for disorder, accessability etc')
        http://www.modpred.org/
        other predictors..

        '''
        seq = self.seq
        length = self.length
        'Potential N Glycosylation sites'
        # http://prosite.expasy.org/PS00001
        NSites = len(re.findall(r'N[^P][ST][^P]', seq))
        ' Aspartic acid, asparagine hydroxylation sites'
        hydroxSites = len(re.findall('CC.{13}C.{2}[GN].{12}C.C.{2, 4}C',seq))
        'C-tail ER targetting'
        # http://prosite.expasy.org/unirule/PRU10138/
        NSites_freq = (NSites / length)*100
        HydroxSites_freq = (hydroxSites / length)*100
        'TODO: pot. phosphorylation sites? glycation sites? '

        res = {'Potential N-Glycosylation Sites:Frequency':round(NSites_freq,2),'Potential Hydroxylation Sites:Frequency':round(HydroxSites_freq,3)}
        return res

    'TODO: Dibasic Cleavage sites feature - expand. Length normalize? (By log maybe?)'
    'TODO: Furin 20 AA surrounding region model'
    def GetCleavageCounts (self,AltLengthNormalize=False):
        '''
        Get various Potential sites (and/or regions/products),
         associated with cleavage of a precursor sequence/pre-pro-protein into
         products.
         EG: Neuropeptide Precursor processing.
         Many different potential cleavage motifs, dep. on organism, enzymes..!
         Default is Known motif model (Neuropeptides Precursor cleavage via flanking dibasic).:
         Southey et al. (2006b)  Known Motif model comprised of several prevalent motifs
         associated with neuropeptide precursor cleavage:
                Xxx-Xxx-Lys-Lys#,
                Xxx-Xxx-Lys-Arg#,
                Xxx-Xxx-Arg-Arg#,
                Arg-Xxx-Xxx-Lys#,
                Arg-Xxx-Xxx-Arg#.
          (Xxx = Any AA. # - Cleavage start site. Lys = K. Arg = R)
            Note - In many cases, Proline should NOT be present prior or adjacent to the cl)

         TODO: "Extract" all putative cleaved products/regions, and analyze them
         using AA propensity scales, etc'  (hydrophobicity, flexibility, products...).
         TODO: More sophisticated Furin Cleavage model suggested in:
         "A 20 Residues Motif Delineates the Furin Cleavage Site[...]"
        http://www.nuolan.net/motif.html

            TODO: Use external predictor to filter false positives.
                Eg: ELM. (Disorder,SS..)
            TODO: Option to input list of RegEx motifs. For other PC (Pre-convertases. not just furin..)

        NOTE: Currently, doesn't normalize for length..
        '''
        # from math import log2,sqrt
        seq_length = self.length
        seq = self.seq
        'If using alternative length for normalization, eg /log2, /sqrt, /fixed num...'
        'TODO - implement?'
        if AltLengthNormalize:
            seq_length=seq_length
        else:
            seq_length=int(seq_length)//4

        #Cleavage sites Based on to Known motif dibasic cleavage model. (+Strong likelihood of NO Proline adjacent)
        #Pay attention when building RegEx to avoid "double counting"!
        Motif_count1 = len(re.findall('R.[^P][RK]', seq)) #Known motif model. Arg-Xxx-Xxx-[Arg|Lys]
        Motif_count2 = len(re.findall('[^R][^P]K[RK]', seq)) #Known motif model. Not R at P1, to avoid overlap with other motif!
        Motif_count3 = len(re.findall('[^R][^P]R{2}', seq)) #Known motif model - Xxx-Xxx-Arg-Arg

        #Potential cleavage sites not from canonical known motif model. (EG - single, Triple..):
        # count_1 = len(re.findall('R.[RK]R.', seq))
        # count_4 = len(re.findall('[^P][RK]{3,4}[^P]', seq)) #3-4 dibasic site(s)

        count_putativepeps = len(re.findall('[^P][RK][^C]]{2,12}[][KR][KR]', seq)) #ad-hoc: search for potential pep. products

        known_motif_counts = Motif_count1+Motif_count2+Motif_count3
        known_motif_freq = (known_motif_counts/(seq_length-3))*100
        # adjacent_known_motif_counts= Motif_count2+Motif_count3
        # other_counts = count_1+count_4

        res = {"Dibasic Cleavage: Total Known Motif counts":known_motif_counts,
        "Dibasic Cleavage: Known Motif frequency":known_motif_freq,
        # "Adjacent_KnownMotif_counts":adjacent_known_motif_counts,
        # "Other_DBMotif_Counts":other_counts,
        "Dibasic Cleavage: Putative Cleavage Products":count_putativepeps,
        }
        return res



    def Dict_Keys_prefix(self,multilevelDict,PrefixStr):
        '''
        given a dict, Returns a new dict with prefix_str added before each key.
        Also used when using "alternative" alphs (so feature names will be unique,
            i.e "AA_14_KR Freq" vs "KR Freq" (for self.alph==AA_20))
        i.e:

        >>> aa_dict = {'fruit':'orange','K':0.23}
        >>>  print(cds.transform(aa_dict,"entropy"))
        '''
        return {PrefixStr+": "+str(key): (self.Dict_Keys_prefix(value,PrefixStr) if isinstance(value, dict) else value) for key, value in multilevelDict.items()}



    # 'TODO: Multiple windows/Windows per scale? '
    'TODO: Speed! (Store scales, rather than recalcing for subsegs or diff. window sizes). Also, not asnparray?'
    # 'TODO - Use own implementation of https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py - protein_scale '
    def Get_ParamScales_Features(self,window=7,edge=1.0,PickScales = None,trim=5): #,PP_scales):
        '''
        Extract features for given scale/propensity represention of a protein sequence.
        Calls on def Get_ParamScales.

        Default is to get all "preloaded" scales.
        User can choose to use only a limited number of scales, using PickScales

        Get 3: max, min, average for each scale.

        trim : get average of top/bottom "trim"'%' of values, and av. of in between values.
        (i.e. _AverageTrimMAX = average of top 20%,AverageTrimMIN - bottom 20%, and trimAv = Average of remaining 60%)
        '''
        #Get list of scales to use.
        PP_scales=self.Get_ParamScales(window,edge,PickScales) # 0.02 micro

        res = {}
        window_prefix = (' aaScale - Window size:'+str(window))
        for scale in PP_scales: # 0.06 micro

                name=scale[0]
                'Original:'
                # arr = np.asarray(scale[1])
                arr = scale[1]

                arr.sort() #Sort the array
                length = int(len(arr)/trim) #Get top/bottom X% (acc. to "trim" param)
                if length==0:
                    outlierLength = length+1  #min. 1
                else:
                    outlierLength = length

                scale_TrimMax = arr[:-outlierLength]
                scale_TrimMin = arr[0:outlierLength]
                scale_trimMean = arr[(outlierLength-1):-(outlierLength-1)] #Av of trimmed 80# +-1
        #Get max/min from arrays already sorted for top max/min (instead of searching orig):
                #scale_max=np.amax(arr)
                #scale_min=np.amin(arr)
                scale_max=scale_TrimMax[-1]
                scale_min=scale_TrimMin[0]
                #Get max and min + bottom/top of range for trimmed min/max .
                #MinScale_min=np.amax(scale_max)
                #MaxScale_max=np.amax(scale_max)

                # print('arr',len(arr))
                # print('arr[(outlierLength-1):-(outlierLength-1)]',len(arr[(outlierLength-1):-(outlierLength-1)]))
                # median = np.median(arr)

                res[name+window_prefix+' MAXIMUM ']=(scale_max)
                res[name+window_prefix+' MINIMUM ']=(scale_min)
                # res[name+window_prefix+'_AverageTrimMAX']=round((np.mean(scale_TrimMax)),2)
                # res[name+window_prefix+'_AverageTrimMIN']=round(np.mean(scale_TrimMin),2)
                # res[name+window_prefix+'_AverageTrimmed']=round((np.mean(scale_trimMean)),2)
                res[name+window_prefix+' Average-Trimmed MAX']=np.mean(scale_TrimMax)
                res[name+window_prefix+' Average-Trimmed MIN']=np.mean(scale_TrimMin)
                res[name+window_prefix+' Trimmed Average']=np.mean(scale_trimMean)
                #res[name+window_prefix+'_MEDIAN']=round(median,2)
                # scale_mode = 0  #Mode = Most common val. MICHAEL: I REMOVED MODE FOR BETTER PERFORMANCE
                # res[name+window_prefix+'_MODE']=np.around(scale_mode,decimals=2)

        return res

    def Get_SubSeqParamScales_Features(self,window=5,edge=1.0,PickScales = MinScales_Dict,segs=3):
        '''
        Similar to "Get_ParamScales_Features", but gets a small(er) set of features (per scale),
        for a smaller set (by efault) of AA scales, extracted from multiple segments of the sequence.

        Extract features for given scale/propensity represention of a protein sequence.
        By default, gets less features (min, max, average), with a minSet of scales. (Includes Atchley).
        An alternative, novel implementation is to divide into unequal segments -
        first,last 20%, middle 60% (or middle 30% *2). Dan Ofer.
        '''
       #Get list of scales to use.
        PP_scales=self.Get_ParamScales(window,edge,PickScales=MinScales_Dict)
        # res = defaultdict(float)
        res = {}
        window_prefix = (' aaScale:'+str(window))
        for scale in PP_scales:
                name=scale[0]
                arr = np.asarray(scale[1])
                # arr = scale[1]

                seg_size = int(len(arr) / segs)  # window size of each segment of the protein
                pos = 0

                for i in range(segs-1):
                    FeatPrefix=str('SubSequence(segment-'+str(i)+') '+window_prefix+str(name))
                    subArr = arr[pos:pos + seg_size]
                    res[FeatPrefix+' Average'] = np.mean(subArr)
                    res[FeatPrefix+' Maximum'] = np.amax(subArr)
                    res[FeatPrefix+' Minimum'] = np.amin(subArr)
        return res

    'TODO: Make it use/support a list of provided scales, or scales from a dict of scales. '
    'TODO: Copy BioPy protein scales data over to AAlphabets.py. Add Disorder affinity, Kidera factors..'
    def Get_ParamScales(self,window=7,edge=1.0,PickScales = None): ##,SA_window = 9, TMD_window = 19):
        '''
        Gets numerical represention of sequence,
        for each amino acid propensity scale (and window size). This is then
        used (later) to Get values for sequences using different amino propensities,
        via def Get_paramScales.

        Default is "built in" scales, but can be expanded easily.

        Default window size from literature is !17-19 for detecting TMDs' using hydrophobicity.

        Returns a list of (string:list) tuples.
        Uses:  Bio.SeqUtils.ProtParam.ProteinAnalysis(protein_scale(self, param_dict, window))
        http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParamData-module.html - builtin scales

        TODO: input "scales to use", "select which".
        TODO: Copy AA Scales from Bio.SeqUtils.ProtParamData to AAlphabets.py !
        '''
        Bio_PP = self.Bio_ProtParam
        PP_scales = []
        # kd_window_SURFACE = SA_window
        # kd_window_TMD = TMD_window
        # PP_scales = defaultdict(list)

        # Hydrophobicity_kd_TMD = Bio_PP.protein_scale(Bio.SeqUtils.ProtParamData.kd,kd_window_TMD,edge)
        # Hydrophilicity = Bio_PP.protein_scale(Bio.SeqUtils.ProtParamData.hw,window,edge)
        # Surface_access  = Bio_PP.protein_scale(Bio.SeqUtils.ProtParamData.em,kd_window_SURFACE,edge)

        # PP_scales = [
        # ('Hydrophobicity_kd_TMD',Hydrophobicity_kd_TMD),
        # ('Hydrophilicity',Hydrophilicity ),
        # ('Surface_access' ,Surface_access ),
        # ]

                ## #ORIGINAL, as of 10.11. Dan. :
        ## if PickScales is None: #Default is to use all preloaded scales from AAScales.py
        ##     'Scales_Dict is a dict of dicts from AAScales.py'
        ##     for scaleName,v in Scales_Dict.items():
        ##         # PP_scales [str(scaleName)] = Bio_PP.protein_scale(v,window,edge)
        ##         s = str(scaleName)
        ##         v = Bio_PP.protein_scale(Scales_Dict[str(scaleName)],window,edge)
        ##         # print('scaleName:',scaleName)
        ##         # print('v:',np.around(v,decimals=3))
        ##         v=np.around(v,decimals=2)
        ##         t = (s,v)
        ##         PP_scales.append (t)

        ## else:
        ##     for scaleName in PickScales: #Use only listed scales
        ##         # PP_scales [str(scaleName)] = Bio_PP.protein_scale(Scales_Dict[str(scaleName)],window,edge)
        ##         s=str(scaleName)
        ##         v = Bio_PP.protein_scale(Scales_Dict[str(scaleName)],window,edge)
        ##         v=np.around(v,decimals=1)
        ##         t = (s,v)
        ##         PP_scales.append (t)

        if PickScales is None: #Default is to use all preloaded scales from AAScales.py
           aaScales = Scales_Dict
           'Scales_Dict is a dict of dicts from AAScales.py'
        else:
#Use user defined list of scales
            aaScales=PickScales
        for scaleName,v in aaScales.items():
            # PP_scales [str(scaleName)] = Bio_PP.protein_scale(v,window,edge)
            s = str(scaleName)
            v = Bio_PP.protein_scale(aaScales[str(scaleName)],window=window,edge=edge)
            # print('scaleName:',scaleName)
            # print('v:',np.around(v,decimals=3))
            #v=np.around(v,decimals=2)
            t = (s,v)
            PP_scales.append (t)

        return PP_scales


    # 'TODO:FIX!  + Support for more letters! and analysis of periodicity (RQA)?'
    'TODO-  make the map work again (deal with gen)'
    'TODO: Add example of AutoCorr with normal aa scales (not out binary example)'
    def BinaryAutocorrellation(self,signal_letters = ['K','R'],loc=-1,LengthNormalization=None) :
        '''
        Calculate autocorrelation for sequence.
        signal_letters - Convert these letter(s) to "1" , and rest
        of letters in sequence to 0, then calc autocor as in signal analysis.

        http://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
        http://stackoverflow.com/questions/13439718/how-to-interpret-numpy-correlate-and-numpy-corrcoef-values?rq=1

        #How to Normalize - By Sum ("1's") or seq.length?
        lengthNorm = "length" normalization variable.
             By seq.length ("length") or sum of "1"s/signals ("sum").

        Loc = Normalization factor,  By Sum ("1's") or length.?
        '''
        # "str_prefix = 'AutoCor_'+str(signal_letters[0])"
        str_prefix = 'Binary AutoCorrellation ('+','.join(signal_letters)+')'

        sub_seq = self.seq
        for c in signal_letters:
            sub_seq = sub_seq.replace(str(c),'1')
        sub_seq = re.sub("\D", '0', sub_seq)
        # sub_seq = [map(int, sub_seq)] #Map is fastest. But is a generator in Py3.
        # sub_seq = [[map(int,i)] for i in sub_seq] #Map is fastest. But is a generator in Py3.
        sub_seq= [int(i) for i in sub_seq] #Lazy fix.S
        # print("Subseq for correlation: %s" %(sub_seq))
        selfCor = np.correlate(sub_seq, sub_seq, 'full')
        #Avoid divide by zero error:
        if sum(sub_seq) == 0:
            # return 0
            return {str_prefix+"MAX":0}
            # return {'AutoCorrelation: MAX':0}

                #How to Normalize - By Sum ("1's") or seq.length?
        'lengthNorm = "length" normalization variable. By seq.length or sum of "1"s/signals'
        if LengthNormalization=="sum":
            lengthNorm = sum(sub_seq)
        else:
            lengthNorm = float(len(sub_seq))

        '??? How is the autocor arranged, order & values..? TODO FIX'
        selfCor=sorted(selfCor)
        # Max_autoCor = selfCor[loc] / lengthNorm
        #Second highest (NOT "100%" Autocorrelation : loc=-2
        # autoCor1 = selfCor[-1] / lengthNorm
        autoCor2 = selfCor[-2] / lengthNorm
        autoCor3 = selfCor[-3] / lengthNorm
        autoCor4 = selfCor[-4] / lengthNorm
        autoCor5 = selfCor[-5] / lengthNorm
        autoCor6 = selfCor[-6] / lengthNorm
        autoCor8 = selfCor[-8] / lengthNorm
        autoCor9 = selfCor[-9] / lengthNorm
        autoCor12 = selfCor[-12] / lengthNorm
        # 'autoCorrelation_1':autoCor1,

        res = {
        ' Lag:2':autoCor2,
        ' Lag:3':autoCor3,
        ' Lag:4':autoCor4,
        ' Lag:5':autoCor5,
        ' Lag:6':autoCor6,
        ' Lag:8':autoCor8,
        ' Lag:9':autoCor9,
        ' Lag:12':autoCor12
        }
                # 'autoCorrelation_MAX':Max_autoCor
        return self.Dict_Keys_prefix(res,str(str_prefix))

    '''
    IDEAs: for more Entropy/complexity based measures: http://www.nbb.cornell.edu/neurobio/land/PROJECTS/Complexity/'
    Another idea - Permutation Entropy - http://tocsy.pik-potsdam.de/petropy.php'
    '''

    'TODO Check values, magnitudes - should be maximum ~4.xx bits max for entropy? '
    def GetEntropy (self,normalizeTotEntropy=False, getLettersEntropy=True):
        '''
        http://bugra.github.io/work/notes/2014-05-16/entropy-perplexity-image-text/
        VS Entropy_Kap_AA , Entropy_seq..

        normalizeTotEntropy : If True, then divide total sequence entropy by log2(length).
                                (Check if correct!!)
        getLettersEntropy : If True, also return  entropy per letter; If false, then
        return only the "total entropy".
        '''
        # word_count_information = []
        AA_information = {}

        seq = self.seq
        length = float(len(seq))
        wordset = set(self.alph)
        # freq = self.AA_Counts
        "Don't count entropy for letters which don't appear:"
        freq = {k:v for (k, v) in self.AA_Counts.items() if v != 0}
        entropy = 0
        # for word in wordset:
        for word in freq.keys():
            probability = freq[word] / (1.0 * length)  #Could be replaced with copied use of AA_Freq
            self_information = np.log2(1.0/probability)
            entropy += (probability * self_information)
            # word_count_information.append([word, freq[word], self_information])
            if getLettersEntropy==True:
                AA_information[str(word)+' Entropy']=self_information

        if normalizeTotEntropy is True:
            from math import log2
            l = log2(length)
            AA_information['Total Entropy - Normalized By Length']=(float(entropy)/l)
        # else:
        AA_information['Total Entropy']=entropy # Equivalent to old Entropy_Seq.

        return self.alphabet_prefix(AA_information)

    'DEPRECATE in favor of GetEntropy! '
    "IDEA: Maybe get/compare each letter's entropy, VS Whole sequence."
    "AND - Whole seq. Entropy vs (store) max entropy for a sequence of that length and alphabet"
    def Entropy_Seq(self): #,s):
        '''
        DEPRECATE! Use GetEntropy Instead!!

        Calculate the (Shannon) information entropy  of a given input string.
        Entropy is the expected value of the measure of information content in system.
        http://rosettacode.org/wiki/Entropy#Python

        Qu - For later comparing to Max Entropy - use possible AA counts?
        '''
        s=self.seq
        p=self.AA_Counts
        lns = self.length
        def count_e():
            return(-sum( count/lns * log(count/lns, 2) for count in p.values()))
        entropy=count_e()
        return {'Sequence Entropy: ':count_e()}


    def tail_properties (self, tail_end='N',GetFullAAFreq=True,reduced_alph = 'Ofer_N_Tail'):
        '''
        Method that gets a number of N or C end terminii tail/subsequence's
        properties.
        IMPORTANT! This method assumed that the self object/protein constructed
         is JUST the tail!
         (I.E - it is NOT meant to be used on part of an instance's own sequence,
            but to be called in the subsequence as its own, constructd object/).

        Note: Differences in calling netCharge, paramscale...
        PH is currently innacurate in assuming a terminii end.
        Assume N-length as being fixed (~25). Likely reduced in its AA alphabet.
        Consider whether to use AA frequency or counts.

        Signal Peptides in Euk are usually ~15-25 AA long,
        often composed of 3 subsections: "charged+"-"Hydrophobic"-"Polar(uncharged)".
        '''
        seq = self.seq
        length = self.length
        PrefixStr = str(tail_end)+'_Tail '
        res = {}
        'Get AA frequency (for reduced alph if wanted)'
        #TODO - this is a mess..
        if GetFullAAFreq:
            Tail_AAFreq = {'Tail: AA Frequency':(self.GetAA_Freq())}
        else: #AAlphabets - get a translate tail sequence in reduced alphabet
            reduced_tail_sequence = translate_sequence(seq, REDUCED_ALPHABETS_TRANSDICTS[reduced_alph])
            tail_counts = Counter(reduced_tail_sequence)
            tail_counts_freq = self.countToFreq(tail_counts)
            Tail_AAFreq= self.Dict_Keys_prefix(tail_counts_freq,str(reduced_alph)+'Tail_AAFreq')
            # reduced_tail_pp = pp.ProteinAnalysis(str(reduced_tail))
            # Tail_AAFreq = {str(reduced_alph)+'Tail_AAFreq':reduced_tail_pp.get_amino_acids_percent()}
            # Tail_AAFreq= self.Dict_Keys_prefix(reduced_tail_pp.get_amino_acids_percent(),str(reduced_alph)+'Tail_AAFreq')
            # Tail_AAFreq= self.Dict_Keys_prefix(reduced_tail.GetAA_Freq(),str(reduced_alph)+'Tail_AAFreq')

        res.update(Tail_AAFreq)

        PH_range = [4.5,5.5,6.8,7.2,8.1]
        PI = {'PI':self.Get_PI()}
        res.update(PI)
        res.update(self.get_netCharge(PH_ranges = PH_range))

        # res.update(self.Entropy_Seq()) #Consider using this also with reduced alph..
        res.update(self.GetEntropy(getLettersEntropy=True)) #Newer Entropy feature
        res.update(self.GetAliphaticness())
        res.update(self.Get_ParamScales_Features(window=4)) #Check window of size 4

        'Get "some" CTD based features'
        'TODO: Implement a better C getting from CTD Pro..'
        res.update(CTD_pro.CalculateCompositionCharge(seq))
        res.update(CTD_pro.CalculateCompositionHydrophobicity(seq))
        res.update(CTD_pro.CalculateCompositionSecondaryStr(seq))
        res.update(CTD_pro.CalculateCompositionPolarity(seq))
        res.update(CTD_pro.CalculateCompositionPolarizability(seq))

        res.update(CTD_pro.CalculateDistributionCharge(seq))
        res.update(CTD_pro.CalculateDistributionHydrophobicity(seq))
        res.update(CTD_pro.CalculateDistributionSecondaryStr(seq))
        res.update(CTD_pro.CalculateDistributionPolarizability(seq))

        res.update(CTD_pro.CalculateTransitionCharge(seq))
        res.update(CTD_pro.CalculateTransitionHydrophobicity(seq))
        res.update(CTD_pro.CalculateTransitionSecondaryStr(seq))


        return (self.Dict_Keys_prefix(res,PrefixStr))

    'TODO: think what we want as "deafault".. and update docstring'
    'TODO: allow input of variables (param scale window size..)'
    'TODO: Make use of features by a paramter ("GetAliphaticness=True,..." )'
    def GetSimpleFeatures (self,ParamScaleWindow=6,segDivide=1,DisorderSegments=4):
        '''
        Returns A large number of "default" features, good for most cases.
        This is meant to be called when the sequence is
        using the normal 20 letter alphabet. (Otherwise, physical
            parameters returned will be inaccurrate).

        '''
        if self.alph_set != 'AA20' :
            print ('Warning! Reduced Protein sequences should not not be used with this method')
        res = {}
        res.update(self.GetAA_Freq())
        res.update(self.GetEntropy()) #NEwer entropy method
        res.update(self.Get_SimpleProtParam())
        res.update(self.Get_ParamScales_Features(window=ParamScaleWindow))
        res.update(self.GetAliphaticness())

        res.update(self.get_netCharge())
        res.update(self.BinaryAutocorrellation()) #Default is K,R, KK, RR'
        res.update(self.BinaryAutocorrellation(signal_letters=['C']))

        res.update(self.getCysteineMotifs(segDivide=segDivide)) #Tal Arian feature
        res.update(self.getFIDisorder(segments=DisorderSegments)) #Tal Arian feature
        return res

    def cysteineMotif(self,segDivide=3): ##If used as part of regular Prot package, import self.  (Why have it as static submethod? Less consistant..)
        '''
        CHANGED: Threshhold for counts removed, made into variable feature instead of discrete
        Note: We could re-increase the threshhold (min 1 count instead of >1).
        And also then divide into more segs..
        '''
        seq=self.seq
        cysPattern = r'C[^C]{0,3}C'
        length = len(seq)  ##If used as part of regular package,maybe use length=self.length instead; seq=self.seq
        # window_size = length / 5  # window size 20% of the protein length #Orig
        window_size = int(length /segDivide)
        # scores = [0 for _ in range(5)]
        scores = [0 for _ in range(segDivide)]
        prog = re.compile(cysPattern)
        pos = 0
        # for i in range(3): ##This will get from 0 to 3.. Shouldn't it be range(5)? D #ORIG
        # for i in range(4): ##Changed range to 4 from 3. D.
        for i in range(segDivide): #Changed to segDivide
            # print(len(prog.findall(seq[pos:pos + window_size])))
            scores[i] = len(prog.findall(seq[pos:pos + window_size]))
            pos += window_size
        scores[-1] = len(prog.findall(seq[pos:])) #-1? Why not do scores.append? - D
        res = {}
        key = "Cysteine window Motif:"
        for i, score in enumerate(scores):
            #Changed if score > 2 to: if score > 1
            "res[key + str(i)] = 1 if score > 1 else 0   #Why not >= ? 3 means at least 6 cysteins! That's a lot"
            res[key + str(i)] = score
        return res

    #@staticmethod
    def cysteineSpaceMotif(self):
        seq=self.seq
        cysPattern = r'C[^C]{0,3}C[^C]{15,40}C'
        match = re.search(cysPattern, seq)
        score = 1 if match else 0 # if match isn't none, the score is 1 ("space motif" was found)
        return {"Cysteine Spaced Motif:": score}


    def getCysteineMotifs(self,segDivide=3): #(Why have called submethods as static (and not calling to "self")? Less consistant; .)
        '''
        Get Cysteine spacer motif, and counts of frequent CxxC motifs over protein sequence
        '''
        res = {}
        res.update(self.cysteineMotif(segDivide=segDivide))
        res.update(self.cysteineSpaceMotif())
        return res

    def getFIDisorder(self,segments=5):
        '''
        Divide protein sequence into segments, and returns for each seg.
        predicted disorder (Y/N) according to FoldIndex method (Uversky et al).
        Method implemented with help of Tal Arian.
        '''
        seq=self.seq
        return getDisordered(seq=seq,segments=segments)


#####################################################################################################
if __name__=="__main__":

    import timeit

    protein_2="ADGCGVGEGKRKRKTGQGPMCNCMCMKPKKHRHWVYADEDAADLESDSFADEDARRRKKRKKSLESDSFPWSNQRVFCSFADEDASAGCVDEDASLESDSFPWSNQRVFCSFADEDASADGFPWSNQRVFCSFADEDAS"
    protein="CACACAACAKKKKKKRMQMQQCKRAQRCARRCAKFDKQRYKRYYQRMCAQMNEEDACLEAGYCLGTTLSSWPKRKRLHFMEEQSQSTMLMGIGIGALLTLAFVGIFAFCCFVYRRVRRLRRAEPTPQYRFRKRDKVMFYGRKIMRKVTTLPHTLVGNTSAPRQRVRKLGLTTELFNPESQAIPLLSVA"
    cds=ProtFeat(protein)

    # print(cds.GetAA_Freq())

    print(cds.getFIDisorder())
    print(cds.getCysteineMotifs())

    # print (cds.GetSimpleFeatures())
    # print(cds.Get_SimpleProtParam())

    print (cds.GetEntropy(normalizeTotEntropy=True))
    print (cds.BinaryAutocorrellation())

    # print (cds.GetAA_Freq())
    # print ('tail \n',cds.tail_properties(GetFullAAFreq=True,tail_end='C'))

    # print(cds.Get_ParamScales_Features(6))
    # print(cds.Get_SubSeqParamScales_Features())

    # print(timeit.timeit(stmt=cds.Get_ParamScales_Features, number=40))
    # print(timeit.timeit(stmt=cds.Get_SubSeqParamScales_Features, number=40))

#    print(cds.GetCleavageCounts())
    # print(cds.GetPTMMotifs())

    # print(NGram_combinations(gbm4,3))
    # print(cds.calculateProteinCharge(7.2))
    # print(cds.get_netCharge())
    # print('length:',cds.length)
    print("\n \n")
    print (cds.GetCTD())

    # print ('KGram Frequencies: %s' %(cds.GetkgramFreq(2)))
    # b=cds.GetKMirrorsFreq(2)
    # d=cds.KMirrors(b)
    # print('b-fun-MIRROR Freq: %s' %(d))


