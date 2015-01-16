import os
import sys
import itertools
import math
from functools import wraps

import numpy

from biopy import file_io



###############################################################################
# Helper functions / decorators
###############################################################################


def memoize(func):
    '''
    This function can be used as decorator to store results of computationally
    expensive functions in a cache store (dict)
    '''
    cache = {}

    @wraps(func)
    def wrapper(*args):
        if args not in cache.keys():
            cache[args] = func(*args)
        return cache[args]

    return wrapper


###############################################################################
# Paths to data files
###############################################################################

# path to data directory
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# amino acid scale files
AA_SCALE_DIR = os.path.join(DATA_DIR, 'amino_acid_scales')
AA_SCALE_GEORGIEV_F = os.path.join(AA_SCALE_DIR, 'georgiev.txt')
AA_SCALE_AAINDEX_F = os.path.join(AA_SCALE_DIR, 'aaindex1')
AA_SCALE_PSEAAC_F = os.path.join(AA_SCALE_DIR, 'pseaac.txt')

# amino acid matrix files
AA_MATRIX_DIR = os.path.join(DATA_DIR, 'amino_acid_matrices')
AA_MATRIX_AAINDEX_F = os.path.join(AA_MATRIX_DIR, 'aaindex2')
AA_MATRIX_SW_F = os.path.join(AA_MATRIX_DIR, 'schneider_wrede.txt')

###############################################################################
#
# AMINO ACID ALPHABET
#
# This section defines the amino acid alphabet and corresponding amino acid
# short and long names. The data is derived from:
# www.ebi.ac.uk/2can/biology/molecules_small_aatable.html
#
#
# IMPORTANT: order of the alphabet letters is important! do not change this.
#
###############################################################################

# unambiguos amino acid alphabet (alphabitic order full names)
aa_unambiguous_alph = 'ARNDCEQGHILKMFPSTWYV'
aa_unambiguous_short = [
    'ala', 'arg', 'asn', 'asp', 'cys', 'glu', 'gln', 'gly', 'his', 'ile',
    'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val'
]
aa_unambiguous_name = [
    'alanine', 'arginine', 'asparagine', 'aspartic acid', 'cysteine',
    'glutamic acid', 'glutamine', 'glycine', 'histidine', 'isoleucine',
    'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine',
    'threonine', 'tryptophan', 'tyrosine', 'valine'
]

# ambiguous amina acids
aa_ambiguous_alph = 'BJZX'
aa_ambiguous_short = ['asx', 'xle', 'xaa', 'glx']
aa_ambiguous_name = [
    'aspartic acid or asparagine', 'leucine or isoleucine',
    'unknown amino acid', 'glutamic acid or glutamine'
]

# special amino acids
aa_special_alph = 'UO'
aa_special_short = ['sec', 'pyl']
aa_special_name = ['selenocysteine', 'pyrralysine']

# stop codon translation
aa_ter_alph = '*'
aa_ter_short = ['ter']
aa_ter_name = ['terminal']

# full alphabet
aa_alph = aa_unambiguous_alph +\
    aa_ambiguous_alph + aa_special_alph + aa_ter_alph
aa_short = list(itertools.chain.from_iterable([
    aa_unambiguous_short, aa_ambiguous_short, aa_special_short, aa_ter_short
]))
aa_name = list(itertools.chain.from_iterable([
    aa_unambiguous_name, aa_ambiguous_name, aa_special_name, aa_ter_name
]))

# mapping from list of (sorted) unambiguous aas to ambiguous aa letter
map_to_ambiguous_aa = dict(zip(list(aa_unambiguous_alph),
                           list(aa_unambiguous_alph)))
map_to_ambiguous_aa['*'] = '*'
map_to_ambiguous_aa['DN'] = 'B'
map_to_ambiguous_aa['IL'] = 'J'
map_to_ambiguous_aa['GQ'] = 'Z'


###############################################################################
#
# AMINO ACID SCALES
#
# Amino acid scales are mappings from the unambiguous amino acids to a value
# that describes some kind of property of the amino acid, such as the size or
# the hydrophobicity.
#
# Each of the rows in the following arrays is a scale, and each column contains
# the values for one amino acid. The order of the columns (amino acids) is the
# same as the order of the amino acids in aa_unambiguous_alph.
#
###############################################################################


def _standardized_scale(scale):
    letters, values = zip(*scale.iteritems())
    mean = numpy.mean(values)
    std = numpy.std(values)
    st_values = [(v - mean) / std for v in values]
    return dict(zip(letters, st_values))


' 19 varimax Georgiev scales'
georgiev_scales, georgiev_scale_ids = file_io.read_scales(AA_SCALE_GEORGIEV_F)
georgiev_st_scales = [_standardized_scale(s) for s in georgiev_scales]

# 544 amino acid scales from AAindex ver.9.1
aas_tuples = file_io.read_scales_db(AA_SCALE_AAINDEX_F)
aaindex_scale_ids, aaindex_scale_descr, aaindex_scales = zip(*aas_tuples)
aaindex_st_scales = [_standardized_scale(s) for s in aaindex_scales]

# 6 amino acid scales used by PseAAC method
pseaac_scales, pseaac_scale_ids = file_io.read_scales(AA_SCALE_PSEAAC_F)


# TODO index_id? descr?
def get_aaindex_scale(index, ambiguous=True, standardized=True):
    '''
    If standardize, values (of the unambiguous amino acids) are standardized
    to mean 0.0 and standard deviation 1.0.
    '''
    if(standardized):
        scale = aaindex_st_scales[index].copy()
    else:
        scale = aaindex_scales[index].copy()

    return _get_scale(scale, ambiguous)


def get_georgiev_scale(index, ambiguous=True, standardized=True):
    '''
    If standardize, values (of the unambiguous amino acids) are standardized
    to mean 0.0 and standard deviation 1.0.
    '''
    if(standardized):
        scale = georgiev_st_scales[index].copy()
    else:
        scale = georgiev_scales[index].copy()

    return _get_scale(scale, ambiguous)


def get_georgiev_scales(ambiguous=True, standardized=True):
    return [get_georgiev_scale(i, ambiguous, standardized) for i in xrange(19)]


def get_pseaac_scale(index, ambiguous=True):
    scale = pseaac_scales[index].copy()
    return _get_scale(scale, ambiguous)


def _get_scale(scale, ambiguous):
    '''
    '''
    if(ambiguous):
        scale.update(_get_non_aa_letter_dict())
    return scale


def _get_non_aa_letter_dict():
    '''
    '''
    other_letters = aa_ambiguous_alph + aa_special_alph + aa_ter_alph
    return dict(zip(other_letters, len(other_letters) * [0.0]))


@memoize
def normalized_pseaac_scale(pseaac_scale_index):

    if not(pseaac_scale_index in xrange(len(pseaac_scales))):
        raise ValueError('PseAAC scale index out of range.')

    letters, values = zip(*[i for i in
                            pseaac_scales[pseaac_scale_index].iteritems()])

    # list with scale values
    values = numpy.array(values)

    # mean scale value
    mean = values.mean()

    # scale values with the mean subtracted
    values_minus_mean = values - mean

    # squared distances between original and mean
    summed_sq_dists_to_mean = sum(numpy.power(values_minus_mean, 2))

    # the denominator
    denominator = math.sqrt(summed_sq_dists_to_mean / 20)

    # normalized scale values
    norm_values = values_minus_mean / denominator

    return dict(zip(letters, norm_values))


@memoize
def pseaac_correlation_matrix(pseaac_scale_indices):
    '''
    pseaac_scale_indices should be sorted tuple with scale indices...!
    don't forget, tuple with one item is (0,)!
    '''

    scales = [normalized_pseaac_scale(i) for i in pseaac_scale_indices]
    alph = aa_unambiguous_alph

    matrix = {}
    for index, aa1 in enumerate(alph):
        for aa2 in alph[index:]:
            sqr_dists = [math.pow(s[aa1] - s[aa2], 2) for s in scales]
            avg_dists = sum(sqr_dists) / len(scales)
            matrix[aa1 + aa2] = avg_dists
            matrix[aa2 + aa1] = avg_dists

    return matrix


@memoize
def pseaac_type2_correlation_matrix(pseaac_scale_index):

    scale = normalized_pseaac_scale(pseaac_scale_index)
    alph = aa_unambiguous_alph

    matrix = {}
    for index, aa1 in enumerate(alph):
        for aa2 in alph[index:]:
            prod = scale[aa1] * scale[aa2]
            matrix[aa1 + aa2] = prod
            matrix[aa2 + aa1] = prod

    return matrix

###############################################################################
#
# AMINO ACID DISTANCE MATRICES
#
###############################################################################


aa_matrix_ids, aa_matrix_descr, aa_matrices =\
    zip(*file_io.read_aa_matrix_db(AA_MATRIX_AAINDEX_F))

aa_matrix_sw = file_io.read_aa_matrix(AA_MATRIX_SW_F)


###############################################################################
#
# AMINO ACID CLUSTERS
#
# 4 (www.ebi.ac.uk/2can/biology/molecules_small_aatable.html)
# 7 source: taylor85 adjusted version on the url above
# 2 wikipedia aa propensities
#
# #############################################################################


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

###############################################################################
#
# AMINO ACID STATE SEQUENCES
#
###############################################################################


# secondary structure alphabet
ss_alph = 'CHE'
ss_short = ['col', 'hel', 'str']
ss_name = ['random coil', 'helix', 'strand']

# solvent accessibility alphabet
sa_alph = 'BE'
sa_short = ['bur', 'exp']
sa_name = ['buried', 'exposed']


###############################################################################
#
# NUCLEOTIDE ALPHABET
#
# IMPORTANT: order of the alphabet letters is important! do not change this.
#
###############################################################################


# ambiguous nucleotide alphabet
nucleotide_unambiguous_alph = 'TCAG'

# unambiguous nucleotide alphabet
nucleotide_ambiguous_alph = 'MRWSYKVHDBN'

# full nucleotide alphabet
nucleotide_alph = nucleotide_unambiguous_alph + nucleotide_ambiguous_alph

# mapping from ambiguous alphabet to list of unambiguous nucleotides
map_to_unambiguous_nucleotides = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'W': 'AT', 'S': 'GC', 'M': 'AC', 'K': 'GT', 'R': 'AG', 'Y': 'CT',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT'
}


###############################################################################
#
# CODON ALPHABET
#
# IMPORTANT: order of the alphabet is important! do not change this.
#
###############################################################################


# unambiguous codon alphabet
codons_unambiguous = [a+b+c for a in nucleotide_unambiguous_alph
                      for b in nucleotide_unambiguous_alph
                      for c in nucleotide_unambiguous_alph]

# full codon alphabet
codons = [a+b+c for a in nucleotide_alph for b in nucleotide_alph
          for c in nucleotide_alph]


###############################################################################
#
# CODON TRANSLATION UTILS
#
###############################################################################


# amino acids corresponding to codons in codons_unambiguous
codon_aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

# unambiguous codon table mapping (codon --> amino acid)
codon_table_unambiguous = dict(zip(codons_unambiguous, codon_aas))


def unambiguous_codons(codon):
    '''
    This function returns all possible unambiguous codons for the provided,
    possibly ambiguous, codon.

    Args:
        codon (str): A possibly ambiguous codon.
    Raises:
        TODO

    If codon is an unambiguous codon, the list with only this codon is
    returned.

    >>> unambiguous_codons('GCG')
    ['GCG']

    If codon is an ambiguous codon, the list with all possible unambiguous
    codons for this particular codon is returned.

    >>> unambiguous_codons('ATN')
    ['ATA', 'ATC', 'ATG', 'ATT']
    '''
    return [x + y + z for x in map_to_unambiguous_nucleotides[codon[0]]
            for y in map_to_unambiguous_nucleotides[codon[1]]
            for z in map_to_unambiguous_nucleotides[codon[2]]]


def unambiguous_aas(codon):
    '''
    This function returns the amino acids for which the provided (possibly
    ambiguous) codon encodes.

    Args:
        codon (str): A possibly ambiguous codon.
    Raises:
        TODO

    If the provided codon is ambiguous, the amino acid letter for which this
    codon encodes is returned.

    >>> unambiguous_aas('ATG')
    'M'

    If the provided codon is unambiguous, a string with all possible amino
    acids that can be encoded by this codon is returned. Thes returned string
    is sorted alphabetically.

    >>> unambiguous_aas('ATN')
    'IM'
    >>> unambiguous_aas('ANN')
    'IKMNRST'
    >>> unambiguous_aas('NNN')
    '*ACDEFGHIKLMNPQRSTVWY'

    '''
    return ''.join(sorted(set([codon_table_unambiguous[c]
                   for c in unambiguous_codons(codon)])))


def ambiguous_aa(codon):
    '''
    This function returns the posibbly ambiguous amino acid that is encoded by
    the provided codon. In all cases, this method returns only one letter.

    Args:
        codon (str):
    Raises:
        TODO

    If the codon is unambiguous, the corresponding unambiguous amino acid is
    returned.

    >>> ambiguous_aa('ATG')
    'M'

    If the codon is ambiguous but encodes for a single amino acid, this amino
    acid is returned.

    >>> ambiguous_aa('TCN')
    'S'

    If the codon is ambiguous and does not encode for a single amino acid, the
    ambiguous amino acid that best represents the set of unambiguous amino
    acids that can be encoded by the ambiguous codon is returned.

    >>> ambiguous_aa('MTT')
    'J'

    In most cases this will be an X, which encodes for any amino acid.

    >>> ambiguous_aa('ATN')
    'X'
    '''
    sorted_unamb_aas = unambiguous_aas(codon)
    return map_to_ambiguous_aa.get(sorted_unamb_aas, 'X')

# ambiguous codon table mapping
codon_table = dict(zip(codons, [ambiguous_aa(c) for c in codons]))


def amino_acid_codons(aa):
    '''
    Returns the list of codons that encode for the provided amino acid aa.
    '''
    codon_is = [i for i, a in enumerate(codon_aas) if a == aa]
    return [codons_unambiguous[i] for i in codon_is]

'''
# number of codons per amino acid
ncodon_per_aa = numpy.array([codon_aas.count(l) for l in codon_aas])
'''


###############################################################################
#
# SEQUENCE OPERATIONS
#
###############################################################################


def translate(orf):
    '''
    This function translates a (possibly ambiguous) ORF nucleotide sequence
    into an amino acid protein sequence.

    Args:
        orf (str): The open reading frame (nucleotide) sequence.

    Returns:
        str The translated protein (amino acid) sequence.

    Raises:
        ValueError: if the orf length is not a multiple of 3.
        ValueError: if orf is not an (possibly ambiguous) nucleotide sequence.

    Translation of a random ORF sequence part can be done with:

    >>> translate('ATGTTTAGTAACAGACTACCACCTCCAAAA')
    'MFSNRLPPPK'

    Ambiguous codons will be translated to corresponding (possibly ambiguous)
    amino acids.

    >>> translate('ATGMTTAGTAACAGACTACCACCTCCAAAA')
    'MJSNRLPPPK'
    '''
    if not(len(orf) % 3 == 0):
        raise ValueError('ORF sequence length is not a multiple of 3.')
    if not(is_nucleotide_sequence(orf)):
        raise ValueError('ORF sequence is not a nucleotide sequence.')

    return ''.join([codon_table[orf[i:i+3]] for i in xrange(0, len(orf), 3)])


def segment(seq, num_segments):
    '''
    This methods chops seq in num_segments (about) equals sized segments and
    returns the list of sequence segments.

    >>> s = 'AAAABBBBCCCC'
    >>> segment(s, 4)
    ['AAA', 'ABB', 'BBC', 'CCC']
    >>> segment(s, 5)
    ['AA', 'AA', 'BB', 'BBC', 'CCC']

    '''

    if(num_segments <= 0):
        raise ValueError('Number of segments must be a positive integer.')

    stepsize = len(seq) / num_segments
    rest = len(seq) % num_segments

    # number of segments larger than sequence...
    if(stepsize == 0):
        segments = list(seq)
        segments.extend([''] * (num_segments - len(seq)))

    else:
        chop_indices = range(0, len(seq), stepsize)

        add_rest = [0] * (num_segments - (rest - 1))
        add_rest.extend(range(1, rest + 1))

        chop_indices = [sum(i) for i in zip(chop_indices, add_rest)]
        if(rest == 0):
            chop_indices.append(len(seq))

        segments = []
        for i in range(len(chop_indices) - 1):
            segments.append(seq[chop_indices[i]:chop_indices[i + 1]])

    return segments


###############################################################################
#
# HELPER FUNCTIONS FOR SEQUENCE FEATURE CALCULATION
#
###############################################################################


def ordered_alph_pairs(alph):
    '''
    This function returns all possible ordered letter pairs for the provided
    alphabet alph.

    For the intended use, alph should be a non-empty string without duplicate
    characters.

    Args:
        alph (str): The letter alphabet

    Returns:
        [str] List with length-2 strings.

    >>> ordered_alph_pairs('AB')
    ['AA', 'AB', 'BA', 'BB']
    '''
    return [''.join((a, b)) for a in alph for b in alph]


def ordered_seq_pairs(seq, distance):
    '''
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
    if(distance <= 0):
        raise ValueError('Distance should be larger than 0.')

    return [''.join(t) for t in zip(seq[:-distance], seq[distance:])]


###############################################################################
#
# SEQUENCE FEATURE CALCULATION FUNCTIONS
#
###############################################################################


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


def diletter_count(seq, alph, distance):
    '''
    Returns pair counts that occur at the given distance in seq.

    >>> diletter_count('ABBBA', 'AB', 1)
    array([0, 1, 1, 2])
    >>> diletter_count('ABBBA', 'AB', 2)
    array([0, 1, 1, 1])
    >>> diletter_count('ABBBC', 'AB', 1)
    array([0, 1, 0, 2])
    '''

    # construct empty dictionory to store diletter counts
    alph_pairs = ordered_alph_pairs(alph)
    pair_counts = dict(zip(alph_pairs, len(alph_pairs) * [0]))

    # count pairs while walking over the sequence
    for seq_pair in ordered_seq_pairs(seq, distance):
        if(seq_pair in pair_counts.keys()):
            pair_counts[seq_pair] += 1

    return numpy.array([pair_counts[p] for p in alph_pairs], dtype=int)


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


def diletter_composition(seq, alph, dist):
    '''
    This function returns the ordered pair compositions of pairs that occur
    at distance dist in seqeunce seq.

    >>> diletter_composition('ABBBA', 'AB', 1)
    array([ 0.  ,  0.25,  0.25,  0.5 ])
    '''
    if(len(seq) == 0):
        raise ValueError('Cannot calculate composition of empty sequence.')

    return diletter_count(seq, alph, dist) / float(len(seq) - dist)


def autocorrelation(ac_type, sequence, scale, lag):
    '''
    Based on the autocorrelation type (ac_type), this function calls that
    autocorrelation function and returns the result.

    ValueError is raised if wrong ac_type is provided.
    '''

    try:
        return getattr(sys.modules[__name__],
                       'autocorrelation_%s' % (ac_type))(sequence, scale, lag)
    except AttributeError:
        raise ValueError('Wrong autocorrelation type provided.')


def autocorrelation_mb(sequence, scale, lag):
    '''
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


def autocorrelation_moran(sequence, scale, lag):
    '''
    This function returns uses the provided scale to transform the sequence seq
    into a signal and returns the autocorrelation of this signal for the given
    lag.

    Moran autocorrelation as given in Li (2006) PROFEAT.

    TODO formula

    >>> autocorrelation_moran('ABACABAC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    1.0
    >>> autocorrelation_moran('ABACABAC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 2)
    -1.0
    >>> autocorrelation_moran('BBBBCCCC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    -1.0
    >>> autocorrelation_moran('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    0.0
    >>> autocorrelation_moran('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 8)
    0.0
    >>> autocorrelation_moran('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 12)
    0.0
    '''

    if(lag < 1):
        raise ValueError('The provided lag should be a positive integer.')

    # TODO check if this is the correct approach
    if(len(sequence) <= lag):
        return 0.0

    # transform sequence to signal using the provided scale
    signal = numpy.array(seq_signal_raw(sequence, scale))

    # calculate average signal value
    avg_signal = numpy.mean(signal)

    # subtracted average from signal
    signal = signal - avg_signal

    # calculate autocorrelation
    autocorr = sum(signal[:-lag] * signal[lag:])

    # normalize
    norm_autocorr = autocorr / float(len(sequence) - lag)

    # second normalization
    sec_norm = sum(numpy.square(signal)) / float(len(sequence))

    if(sec_norm == 0):
        return 0.0  # TODO should this be zero???
    else:
        return norm_autocorr / sec_norm


def autocorrelation_geary(sequence, scale, lag):
    '''
    This function returns uses the provided scale to transform the sequence seq
    into a signal and returns the autocorrelation of this signal for the given
    lag.

    Geary autocorrelation as given in Li (2006) PROFEAT.

    TODO formula

    >>> autocorrelation_geary('ABACABAC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    0.0
    >>> autocorrelation_geary('ABACABAC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 2)
    1.75
    >>> autocorrelation_geary('BBBBCCCC', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    1.75
    >>> autocorrelation_geary('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 4)
    0.0
    >>> autocorrelation_geary('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 8)
    0.0
    >>> autocorrelation_geary('BBBBBBBB', {'A': 0.0, 'B': 1.0, 'C': -1.0}, 12)
    0.0
    '''

    if(lag < 1):
        raise ValueError('The provided lag should be a positive integer.')

    # TODO check if this is the correct approach
    if(len(sequence) <= lag):
        return 0.0

    # transform sequence to signal using the provided scale
    signal = numpy.array(seq_signal_raw(sequence, scale))

    # calculate average signal value
    avg_signal = numpy.mean(signal)

    #
    denom = sum(numpy.square(signal[:-lag] - signal[lag:]))

    norm_denom = denom / (2.0 * (len(sequence) - lag))

    #
    enume = sum(numpy.square(signal - avg_signal))

    norm_enume = enume / (len(sequence) - 1.0)

    if(norm_enume == 0):
        return 0.0  # TODO should this be zero???
    else:
        return norm_denom / norm_enume


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

    'Dan: First, last occurence ? '

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

                for i in xrange(letter_positions.count(letter_count)):
                    seq_fractions.append((index + 1.0) / len(seq))

    return seq_fractions


def seq_order_coupling_number(seq, r, aa_dists):
    '''
    This function returns the r-th rank sequence order coupling number for
    sequence seq using aa_dist_matrix as amino acid distance function.

    >>> dist_m = {'AA': 0.3, 'AB': 0.6, 'BA': 0.9, 'BB': 0.5}
    >>> d = 2
    >>> s = 'AABABBAAAB'
    >>> so = seq_order_coupling_number(s, d, dist_m)
    >>> round(so, 2)
    3.13

    '''
    return sum([math.pow(aa_dists[p], 2) for p in ordered_seq_pairs(seq, r)])


def quasi_sequence_order_descriptors(seq, aa_dists, rank, weight=0.1):
    '''
    This

    Number returned values is rank + 20
    '''

    if(rank <= 0):
        raise ValueError('The max rank should be larger than 0.')

    # sequence order coupling numbers 1 ... rank
    socns = [seq_order_coupling_number(seq, r, aa_dists)
             for r in xrange(1, rank + 1)]

    # denominator terms
    sum_aa_freqs = 1.0
    sum_socns = sum(socns)

    # denominator
    denominator = sum_aa_freqs + weight * sum_socns

    # type-1 quasi-sequence-order (amino acid composition / denominator)
    aacomp = letter_composition(seq, aa_unambiguous_alph)
    type1_qsod = aacomp / denominator

    # type-2 quasi-sequence-order
    type2_qsod = (weight * numpy.array(socns)) / denominator

    return numpy.concatenate((type1_qsod, type2_qsod))


def sequence_order_correlated_factors(seq, r, aa_corr):
    '''
    NOTE: undesired result if len(seq) <= r. For now zero is returned in this
    case. Need to check if this is the correct behaviour.

    NOTE: Correlation value zero is used as default, in case of ambiguous amino
    acids occur.
    '''
    if(len(seq) < r):
        return 0.0
    else:
        corr_values = [aa_corr.get(p, 0.0) for p in ordered_seq_pairs(seq, r)]
        return sum(corr_values) / (len(seq) - r)


def pseaac_type1(seq, aa_scales, lambda_, weight=0.05):
    '''
    This function calculates the pseudo amino acid composition for amino acid
    sequence seq using the provided liste of amino acid scales and lambda as
    maximal distance parameter.

    Raises:
        ValueError: if one of the provided aa_scales is not in the list of
                    pseaac scales
        ValueError: if lambda is zero or smaller.
    '''

    if(lambda_ <= 0):
        raise ValueError('The max rank should be larger than 0.')

    # TODO move to correlation matrix function???
    if(all(type(a) == str for a in aa_scales)):
        aa_scales = [pseaac_scale_ids.index(a) for a in aa_scales]
    elif not(all([type(a) == int and a in range(len(pseaac_scale_ids))
                  for a in aa_scales])):
        raise ValueError('Wrong amino acid scale index provided.')

    aa_scales = tuple(sorted(aa_scales))

    aa_corr = pseaac_correlation_matrix(aa_scales)

    # sequence order correlated factors
    socfs = [sequence_order_correlated_factors(seq, l, aa_corr)
             for l in xrange(1, lambda_ + 1)]

    # denominator terms
    sum_aa_freqs = 1.0
    sum_socfs = sum(socfs)

    # denominator
    denominator = sum_aa_freqs + weight * sum_socfs

    # first 20 features, aa composition
    aacomp = letter_composition(seq, aa_unambiguous_alph)
    first_20_aacomp = aacomp / denominator

    # order stuff
    rest = (weight * numpy.array(socfs)) / denominator

    return numpy.concatenate((first_20_aacomp, rest))


def pseaac_type2(seq, aa_scales, lambda_, weight=0.05):

    if(lambda_ <= 0):
        raise ValueError('The max rank should be larger than 0.')

    # TODO move to correlation matrix function???
    if(all(type(a) == str for a in aa_scales)):
        aa_scales = [pseaac_scale_ids.index(a) for a in aa_scales]
    elif not(all([type(a) == int and a in range(len(pseaac_scale_ids))
                  for a in aa_scales])):
        raise ValueError('Wrong amino acid scale index provided.')

    aa_cors = [pseaac_type2_correlation_matrix(i) for i in aa_scales]

    # sequence order correlated factors
    socfs = [sequence_order_correlated_factors(seq, l, aa_cor)
             for l in xrange(1, lambda_ + 1)
             for aa_cor in aa_cors]

    # denominator terms
    sum_aa_freqs = 1.0
    sum_socfs = sum(socfs)

    # denominator
    denominator = sum_aa_freqs + weight * sum_socfs

    # first 20 features, aa composition
    aacomp = letter_composition(seq, aa_unambiguous_alph)
    first_20_aacomp = aacomp / denominator

    # order stuff
    rest = (weight * numpy.array(socfs)) / denominator

    return numpy.concatenate((first_20_aacomp, rest))


def state_subseq(seq, state_seq, state_letter):
    '''
    This function returns the returns those parts of seq where the
    state_seq letter is equal to state_letter. The subparts are glued together
    and returned as a single string.

    >>> state_subseq('ABCDEFGHIJKLMNO', 'AAABBBCCCAAABBB', 'B')
    'DEFMNO'
    >>> state_subseq('ABCDEFGHIJKLMNO', 'AAABBBCCCAAABBB', 'D')
    ''
    >>> state_subseq('ABCDEFGHIJKLMNO', 'AAAAAACCCAAACCC', 'B')
    ''
    >>> state_subseq('ABCDEFGHIJKLMNO', 'AAAAAACCCAAACCC', '')
    ''
    >>> state_subseq('', '', 'A')
    ''
    '''
    if not(len(seq) == len(state_seq)):
        raise ValueError('The state_seq should have the same length as seq.')

    return ''.join([i for i, j in zip(seq, state_seq) if j == state_letter])


def state_subseq_composition(seq, state_seq, seq_alph, state_alph):
    '''
    This function returns the seq_alph composition of seq, but only for the
    parts where the state_seq has a letter that is in state_alph.

    >>> s = 'SSSSSTTTTTSSSS'
    >>> t = 'AAABBCCAAAAAAA'
    >>> comp = state_subseq_composition(s, t, 'ST', 'A')
    >>> print(round(comp[0], 1))
    0.7
    >>> print(round(comp[1], 1))
    0.3
    '''
    result = []
    # TODO is extend correct??? Now all composition are added to the same list
    for l in state_alph:
        comp = letter_composition(state_subseq(seq, state_seq, l), seq_alph)
        result.extend(comp)
    return result


# TODO add cluster(s) as parameter... instead of aa_subsets
def aa_cluster_composition(protein):
    '''
    This function returns the protein sequence composition of the different
    defined amino acid clusters.
    '''
    comp = dict(zip(aa_unambiguous_alph,
                    letter_composition(protein, aa_unambiguous_alph)))

    return numpy.array([sum([comp[l] for l in aa_subset_dict[subset]])
                       for subset in aa_subsets])


def window_seq(seq, window_size, overlapping=False):
    '''
    This function returns a chopped version of seq, in which it is chopped in
    subsequences of length window_size. By default, non-overlapping
    subsequences are returned, if overlapping is set to True, overlapping sub-
    sequences are returned.

    IMPORTANT: If the length of the sequence is not a multiple of the window
    size, the last letters are NOT returned.

    >>> s = 'ACCACCAAAA'
    >>> window_seq(s, 3)
    ['ACC', 'ACC', 'AAA']
    >>> window_seq(s, 3, overlapping=True)
    ['ACC', 'CCA', 'CAC', 'ACC', 'CCA', 'CAA', 'AAA', 'AAA']
    >>> window_seq(s, 1)
    ['A', 'C', 'C', 'A', 'C', 'C', 'A', 'A', 'A', 'A']
    >>> window_seq(s, 10)
    ['ACCACCAAAA']
    >>> window_seq(s, 11)
    []
    '''

    if(window_size < 2):
        return list(seq)
    else:
        start = 0
        stop = len(seq) - window_size + 1
        step = window_size
        if(overlapping):
            step = 1
        return [seq[i:i + window_size] for i in range(start, stop, step)]


@memoize
def convolution_filter(window, edge):
    '''
    This function returns a triangular convolution filter. The filter values
    add up to 1.0.

    Args:
        window (int): The width of the filter
        edge (float): The weight of the edges of the window [0.0, 1.0]
        NOTE edge is now 0..100 (in percentage)
    Raises:
        ValueError: if the window is not an uneven number.
        ValueError: if the window is too small, smaller than 3.
        ValueError: if the edge parameter is out of range [0.0, 1.0].

    >>> convolution_filter(9, 0.0)
    array([ 0.    ,  0.0625,  0.125 ,  0.1875,  0.25  ,  0.1875,  0.125 ,
            0.0625,  0.    ])
    >>> convolution_filter(3, 33.333333)
    array([ 0.2,  0.6,  0.2])
    '''

    if(window % 2 == 0):
        raise ValueError('Window must be an uneven number.')
    if(window < 3):
        raise ValueError('Window must be 3 or larger.')
    if(edge < 0 or edge > 100):
        raise ValueError('The edge parameter must be in the range 0 to 100.')

    if(edge == 100):
        result = numpy.ones(window) / window
        return result
    else:
        result = numpy.ones(window)
        num = window / 2
        forw = numpy.linspace(edge / 100.0, 1.0, num, endpoint=False)
        result[:num] = forw
        result[-num:] = forw[::-1]
        result = result / result.sum()
        return result


def seq_signal_raw(sequence, scale):
    '''
    This function maps the sequence to a raw value signal using the provided
    alphabet scale.

    Args:
        sequence (str): The sequence to be mapped.
        scale (dict(str --> float)): A mapping from each letter in the sequence
                                     alphabet to a float value.
    Raises:
        KeyError: if sequence contains letters that have no scale value.

    >>> s = 'AABBAA'
    >>> sc = {'A': -1.0, 'B': 1.0}
    >>> seq_signal_raw(s, sc)
    [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0]
    '''
    return [scale[letter] for letter in sequence]


def seq_signal(sequence, scale, window=9, edge=0):
    '''
    This function returns a smoothed sequence signal using the provided letter
    to value scale and a triangular smoothing filter with the provided window
    width and edge weights.

    Args:
        sequence (str):
        scale (float):
        window (int):
        edge (float):
    Raises:
        ValueError: if window is too small or large.

    >>> s = 'AABBAA'
    >>> sc = {'A': -1.0, 'B': 1.0}
    >>> seq_signal(s, sc, window=3, edge=100)
    array([-0.33333333,  0.33333333,  0.33333333, -0.33333333])
    '''

    if(window > len(sequence) or window < 1):
        raise ValueError('1 <= window <= sequence length.')

    # obtain raw signal
    signal = seq_signal_raw(sequence, scale)

    # return the raw signal if window size is one
    if(window == 1):
        return signal

    # otherwise return convolved signal
    else:
        conv = convolution_filter(window, edge)
        return numpy.convolve(signal, conv, 'valid')


def avg_seq_signal(sequence, scale, window=9, edge=0):
    '''
    This function returns the average value of the smoothed sequence signal
    that is constructed using scale and a triangular filter with width window
    and edge weights.
    '''
    sig = seq_signal(sequence, scale, window, edge)
    return sum(sig) / len(sequence)


def auc_seq_signal(sequence, scale, window=9, edge=0, threshold=10):
    '''
    This function returns sequence signal area above and underneeth the
    specified threshold, normalized by the sequence length.

    The most basic area estimation is used.

    >>> s = 'AAABBBCCCAAA'
    >>> sc = {'A': 0.0, 'B': 1.0, 'C': -1.0}
    >>> auc_seq_signal(s, sc, window=3, edge=0, threshold=5)
    (0.125, 0.125)
    '''

    # HACK for spice website...
    threshold = threshold / 10.0

    area_above = 0.0
    area_below = 0.0

    sig = seq_signal(sequence, scale, window, edge)
    for value in sig:
        if value > threshold:
            area_above += value - threshold
        elif value < -1.0 * threshold:
            area_below += -1.0 * value - threshold
        else:
            pass

    return (area_above / len(sequence), area_below / len(sequence))


def codon_count(orf):
    '''
    This function returns the codon counts in the orf sequence.

    Args:
        orf (str): Open reading frame (nucleotide) sequence

    >>> cc = codon_count('TTTTTCTTTTTC')
    >>> cc[0]
    2
    >>> cc[1]
    2
    >>> all([c == 0 for c in cc[2:]])
    True
    '''

    if not(len(orf) % 3 == 0):
        raise ValueError('ORF sequence length is not a multiple of 3.')
    if not(is_nucleotide_sequence(orf)):
        raise ValueError('ORF sequence is not a nucleotide sequence.')

    wseq = window_seq(orf, 3, overlapping=False)
    return numpy.array([wseq.count(c) for c in codons_unambiguous])


def codon_composition(orf):
    '''
    This function returns the codon composition of the provided orf sequence.

    Args:
        orf (str): Open reading frame (nucleotide) sequence

    >>> cc = codon_composition('TTTTTCTTTTTC')
    >>> cc[0]
    0.5
    >>> cc[1]
    0.5
    >>> all([c == 0.0 for c in cc[2:]])
    True
    '''
    return codon_count(orf) / float(len(orf) / 3)


def codon_usage(orf):
    '''

    >>> cu = codon_usage('TTTTTTTTTTTC')
    >>> cu[0]
    0.75
    >>> cu[1]
    0.25
    >>> all([c == 0.0 for c in cu[2:]])
    True
    '''

    # TODO leave out codons that encode for only one amino acid?
    # What about start and stop codon?

    # count amino acids for translated orf sequence (dict: aa --> count)
    aa_c_dict = dict(zip(aa_unambiguous_alph,
                         letter_count(translate(orf), aa_unambiguous_alph)))

    # turn into list with amino acid count per codon
    # change 0 to 1, to prevent / 0 (answer will still always be 0 (0/1))
    aa_c = numpy.array([float(aa_c_dict.get(a, 0)) if aa_c_dict.get(a, 0) > 0
                       else 1.0 for a in codon_aas])

    # get the codon counts
    codon_c = codon_count(orf)

    # divide codon count by corresponding amino acid count
    return codon_c / aa_c


###############################################################################
#
# SEQUENCE CHECKS
#
###############################################################################


def is_amino_acid_sequence(sequence):
    return set(sequence).issubset(set(aa_alph))


def is_not_an_amino_acid_sequence(sequence):
    return not is_amino_acid_sequence(sequence)


def is_unambiguous_amino_acid_sequence(sequence):
    return set(sequence).issubset(set(aa_unambiguous_alph))


def is_not_an_unambiguous_amino_acid_sequence(sequence):
    return not is_unambiguous_amino_acid_sequence(sequence)


def is_nucleotide_sequence(sequence):
    return set(sequence).issubset(set(nucleotide_alph))


def is_not_a_nucleotide_sequence(sequence):
    return not is_nucleotide_sequence(sequence)


def is_unambiguous_nucleotide_sequence(sequence):
    return set(sequence).issubset(set(nucleotide_unambiguous_alph))


def is_not_an_unambiguous_nucleotide_sequence(sequence):
    return not is_unambiguous_nucleotide_sequence(sequence)


def is_sec_struct_sequence(sequence):
    return set(sequence).issubset(set(ss_alph))


def is_not_a_sec_struct_sequence(sequence):
    return not is_sec_struct_sequence(sequence)


def is_solv_access_sequence(sequence):
    return set(sequence).issubset(set(sa_alph))


def is_not_a_solv_access_sequence(sequence):
    return not is_solv_access_sequence(sequence)


def probably_nucleotide(sequence):
    pass


def is_empty(sequence):
    return len(sequence) == 0


###############################################################################
#
# Sequence properties TODO order and document this...
#
###############################################################################


def hamming_distance(s0, s1):
    '''
    Returns the Hamming distance between the two equal lengths sequences s0 and
    s1. It returns the sum of the pairwise character distances. Distance
    between two characters is zero if the two characters are the same, one
    otherwise.

    Args:
        s0 (str):
        s1 (str):

    >>> hamming_distance('AAA', 'AAA')
    0
    >>> hamming_distance('AAA', 'AAB')
    1
    >>> hamming_distance('AAA', 'BBB')
    3
    '''
    assert(len(s0) == len(s1))
    return sum([not s0[i] == s1[i] for i in range(len(s0))])


def dist_one_codons(codon):
    '''
    This function returns all (unambiguous) codons that have a Hamming distance
    of one to the provided (unambiguous) codon.

    Args:
        codon (str): An unambiguous codon.
    Raises:
        ValueError: TODO

    >>> sorted(dist_one_codons('AAA'))
    ['AAC', 'AAG', 'AAT', 'ACA', 'AGA', 'ATA', 'CAA', 'GAA', 'TAA']
    '''
    return [c for c in codons_unambiguous if hamming_distance(codon, c) == 1]


def dist_one_amino_acids(codon):
    '''
    This function returns all amino acids for which one of their codons is one
    Hamming distance away of the provided codon. This means that all amino
    acids are returned that can be obtained by a single mutation of the codon.

    Args:
        codon (str):
    Raises:
        TODO

    >>> dist_one_amino_acids('AAA')
    ['*', 'E', 'I', 'K', 'N', 'Q', 'R', 'T']
    '''
    codons = dist_one_codons(codon)
    return sorted(set([codon_table_unambiguous[c] for c in codons]))


def aa_substitutions():
    tuples = [(l0, l1) for l0 in aa_unambiguous_alph
              for l1 in aa_unambiguous_alph if not l0 == l1]
    return tuples


def mutations():
    return [(l0, l1) for l0 in nucleotide_unambiguous_alph
            for l1 in nucleotide_unambiguous_alph if not l0 == l1]


def codon_mutations():
    result = []
    for codon in codons_unambiguous:
        result.extend(codon_mutation(codon))
    return result


def codon_mutation(codon):
    assert(codon in codons_unambiguous)
    return [(codon, c1) for c1 in codons_unambiguous
            if hamming_distance(codon, c1) == 1]


def single_mutation_aa_substitutions():
    return [(codon_table_unambiguous[c[0]], codon_table_unambiguous[c[1]])
            for c in codon_mutations()]


def possible_single_mutation_aa_substitutions():
    return sorted(set([s for s in single_mutation_aa_substitutions()
                  if not(s[0] == s[1] or (s[0] == '*' or s[1] == '*'))]))


def impossible_single_mutation_aa_substitutions():
    return sorted(set(aa_substitutions()) -
                  set(possible_single_mutation_aa_substitutions()))


def single_mutation_aa_substitution_dict():
    result = {}
    for fr, to in possible_single_mutation_aa_substitutions():
        result[fr] = result.setdefault(fr, '') + to
    return result


def single_mutation_aa_substitution_stats():
    subs = single_mutation_aa_substitutions()
    no_sub = 0
    stop_sub = 0
    aa_subs = []
    for s in subs:
        if(s[0] == s[1]):
            no_sub += 1
        elif(s[0] == '*' or s[1] == '*'):
            stop_sub += 1
        else:
            aa_subs.append(s)

    print ('No substitution: %i') % (no_sub)
    print ('Stop codon substitution: %i') % (stop_sub)
    print ('Amino acid substitution: %i') % (len(aa_subs))
    print ('TOTAL: %i' % (len(subs)))

    possible_subs = set(aa_subs)
    impossible_subs = set(aa_substitutions()) - possible_subs
    print ('Possible amino acid substitutions: %i') % (len(possible_subs))
    print ('Impossible amino acid substitutions: %i') % (len(impossible_subs))
    print ('TOTAL: %i') % (len(aa_substitutions()))
    #print impossible_subs
    print ('')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
