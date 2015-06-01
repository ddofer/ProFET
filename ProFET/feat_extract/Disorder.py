__author__ = 'DanaLab'
'''
Look at using:
Get_ParamScales()  from protfeat - and import scales from AAScales.py?
(Also, calc. normalized KD scale, and save it (modify file) to AAScales.py and import from there
= performance.
Also, AAScales should hold (import) the TDP-IDP scale.  - Dan. )

netCharge; calculateAminoAcidCharge - why not use built in ones from main ProtFeat?
 (also - names of methods here/there are SAME!
    => asking for bugs when importing, calling methods..) | Change method names.

Look at using different PHs for netcharge calcing. (This would be a Different feature of
course; i.e diff key name in res-dict)
'''
from collections import Counter
# from ProtFeat import pKa

pKa     = {'D':3.9, 'E':4.3, 'H':6.1, 'C':8.3, 'Y':10.1, 'K':10.5, 'R':12, 'N-term':8, 'C-term':3.1}
charges = {'D':-1,  'E':-1,  'H':+1,  'C':-1,  'Y':-1,   'K':1,    'R':1,  'N-term':1, 'C-term':-1}


# @staticmethod
def netCharge(seq,pH = 7.2): #maybe ReName, to "subseq_ .." , to avoid confusion with "calculateProteinCharge", get_netCharge, From ProtFeat.py ? (OR use them directly) - D
    """

    :param seq:
    :return:
    """
    aa_counts = Counter(seq)
    # pH = 7.2
    res = 0.0

    def calculateAminoAcidCharge(amino_acid, pH):
        ratio = 1 / (1 + 10 ** (pH - pKa[amino_acid]))
        if charges[amino_acid] == 1:
            return ratio
        else:
            return ratio - 1

    for amino_acid in pKa:
        res += aa_counts[amino_acid] * calculateAminoAcidCharge(amino_acid, pH)
    return res


# @staticmethod
def hydrophobicity(seq):
    hydropathy = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
          'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'U': 0.0, 'V': 4.2, 'W': -0.9,
          'Y': -1.3, 'B': -3.5, 'X': -0.49, 'Z': -3.5}
    WINDOW_SIZE = 5
    NORME = True  # NORME = ?

    def normalizeHydropathy():
        """


        """
        minimum = min(hydropathy.values())
        maximum = max(hydropathy.values())
        for key in hydropathy.keys():
            oldVal = hydropathy[key]
            hydropathy[key] = (oldVal - minimum) / (maximum - minimum)

    def kyteDoolittle(seq, windowSize, normaliz):
        """

        :param seq:
        """
        seq = seq.strip()
        if normaliz:
            normalizeHydropathy()
        maxJump = int((windowSize - 1) / 2) #MOD
        result = 0
        subResultsArr = []
        for i in range(maxJump):
            subResultsArr.append(0)
        for i in range(maxJump, len(seq) - maxJump):
            summ = 0
            for j in range(-maxJump, maxJump + 1):
                key = seq[i + j]
                if key in hydropathy.keys():
                    summ += hydropathy[key]
                else:
                    print(key)
            subResultsArr.append(summ / windowSize)
            result += summ / windowSize

        result /= len(seq)
        return result, subResultsArr

    return kyteDoolittle(seq, WINDOW_SIZE, NORME)[0]


def uversky(seq): #As implemented, this gets the foldindex for WHOLE seq; vs segments/window.
    '''
    FoldIndex method prediction of disorder.
    '''
    #Why use sep. Seq? Use seq=self.seq for consistancy/less confusion, ne?
    R = netCharge(seq) #Maybe have this for different PHs.
    H = hydrophobicity(seq)
    uScore = (2.785 * float(H) - 1.151 - float(R))
    return uScore

def getDisordered(seq,segments=5):
    '''
    Get predicted disorder for protein, divided into segments (default=5).,
    predicted individually for entirety of each segment; using:
    A) FoldIndex (Uversky) method.

    #I'd add other method(s) for getting predicted disordered here also. (EG, TDP-IDP scale, whether seperate or "joint"feature) - D
    '''

    # seq = self.seq
    length = len(seq)
    window_size = int(length / segments)  # window size 20% of the protein length
    pos = 0
    scores = [0 for _ in range(segments)]

    for i in range(segments-1):
        scores[i] = uversky(seq[pos:pos + window_size])
        pos += window_size
    scores[-1] = uversky(seq[pos:])
    res = {}
    key = "disordered_window_"
    for i, uscore in enumerate(scores):
        res[key + str(i)] = 1 if uscore < 0 else 0 #ORIG
        # res[key + str(i)] = 1 if uscore != 0 else 0
        'Binary feature - presence of ANY disordered window:'
        if uscore < (-0.1):
            res['anyDISORDER_'+str(segments)]=1

    return res