# Copyright Yair Benita Y.Benita@pharm.uu.nl
# Biopython (http://biopython.org) license applies

"""Calculate isoelectric points of polypeptides using methods of Bjellqvist.

pK values and the methos are taken from:

* Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F., Sanchez,
J.-Ch., Frutiger, S. & Hochstrasser, D.F.
The focusing positions of polypeptides in immobilized pH gradients can be predicted
from their amino acid sequences. Electrophoresis 1993, 14, 1023-1031.

* Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
Reference points for comparisons of two-dimensional maps of proteins from
different human cell types defined in a pH scale where isoelectric points correlate
with polypeptide compositions. Electrophoresis 1994, 15, 529-539.

I designed the algorithm according to a note by David L. Tabb, available at:
http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf

"""

positive_pKs = {'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H': 5.98}
negative_pKs = {'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0}
pKcterminal = {'D': 4.55, 'E': 4.75}
pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.7}
charged_aas = ('K', 'R', 'H', 'D', 'E', 'C', 'Y')


# access this module through ProtParam.ProteinAnalysis class.
# first make a ProteinAnalysis object and then call its isoelectric_point method.
class IsoelectricPoint(object):
    def __init__(self, ProteinSequence, AminoAcidsContent):
        self.sequence = ProteinSequence
        self.charged_aas_content = self._select_charged(AminoAcidsContent)

    # This function creates a dictionary with the contents of each charged aa,
    # plus Cterm and Nterm.
    def _select_charged(self, AminoAcidsContent):
        charged = {}
        for aa in charged_aas:
            charged[aa] = float(AminoAcidsContent[aa])
        charged['Nterm'] = 1.0
        charged['Cterm'] = 1.0
        return charged

    #This function calculates the total charge of the protein at a given pH.
    def chargeR(self, pH, pos_pKs, neg_pKs):
        PositiveCharge = 0.0
        for aa, pK in pos_pKs.items():
            CR = 10**(pK-pH)
            partial_charge = CR/(CR+1.0)
            PositiveCharge += self.charged_aas_content[aa] * partial_charge

        NegativeCharge = 0.0
        for aa, pK in neg_pKs.items():
            CR = 10**(pH-pK)
            partial_charge = CR/(CR+1.0)
            NegativeCharge += self.charged_aas_content[aa] * partial_charge

        return PositiveCharge - NegativeCharge

    # This is the action function, it tries different pH until the charge of the protein is 0 (or close).
    def pi(self):
        pos_pKs = dict(positive_pKs)
        neg_pKs = dict(negative_pKs)
        nterm = self.sequence[0]
        cterm = self.sequence[-1]
        if nterm in pKnterminal:
            pos_pKs['Nterm'] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs['Cterm'] = pKcterminal[cterm]

        # Bracket between pH1 and pH2
        pH = 7.0
        Charge = self.chargeR(pH, pos_pKs, neg_pKs)
        if Charge > 0.0:
            pH1 = pH
            Charge1 = Charge
            while Charge1 > 0.0:
                pH = pH1 + 1.0
                Charge = self.chargeR(pH, pos_pKs, neg_pKs)
                if Charge > 0.0:
                    pH1 = pH
                    Charge1 = Charge
                else:
                    pH2 = pH
                    Charge2 = Charge
                    break
        else:
            pH2 = pH
            Charge2 = Charge
            while Charge2 < 0.0:
                pH = pH2 - 1.0
                Charge = self.chargeR(pH, pos_pKs, neg_pKs)
                if Charge < 0.0:
                    pH2 = pH
                    Charge2 = Charge
                else:
                    pH1 = pH
                    Charge1 = Charge
                    break

        # Bisection
        while pH2 - pH1 > 0.0001 and Charge != 0.0:
            pH = (pH1 + pH2) / 2.0
            Charge = self.chargeR(pH, pos_pKs, neg_pKs)
            if Charge > 0.0:
                pH1 = pH
                Charge1 = Charge
            else:
                pH2 = pH
                Charge2 = Charge

        return pH
#############################################################################

'Pyteomics package:'

# from . import parser #Pyteomics

pK_lehninger = {
    'E': [(4.25, -1),],
    'R': [(12.48, +1),],
    'Y': [(10.07, -1),],
    'D': [(3.65, -1),],
    'H': [(6.00, +1),],
    'K': [(10.53, +1),],
    'C': [(8.18, -1),],
    'H-': [(9.69, +1),],
    '-OH':  [(2.34, -1),],
    }

def charge(sequence, pH, **kwargs):
    """Calculate the charge of a polypeptide in given pH or list of pHs using
    a given list of amino acid electrochemical properties.

    .. warning::

        Be cafeful when supplying a list with a parsed sequence or a dict with
        amino acid composition as `sequence`. Such values must be obtained
        with enabled `show_unmodified_termini` option.

    Parameters
    ----------
    sequence : str or list or dict
        A string with a polypeptide sequence, a list with a parsed
        sequence or a dict of amino acid composition.
    pH : float or list of floats
        pH or list of pHs for which the charge is calculated.
    pK : dict {str: [(float, int),]}, optional
        A set of pK of amino acids' ionizable groups. It is a dict, where keys
        are amino acid labels and the values are lists of tuples (pK,
        charge_in_ionized_state), a tuple per ionizable group. The default
        value is `pK_lehninger`.

    Returns
    -------
    out : float or list of floats or None
        A single value of charge or a list of charges. Returns None if
        `sequence` is not of supported type.
    """

    # Get the list of valid modX labels.
    pK = kwargs.get('pK', pK_lehninger)
    labels = list(parser.std_labels)
    for label in pK:
        if label not in labels:
            labels.append(label)

    # Parse the sequence.
    if isinstance(sequence, str) or isinstance(sequence, list):
        peptide_dict = parser.amino_acid_composition(sequence, True, False,
                                                     labels=labels)
    elif isinstance(sequence, dict):
        peptide_dict = sequence
    else:
        raise PyteomicsError('Unsupported type of sequence: %s' % type(sequence))

    # Check if a sequence was parsed with `show_unmodified_termini` enabled.
    num_term_mod = 0
    for aa in peptide_dict:
        if parser.is_term_mod(aa):
            num_term_mod += 1
    if num_term_mod != 2:
        raise PyteomicsError('Parsed sequences must contain unmodified termini.')

    # Process the case when pH is a single float.
    pH_list = pH if isinstance(pH, list) else [pH,]

    # Calculate the charge for each value of pH.
    charge_list = []
    for pH_value in pH_list:
        charge = 0
        for aa in peptide_dict:
            for ionizable_group in pK.get(aa, []):
                charge += peptide_dict[aa] * ionizable_group[1] * (
                    1.0
                    / (1.0 + 10 ** (ionizable_group[1]
                                  * (pH_value - ionizable_group[0]))))
        charge_list.append(charge)

    return charge_list[0] if len(charge_list) == 1 else charge_list

protein="MQNEEDACLEAGYCLGTTLSSWRLHFMEEQSQSTMLMGIGIGALLTLAFVGIFFFVYRRVRRLRRAEPTPQYRFRKRDKVMFYGRKIMRKVTTLPHTLVGNTSAPRQRVRKRTKVLSLAKRILRFKKEYPTLQPKEPPPSLLEADLTEFDVKNSHLPSEVLYMLKNVRVLGHFEKPLFLELCKHMVFVQLQEGEHVFQPGEPDISIYVVQDGRLEVCIQDADGTEVVVKEVLPGDSVHSLLSILDVITGHTAPYKTVSARAAVSSTVLWLPAAAFQGVFEKYPETLVRVVQIIMVRLQRVTFLALHNYLGLTTELFNPESQAIPLLSVASVAGRAKRQMSYGPEEQLERSLRPSEFSSSDHGSSCVTVSGPLLKRSCSVPLPSNHGEVDELRQSQGSGSNTSAFQESHEGATSDLGMAYNRARILPHSDEQLGNSLASKSKKSVVAETPSAIFHYSENFRDETGACGKTDAIFRAATKDLLTLMKLDDPSLLDGRVAFLHVPAGTLVSKQGDQDVNILFVVSGMLHVYQQKIDSLEDTCLFLTHPGEMVGQLAVLTGEPLMFTIRANRDCSFLSISKAHFYEIMRKRPDVVLGVAHTVVKRMSSFVRQIDFALDWMEVEAGRAIYRQGDKSDCTYIVLSGRLRSVIRKDDGKKRLAGEYGRGDLVGVVETLTHQARATTVHAVRDSELAKLPAGALTSIKRRYPQVVTRLIHLLGEKILGSLQQGSATGHQLGFNTASSKWDLGNPPGNLSTVAALPASEDVPLTAFALELQHALSAIGPVLLLTSDNIKQRLGSAALDSIHEYRLSSWLGQQEDIHRIVLYQADGTLTPWTQRCIRQADCILIVGLGEQEPAVGELEQMLESTAVRAQKQLILLHKEDGPVPSRTVEWLNMRSWCSGHLHLCCPRRVFSKRSLPKLVEMYTRVFQRPPDRHSDFSRLARMLTGNAIALVLGGGGARGCAQVGILRALAECGVPVDIIGGTSIGAFMGALFAEERSYSQTRIRAKQWAEGMTSMMKTILDLTYPITSMFSGTGFNSSISNIFKDRQIEDLWLPYFAITTDITASAMRVHTDGSLWRYVRASMSLSGYMPPLCDPKDGHLLMDGGYINNLPADVARSMGAKVVIAIDVGSRDETDLTNYGDALSGWWLLWKRWNPLATKVKVLNMAEIQTRLAYVCCVRQLEMVKNSDYCEYLRPPIDSYRTLDFGKFDEICEVGYQHGRTVFDIWVRSGVLEKMLQDQQGTSKRKDCGVFTCPNSSFTDLAEIVSRIEPAKVAAVDDESDYQTEYEEELPAIPKETYADFQSTGIELDSDSEYEPSMLQGPPSLTSPEQSQDSFPWLPNQDDQGPRLEHPS"
print(charge(protein,7.2))
