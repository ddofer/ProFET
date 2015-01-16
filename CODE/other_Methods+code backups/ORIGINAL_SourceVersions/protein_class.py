import math

import numpy

from biopy import sequtil


class Protein(object):

    def __init__(self, pid):

        self.pid = pid

        self.missense_mutations = []

        self.orf_sequence = None
        self.protein_sequence = None
        self.ss_sequence = None
        self.sa_sequence = None
        self.protein_structure = None

        self.rasa = None

        # TODO depricate
        # rank score per residue
        self.msa_residue_rank = None
        #self.msa_variability = None
        self.msa_coverage = None

        # updated own hhblits MSA, list of aligned sequences, first is this seq
        self.msa = None

        self.pfam_annotations = None
        self.backbone_dynamics = None

        # protein interaction counts (maybe combine in tuple, index as consts)
        self.ppi_count = None
        self.metabolic_count = None
        self.genetic_count = None
        self.phosphorylation_count = None
        self.regulatory_count = None
        self.signaling_count = None


    # set attributes, sequence data
    # TODO turn this into proper setters...


    def set_protein_sequence(self, seq):
        self.protein_sequence = seq

    def set_protein_structure(self, struct):
        self.protein_structure = struct

    def set_ss_sequence(self, seq):
        self.ss_sequence = seq

    def set_sa_sequence(self, seq):
        self.sa_sequence = seq

    # TODO depricate
    def set_msa_data(self, msa_data):

        # 'unzip' lists
        if not(msa_data is None):
            i0, i1, r, cov, v0, v1, rank = zip(*msa_data)

            # check if sequence corresponds to protein sequence
            #assert(''.join(r) == self.protein_sequence)
            # I currently allow for at most 5 residue differences...
            assert(len(r) == len(self.protein_sequence))
        else:
            # TODO default values!!! check this!
            cov = [0.0] * len(self.protein_sequence)
            v1 = [[]] * len(self.protein_sequence)  # not sure about this...
            rank = [0.0] * len(self.protein_sequence)

        # store coverage, variability, and rank score
        self.msa_coverage = cov
        self.msa_variability = v1
        self.msa_residue_rank = rank

    def set_msa(self, msa):
        '''
        msa is list with (equal length) aligned sequences. The first sequence
        is the sequence of the protein (without gaps).
        '''

        if(msa is None):
            # if not available, use own sequence as only sequence
            msa = [self.protein_sequence]

        else:
            # checks
            if not(msa[0] == self.protein_sequence):
                raise ValueError('First sequence in MSA does not correspond ' +
                                 'to this protein sequence')
            if not(all([len(msa[0]) == len(m) for m in msa])):
                raise ValueError('Not all sequences in MSA have the same ' +
                                 'length')

        # store data
        self.msa = msa

    def set_rasa(self, rasa):
        assert(rasa is None or type(rasa) == list)
        self.rasa = rasa

    def set_pfam_annotations(self, pfam_annotations):
        self.pfam_annotations = [Pfam(a[0], a[1], a[2], a[3], a[4], a[5], a[6],
                                 a[7], a[8]) for a in pfam_annotations]

    def set_backbone_dynamics(self, backbone_dynamics):
        assert(type(backbone_dynamics) == list)
        assert(len(backbone_dynamics) == len(self.protein_sequence))
        self.backbone_dynamics = backbone_dynamics

    def set_interaction_counts(self, interaction_counts):
        self.ppi_count = interaction_counts[0]
        self.metabolic_count = interaction_counts[1]
        self.genetic_count = interaction_counts[2]
        self.phosphorylation_count = interaction_counts[3]
        self.regulatory_count = interaction_counts[4]
        self.signaling_count = interaction_counts[5]

    ###########################################################################
    # feature calculation functions
    ###########################################################################

    def amino_acid_composition(self, num_segments, feature_ids=False):

        alph = sequtil.aa_unambiguous_alph

        if not(feature_ids):

            seq = self.protein_sequence

            if(num_segments == 1):
                return sequtil.letter_composition(seq, alph)
            else:
                seqs = sequtil.segment(seq, num_segments)
                return numpy.concatenate(
                    [sequtil.letter_composition(s, alph) for s in seqs])
        else:

            feat_ids = []
            feat_names = []

            for si in xrange(1, num_segments + 1):
                for aa in alph:
                    feat_ids.append('%s%i' % (aa, si))
                    feat_names.append('%s, segment %i' % (aa, si))

            return (feat_ids, feat_names)

    def dipeptide_composition(self, num_segments, feature_ids=False):

        alph = sequtil.aa_unambiguous_alph

        if not(feature_ids):

            seq = self.protein_sequence

            if(num_segments == 1):
                return sequtil.diletter_composition(seq, alph, 1)
            else:
                seqs = sequtil.segment(seq, num_segments)
                return numpy.concatenate(
                    [sequtil.diletter_composition(s, alph, 1) for s in seqs])

        else:

            pairs = sequtil.ordered_alph_pairs(alph)

            feat_ids = []
            feat_names = []

            for si in xrange(1, num_segments + 1):
                for p in pairs:
                    feat_ids.append('%s%i' % (p, si))
                    feat_names.append('%s, segment %i' % (p, si))

            return (feat_ids, feat_names)

    def terminal_end_amino_acid_count(self, terminal_end, length,
                                      feature_ids=False):

        if not(feature_ids):
            alph = sequtil.aa_unambiguous_alph
            seq = self.terminal_end_seq(terminal_end, length)
            return sequtil.letter_count(seq, alph)
        else:
            aa_alph = sequtil.aa_unambiguous_alph
            feat_ids = ['%s%s' % (terminal_end, aa) for aa in aa_alph]
            feat_names = ["%s'-terminal end amino acid count %s" %
                          (terminal_end, aa) for aa in aa_alph]
            return (feat_ids, feat_names)

    def _parse_scales(self, scales):

        if(type(scales) == str):
            try:
                scales = int(scales)
            except ValueError:
                pass

        # retrieve the set of Georgiev aa scales
        if(scales == 'gg'):
            scale_list = sequtil.get_georgiev_scales()
            scale_ids = ['gg%i' % (i) for i in xrange(1, len(scale_list) + 1)]
            scale_names = ['Georgiev scale %i' % (i)
                           for i in xrange(1, len(scale_list) + 1)]

        # retrieve list of pseaac scale
        elif(scales[0] == 'p' and len(scales) > 2):
            scale_indices = [int(i) for i in scales.split('p')[1:]]
            scale_list = scale_indices
            scale_ids = ['pseaac%i' % (i + 1) for i in scale_indices]
            scale_names = ['PseAAC scale %i' % (i + 1) for i in scale_indices]

        # retrieve pseaac scale
        elif(scales[0] == 'p'):
            scale_index = int(scales[1:])
            scale_list = [int(scales[1:])]
            scale_ids = ['pseaac%i' % (scale_index + 1)]
            scale_names = ['PseAAC scale %i' % (scale_index + 1)]

        # retrieve AAIndex scale with index scales
        elif(type(scales) == int):
            scale_list = [sequtil.get_aaindex_scale(scales)]
            scale_ids = ['aai%i' % (scales)]
            scale_names = ['amino acid index %i' % (scales)]

        # retrieve list of AAIndex scales... (still used somewhere?)
        elif(type(scales) == list and all([type[i] == int for i in scales])):
            scale_list = [sequtil.get_aaindex_scale(i) for i in scales]
            scale_ids = ['aai%i' % (i) for i in scales]
            scale_names = ['amino acid index %i' % (i) for i in scales]

        else:
            raise ValueError('Incorrect scale provided: %s\n' % (str(scales)))

        return (scale_list, scale_ids, scale_names)

    def _parse_aa_matrix(self, aa_matrix_id):

        try:
            aa_matrix_id = int(aa_matrix_id)
        except ValueError:
            pass

        if(aa_matrix_id == 'sw'):
            aa_matrix = sequtil.aa_matrix_sw
            aa_matrix_id = 'sw'
            aa_matrix_name = 'schneider-wrede'
        elif(type(sc) == int):
            aa_matrix = sequtil.aa_matrices[aa_matrix_id]
            aa_matrix_id = ['aam%i' % (aa_matrix_id)]
            aa_matrix_name = ['amino acid index matrix %i' % (aa_matrix_id)]
        else:
            raise ValueError('Incorrect matrix provided: %s\n' % (str(scales)))

        return (aa_matrix, aa_matrix_id, aa_matrix_name)

    def average_signal(self, scales, window, edge, feature_ids=False):
        '''
        scales: 'gg',1 ,2 ,3, ..., '1', '2', ...
        window: 5, 7, ...,
        edge: 0...100
        '''

        # fetch scales for provided scales param
        scales, scale_ids, scale_names = self._parse_scales(scales)

        if not(feature_ids):

            result = []

            for scale in scales:
                seq = self.protein_sequence
                result.append(sequtil.avg_seq_signal(seq, scale, window, edge))

            return result

        else:
            return (scale_ids, scale_names)

    def signal_peaks_area(self, scales, window, edge, threshold,
                          feature_ids=False):

        # fetch scales for provided scales param
        scales, scale_ids, scale_names = self._parse_scales(scales)

        if not(feature_ids):

            result = []

            for scale in scales:
                seq = self.protein_sequence
                top, bot = sequtil.auc_seq_signal(seq, scale, window, edge,
                                                  threshold)
                result.append(top)
                result.append(bot)

            return result
        else:

            feat_ids = []
            feat_names = []

            for sid, sname in zip(scale_ids, scale_names):
                feat_ids.append('%stop' % (sid))
                feat_ids.append('%sbot' % (sid))
                feat_names.append('%stop' % (sname))
                feat_names.append('%sbot' % (sname))

            return (feat_ids, feat_names)

    def autocorrelation(self, ac_type, scales, lag, feature_ids=False):

        scale_list, scale_ids, scale_names = self._parse_scales(scales)

        # calculatie features
        if not(feature_ids):

            #num_feat = len(scales) * len(lags)
            result = []

            for scale in scale_list:
                seq = self.protein_sequence
                result.append(sequtil.autocorrelation(ac_type, seq, scale,
                              lag))

            return result
        # or return feature ids and names
        else:
            return (scale_ids, scale_names)

    def property_ctd(self, property, feature_ids=False):

        property_map = {
            'hyd': 'hydrophobicity',
            'vdw': 'normvdw',
            'plr': 'polarity',
            'plz': 'polarizability',
            'chr': 'charge',
            'ss': 'ss',
            'sa': 'sa'
        }

        if not(property in property_map.keys()):
            raise ValueError('Invalid property id provided.')

        # calculatie features
        if not(feature_ids):

            property_name = property_map[property]
            seq = self.protein_sequence

            return sequtil.property_ctd(seq, property_name)

        # or return feature ids and names
        else:
            feat_ids = ['c1', 'c2', 'c3', 't1', 't2', 't3',
                        'd10', 'd11', 'd12', 'd13', 'd14',
                        'd20', 'd21', 'd22', 'd23', 'd24',
                        'd30', 'd31', 'd32', 'd33', 'd34']
            feat_names = ['composition letter %i' % (i) for i in xrange(1, 4)]
            feat_names.extend([
                'transitions 1-2 2-1',
                'transitions 1-3 3-1',
                'transitions 2-3 3-2',
                'distribution letter 1 frac to first letter',
                'distribution letter 1 frac to 25%',
                'distribution letter 1 frac to 50%',
                'distribution letter 1 frac to 75%',
                'distribution letter 1 frac to 100%',
                'distribution letter 2 frac to first letter',
                'distribution letter 2 frac to 25%',
                'distribution letter 2 frac to 50%',
                'distribution letter 2 frac to 75%',
                'distribution letter 2 frac to 100%',
                'distribution letter 3 frac to first letter',
                'distribution letter 3 frac to 25%',
                'distribution letter 3 frac to 50%',
                'distribution letter 3 frac to 75%',
                'distribution letter 3 frac to 100%'])

            return (feat_ids, feat_names)

    def quasi_sequence_order_descriptors(self, aa_matrix, rank, weight=0.1,
                                         feature_ids=False):

        # fetch aa distance matrix
        aam, aam_id, aam_name = self._parse_aa_matrix(aa_matrix)

        alph = sequtil.aa_unambiguous_alph

        if not (feature_ids):

            seq = self.protein_sequence
            return sequtil.quasi_sequence_order_descriptors(seq, aam, rank,
                                                            weight)

        else:
            feat_ids = ['%s' % (l) for l in alph]
            feat_ids.extend(['r%i' % (i) for i in xrange(1, rank + 1)])
            feat_names = ['amino acid %s' % (aa) for aa in alph]
            feat_names.extend(['rank %i' % (i) for i in xrange(1, rank + 1)])
            return (feat_ids, feat_names)

    def pseaac_type1(self, aa_scales, lambda_, weight=0.05, feature_ids=False):
        '''
        aa_scale should be paac0, paac1, ..., paac5
        '''

        alph = sequtil.aa_unambiguous_alph

        pseaac_scale_indices, _, _ = self._parse_scales(aa_scales)

        if not (feature_ids):

            seq = self.protein_sequence
            return sequtil.pseaac_type1(seq, pseaac_scale_indices, lambda_,
                                        weight)

        else:

            feat_ids = ['%s' % (l) for l in alph]
            feat_ids.extend(['l%i' % (i) for i in xrange(lambda_)])
            feat_names = ['amino acid %s' % (aa) for aa in alph]
            feat_names.extend(['lambda %i' % (i)
                               for i in xrange(lambda_)])

            return (feat_ids, feat_names)

    def pseaac_type2(self, aa_scales, lambda_, weight=0.05, feature_ids=False):
        '''
        aa_scale should be paac0, paac1, ..., paac5
        '''

        alph = sequtil.aa_unambiguous_alph

        pseaac_scale_indices, _, _ = self._parse_scales(aa_scales)

        if not (feature_ids):

            seq = self.protein_sequence
            return sequtil.pseaac_type2(seq, pseaac_scale_indices, lambda_,
                                        weight)

        else:

            feat_ids = ['%s' % (l) for l in alph]
            feat_ids.extend(['l%i-s%i' % (i, s)
                             for i in xrange(lambda_)
                             for s in xrange(len(pseaac_scale_indices))])
            feat_names = ['amino acid %s' % (aa) for aa in alph]
            feat_names.extend(['lambda %i scale %i' % (i, s)
                               for i in xrange(lambda_)
                               for s in xrange(len(pseaac_scale_indices))])

            return (feat_ids, feat_names)

    def length(self, feature_ids=False):
        if not(feature_ids):
            return [len(self.protein_sequence)]
        else:
            return (['len'], ['Protein length'])

    def ss_composition(self, num_segments, feature_ids=False):

        alph = sequtil.ss_alph

        if not(feature_ids):

            seq = self.ss_sequence

            if(num_segments == 1):
                return sequtil.letter_composition(seq, alph)
            else:
                seqs = sequtil.segment(seq, num_segments)
                return numpy.concatenate(
                    [sequtil.letter_composition(s, alph) for s in seqs])
        else:

            feat_ids = []
            feat_names = []

            for si in xrange(1, num_segments + 1):
                for ss_id, ss_name in zip(alph, sequtil.ss_name):
                    feat_ids.append('%s%i' % (ss_id, si))
                    feat_names.append('%s, segment %i' % (ss_name, si))

            return (feat_ids, feat_names)

    def sa_composition(self, num_segments, feature_ids=False):

        alph = sequtil.sa_alph

        if not(feature_ids):

            seq = self.sa_sequence

            if(num_segments == 1):
                return sequtil.letter_composition(seq, alph)
            else:
                seqs = sequtil.segment(seq, num_segments)
                return numpy.concatenate(
                    [sequtil.letter_composition(s, alph) for s in seqs])
        else:

            feat_ids = []
            feat_names = []

            for si in xrange(1, num_segments + 1):
                for sa_id, sa_name in zip(alph, sequtil.sa_name):
                    feat_ids.append('%s%i' % (sa_id, si))
                    feat_names.append('%s, segment %i' % (sa_name, si))

            return (feat_ids, feat_names)

    def ss_aa_composition(self, feature_ids=False):

        alph = sequtil.aa_unambiguous_alph
        st_alph = sequtil.ss_alph

        if not(feature_ids):

            seq = self.protein_sequence
            st_seq = self.ss_sequence

            return sequtil.state_subseq_composition(seq, st_seq, alph, st_alph)
        else:
            ids = ['%s%s' % (s, a) for s in st_alph for a in alph]
            names = ['%s-%s' % (s, a)
                     for s in sequtil.ss_name
                     for a in sequtil.aa_unambiguous_name]
            return (ids, names)

    def sa_aa_composition(self, feature_ids=False):

        alph = sequtil.aa_unambiguous_alph
        st_alph = sequtil.sa_alph

        if not(feature_ids):

            seq = self.protein_sequence
            st_seq = self.sa_sequence

            return sequtil.state_subseq_composition(seq, st_seq, alph, st_alph)
        else:
            ids = ['%s%s' % (s, a) for s in st_alph for a in alph]
            names = ['%s-%s' % (s, a)
                     for s in sequtil.sa_name
                     for a in sequtil.aa_unambiguous_name]
            return (ids, names)

    def cluster_composition(self, feature_ids=False):
        if not(feature_ids):
            return sequtil.aa_cluster_composition(self.protein_sequence)
        else:
            return (sequtil.aa_subsets, sequtil.aa_subsets)

    def codon_composition(self, feature_ids=False):
        if not(feature_ids):
            return sequtil.codon_composition(self.orf_sequence)
        else:
            names = ['%s (%s)' % (c, sequtil.translate(c))
                     for c in sequtil.codons_unambiguous]
            return (sequtil.codons_unambiguous, names)

    def codon_usage(self, feature_ids=False):
        if not(feature_ids):
            return sequtil.codon_usage(self.orf_sequence)
        else:
            names = ['%s (%s)' % (c, sequtil.translate(c))
                     for c in sequtil.codons_unambiguous]
            return (sequtil.codons_unambiguous, names)

    # feature calculation help functions

    def terminal_end_seq(self, terminal_end, length):

        '''
        This function returns the C- or N-terminal end side of the protein
        amino acid sequence. The parameter terminal_end must be either 'N' or
        'C' and indicates if the N- or C-terminal end should be returned. The
        length indicates how long the returned sequence will be. If length is
        greater than or equal to the full protein sequence length, than the
        full protein sequence will be returned.

        >>> p = Protein('test')
        >>> p.set_protein_sequence('AAAAACCCCC')
        >>> p.terminal_end_seq('N', 2)
        'AA'
        >>> p.terminal_end_seq('C', 9)
        'AAAACCCCC'
        >>> p.terminal_end_seq('C', 10)
        'AAAAACCCCC'
        >>> p.terminal_end_seq('C', 11)
        'AAAAACCCCC'
        '''
        if not(terminal_end in ['N', 'C']):
            raise ValueError('terminal_end must be either N or C')
        if(length < 1):
            raise ValueError('length must be a positive integer.')

        if(terminal_end == 'N'):
            return self.protein_sequence[:length]
        else:  # terminal_end == 'C'
            return self.protein_sequence[-length:]

    def sequence_signal(self, scale, window, edge):
        return sequtil.seq_signal(self.protein_sequence, scale, window, edge)

    def pfam_family(self, position):
        return self.pfam_hmm_acc(position, 'Family')

    def pfam_domain(self, position):
        return self.pfam_hmm_acc(position, 'Domain')

    def pfam_repeat(self, position):
        return self.pfam_hmm_acc(position, 'Repeat')

    def pfam_hmm_acc(self, position, type):
        if not(self.pfam_annotations is None):
            for annotation in self.pfam_annotations:
                if(annotation.type_ == type):
                    if(position >= annotation.start_pos and
                            position <= annotation.end_pos):
                        # assuming no overlap, return the first one found
                        return annotation.hmm_acc
        return None

    def pfam_clan(self, position):
        if not(self.pfam_annotations is None):
            for annotation in self.pfam_annotations:
                if not(annotation.clan is None):
                    if(position >= annotation.start_pos and
                            position <= annotation.end_pos):
                        # assuming no overlap, return the first one found
                        return annotation.clan
        return None

    def pfam_clan_index(self, position):
        if not(self.pfam_annotations is None):
            for annotation in self.pfam_annotations:
                if not(annotation.clan is None):
                    if(position >= annotation.start_pos and
                            position <= annotation.end_pos):
                        # assuming no overlap, return the first one found
                        return annotation.clan_index
        return None

    def pfam_active_residue(self, position):
        if not(self.pfam_annotations is None):
            for annotation in self.pfam_annotations:
                for r in annotation.active_residues:
                    if(r == position):
                        return True
        return False

    def msa_column(self, position, with_gaps=True):
        '''
        Returns the aligned amino acids of the give positions, without the
        gaps (-).
        '''
        index = position - 1

        if(with_gaps):
            return [s[index] for s in self.msa]
        else:
            return [s[index] for s in self.msa if not s[index] == '-']

    def msa_num_ali_seq(self, position):
        return len(self.msa_column(position, with_gaps=True))

    def msa_num_ali_let(self, position):
        return len(self.msa_column(position, with_gaps=False))

    def msa_variability(self, position, with_gaps=False):
        '''
        Returns the set of (unambiguous) amino acid letters found on the given
        position in the multiple sequence alignment. All other letters (exept
        the gap character '-') are disregarded.

        with_gaps: If set to True, a gap is also part of the column variability
        '''

        # amino acids + gap
        aas = sequtil.aa_unambiguous_alph + '-'
        column_letters_set = set(self.msa_column(position, with_gaps))
        return [l for l in column_letters_set if l in aas]

    def msa_fraction(self, position, letter, with_gaps):
        '''
        TODO: what to do if no aligned seqs, or only few...
        !!! with or without gaps...
        '''
        # obtain all letters on this position in the MSA
        col = self.msa_column(position, with_gaps=with_gaps)

        if(len(col) <= 1):
            assert(len(col) == 1)
            return 0.5
        else:
            # return the fraction of letter
            return float(col.count(letter)) / len(col)

    def msa_conservation_index(self, position):
        '''
        SNPs&GO conservation index, I don't understand the formula...
        '''
        #col = self.msa_column(position, with_gaps=True)

    def msa_entropy21(self, position, with_gaps):

        # amino acids + gap
        aas = sequtil.aa_unambiguous_alph + '-'

        # obtain MSA column letters
        col = self.msa_column(position, with_gaps=with_gaps)

        # is not always the case...
        #assert(all([c in aas for c in col]))

        if(len(col) <= 1):

            assert(len(col) == 1)

            # default entropy in case of no aligned sequences
            entropy = 0.5

            # TODO num seqs < some threshold? Do some other default thing?
        else:

            n = len(col)
            k = len(aas)  # should be 21, 20 amino acids + 1 gap

            # fraction per letter
            na_list = [col.count(l) for l in set(col)]
            pa_list = [float(na) / n for na in na_list]
            na_log_sum = sum([pa * math.log(pa, 2) for pa in pa_list])

            # calculate entropy and return that
            entropy = (-1.0 * na_log_sum) / math.log(min(n, k), 2)

        return entropy

    # check attribute availability functions (simple getters)

    def get_protein_sequence(self):
        return self.protein_sequence

    def get_orf_sequence(self):
        return self.orf_sequence

    def get_secondary_structure_sequence(self):
        return self.ss_sequence

    def get_solvent_accessibility_sequence(self):
        return self.sa_sequence

    def get_msa(self):
        return self.msa

    def get_structure(self):
        return self.protein_structure

    def get_missense_mutations(self):
        return self.missense_mutations

    def get_rasa(self):
        return self.rasa


class Pfam(object):
    '''
    Class that contains Pfam annotation for a protein. Nothing more than a
    data store currently.
    '''

    def __init__(self, start_pos, end_pos, hmm_acc, hmm_name, type_,
                 bit_score, e_value, clan, active_residues):
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.hmm_acc = hmm_acc
        self.hmm_name = hmm_name
        self.type_ = type_
        self.bit_score = bit_score
        self.e_value = e_value
        self.clan = clan
        self.active_residues = active_residues
    '''
    def single_line_str(self):
        return '%i\t%i\t%s\t%s\t%s\t%.1f\t%e\t%s\t%s' % (
            self.start_pos, self.end_pos, self.hmm_acc, self.hmm_name,
            self.type_, self.bit_score, self.e_value, self.clan,
            self.active_residues)

    @classmethod
    def parse(self, s):
        tokens = s.split()
        start_pos = int(tokens[0])
        end_pos = int(tokens[1])
        hmm_acc = tokens[2]
        hmm_name = tokens[3]
        type_ = tokens[4]
        bit_score = float(tokens[5])
        e_value = float(tokens[6])
        clan = None if tokens[7] == 'None' else tokens[7]
        active_residues = eval(' '.join(tokens[8:]))

        return self(start_pos, end_pos, hmm_acc, hmm_name, type_, bit_score,
                    e_value, clan, active_residues)
    '''

if __name__ == "__main__":
    import doctest
    doctest.testmod()
