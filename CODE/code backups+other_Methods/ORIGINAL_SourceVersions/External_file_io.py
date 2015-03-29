import os
from functools import wraps

'''
Created on Sep 10, 2011

@author: Bastiaan van den Berg

- Includes file reading functions, and to read AA matrices, AA scales from file.
'''


def file_reader(func):
    '''
    This function can be used as a decorator around functions read data from a
    file. It assumes that the first argument if a reference to a file, either
    by a path to the or as a file object. If a path is provided, the wrapper
    turns this into a file object.
    '''

    @wraps(func)
    def wrapper(*args, **kwargs):

        # extract file argument (should always be the first argument)
        f, args_rest = args[0], args[1:]

        if(isinstance(f, str) or isinstance(f, unicode)):
            with open(f, 'r') as fin:
                for item in func(fin, *args_rest, **kwargs):
                    yield item
        else:
            for item in func(*args, **kwargs):
                yield item

    return wrapper


@file_reader
def read_ids(f):

    for line in f:
        yield(line.split()[0])


def write_ids(f, ids):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    for uid in ids:
        handle.write('%s\n' % (uid))

    # close file if we opened it
    if not(type(f) == file):
        handle.close()


@file_reader
def read_fasta(f, filter_ids=None):
    '''
    '''

    # initialize sequence id and string to an empty string
    seq_id = ''
    seq_str = ''

    # iterate over each line in the fasta file
    for line in f:

        if(seq_id == '' and seq_str == ''):
            if(line[0] == ">"):
                seq_id = line.split()[0][1:]
                if(seq_id == ''):
                    raise Exception('FASTA file error: Empty id encountered.')
            elif(line[0] == '#'):
                pass
            elif(line.strip()):
                # non-empty line...
                #print(line.strip())
                raise Exception('Error in FASTA file.')
        else:
            if((line.strip() == '' or line[0] == '>') or line[0] == '#'):
                if(filter_ids is None or seq_id in filter_ids):
                    yield (seq_id, seq_str)
                seq_str = ''
                if(line[0] == '>'):
                    seq_id = line.split()[0][1:]
                else:
                    seq_id = ''
            else:
                seq_str += line.strip()

    # return the last sequence (not if the file was empty)
    if not(seq_id == ''):
        if(filter_ids is None or seq_id in filter_ids):
            yield (seq_id, seq_str)


def write_fasta(f, seqs):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    for s in seqs:
        handle.write('>' + s[0] + '\n')
        for i in range(0, len(s[1]), 70):
            handle.write(s[1][i:i + 70] + '\n')
        handle.write('\n')
        handle.flush()

    # close file if we opened it
    if not(type(f) == file):
        handle.close()


@file_reader
def read_ensembl_fasta(f, filter_ids=None):
    '''
    '''

    # initialize sequence id and string to an empty string
    seq_id = ""
    seq_str = ""

    # iterate over each line in the fasta file
    for line in f:

        if(seq_id == "" and seq_str == ""):
            if(line[0] == ">"):
                seq_id = line.split()[0][1:]
                # added for ensembl stuff
                ensembl_tokens = line.split()[1:]
                ensembl_dict = {}
                for token in ensembl_tokens:
                    items = token.split(':')
                    ensembl_dict[items[0]] = ':'.join(items[1:])
            elif(line[0] == '#'):
                pass
            elif(line.strip()):
                # non-empty line...
                print(line.strip())
                raise(Exception, "Error in fasta file")
        else:
            if(line.strip() == "" or line[0] == ">"):
                if(filter_ids is None or seq_id in filter_ids):
                    yield (seq_id, seq_str, ensembl_dict)
                seq_str = ""
                if(line[0] == ">"):
                    seq_id = line.split()[0][1:]
                    # added for ensembl stuff
                    ensembl_tokens = line.split()[1:]
                    ensembl_dict = {}
                    for token in ensembl_tokens:
                        items = token.split(':')
                        ensembl_dict[items[0]] = ':'.join(items[1:])
                else:
                    seq_id = ""
            else:
                seq_str += line.strip()

    # return the last sequence (not if the file was empty)
    if not(seq_id == ""):
        if(filter_ids is None or seq_id in filter_ids):
            yield (seq_id, seq_str, ensembl_dict)


def read_flex(f):
    ids, seqs = zip(*read_fasta(f))
    flex_strings = [s.split(',') for s in seqs]
    flex_floats = []
    for f in flex_strings:
        flex_floats.append([float(i) for i in f if not i.strip() == ''])
    return zip(ids, flex_floats)


def write_flex(f, flex_data):
    ids, flex_floats = zip(*flex_data)
    flex_strings = []
    for flex in flex_floats:
        flex_strings.append(','.join(['%.3f' % (fl) for fl in flex]))
    write_fasta(f, zip(ids, flex_strings))


def read_interaction_counts(f):
    data = read_tuple_list(f, (str, int, int, int, int, int, int))
    return [(item[0], item[1:]) for item in data]


def write_interaction_counts(f, interaction_counts_data):
    tuples = [(i[0], i[1][0], i[1][1], i[1][2], i[1][3], i[1][4], i[1][5])
              for i in interaction_counts_data]
    write_tuple_list(f, tuples)


@file_reader
def read_pfam(f):

    # initialize sequence id and list of annotations
    current_id = ""
    current_annotations = []

    # iterate over lines
    for line in f:

        if(current_id == "" and len(current_annotations) == 0):
            if(line[0] == ">"):
                current_id = line.split()[0][1:]
            elif(line[0] == '#'):
                pass
            elif(line.strip()):
                # non-empty line...
                print(line.strip())
                raise(Exception, "Error in pfam file")
        else:
            if(line.strip() == "" or line[0] == ">"):
                yield (current_id, current_annotations)
                current_annotations = []
                if(line[0] == ">"):
                    current_id = line.split()[0][1:]
                else:
                    current_id = ""
            else:
                #annotation = protein.Pfam.parse(line.strip())
                tokens = line.strip().split()
                start_pos = int(tokens[0])
                end_pos = int(tokens[1])
                hmm_acc = tokens[2]
                hmm_name = tokens[3]
                type_ = tokens[4]
                bit_score = float(tokens[5])
                e_value = float(tokens[6])
                clan = None if tokens[7] == 'None' else tokens[7]
                active_residues = eval(' '.join(tokens[8:]))
                annotation = (start_pos, end_pos, hmm_acc, hmm_name, type_,
                              bit_score, e_value, clan, active_residues)
                current_annotations.append(annotation)

    # return annotations for the last item (not if the file was empty)
    if not(current_id == ""):
        yield (current_id, current_annotations)


def write_pfam(f, pfam_data):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    for uni, pfam in pfam_data:
        handle.write('>%s\n' % (uni))
        for annotation in pfam:
            handle.write('%i\t%i\t%s\t%s\t%s\t%.1f\t%e\t%s\t%s\n' %
                         (annotation))
        handle.write('\n')
        handle.flush()

    # close file if we opened it
    if not(type(f) == file):
        handle.close()


@file_reader
def read_mutation(f):

    for line in f:
        tokens = line.split()

        uni_id = tokens[0]
        pos = int(tokens[1])
        fr = tokens[2]
        to = tokens[3]
        label = int(tokens[4])
        pep = tokens[5]
        pep_i = int(tokens[6])
        codons = tokens[7]
        fr_codon = tokens[8]
        to_codons = tokens[9]
        pdb_id = tokens[10]
        pdb_resnum = int(tokens[11])

        if(pdb_id == 'None'):
            pdb_id = None

        yield((uni_id, pos, fr, to, label, pep, pep_i, codons, fr_codon,
               to_codons, pdb_id, pdb_resnum))


def write_mutation(f, mutations):
    assert(all([len(m) == 12 for m in mutations]))
    write_tuple_list(f, mutations)


@file_reader
def read_mut(f):

    types = (str, int, str, str)

    for line in f:
        tokens = line.split()
        row = []
        for index, t in enumerate(types):
            row.append(t(tokens[index]))
        yield(tuple(row))


def write_mut(f, mutations):
    assert(all([len(m) == 4 for m in mutations]))
    write_tuple_list(f, mutations)


def read_pdb_dir(pdb_fs, pdb_dir):
    '''
    Only .ent files are read...
    Returns AtomGroup object (prody)
    '''

    # import prody for pdb read/write
    import prody
    prody.confProDy(verbosity='none')

    struct_data = []

    for pdb_f in pdb_fs:

        pdb_f = os.path.join(pdb_dir, pdb_f)

        if(os.path.exists(pdb_f)):
            struct = prody.parsePDB(pdb_f)
        else:
            struct = None

        struct_data.append((pdb_f, struct))

    return struct_data


def write_pdb_dir(pdb_dir, struct_data):

    # import prody for pdb read/write
    import prody

    # create output directory, if not yet present
    if not(os.path.exists(pdb_dir)):
        os.makedirs(pdb_dir)

    # write protein chain pdb files to output directory
    for (pdb_f, struct) in struct_data:
        if not(struct is None):
            out_f = os.path.join(pdb_dir, pdb_f)
            prody.writePDB(out_f, struct)


def read_rasa_dir(rasa_fs, rasa_dir):
    rasa_data = []
    for rasa_f in rasa_fs:
        rasa_f = os.path.join(rasa_dir, rasa_f)
        if(os.path.exists(rasa_f)):
            rasa = read_python_list(rasa_f)
        else:
            rasa = None
        rasa_data.append((rasa_f, rasa))
    return rasa_data


def write_rasa_dir(rasa_dir, rasa_data):
    # create output directory, if not yet present
    if not(os.path.exists(rasa_dir)):
        os.makedirs(rasa_dir)

    # write protein chain pdb files to output directory
    for (rasa_f, rasa) in rasa_data:
        if not(rasa is None):
            out_f = os.path.join(rasa_dir, rasa_f)
            write_python_list(out_f, rasa)


def read_python_list(f):

    # read list (not very neat, but whatever :)
    result = eval(f.read())
    assert(type(result) == list)
    assert(all([type(item) == float for item in result]))
    return result


def write_python_list(f, l):

    assert(type(l) == list)
    assert(all([type(item) == float for item in l]))

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')
    handle.write('%s\n' % (str(l)))
    if not(type(f) == file):
        handle.close()


def read_residue_rank_dir(rank_fs, rank_dir):

    rank_data = []

    for rank_f in rank_fs:

        rank_f = os.path.join(rank_dir, rank_f)

        if(os.path.exists(rank_f)):
            rank = read_residue_rank(rank_f)
        else:
            rank = None

        rank_data.append((rank_f, rank))

    return rank_data


def write_residue_rank_dir(rank_dir, rank_data):

    if not(os.path.exists(rank_dir)):
        os.makedirs(rank_dir)

    for (rank_f, rank) in rank_data:
        if not(rank is None):
            out_f = os.path.join(rank_dir, rank_f)
            write_residue_rank(out_f, rank)


@file_reader
def read_residue_rank(f):

    # iterate over lines in file
    for line in f:

        # ignore empty lines and comments
        if(line.strip() and not line[0] == '%'):

            # each other line should be composed of 7 items
            tokens = line.split()

            # read the 7 items
            ali_pos = int(tokens[0])
            seq_pos = int(tokens[1])
            aa = tokens[2]
            coverage = float(tokens[3])
            var_count = int(tokens[4])
            var_letters = tokens[5]
            rvet_score = float(tokens[6])

            # store ass tuple and add to result
            yield((ali_pos, seq_pos, aa, coverage, var_count, var_letters,
                   rvet_score))


def write_residue_rank(f, rank_data):
    assert(all([len(r) == 7 for r in rank_data]))
    write_tuple_list(f, rank_data)


def read_msa_dir(msa_fs, msa_dir):

    msa_data = []

    for msa_f in msa_fs:

        msa_f = os.path.join(msa_dir, msa_f)

        if(os.path.exists(msa_f)):
            msa = read_msa(msa_f)
        else:
            msa = None

        msa_data.append((msa_f, msa))

    return msa_data


def write_msa_dir(msa_dir, msa_data):

    if not(os.path.exists(msa_dir)):
        os.makedirs(msa_dir)

    for (msa_f, msa) in msa_data:
        if not(msa is None):
            out_f = os.path.join(msa_dir, msa_f)
            write_msa(out_f, msa)


def read_msa(f):
    return [seq for seq in read_ids(f)]


def write_msa(f, msa):
    assert(len(msa) > 0)
    assert(all([type(m) == str for m in msa]))
    assert(all([len(msa[0]) == len(m) for m in msa]))
    write_ids(f, msa)


def read_classification_result(f, score=None):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')

    score_names = handle.readline().split(',')

    score_lists = []
    for sn in score_names:
        score_lists.append(eval(handle.readline()))

    if(score):
        if not(score in score_names):
            raise ValueError('Requested score not available: %s' % (score))
        else:
            result = (score, score_lists[score_names.index(score)])
    else:
        result = (score_names, score_lists)

    if not(type(f) == file):
        handle.close()

    return result


def write_classification_result(f, result):
    # TODO
    pass


def read_labeling(f):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')

    class_names = handle.readline().split()

    label_dict = {}
    for line in handle:
        tokens = line.split()
        label = int(tokens[1])
        assert(label < len(class_names))
        label_dict[tokens[0]] = label

    if not(type(f) == file):
        handle.close()

    return (label_dict, class_names)


def write_labeling(f, object_ids, labels, class_names):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    handle.write('%s\n' % ('\t'.join(class_names)))

    for (obj, lab) in zip(object_ids, labels):
        handle.write('%s\t%s\n' % (obj, lab))

    # close file if we opened it
    if not(type(f) == file):
        handle.close()


def read_propka30(filename):

    # values to be read
    feph = []           # 15 times free energy per pH
    chphf = []          # 15 charge per pH folded
    chphu = []          # 15 charge per pH unfolded
    femin = 100000.0
    feminph = -1.0
    pif = -1.0
    piu = -1.0

    # parse status
    first = True
    fe = False
    fe_count = 0
    ch = False
    ch_count = 0
    max_count = 15

    with open(filename, 'r') as fin:

        for line in fin:
            tokens = line.split()
            if(first):
                assert(tokens[0] == 'propka3.0,' and tokens[2] == '182')
                first = False
            if(tokens):
                if(fe and fe_count < max_count):
                    feph.append(float(tokens[1]))
                    fe_count += 1
                if(ch and ch_count < max_count):
                    chphf.append(float(tokens[1]))
                    chphu.append(float(tokens[2]))
                    ch_count += 1
                if(tokens[0] == 'Free'):
                    fe = True
                elif(tokens[0] == 'pH'):
                    ch = True
                if(tokens[0] == 'The' and tokens[1] == 'pI'):
                    pif = float(tokens[3])
                    piu = float(tokens[6])
                if(tokens[0] == 'The' and tokens[1] == 'pH'):
                    feminph = float(tokens[6])
                    femin = float(tokens[13])

    assert(len(feph) == max_count)
    assert(len(chphf) == max_count)
    assert(len(chphu) == max_count)
    assert(not femin == 100000.0)
    assert(not feminph == -1.0)
    assert(not pif == -1.0)
    assert(not piu == -1.0)

    # bit tricky... return as dict???
    return((feph, chphf, chphu, femin, feminph, pif, piu))


@file_reader
def read_names(f):

    for line in f:
        name = line.strip()
        yield(name)


def write_names(f, names):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    for name in names:
        handle.write('%s\n' % (name.strip()))

    # close file if we opened it
    if not(type(f) == file):
        handle.close()


def read_dict(handle, value_type, num_cols=1):
    result = {}
    for line in handle:
        tokens = line.split()
        if(num_cols > 1):
            end = num_cols + 2
            result[tokens[0]] = tuple([value_type(i) for i in tokens[1:end]])
        else:
            result[tokens[0]] = value_type(tokens[1])
    return result


def write_dict(handle, d):
    for key in d.keys():
        handle.write('%s\t%s\n' % (key, str(d[key])))


def read_settings_dict(f):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')

    settings_dict = {}

    keys = [k.strip() for k in handle.readline().split(',')]
    for key in keys:
        line = handle.readline()
        try:
            settings_dict[key] = eval(line)
        except(NameError):
            settings_dict[key] = line.strip()

    # close file if we opened it
    if not(type(f) == file):
        handle.close()

    return settings_dict


@file_reader
def read_tuple_list(f, types):

    for line in f:
        tokens = line.split()
        row = []
        for index in xrange(len(types)):
            row.append(types[index](tokens[index]))
        yield(tuple(row))


def write_tuple_list(f, tuple_list):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'w')

    for tup in tuple_list:
        for item in tup:
            handle.write('%s\t' % (str(item)))
        handle.write('\n')

    # close file if we opened it
    if not(type(f) == file):
        handle.close()


def read_cross_validation(cv_file):
    '''
    Each line contains a list of test set indices, all other indices are
    assumed to be train set.

    returns list of tuples, each tuple containing a list of train indices
    and a list of test indices, these can be used as cv parameter.
    '''

    # store test indices per CV-fold
    tst_is = []

    # featch the cv-fold test indices from file
    with open(cv_file, 'r') as fin:
        for line in fin:
            tst_is.append([int(t) for t in line.split()])

    # get all indices
    all_is = []
    for i_list in tst_is:
        all_is.extend(i_list)
    all_is_set = set(all_is)

    # test if none of the indices is in multiple sets and make a full range
    assert(len(all_is) == len(all_is_set))
    assert(range(len(all_is)) == sorted(all_is))

    # create cv indices
    cv = []
    for tst_i_list in tst_is:
        trn_i_list = all_is_set - set(tst_i_list)
        cv.append((sorted(trn_i_list), sorted(tst_i_list)))

    return cv


def read_scales(f):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')

    scales = []

    first = True
    second = False
    for line in handle:
        if(first):
            letters = line.split(',')
            letters = [l.strip() for l in letters]
            first = False
            second = True
        elif(second):
            scale_ids = line.split(',')
            scale_ids = [s.strip() for s in scale_ids]
            second = False
        else:
            scale = [float(v) for v in line.split(',')]
            assert(len(scale) == len(letters))
            scales.append(dict(zip(letters, scale)))

    # close file if we opened it
    if not(type(f) == file):
        handle.close()

    assert(len(scales) == len(scale_ids))

    return (scales, scale_ids)


@file_reader
def read_scales_db(f):

    scale_id = None
    scale_descr = None
    scale = []
    letters = None

    for line in f:
        tokens = line.split()

        if(len(tokens) > 0):

            if(tokens[0] == '//'):

                # yield scale
                assert(len(letters) == len(scale))
                yield((scale_id, scale_descr, dict(zip(letters, scale))))

                # reset variables for next scale
                scale_id = None
                scale_descr = None
                scale = []
                letters = None
            elif(tokens[0] == 'H'):
                scale_id = tokens[1]
            elif(tokens[0] == 'D'):
                scale_descr = ' '.join(tokens[1:])
            elif(tokens[0] == 'I'):
                firstrow = [t[0] for t in tokens[1:]]
                lastrow = [t[-1] for t in tokens[1:]]
                letters = firstrow + lastrow
            # reading the scale...
            elif not(letters is None):
                for token in tokens:
                    try:
                        scale.append(float(token))
                    except(ValueError):
                        # TODO now append 0.0 if NA value...
                        scale.append(0.0)
            else:
                pass


def read_aa_matrix(f):

    # open file if path is provided instead of file
    if(type(f) == file):
        handle = f
    else:
        handle = open(f, 'r')

    aa_matrix = {}

    row_letters = handle.readline().strip()
    col_letters = handle.readline().strip()

    for row_i, row_letter in enumerate(row_letters):

        values = handle.readline().strip().split(',')

        for col_i, (col_letter, value) in enumerate(zip(col_letters, values)):

            aa_matrix['%s%s' % (row_letter, col_letter)] = float(value)

    # close file if we opened it
    if not(type(f) == file):
        handle.close()

    return aa_matrix


@file_reader
def read_aa_matrix_db(f):

    mat_id = None
    mat_descr = None
    aa_matrix = {}
    row_letters = None
    col_letters = None

    for line in f:

        tokens = line.split()

        if(len(tokens) > 0):

            if(tokens[0] == '//'):

                # yield matrix
                yield((mat_id, mat_descr, aa_matrix))

                # reset variables for next scale
                mat_id = None
                mat_descr = None
                aa_matrix = {}
                row_letters = None
                col_letters = None

            elif(tokens[0] == 'H'):
                mat_id = tokens[1]
            elif(tokens[0] == 'D'):
                mat_descr = ' '.join(tokens[1:])
            elif(tokens[0] == 'M'):
                row_letters = tokens[3][:-1]
                col_letters = tokens[-1]
                row_index = 0

            # reading the scale...
            elif not(row_letters is None):
                row_aa = row_letters[row_index]
                for col_index, token in enumerate(tokens):
                    col_aa = col_letters[col_index]
                    if(token == '-'):  # this is not nice... what to do...?
                        aa_matrix[row_aa + col_aa] = -1
                        aa_matrix[col_aa + row_aa] = -1
                    else:
                        aa_matrix[row_aa + col_aa] = float(token)
                        aa_matrix[col_aa + row_aa] = float(token)
                row_index += 1
            else:
                pass
