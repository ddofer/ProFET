#!/sw/bin/python3.3
#! E:\Python33\python
#Read FASTA files from given directory, generate output csv file with values for features.
'Modified from Michael doron - 7.8.2014'

'4.4.2014 - New version to use in aa molec.info proj. '

'''This script can:
1. a) Parses either all .fasta files in a given directory ,
 OR Parses a single .Fasta file by path location.
 Uses BioPython's SimpleFastaParser, which yields seqid + sequence).
(Each fasta is assumed to hold multiple protein sequences),
1. b) For each protein read, an individual Protein-class object is made, according
to the sequence,  using the helper class "protein_class.py".
2.  a) Features are extracted for each sequence, using the methods in "protein_class.py".
    b) Features are stored as a DICTIONARY (ordered?); with key:value pairs,
    - "key" = Feature NAME (string) , and Value = that feature's value. (For that seq).
    c) Return a multilevel dictionary
3. SeqID:
    c) Output for each protein is written to output file.

TODO
TODO: Remove (extra) index of protein names if unused. (Keeping just file name as "class"/label)
NOTE - LATER methods - will want to add graphing,output ... Output to fasta, prediction of new test samples..
Also - removal of duplicate features. (E.G - double counting of amino acid frequencies in different reduced alphabets.
Also: add multiprocessing.
TODO - multi label capacity.
TODO - Handling of same protein sequence, in different (multiple) fast files. (I.E - multi label case).


NOTE: It MAY be worthwhile to try implementing the tool/package FeatureForge:
http://www.machinalis.com/blog/machine-learning-feature-forge/
http://feature-forge.readthedocs.org/en/latest/


Process repeats for each sequence read.
(Note - Additional filter & transformation steps exist: unwanted aa residue filtering, duplicate proteins, etc').

# Handling multilabel files (from uniprot) should be done differently than handling multiclass .fasta files.

# Handling of group class-labels (0,1,2..) - could be done with string and multiindex,
http://stackoverflow.com/questions/17677154/pandas-normalize-multiindex-dataframe?rq=1
(And later - replace/store. Use split of filename as class/group label. LAter, use labelbinarizer)


http://stackoverflow.com/questions/13226029/benefits-of-pandas-multiindex
"Hierarchical indexing"
http://stackoverflow.com/questions/18262962/convert-dataframe-columns-to-multiindex?rq=1

For simple multiclass (instead of hierarchical labelling/index):
 http://pandas.pydata.org/pandas-docs/stable/basics.html#renaming-mapping-labels

Look at:
http://pandas.pydata.org/pandas-docs/dev/groupby.html#transformation
(Groupby, key.. )

# Handling of sklearn-pandas - options: df.values, or an existing (confusing) external module,
http://www.markhneedham.com/blog/2013/11/09/python-making-scikit-learn-and-pandas-play-nice/

'''

'AAlphabets Holds alt.alphabet, trans.dicts+methods and alph letters.'
import os
from multiprocessing import Pool #Currently unused in multiclass
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import numpy as np
# from Feature_Extract.AAlphabets import *
# import Feature_Extract.ProtFeat
from AAlphabets import *
import ProtFeat


FILE_LOC = "./test_seq/Organellas"
FILE_LOC_TRY ="./test_seq"

'ILLEGALS  - aa stored in Aalphabet'
N_TAIL_REGION = 26
C_TAIL_REGION = 23
MIN_PROTEIN_LENGTH = 20 #Minimum length for protein sequences to be processed by FeatureGen().
MIN_PROTEIN_LENGTH_TAILED = 47 #If Terminal Tail split (27) Minimum length for protein sequences to be processed by FeatureGen().
'Note: Re min length - this MUST take into account the length of the sequence '
' post "cleaving" for tail(s) AND the protparam windows! (minimum length ~17+)'



def remove_unknown_AA(s):
  '''
  Checks for and Removes unknown-Illegal AA in string/protein sequence.
  Does NOT deal with any ambigous or nonstandard AA in the sequence.
  Code for 'unknown amino acid' = Z
  UNKNOWN_AA = "Z"
  What best to do? (Currently whole protein ignored).
  Possible: replace modified with standard, and count presence of
  nonstandard AA as a feature (add to dict). Remove or replace others?
  '''
  if UNKNOWN_AA in s.upper():
    s.replace(UNKNOWN_AA,"")
  return s

def contain_illegals(seq, illegals=ILLEGALS):
  '''
  Used tofilter  U,B,Z,X.. (non standard AA) from a given sequence.
   Returns True if illegals in sequence.
  >>>print(contain_illegals('XSC'))
  True
  '''
  for c in illegals:
    if c in seq:
      return True
  return False

def dict2df (seq_feature_dict, orientation='index'):
  '''
  Takes a dict of dicts as input, returns a pandas DataFrame.
  NOTE: We also replace NaN with 0 ..
  '''
  Dframe = pd.DataFrame.from_dict(seq_feature_dict, orient=orientation)
  Dframe.fillna(0, inplace=True) #Fill missing data
  return(Dframe)

'Remove scale factor maybe... '
def col_MAD (col):
  '''
  Could be replaced by :
  http://statsmodels.sourceforge.net/stable/_modules/statsmodels/robust/scale.html#stand_mad
  http://statsmodels.sourceforge.net/devel/generated/statsmodels.robust.scale.mad.html
  Calc Median Absolute Deviation (MAD) of a given list/series. (col)
  returns scaled (by 0.75th percentile) MAD.
  '''

  '''
  # This if we want to use a custom, consistent, non gaussian scaling factor:
  MAD_unscaled = calc_MAD(col) #,listMedian=col.median)
  "scale factor, If we don't assume an underlying normal distribution:"
  q_75 = col.quantile(0.75) #75th quantile
  MAD_scaled = MAD_unscaled/q_75

  We can also use statsmodels directly; statsMAD
  '''
  return (calc_MAD(col))
#==============================================================================
#
# #main CODE:
#
#==============================================================================

"http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/11%20-%20Lesson.ipynb"
def Get_Dirr_All_Fasta (classType, Dirr = '.'):
    '''
    Get all FASTA (*.fasta) files from current working directory,
    returns a list of files.
    If not additional param given, default is in current dir
    CURRENTLY - Does not get files from subdirectories.
    '''
    files_dict = {}
    '''old -   for f in os.listdir(sys.argv[1]) :
    We could also do:
    for file in glob("*.fasta"):
    '''
    if Dirr != '.':
        os.chdir(str(Dirr))
        print ("dirr change to: %s" %(Dirr))

    for root, dirs, files in os.walk(Dirr) :
        for name in files:
            if (name.endswith('.fasta')):
                className = ''
                if classType == 'dir':
                    className = os.path.basename(root)
                elif classType == 'file':
                    filename, ext = os.path.splitext(name)
                    className = filename
                if name.endswith('.fasta'):
                    files_dict[os.path.join(root, name)] = className
    return files_dict


# print(Get_Dirr_All_Fasta())

def get_reduced_prot (seq,alph_name="ofer14"):
  '''
  Given a sequence and name of a reduced a.a set,
  gets the translated, reduced AA version of the sequence, and
  returns a protFeat object, constructed from the "reduced" seq.
  Uses the alphabets in AAAlph. (Custom alphabets should be implemented there).

  '''
  red_seq = translate_sequence(seq, REDUCED_ALPHABETS_TRANSDICTS[alph_name])
  # reduced_protein = Feature_Extract.ProtFeat.ProtFeat(ProteinSequence=red_seq,alph_letters=alph_name)
  reduced_protein = ProtFeat.ProtFeat(ProteinSequence=red_seq,alph_letters=alph_name)
  return(reduced_protein)

'TODO: Make it capable of accepting a list/dict of the various functions to call per prot or subseq..'
def Get_Protein_Feat(seq,body_reduced_alph='ofer14',split_N=False,split_C=False,N_TAIL_REGION=27,C_TAIL_REGION=27):
    '''
    Get most/default protein features.
    This includes Features for the N-tail end (subseq), seperately from the
    remaining sequence, and 14_letter (default) 2-mer composition.

    Can be expanded, to pass on alt. arguments (= which ProtFeat functions called).
    Reduced AA representation can be done here or seperately.

    N_TAIL_REGION, C_TAIL_REGION - Rough guess.
     Av. Metazoan tail length for sub 300 length proteins is 67 for both tails together. (?)
     From: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002364#s2

     N tail is manually based on SP av length..
    '''
    features_dict = {}

    main_seq = seq
    if split_N is True:
        N_seq = seq[0:N_TAIL_REGION]
        main_seq = main_seq[N_TAIL_REGION:]
        # N_protein = Feature_Extract.ProtFeat.ProtFeat(N_seq,HAS_N=True,HAS_C=False)
        N_protein = ProtFeat(N_seq, HAS_N=True, HAS_C=False)
    if split_C is True: #Slice out end of seq
        C_seq = seq[-C_TAIL_REGION:]
        main_seq = main_seq[:C_TAIL_REGION]
        # C_protein = Feature_Extract.ProtFeat.ProtFeat(C_seq,HAS_N=False,HAS_C=True)
        C_protein = ProtFeat.ProtFeat(C_seq, HAS_N=False, HAS_C=True)
    # main_protein = Feature_Extract.ProtFeat.ProtFeat(main_seq,HAS_N=(not(split_N)),HAS_C=(not(split_C))) #new object of the proteinFeat class
    main_protein = ProtFeat.ProtFeat(main_seq, HAS_N=(not (split_N)),HAS_C=(not (split_C)))

    'Get simple protein features'
    features_dict.update (main_protein.GetSimpleFeatures(ParamScaleWindow=9)) #call on class method
    'N-Tail properties:'
    if split_N is True:
        features_dict.update (N_protein.tail_properties(tail_end='N',reduced_alph = 'Ofer_N_Tail'))
    'N-Tail properties:'
    if split_C is True:
        features_dict.update (N_protein.tail_properties(tail_end='C',reduced_alph = 'Ofer_N_Tail'))

    'Get features after any potential splitting of tails'
    features_dict.update (main_protein.GetCTD('CTD')) # CTD = 'CTD'...
    features_dict.update (main_protein.GetPTMMotifs())
    features_dict.update (main_protein.GetCleavageCounts())

    'Reduced alphabet representation(s) of main protein '

    ofer13_seq = get_reduced_prot (seq=main_seq,alph_name='ofer13')
    features_dict.update(ofer13_seq.GetkgramFreq(k=2))
    # features_dict.update(ofer13_seq.GetKMirrorsFreq(k=2,getFreq=False))

    # ofer_Tail8_seq = get_reduced_prot (seq=main_seq,alph_name='Ofer_N_Tail')
    # features_dict.update(ofer_Tail8_seq.GetKMirrorsFreq(k=3,getFreq=False))
    # features_dict.update (ofer_Tail8_seq.GetEntropy(normalizeTotEntropy=True))

    # alex6_seq = get_reduced_prot (seq=main_seq,alph_name="alex6")
    # features_dict.update (alex6_seq.GetAA_Freq())
    # features_dict.update (alex6_seq.GetEntropy())
    # features_dict.update(alex6_seq.GetkgramFreq(k=2))

# UseVariableCombinations not Tested if working yet!!
    # hp3_seq = get_reduced_prot (seq=main_seq,alph_name="hp3_Plus")
    # features_dict.update (hp3_seq.GetEntropy())
    # features_dict.update (hp3_seq.GetkgramCounts(k=4))

    hp2_seq = get_reduced_prot(seq=main_seq,alph_name="hp2")
#     features_dict.update (hp2_seq.GetEntropy())
#    features_dict.update (hp2_seq.GetkgramCounts(k=5,UseVariableCombinations=True)) # name buggy?
    features_dict.update (hp2_seq.GetkgramFreq(k=5))

    # sdm12_seq = get_reduced_prot (seq=main_seq,alph_name="sdm12")
    # features_dict.update (sdm12_seq.GetEntropy())
    # features_dict.update(sdm12_seq.GetkgramFreq(k=2))

    selected3alph_3 = ['Charge_3','SolventA_3','NormVDWV_3'] #],'Hydrophobicity_3' ]
    selected3alph_5 = ['Disorder_3','SecondaryStr_3']
    ' Get Profeat three letter group properties on whole sequence, (can be just main-seq)'
###    for alph_3 in THREE_LETTER_ALPH_NAMES: #THREE_LETTER_ALPH_NAMES - from Aalphabets, ProFEAT groups
    for alph_3 in selected3alph_3:
        seq_3 = get_reduced_prot (seq=seq,alph_name=alph_3)
        # features_dict.update (seq_3.GetKMirrorsFreq(k=3))
        features_dict.update (seq_3.GetkgramFreq(k=3))
    for alph_3 in selected3alph_5:
        seq_3 = get_reduced_prot (seq=seq,alph_name=alph_3)
        features_dict.update (seq_3.GetKMirrorsFreq(k=5))
        # features_dict.update (seq_3.GetkgramFreq(k=5))


    return features_dict


'''
Go Over each fasta file's sequences individually.
# http://biopython.org/DIST/docs/api/Bio.SeqIO.FastaIO-module.html

We output it as pandas dataframe -> csv , from nested dict of dicts
http://pandas.pydata.org/pandas-docs/stable/dsintro.html#from-a-list-of-dicts

http://stackoverflow.com/questions/19436789/biopython-seqio-to-pandas-dataframe/19452991#19452991
Pandas dataframe from dict of dicts - pd.DataFrame.from_dict(prot, orient='index')
http://stackoverflow.com/questions/13575090/construct-pandas-dataframe-from-items-in-nested-dictionary
'''

'TODO: Make multicore-map work. (Currently it never gets called for multifile multiclass fasta..)'
'TODO: Make Get_Protein_Feat flexible in terms of called features, add sending of params.. (EG - getNTail)'
def GenFeaturesFromFasta(files=None, Dirr = '.'): #, Multiclass=True):
    '''
    Gets protein features for each protein in the given list of file(s),
    or gets list of fasta files from provided Dirr.

    Returns a dict of dicts.,
    {protid: {'feature':feature value}}.

    TODO: Make it return by default the dict2df() -> dataframe.
    (and return the dict only if requested).
    This is more conveinet and avoids forgetting to convert.
    '''
    all_feature_dict = {}
    def GetFeatFasta(filename, classname):
        with open(filename) as fasta:
            #for seq_record in SeqIO.parse(fasta_file, "fasta"):
            # id_list =[] #Will hold the protein/sequence fastaheader/name, in order. (For use w' pandas)
            # Feature_List = []
            for record in SimpleFastaParser(fasta): #SeqIO.SimpleFastaParser(fasta):
                sequence = record[1]
                title = record[0]
                seq_id = title.split(None, 1)[0]
              # for record in SeqIO.FastaIO.FastaIterator(fasta):
              #   title = record.id
              #   sequence = str(record.seq)
                'Filter out nonstandard sequences: Illegal AA, short, or already processed: '
                # if ((contain_illegals(sequence,ILLEGALS))==False) & (len(sequence)>MIN_PROTEIN_LENGTH) & (all_feature_dict.get(seq_id) is None):
                if (((contain_illegals(sequence,ILLEGALS))==False) & (len(sequence)>45) & (all_feature_dict.get(seq_id) is None)):
                    all_feature_dict[seq_id]=Get_Protein_Feat(sequence)
    # print ('type(files) %s' %(type(files)))

    for filename, classname in files.items():
        print('Getting features from a single fasta file- %s' %(files.keys()))
        GetFeatFasta(filename,classname)

    """
    'TODO: Multiprocessing here'
    if isinstance(files,str): #Single file to open
        print('Getting features from a single fasta file- %s' %(files.keys()))
        GetFeatFasta(files)

    else: #List of files
        print ('getting multiple fasta files')
        # for fasta_file in files:
        #     print (fasta_file)
        #     GetFeatFasta(fasta_file)
        #Multiprocess map
        multiproc = Pool()
        multiproc.map(GetFeatFasta,files)
    """
    # print ('returning all_feature_dict')
    return all_feature_dict


'We will later want to store the MAD, median for each feature externally.'
'''
Alt calcMAD ?:
Should use http://statsmodels.sourceforge.net/stable/_modules/statsmodels/robust/scale.html#stand_mad
CHECK!
'''
def calc_MAD(a, axis=0): #, c=Gaussian.ppf(3/4.)):  # c \approx .6745
    """
    The Median Absolute Deviation along given axis of an array. From statsmodels:
    http://statsmodels.sourceforge.net/stable/generated/statsmodels.robust.scale.mad.html

    REPLACE!  with -
     http://statsmodels.sourceforge.net/stable/_modules/statsmodels/robust/scale.html#stand_mad
    Parameters:
    ----------
    a : array-like
        Input array.
    c : float, optional
        The normalization, scaling constant.
        For a normal distribution, defined as scipy.stats.norm.ppf(3/4.), ~ .6745.
    axis : int, optional.         The default is 0.
    Returns
    -------
    mad : float
        `mad` = median(abs(`a`))/`c`
    """
    c = 0.6745
    'MAD, not standardized.. ver: '
    # a = np.asarray(a)
    # c = np.percentile(a,0.75) #75th quantile
    # return np.median((np.fabs(a))/c, axis=axis)
    'MAD_Stand ver:'
    a = np.asarray(a)
    d = np.median(a, axis = axis)
    ### d = tools.unsqueeze(d, axis, a.shape)
    return np.median(np.fabs(a - d)/c, axis = axis)

'Check file vs files... ? '
'TODO: Memory bloat likely here!'
def get_features(trainingSetFlag, classType, files=None,returndf = True,saveCSV = True,Dirr = '.',\
  saveName = "Features_seq", multiClass=True,  multiClassLabels = ["one","two","three","four"]):
    '''
    Get file or dir location,
    call other def to extract data from fasta files/seq there,
    convert to data-frame, then save frame to .csv.

    Multiclass - is Similar to get_features, but adds a label to each dataframe seperately.
    param:   Multiclass: assumes each files represents  a different class.
    For each file/class, a 'labels' column will be added to the returned dataframe.
    'labels' can be numeric or string , if label_binarizer from sklearn is implemented.

    TODO: MAke multiClassLabels use filename split.

    Also possible:
     http://pandas.pydata.org/pandas-docs/version/0.13.1/basics.html#renaming-mapping-labels
        http://stackoverflow.com/questions/19851005/rename-pandas-dataframe-index?rq=1
        dataframe.index.names = ['samples']

     Current, best ?:
     http://www.bearrelroll.com/2013/05/python-pandas-tutorial/
     OR
     http://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex

    '''
    ## FIX ?
    # print(multiClassLabels)
    def get_MultiClass_features(trainingSetFlag, classType):
        fasta_files_dict = Get_Dirr_All_Fasta (classType,Dirr)

        print ('Multiclass fasta_files list found: %s' %(fasta_files_dict) )
        # i=0
        df = pd.DataFrame()
        'TODO: Modify this bit, so Mult.files can be sent for processing at once/multiprocessing'
        for key, value in fasta_files_dict.items():
            'Get seperate DataFrame of features for all sequences in each seperate file:'
            # dataframe = dict2df(GenFeaturesFromFasta(files)) # ORIG?
            features = GenFeaturesFromFasta(files={key:value}) #,Dirr=Dirr)
            dataframe = dict2df(features)
            # dataframe['labels']=pd.Series(multiClassLabels[i],index=dataframe.index)
            """
            Set new index by "labels"
            http://www.bearrelroll.com/2013/05/python-pandas-tutorial/
            (We can also column-wise (axis=1) groupby - ['labels'] ) ;
            or:         df.index = ['two' for x in df.index]
            multilevel index:
            http://stackoverflow.com/questions/20085308/add-multi-index-to-pandas-dataframe-and-keep-current-index?rq=1
            http://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex
            """
            # dataframe['labels'] = str(multiClassLabels[i])
            # print(file.split())
            filename = os.path.basename(key)
            classname = value
            dataframe['proteinname'] = [name for name in dataframe.index]
            if trainingSetFlag == True:
                dataframe['classname'] = classname


            #dataframe.reindex(index, level=0)

            'Change this bit for multilevel / hierarchical index..'
            # dataframe.set_index('labels', append=False, inplace=True,drop=True) #ORIGINAL #WORKS
            # dataframe.set_index('labels', append=True, inplace=True,drop=True) #Try?

            ## dataframe.reindex('labels', append=True, inplace=True,drop=True)
            # dataframe.reorder_levels(['labels',dataframe.index.name])
            ### i += 1
            df = df.append(dataframe) #,inplace=True) #,ignore_index = True)
        # print(df)
        if trainingSetFlag == True:
            df.set_index(['proteinname','classname'], append=True, inplace=True,drop=True)
        else:
            df.set_index(['proteinname'], append=True, inplace=True,drop=True)
        return df

    if multiClass == True:
      print ("Getting multiclass dataframes")
      features = get_MultiClass_features(trainingSetFlag, classType)
      dataframe = features
    else:
      print ('Getting & merging features from all subfiles')
      features = GenFeaturesFromFasta(files,Dirr)
      dataframe = dict2df(features)

    print('Features generated')
    # print('Data converted to dataframe, NaN cleaned')
    if saveCSV == True:
        print('Now Saving')
        dataframe.to_csv(saveName+'.csv')  #TODO: 'Save to output folder?'
    if returndf== True:
        return dataframe


def UFF_features(files=None):
  '''
  Like save_features, but pivoted and gets just feature names - for
  unsupervised feature selection (UFFizi).
  Max 500 samples.

  >>> uf2 = features.iloc[:500]
  >>> uf2 = uf2.T
  >>> uf2 = uf2.iloc[1:]
  >>> uf3 = uf2.drop_duplicates()
  >>> uf3.to_csv('transposed_uff.csv',header=False,index_label=None)

  http://adios.tau.ac.il/UFFizi/
  '''
  features = GenFeaturesFromFasta(files)
  print('Features generated')
  "Could be replaced with call to (get_features(returndf = True ,saveCSV = False')).T "
  Dframe = pd.DataFrame.from_dict(features, orient='columns') #We could also just use .T Transpose
  Dframe.fillna(0, inplace=True) #Fill missing data
  print('UFF Data converted to dataframe, NaN cleaned. Now saving')
  Dframe.to_csv('UFF_Feat.csv')

'TODO: Pickle Medians & MAD from population, or from training data'
def normalize_median(df,saveCSV = True,filename='Feat_MAD.csv'):
    '''
    This functions normalizes a Dataframe's values, Col by Column (default),
    using the median and MAD (Median Absolute Deviation)!
    This is a robust statistics, as opposed to Z-score normalizing.

    IF We do NOT assume a Normal distribution, MAD Scale factor
    is set to  1/(0.75th quantile) of the sampled distribution.
    (In a normal distribution, it would be ~ 1/0.675 = ~1.45).

    Note - Scaling [0,1] is not included currently.

    IMPORTANT: We will want to store the median and MAD for each feature,
    for use with test sets. (Via pickling).

    TODO: We calc the MAD using an external function -  CHECK if CORRECT!!??

    "WARNING - assumes no name on col holding names. might change!!"
    '''
    "'Applying a applying a function on 1D arrays to each column or row."
    'or groupBy  - p 252'

    #Put here a default dict for median, MAD per col.  For storage/pickle
    ##
    'groupby to avoid error by dividing strings. ? '

    'http://stackoverflow.com/questions/12725417/drop-non-numeric-columns-from-a-pandas-dataframe'
    'http://stackoverflow.com/questions/18137341/applying-functions-to-groups-in-pandas-dataframe'

    'x is a column (series) in the dataframe. p 133 in pd  book'
    f = lambda x: ((x-x.median( numeric_only=True))/col_MAD(x))
    # number_df = df._get_numeric_data()
    MAD_df =df.apply(f)

    MAD_df.replace([np.inf, -np.inf], 0)
    MAD_df.fillna(0, inplace=True)

    'NOTE - we could also do something like: math.fabs((df-df.median)/df.MAD_scaled)) ? '
    'save results'

    if saveCSV:
        MAD_df.to_csv(filename)
    return MAD_df

def normalize_zscore(df,saveCSV = True,filename='Feat_normalized.csv', normParams = pd.DataFrame.from_records([{'a':1}])):
  '''
  Normalizes (column wise) the dataframe .
  We will now have  now have mean 0 and standard deviation 1.
  http://pandas.pydata.org/pandas-docs/stable/groupby.html#transformation

  Note that this operation is not robust to outliers, (unlike MAD),
   and that future data MUST be normalized with the same params used here.
  Alternative: Sci-kit learn Normalizer.
  '''
  if (len(normParams.index) == 1):
      f = lambda x: (x - x.mean()) / x.std()
      meanNormParams = df.copy().apply(np.mean)
      stdNormParams = df.copy().apply(np.std)
      normParams = pd.concat([meanNormParams, stdNormParams], axis=1)
      normParams.columns = ['mean','std']
      z_df = df.apply(f)
      z_df.replace([np.inf, -np.inf], 0)
      z_df.fillna(0, inplace=True)
      if saveCSV:
          z_df.to_csv(filename) #,index=False)
  else:
      z_df = df.copy()
      for index, row in normParams.iterrows():
          z_df[index] = (z_df[index] - row['mean']) / row['std']
      z_df.replace([np.inf, -np.inf], 0)
      z_df.fillna(0, inplace=True)
  return z_df, normParams

'TODO: NOT Working'
def filterDF (df,RemoveZeroes=True,RemoveDuplicateCol=False,RemoveNoVarCol=True):
    '''
    Filter "noisy" data columns, of all zeroes, or duplicate values. ()
    Creates modified, "cleaned" version of the input dataframe.
    (Modify inplace? Not sure...)

        http://stackoverflow.com/questions/12411649/filter-columns-of-only-zeros-from-a-pandas-data-frame)
    http://pandas.pydata.org/pandas-docs/dev/groupby.html#filtration
    http://wesmckinney.com/blog/?p=340
    '''

    'Drop duplicate columns/features'
    'Currently not working - due to duplicates in labels column? '
    # Training_Df.T.drop_duplicates(inplace=True).T

    if removeZeroes==True:
        # mod_DF = df.filter(lambda x: x != 0)
##        mod_df = df.groupby(df).filter(lambda x: x != 0)
        df_cleaned = df[f.columns[(df != 0).any()]]
        return df_cleaned


# import click
#
# @click.command()
# @click.option('--dir','-d','directory',default='.',help='The path to the fasta files', type = str)
# @click.option('--training_set','-r','trainingSetFlag',default=True,help='Whether this is a training set', type = bool)
# @click.option('--type_of_class','-c','classType',default='dir',help='Defines the class to each protein by its \
#                                                                         directory "dir" or by its file name "file"', type = str)
# @click.option('--normalize','-n','normFlag',default=True,help='Whether to normalize the data', type = bool)
# @click.option('--normParams','-np','normParams',default='.',help='The path to the csv containing the normalization \
#                                                                         parameters for the testing data (irrelevant \
#                                                                         in training data)', type = str)
def featExt(directory, trainingSetFlag, classType, normFlag, normParams):
    #file_location = FILE_LOC_TRY
    # file_location = FILE_LOC
    # file_location ="./test_seq/Extracellular"
    # file_location ="./test_seq/Viri"
    # file_location ="./test_seq/Organellas"
    # file_location ="./test_seq/Viri/Viri-Capsids"
    # file_location ="./test_seq/NP"
    file_location = "./test_seq/NP/NP+SP+Negs_SmallBalanced"

    df = get_features(trainingSetFlag, classType, returndf = True,
      saveCSV = False, saveName = "raw_Train_Feat",
     multiClass=True, Dirr = directory)

    'Removal of all zero columns - should not be done this way, disabled. Dan'

    # print("Removing all zero features")

    # Training_Df =filterDF (Training_Df,removeZeroes=True,RemoveDuplicateCol=False) #OLD
    # Training_Df = Training_Df.groupby(Training_Df).filter(lambda x: x != 0) #OLD

    # df_cleaned = df[df.columns[(df != 0).any()]] #All zero columns removed #ORIG
    df_cleaned = df

    # z_df = normalize_zscore(Training_Df,saveCSV = True,filename='Feat.csv')
    # MAD_df = normalize_median(Training_Df,saveCSV = True,filename='MAD_Feat.csv')

    'MAD then Z-score'
    ##    MAD_df = normalize_median(Training_Df,saveCSV = False,filename='MAD_Feat.csv')
    ##    df_cleaned = normalize_median(df_cleaned,saveCSV = False)

    if normFlag == True:
        if (trainingSetFlag == True):
            z_df, normParams = normalize_zscore(df_cleaned,saveCSV = False,filename='Z+MAD_Feat.csv')
            normParams.to_csv('trainingSetNormParams.csv')
        else:
            normParamDf = pd.read_csv(normParams, index_col=[0])
            z_df, normParams = normalize_zscore(df_cleaned,saveCSV = False,filename='Z+MAD_Feat.csv',normParams = normParamDf)
    else:
        z_df = df_cleaned

    'Removal of all zero columns - should not be done this way, disabled. Dan'
    # z_df_2 = z_df[z_df.columns[((abs(z_df) > 0.00001)).any()]] #All zero columns removed #ORIG
    z_df_2 = z_df #NEW

    'Remove sample/sequence names, keep only file names as index'
    'http://stackoverflow.com/questions/18624039/pandas-reset-index-on-series-to-remove-multiindex'
    z_df_2.reset_index(0,drop=True,inplace=True) #Keep only class labels
    if (trainingSetFlag == True):
        z_df_2.to_csv('trainingSetFeatures.csv')
    else:
        z_df_2.to_csv('testingSetFeatures.csv')

    # print (('ORIGINAL Df: %s') %(Df))
    # print (('MAD normalized Df: %s') %(MAD_df))
    print("Done")

    "http://stackoverflow.com/questions/22485375/efficiently-select-rows-that-match-one-of-several-values-in-pandas-dataframe?rq=1"
    "df[df.Name.isin(['cytosol', 'golgi'])]"
    return z_df, normParams

if __name__=="__main__":
    featExt()

