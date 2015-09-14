#!/sw/bin/python3.3
#! E:\Python33\python
#Read FASTA files from given directory, generate output csv file with values for features.

"BUG: Handling of features with STD, mean =0. (df_cleaned... and normalize_zscore ) - Dan. 25.1"

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
NOTE - LATER methods - will want to add graphing,output ... Output to fasta, prediction of new test samples..
Also - removal of duplicate features. (E.G - double counting of amino acid frequencies in different reduced alphabets.
Also: add multiprocessing.
TODO - multi label capacity.
TODO - Handling of same protein sequence, in different (multiple) fast files. (I.E - multi label case).

NOTE: It MAY be worthwhile to try implementing the tool/package FeatureForge:
http://www.machinalis.com/blog/machine-learning-feature-forge/

'''
'AAlphabets Holds alt.alphabet, trans.dicts+methods and alph letters.'
import os
from multiprocessing import Pool #Currently unused in multiclass
from Bio.SeqIO.FastaIO import SimpleFastaParser, FastaIterator, FastaWriter
import pandas as pd
import numpy as np
from collections import OrderedDict
# from Feature_Extract.AAlphabets import *
# import Feature_Extract.ProtFeat
from AAlphabets import *
from AAScales import MinScales_Dict
import ProtFeat
#import time

'ILLEGALS  - aa stored in Aalphabet'
N_TAIL_REGION = 26
C_TAIL_REGION = 24
MIN_PROTEIN_LENGTH = 27 #Minimum length for protein sequences to be processed by FeatureGen().
MIN_PROTEIN_LENGTH_TAILED = 45 #If Terminal Tail split (27) Minimum length for protein sequences to be processed by FeatureGen().
'Note: Re min length - this MUST take into account the length of the sequence '
' post "cleaving" for tail(s) AND the protparam windows! (minimum length ~17+)'


def Get_Protein_Feat(seq,
  SeqReducedAlph='ofer14',ReducedK=2,
  GetSimpleFeatSet=True,
  GetExtraScaleFeatSet=True,
  aaParamScaleWindow=7,
  ExtraScaleWindow=17,
   GetSubSeqSegs=True,SubSeqSegs=3,
   GetTriLetterGroupAlphKFreq=True,TriLetterGroupAlphK=5,
   GetSeqReducedGroups=True,SeqReducedGroups='ofer_w8',
   GetSeqReducedAlph=True,GetCTDFeatSet = True,
   GetPTMFeatSet= True, GetDBCleavageFeatSet=True,
   split_N=False,split_C=False,N_TAIL_REGION=30,C_TAIL_REGION=30
   ):
    '''
    Get most/default protein features.
    This includes Features for the N-tail end (subseq), seperately from the
    remaining sequence, and ofer14 (default) k-mer composition.

    Can be expanded, to pass on alt. arguments (= which ProtFeat functions called).
    Reduced AA representation can be done here or seperately.

    N_TAIL_REGION, C_TAIL_REGION : Av. Metazoan tail length for sub 300 length proteins is 67 for both tails together.
     From: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002364#s2

     N tail is based on SP average length..

     SeqReducedAlph : Reduced alphabet representation to use with extracting K-Mers.
     Can be any saved Reduced alphabet - e.g. shen7, ofer_w8.
     '''
    features_dict = {}
    main_seq = seq
    if split_N is True:
        N_seq = seq[0:N_TAIL_REGION]
        main_seq = main_seq[N_TAIL_REGION:]
        N_protein = ProtFeat(N_seq, HAS_N=True, HAS_C=False)
    if split_C is True: #Slice out end of seq
        C_seq = seq[-C_TAIL_REGION:]
        main_seq = main_seq[:C_TAIL_REGION]
        # C_protein = Feature_Extract.ProtFeat.ProtFeat(C_seq,HAS_N=False,HAS_C=True)
        C_protein = ProtFeat.ProtFeat(C_seq, HAS_N=False, HAS_C=True)
    # main_protein = Feature_Extract.ProtFeat.ProtFeat(main_seq,HAS_N=(not(split_N)),HAS_C=(not(split_C))) #new object of the proteinFeat class
    main_protein = ProtFeat.ProtFeat(main_seq, HAS_N=(not (split_N)),HAS_C=(not (split_C)))

    if GetSimpleFeatSet==True:
        'Get simple protein features:'
        features_dict.update (main_protein.GetSimpleFeatures(ParamScaleWindow=aaParamScaleWindow)) #call on class method

    if GetExtraScaleFeatSet==True:
      features_dict.update (main_protein.Get_ParamScales_Features(window=ExtraScaleWindow,edge=0.8))

    'N-Tail properties:'
    if split_N is True:
        features_dict.update (N_protein.tail_properties(tail_end='N',reduced_alph = 'Ofer_N_Tail'))
    'N-Tail properties:'
    if split_C is True:
        features_dict.update (N_protein.tail_properties(tail_end='C',reduced_alph = 'Ofer_N_Tail'))

    'Get features (after any potential tail splitting)'
    if GetCTDFeatSet == True:
        features_dict.update (main_protein.GetCTD('CTD')) # CTD = 'CTD'...
    if GetPTMFeatSet == True:
        features_dict.update (main_protein.GetPTMMotifs())
    if GetDBCleavageFeatSet == True:
        features_dict.update (main_protein.GetCleavageCounts())

    if GetSubSeqSegs==True:
      features_dict.update(main_protein.Get_SubSeqParamScales_Features(window=4,
                                                                       edge=1.0,
                                                                       PickScales = MinScales_Dict,
                                                                       segs=SubSeqSegs))
    'Reduced alphabet representation(s) of protein + Kmer Freqs'
    if GetSeqReducedAlph==True:
      ReducedAlph_seq = get_reduced_prot (seq=main_seq,alph_name=SeqReducedAlph)
      # features_dict.update(ReducedAlph_seq.GetkgramFreq(k=ReducedK))
      features_dict.update(ReducedAlph_seq.GetKMirrorsFreq(k=ReducedK))


    if GetSeqReducedGroups==True:
      # alph_name='ofer_w8')
      ReducedAlphGroups_seq = get_reduced_prot (seq=main_seq,alph_name=SeqReducedGroups)
      features_dict.update(ReducedAlphGroups_seq.GetAA_Freq())
      features_dict.update (ReducedAlphGroups_seq.GetEntropy())

    # UseVariableCombinations - not Tested if working yet!!
    hp2_seq = get_reduced_prot (seq=main_seq,alph_name="hp2")
    features_dict.update (hp2_seq.GetkgramFreq(k=5))
    #features_dict.update (hp2_seq.GetEntropy())

    if GetTriLetterGroupAlphKFreq==True:
      selected3alph_3 = ['SolventA_3','Disorder_3','SecondaryStr_3','NormVDWV_3'] #]'NormVDWV_3' ,'Charge_3',,'Hydrophobicity_3' ]
      ' Get Profeat three letter group properties on whole sequence, (can be just main-seq)'
  ###    for alph_3 in THREE_LETTER_ALPH_NAMES: #THREE_LETTER_ALPH_NAMES - from Aalphabets, ProFEAT groups

      for alph_3 in selected3alph_3:
          seq_3 = get_reduced_prot (seq=seq,alph_name=alph_3)
          features_dict.update (seq_3.GetKMirrorsFreq(k=TriLetterGroupAlphK))

    return features_dict

'''
'TODO: Performance! '
CHECK is "files.items" indeed a generator object [avoid memory use]
#TODO - Seq length limit needs to be made dynamic!!
'TODO: Make multicore-map work. (Currently it never gets called for multifile multiclass fasta..)'
'TODO: Make Get_Protein_Feat flexible in terms of called features, add sending of params.. (EG - getNTail)'
'''
def GenFeaturesFromFasta(files=None, Dirr = '.', classType = 'dir'): #, Multiclass=True):

  '''
  Given file(s), extract features for each protein sequence in the file.
  Gets protein features for each protein in the given list of file(s),
  or gets list of fasta files from provided Dirr.

  Returns a dict of dicts.,
  {protid: {'feature':feature value}}.
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
#'Filter out nonstandard sequences: Illegal AA, short, or already processed: '
              #TODO - Seq length limit needs to be made dynamic!!
              if ((contain_illegals(sequence,ILLEGALS))==False) and (len(sequence)> MIN_PROTEIN_LENGTH) and (all_feature_dict.get(seq_id) is None):
                  all_feature_dict[seq_id]=Get_Protein_Feat(sequence)
                  if (classType == 'id'):
                      (all_feature_dict[seq_id])['classname'] = title.split()[1]
  # print ('type(files) %s' %(type(files)))

  
  for filename, classname in files.items():  #Is files.items a generator?
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

def get_reduced_prot(seq,alph_name="ofer14"):
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

def dict2df (seq_feature_dict, orientation='index',round=True,roundPrec=3):
  '''
  Takes a dict of dicts as input, returns a pandas DataFrame.
  NOTE: We also replace NaN with 0 ..
  As default, also rounds down data to {3} decimal places. (Saves on file size).

  '''
  from numpy import round
  Dframe = pd.DataFrame.from_dict(seq_feature_dict, orient=orientation)
  Dframe.replace([np.inf, -np.inf], 0) #Remove "infinities"
  Dframe.fillna(0, inplace=True) #Fill missing data
  if round==True:
    #Dframe.apply(npround)
    round(Dframe, decimals=roundPrec)
    #http://stackoverflow.com/questions/25272024/round-each-number-in-a-python-pandas-data-frame-by-2-decimals?rq=1
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
def Get_Dirr_All_Fasta (classType, dir_path='./'):
    '''
    Get all FASTA (*.fasta) files from current working directory,
    returns a list of files.
    If not additional param given, default is in current dir
    CURRENTLY - Does not get files from subdirectories.
    '''
    files_dict = {}
    '''
    We could also do:
    for file in glob("*.fasta"):
    '''
    for root, dirs, files in os.walk(dir_path) :
        for name in files:
            if (name.endswith(('.fasta','.fa'))):
                className = ''
                if classType == 'dir':
                    className = os.path.basename(root)
                elif classType == 'file':
                    filename, ext = os.path.splitext(name)
                    className = filename
                elif classType == 'id':
                    className = 'chooseID'
            files_dict[os.path.join(root, name)] = className
    return files_dict

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


def writeClassifiedFastas(classType,Dirr,resultsDir, df):
    fasta_files_dict = Get_Dirr_All_Fasta (classType,Dirr)
    classDict = {}
    writerDict = {}
    for key, value in fasta_files_dict.items():
        files = {key:value}
        for filename, classname in files.items():
            with open(filename) as fasta:
                for record in FastaIterator(fasta): #SeqIO.SimpleFastaParser(fasta):
                    title = record[0]
                    seq_id = title.split(None, 1)[0]
                    if (record.id in df.index):
                        classname = df[record.id]
                        if (classname not in writerDict):
                            classname = "".join([c for c in classname if c.isalpha() or c.isdigit() or c==' ']).rstrip()
                            file = resultsDir + '\\' + classname + '.fasta'
                            classHandle = open(file, "w")
                            classDict[classname] = classHandle
                            myWriter = FastaWriter(classDict[classname])
                            myWriter.write_header()
                            writerDict[classname] = myWriter
                        writerDict[classname].write_record(record)
    for classname, classHandle in classDict.items():
        writerDict[classname].write_footer()
        classDict[classname].close()



'Check file vs files... ? '
'TODO: Memory bloat likely here!'
'TODO: 1. MultiCore! 2. Save to disk every 10K'
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

    Also possible:
     http://pandas.pydata.org/pandas-docs/version/0.13.1/basics.html#renaming-mapping-labels
        http://stackoverflow.com/questions/19851005/rename-pandas-dataframe-index?rq=1
        dataframe.index.names = ['samples']
    '''
    def get_MultiClass_features(trainingSetFlag, classType):
        fasta_files_dict = Get_Dirr_All_Fasta (classType,Dirr)

        print ('Multiclass fasta_files list found: %s' %(list(fasta_files_dict) )) #Changed New
        'Warning: This way is very slow - the Dataframe structure should have its columns created ONCE! TODO - Fix!!'
        df = pd.DataFrame(columns=['proteinname','classname']) #CHANGED - Dan

        'TODO: Modify this bit, so Mult.files can be sent for processing at once/multiprocessing'
        for key, value in fasta_files_dict.items():
            'Get seperate DataFrame of features for all sequences in each seperate file:'
            # dataframe = dict2df(GenFeaturesFromFasta(files)) # ORIG?
            features = GenFeaturesFromFasta(files={key:value}, classType = classType) #,Dirr=Dirr)
            dataframe = dict2df(features)
            # dataframe['labels']=pd.Series(multiClassLabels[i],index=dataframe.index)
            """
            Set new index by "labels"
            http://www.bearrelroll.com/2013/05/python-pandas-tutorial/
            (We can also column-wise (axis=1) groupby - ['labels'] ) ;
            or:         df.index = ['two' for x in df.index]
            """
            # dataframe['labels'] = str(multiClassLabels[i])
            # print(file.split())
            filename = os.path.basename(key)
            classname = value
            dataframe['proteinname'] = [name for name in dataframe.index]
            if trainingSetFlag == True:
                if (classType != 'id'):
                    dataframe['classname'] = classname

            #dataframe.reindex(index, level=0)

            'Change this bit for multilevel / hierarchical index..'
            # dataframe.set_index('labels', append=False, inplace=True,drop=True) #ORIGINAL #WORKed
            ## dataframe.reindex('labels', append=True, inplace=True,drop=True)
            # dataframe.reorder_levels(['labels',dataframe.index.name])
            'Orig:'
            df = df.append(dataframe) #,inplace=True) #,ignore_index = True) #ORIG
            # df.append(dataframe,inplace=True) #CHANGED  - Dan

        # print(df)
        if trainingSetFlag == True:
            df.set_index(['proteinname','classname'], append=True, inplace=True,drop=True)
        else:
            df.set_index(['proteinname'], append=True, inplace=True,drop=True)
        return df

    if multiClass == True:
      # print ("Getting class dataframes")
      features = get_MultiClass_features(trainingSetFlag, classType)
      dataframe = features
    else:
      print ('Getting & merging features from subfiles')
      features = GenFeaturesFromFasta(files,Dirr)
      dataframe = dict2df(features)

    print('Features generated')
    # print('Data converted to dataframe, NaN cleaned')
    if saveCSV == True:
        print('Saving to disk.')
        dataframe.to_csv(saveName+'.csv')
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
      "TODO - set columns first? Dan"
      normParams = pd.concat([meanNormParams, stdNormParams], axis=1)
      normParams.columns = ['mean','std']
      z_df = df.apply(f)
      z_df.replace([np.inf, -np.inf], 0)
      z_df.fillna(0, inplace=True)
      if saveCSV:
          print("Data normalized and saved to",filename)
          z_df.to_csv(filename) #,index=False)
  else:   #NormParams exist :
      z_df = df.copy()
      for index, row in normParams.iterrows():
          "Check that feature is indeed present in OUR data. (Not absent):"#new - Dan.
          if index in z_df.columns: #Added - Dan
              # print()
              # print("index",index, "row",row)
              # index=str(index) #New. Dan. 22.1
              STD = row['std'] #Hack - avoid divide by 0 error
              if STD==0:
                  print(row,"STD=",row['std'])
                  STD=0.1
              # print("index: ",index)
              # print("z_df [index]: \n",z_df[index])

              # z_df[index] = ((z_df[index] - row['mean'])/row['std'])    # ORIG
              z_df[index] = (z_df[index] - row['mean'])/STD
          else:
            print(index," Feature not present (and removed).")

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

import params  #What is Params ??? #Dan.
def featExt(directory, trainingSetFlag, classType, normParams):

    #print("directory is: ",directory)
    df = get_features(trainingSetFlag, classType, returndf = True,
      saveCSV = False, saveName = "raw_Train_Feat",
     multiClass=True, Dirr = directory)

    'Removal of all zero columns - should not be done this way, disabled. Dan. (must be for TEST as well as Train'
    print("Removing any all zero features")

    # Training_Df =filterDF (Training_Df,removeZeroes=True,RemoveDuplicateCol=False) #OLD
    # Training_Df = Training_Df.groupby(Training_Df).filter(lambda x: x != 0) #OLD

    'Remove all zero features. Dan 2015. New'
    # df_cleaned = df[df.columns[(df != 0).any()]] #All zero columns removed #ORIG . #New - Used.
    df_cleaned = df[[col for col in df.columns if (df[col].std()>0)]] #NEW Dan
    df = df_cleaned  #Changed - was "df_cleaned=df"  mistakenly
    print("df.shape: ",df.shape)
    print("df_cleaned shape: ",df_cleaned.shape)
    # z_df = normalize_zscore(Training_Df,saveCSV = True,filename='Feat.csv')
    # MAD_df = normalize_median(Training_Df,saveCSV = True,filename='MAD_Feat.csv')

    # 'MAD then Z-score'
    ##    MAD_df = normalize_median(Training_Df,saveCSV = False,filename='MAD_Feat.csv')
    ##    df_cleaned = normalize_median(df_cleaned,saveCSV = False)

    if params.normalizeTrainingSetFlag == True:
        if trainingSetFlag == True:
            z_df, normParams = normalize_zscore(df_cleaned,saveCSV = False,filename='Z+MAD_Feat.csv')
            normParams.to_csv('trainingSetNormParams.csv')
        else:
            print("Loading stored NormParams")
            normParamDf = pd.read_csv(normParams, index_col=[0])
            z_df, normParams = normalize_zscore(df_cleaned,saveCSV = False,filename='Z+MAD_Feat.csv',normParams = normParamDf)
    else:
        z_df = df_cleaned

    # 'Removal of all zero columns - should not be done this way, disabled. Dan'
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

    #"http://stackoverflow.com/questions/22485375/efficiently-select-rows-that-match-one-of-several-values-in-pandas-dataframe?rq=1"
    #"df[df.Name.isin(['cytosol', 'golgi'])]"
    return z_df, normParams

def extract_features(trainingset, outputfeatures):
#def featExt(directory, trainingSetFlag, classType, normParams):

    df = get_features(trainingSetFlag, classType, returndf = True,
      saveCSV = False, saveName = "raw_Train_Feat",
        multiClass=True, Dirr = directory)

    'Remove all zero features. Dan 2015. New'
    df_cleaned = df[[col for col in df.columns if (df[col].std()>0)]] #NEW Dan
    df = df_cleaned
    print("df.shape: ",df.shape)
    print("df_cleaned shape: ",df_cleaned.shape)

    if params.normalizeTrainingSetFlag == True:
        if trainingSetFlag == True:
            z_df, normParams = normalize_zscore(df_cleaned,saveCSV=False,filename='Z+MAD_Feat.csv')
            normParams.to_csv('trainingSetNormParams.csv')
        else:
            print("Loading stored NormParams")
            normParamDf = pd.read_csv(normParams, index_col=[0])
            z_df, normParams = normalize_zscore(df_cleaned,saveCSV = False,filename='Z+MAD_Feat.csv',normParams = normParamDf)
    else:
        z_df = df_cleaned

    z_df_2 = z_df #NEW

    'Remove sample/sequence names, keep only file names as index'
    'http://stackoverflow.com/questions/18624039/pandas-reset-index-on-series-to-remove-multiindex'
    z_df_2.reset_index(0,drop=True,inplace=True) #Keep only class labels
    if (trainingSetFlag == True):
        z_df_2.to_csv('trainingSetFeatures.csv')
    else:
        z_df_2.to_csv('testingSetFeatures.csv')
    return z_df, normParams


if __name__=="__main__":
    #featExt()
    f = open('../Mamm_Organellas/Mammal_peroxisome9.fasta')
    record = next(SimpleFastaParser(f))
    sequence = record[1]
    d = Get_Protein_Feat(sequence)
    
