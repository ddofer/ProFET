'''
Go through a given fasta file (later - mod for sets of fastas),
and output a new fasta file, with sequences containg unknown, or non standard AA
removed; too short sequences removed.

Later, can be used to filter sequences whose ID is a classname; (keeping those
 with a minimum amount of examples, e.g. 30+ samples per CATH classs/domain/..)

Look for:
    5.1.2 Iterating over the records in a sequence ﬁle - BioPy Cookbook
    16.1.1 Filtering a sequence ﬁle  (BioPython cookbook)
    Filtering biopython seqio.
    5.5 Writing Sequence Files   - BioPy cookbook

REGEX to remove last letter:
([a-z]+[.][0-9]{1,3}[.][0-9]{1,3})[.][0-9]{1,3}
$1
'''
from sys import argv
import os
from Bio import SeqIO
from collections import Counter
from Bio.Seq import Seq

ILLEGALS = ['B', 'J', 'Z', 'X', 'U', 'O', 'Z']

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

"http://nbviewer.ipython.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/11%20-%20Lesson.ipynb"
def Get_Dirr_All_Fasta (Dirr = '.'):
  '''
  Get all FASTA (*.fasta) files from current working directory,
  returns a list of files.
  If not additional param given, default is in current dir
  CURRENTLY - Does not get files from subdirectories.
  '''

  '''old -   for f in os.listdir(sys.argv[1]) :
   We could also do:
  for file in glob("*.fasta"):
  '''
  if Dirr != '.':
    os.chdir(str(Dirr))
    print ("dirr change to: %s" %(Dirr))

  for f in os.listdir(os.curdir) : #If nothing given, def
    files = [f for f in os.listdir(os.curdir) if (os.path.isfile(f) and f.endswith(".fasta"))]
  return files


def FilterFastaSeqGroups (Dir,fName,minGroupCount = 45,minLength= 27,FilterGroups=True):
    # minGroupCount = 40 #Minimum amount of samples per class.
    # minLength= 29 #Minimum length of a protein sequence.

    KeywordCount = Counter() #Hhold counts of encoutnered classes/KW = second word in fasta description. Used to filter
    os.chdir(str(Dir))
    # testHandle=Dir+fName
    # testHandle=fName
    # handle = open(sys.argv[1])
    # handle = testHandle
    handle=fName
    modHandle = 'FILT'+fName
    '5.1.2 Iterating over the records in a sequence ﬁle - BioPy Cookbook; Modified for Python 3'

    def FiltSequences (handle,minLength=27):
        '''
        Filter and output a set of fasta files,
        creates a new multifasta file containing only sequences of min. length
        and not containing illegal AA.
        Also counts + updates occurences of KW/Group-ID (second word in descriptor).
        '''
        it = SeqIO.parse(handle, "fasta")
        i=0
        filteredRecords = []
        for seq in it:
            i +=1
            if contain_illegals(seq=seq.seq)==False:
                if len(seq.seq)>minLength:
                    filteredRecords.append(seq)
                    KW = seq.description.split()[1]
                    KeywordCount[KW] +=1
        print(i)
        print('KeywordCount', KeywordCount)
        SeqIO.write(filteredRecords, modHandle, "fasta")

    def FiltSeqGroups (handle,Classes):
        '''
        Keeps sequences whose class is in pre-approved list (min. amount of counts);
        Class/KW is second word of seq's description.
        '''

        it = SeqIO.parse(handle, "fasta")
        filteredRecords = []
        i=0
        for seq in it:
            KW = seq.description.split()[1]
            seq.description=KW #remove extra info. LOSSY
            if KW in Classes:
                filteredRecords.append(seq)
                i += 1
        SeqIO.write(filteredRecords, modHandle, "fasta")
        print('Filtered left:', i)


    FiltSequences(handle=handle,minLength=minLength)

    'Keep only those classes with a minimal amount of samples'
    if FilterGroups==True:
        FiltGroup = set()
        for k,v in KeywordCount.items():
            if v>=minGroupCount:
                FiltGroup.add(k)
        print(FiltGroup)

        FiltSeqGroups(handle=modHandle,Classes=FiltGroup)


if __name__=="__main__":
    # Dir = r"D:\SkyDrive\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\TestFilt"

    #minCounts = int(input('Input Min Count'))
    #Dir = str(input('Input target directory with fastas to filter. \n'))

    # print(len(argv))
    print ("Input should be: Dirr with target fasta, MinCountPerClass , minLengthPerSeq")
    Dir = str(argv[1])

    if argv[2] is None:
        minCounts=60
    else:
        minCounts=int(argv[2])

        if argv[3] is None:
          minLength=40
        else:
          minLength=int(argv[3])
# minLength=input('Enter minimum protein length')

    FastaFiles = Get_Dirr_All_Fasta(Dirr=Dir)
    for f in FastaFiles:
        FilterFastaSeqGroups (Dir,f,minGroupCount = minCounts,minLength= minLength,FilterGroups=True)
