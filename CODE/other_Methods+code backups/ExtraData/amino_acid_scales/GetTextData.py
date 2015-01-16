'Extract AA scales from csvs, e.g. Georgiev scales reorganized as csv table (with headers+rows'

# import pandas as pandas
import sys
from sys import argv
import csv
'http://stackoverflow.com/questions/20200433/convert-csv-table-to-dictionary'
'http://stackoverflow.com/questions/17870022/read-two-column-csv-as-dict-with-1st-column-as-key'

'http://courses.cs.washington.edu/courses/cse140/13wi/csv-parsing.html'

def getListShape(list):
    """takes a 2 dimensional list object and returns the numbers of rows and cols"""
    datainfo = dict()
    datainfo["Row-count"]=len(list)
    datainfo["Col-count"]= len(list[0])
    return datainfo


if __name__=="__main__":
    #f='georgievAAScales.csv'
    f = str(argv[1])
    
    input_file = csv.DictReader(open(f), delimiter=',') #'\t' For .tsv files
    # scales={}
    i=1
    for row in input_file:
        #ORIG:
        #print (('gg_{} = {row}').format(i,row=row))
        print (('Atch_%i = %s' %(i,row)))
        
        #OLD, unused
        #print(("'gg_{}' : gg_{},").format(i,i)) #hack to make joint dict easy to save
        i +=1


#Needed to be cleaned with regex  (capture groups)- i) remove spaces per letter
# ii) ['](..[0-9]{1,4})[']   ->  $1   (remove parantheses around numbers using capture groups/s&replace)
# + [']([-]..[0-9]{1,4})[']   -> $1
# [']([0-9]{1,4})[']
