''' Cleans data from a pair of .csv files containing:
1) a list of protein names + sequences. "SP_domains"
2) A list of protein names + Start + Stop locations (Counting from 1)
of domain(s). - "SP_sequences"

We will output a fresh .csv file structured such that each row is a protein,
& columns are the sub-sequences: Tail N, Domain(s), Tail-C, Interdomain.
NOTE! Some of the above fields may be blank!

Additional: We may choose to further split at this stage of later:
2 files - one for single domain proteins (no interdomain region),
second for multidomain proteins.
And also - at the end (may be done "later" in Matlab) - the sum/concatenation
of all sequences in a column/type
'''

import sys      # imports the sys module
import csv     # imports the csv module 
from collections import defaultdict

'''
sp_domains = 'TEST_sp_domains.csv'
sp_sequences = 'test_SP_sequences.csv'
'''
sp_domains = 'SP_domains.csv'
sp_sequences = 'SP_sequences.csv'


P_NAME_COL = 0
DOMAIN_COL = 1
N_TAIL_COL = 2
C_TAIL_COL = 3
INTERDOM_LINKER_COL = 4

#import numpy
'''
http://www.pythonforbeginners.com/systems-programming/using-the-csv-module-in-python/
http://docs.python.org/3.3/library/csv.html#csv.DictReader
'''



def name_loc_join (domain_reader):
    locations_dict = defaultdict(list)
    for r in domain_reader:
        p_name = r[P_NAME_COL] #i-th row, first column
        'Make a corresponding list of (start,stop) tuples. Subtract 1 from locations for later list use' 
        loc = (int(r[1])-1,int(r[2])-1)
    ##        print ("p_name",p_name)
        locations_dict[p_name].append(loc)  
    return locations_dict


def slice_oneD (target_outfile,domain_locations,r):
    ''' Used when there's just one domain + Tails;'''
    Tot_sequence = r[1]
    N_Tail = Tot_sequence[0:domain_locations[0][0]] #if 0 then zero?
    DOM = Tot_sequence[domain_locations[0][0]:domain_locations[0][1]]
    C_Tail = Tot_sequence[domain_locations[0][1]:len(Tot_sequence)] #check +1 len?
    output = ','.join([r[P_NAME_COL],DOM,N_Tail,C_Tail])
    #Write out + extraline char
    target_outfile.write(output+'\n')            


def slice_multD (target_outfile,domain_locations,r):
    ''' Used when there are multiple domains'''
    Tot_sequence = r[1]
    d_count = len(domain_locations)
    #print(d_count)
    N_Tail = Tot_sequence[0:domain_locations[0][0]] #if 0 then zero?
    C_Tail = Tot_sequence[domain_locations[0][1]:len(Tot_sequence)]
    DOMS,LINKS='',''
    for i in range(d_count-2):
        # Stop before final domain
        DOMS += ''.join([Tot_sequence[domain_locations[i][0]:domain_locations[i][1]]])
        LINKS += ''.join([Tot_sequence[domain_locations[i][1]:domain_locations[i+1][0]]])
    DOMS +=Tot_sequence[domain_locations[d_count-1][0]:domain_locations[d_count-1][1]]
    output = ','.join([r[P_NAME_COL],DOMS,N_Tail,C_Tail,LINKS])
    #Write out + extraline char
    target_outfile.write(output+'\n')

def slice_uniD (target_outfile,r):
    '''
    Used when the protein is entirely composed of one domain. (No Tails)
       NOTE! It might be an issue with the downloaded data, such that
    this "type" could also include proteins with no domain locations
    '''
    Tot_sequence = r[1]
    output = ','.join([r[P_NAME_COL],Tot_sequence])
    target_outfile.write(output+'\n')                  

    
'For first line - column headers' 
def write_headers(target_outfile,subregions=4):
    if subregions==4:
        print("PROTEIN,DOMAINS,N_Tail,C_Tail,Linker", file=target_outfile) #Note! Domains vs domain here!
    elif subregions==3:
        print("PROTEIN,DOMAIN,N_Tail,C_Tail", file=target_outfile)
    else: # subregions = 1
        print("PROTEIN,DOMAIN", file=target_outfile)              
                         

'''First we read the start/stops of domains from the domains file
This will be stored in a dict, then used to "split" the corresponding
p.sequences in the protein-seq file. (Which will be read as a dict)
'''

'1. File with the p.name + domain(s):'
with open(sp_domains,"r") as domain_file:
    domain_reader = csv.reader(domain_file) #Open a csv reader on file.
##    if len(domain_locations) != len(domain_file):
##        print ("len(domain_locations) ",len(domain_locations),"different from num rows", len(domain_file))
    loc_dict=name_loc_join(domain_reader)
    
'We now have loc_dict with key=protein name, values = list of (start,stop) tups'
'Close file. (memory). Open next file containing sequences, for parsing'  
#print(loc_dict)

with open(sp_sequences,"r") as seq_file:
    
    ONEDom_output = open('output_singleD.csv', "w") #remem. to close at end!
    MULTDom_output = open('output_multiD.csv', "w")
    UNIDom_output = open('output_uniD.csv', "w")     # proteins with no dom. start/stop.  
    seq_reader = csv.reader(seq_file)
    write_headers(MULTDom_output,4)
    write_headers(ONEDom_output,3)
    write_headers(UNIDom_output,1)
                       
    for r in seq_reader:
# Write to the multidomain or single D file : Target. Check Num(Start,Stop)s - 2d pos' 
        #print("Now:",r)
        domain_locations = loc_dict.get(r[P_NAME_COL]) # get the values = a vector
        #print(r[P_NAME_COL])
        #print(domain_locations,len(domain_locations))
        if  domain_locations is None:
            slice_uniD(UNIDom_output,r)
            
        elif len(domain_locations)>=2:
               slice_multD(MULTDom_output,domain_locations,r)
        else:
               slice_oneD(ONEDom_output,domain_locations,r)
        #print(r[P_NAME_COL],"Done")

ONEDom_output.close()
MULTDom_output.close()
UNIDom_output.close()
