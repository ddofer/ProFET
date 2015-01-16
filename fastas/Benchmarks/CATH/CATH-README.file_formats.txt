----------------------------------------------------------------------------------------------------
FILE:     README.file_formats
AUTHOR:   Mark Dibley
CREATED:  05.05.2006
MODIFIED: ----------
VERSION:  v2.01

This README file lists the file formats of CATH data files and at the bottom has a summary of the changes between format 1.0 and 2.0 for each file format

File Formats
   1. Cath List File (CLF) Format 2.0
   2. Cath Names File (CNF) Format 2.0
   3. Cath Domain Description File (CDDF) Format 2.0
   4. Cath Domall File (CDF) Format 2.0
 
----------------------------------------------------------------------------------------------------

CATH List File (CLF) Format 2.0
-------------------------------
This file format has an entry for each structural entry in CATH.

Column 1:  CATH domain name (seven characters)
Column 2:  Class number
Column 3:  Architecture number
Column 4:  Topology number
Column 5:  Homologous superfamily number
Column 6:  S35 sequence cluster number
Column 7:  S60 sequence cluster number
Column 8:  S95 sequence cluster number
Column 9:  S100 sequence cluster number
Column 10: S100 sequence count number
Column 11: Domain length
Column 12: Structure resolution (Angstroms)
           (999.000 for NMR structures and 1000.000 for obsolete PDB entries)

Comment lines start with a '#' character.

Example:
--------
1oaiA00     1    10     8    10     1     1     1     1     1    59 1.000
1go5A00     1    10     8    10     1     1     1     1     2    69 999.000
1oksA00     1    10     8    10     2     1     1     1     1    51 1.800
1t6oA00     1    10     8    10     2     1     2     1     1    49 2.000
1cuk003     1    10     8    10     3     1     1     1     1    48 1.900
1hjp003     1    10     8    10     3     1     1     2     1    44 2.500
1c7yA03     1    10     8    10     3     1     1     2     2    48 3.100
1p3qQ00     1    10     8    10     4     1     1     1     1    43 1.700
1mn3A00     1    10     8    10     4     1     2     1     1    52 2.300
1nv8B01     1    10     8    10     5     1     1     1     1    71 2.200


CATH Domain Names
-----------------
The domain names have seven characters (e.g. 1oaiA00).

CHARACTERS 1-4: PDB Code
The first 4 characters determine the PDB code e.g. 1oai

CHARACTER 5: Chain Character
This determines which PDB chain is represented.
Chain characters of zero ('0') indicate that the PDB file has no chain field.

CHARACTER 6-7: Domain Number
The domain number is a 2-figure, zero-padded number (e.g. '01', '02' ... '10', '11', '12'). 
Where the domain number is a double ZERO ('00') this indicates that the domain is a whole PDB chain with no domain chopping. 


Hierachy Node Representatives
-----------------------------
Representative structural domains are selected from the CathDomainList based on 
the numbering scheme. For example the S35 sequence family representatives 
for superfamily 1.10.8.10 in the above example are 1oaiA00, 1oksA00, 1cuk003,
1p3qQ00 and 1nv8B01 as these are the first instances in the file with the same
superfamily number i.e. 1.10.8.10 but all have different S35 numbers (1 to 5).


----------------------------------------------------------------------------------------------------

CATH Names File (CNF) Format 2.0
--------------------------------
Description of each node in the CATH hierarchy for class, architecture, topology
and homologous superfamily levels.

Column 1:   Node number
Column 2:   Representative protein domain (CATH seven character domain name)
Free Text:  Node description (starting with ':' to indicate description start)
            (Note: this is not a single column)

Comment lines start with a '#' character.

Example:
--------
1                 1oaiA00       :Mainly Alpha
2                 1h8pA02       :Mainly Beta
1.10              1oaiA00       :Orthogonal Bundle
1.10.10           1i27A00       :Arc Repressor Mutant, subunit A
1.10.10.10        1i27A00       :"winged helix" repressor DNA binding domain
1.10.10.160       1pjr002       :


----------------------------------------------------------------------------------------------------

CATH Domain Description File (CDDF) Format 2.0
----------------------------------------------
Each entry corresponds with a CATH domain for a given release of the CATH database. 
Note: Different releases of CATH may have different domains definitions.
See below for CATH domain and segment naming conventions.

Comment lines start with a '#' character.

MAXIMUM of 80 characters per line (composed of tags that always a maximum of 
10 characters and the rest of line should be no longer than 70 characters).

Tags	    Description
-----------------------
FORMAT      Format definition (CDDF2.0) and first line of each entry
DOMAIN      CATH domain identifier - seven character code (e.g. 1abcA01)
PDBID	    PDB identifier - four character code
            (currently only used in PdbSumData files)
VERSION     CATH version number
VERDATE     CATH version release date
NAME	    PDB entry description
SOURCE      PDB entry organism/source
CATHCODE    CATH superfamily code C.A.T.H e.g. 1.10.10.10
CLASS	    Text description of class level (default: 'void')
ARCH	    Text description of architecture level (default: 'void')
TOPOL	    Text description of topology level (default: 'void')
HOMOL	    Text description of homologous superfam level (default: TOPOL entry)
DLENGTH     Length of the domain sequence
DSEQH	    Domain sequence header in FASTA format (e.g. '>pdb|1abcA01')
DSEQS	    Domain sequence string in FASTA format
NSEGMENTS   Number of segments that comprise the domain (integer)
SEGMENT     Segment identifier (e.g. 1abcA01:1:2)
SRANGE      Start and stop PDB residue identifiers that define range of segment
            (e.g. START=159  STOP=202)
SLENGTH     Length of the segment sequence
SSEQH	    Segment sequence header in FASTA format (e.g. '>pdb|1abcA01:1:2')
SSEQS	    Segment sequence string in FASTA format
ENDSEG      Signifies end of segment entry
COMMENTS    Text (optional line)
//	    Signifies end of entry


Example: Cath Domain Description File (CDDF)
--------------------------------------------
FORMAT    CDDF1.0
DOMAIN    9lprA01
VERSION   3.0.0
VERDATE   04-May-2006
NAME      Alpha-lytic protease (e.c.3.4.21.12) complex with methoxysuccinyl-*ala
NAME      -*ala-*pro-*leucine boronic acid
SOURCE    (Lysobacter enzymogenes 495) cloned and expressed in (escherichia coli
SOURCE    )
CATHCODE  2.40.10.10
CLASS     Mainly Beta
ARCH      Beta Barrel
TOPOL     Thrombin, subunit H
HOMOL     Trypsin-like serine proteases
DLENGTH   87
DSEQH     >pdb|9lprA01
DSEQS     IVGGIEYSINNASLCSVGFSVTRGATKGFVTAGHCGTVNATARIGGAVVGTFAARVFPGNDRAWVSLTSA
DSEQS     QTLLLQPILSQYGLSLV
NSEGMENTS 2
SEGMENT   9lprA01:1:2
SRANGE    START=16  STOP=115
SLENGTH   74
SSEQH     >pdb|9lprA01:1:2
SSEQS     IVGGIEYSINNASLCSVGFSVTRGATKGFVTAGHCGTVNATARIGGAVVGTFAARVFPGNDRAWVSLTSA
SSEQS     QTLL
ENDSEG
SEGMENT   9lprA01:2:2
SRANGE    START=231  STOP=242
SLENGTH   13
SSEQH     >pdb|9lprA01:2:2
SSEQS     LQPILSQYGLSLV
ENDSEG
//

NOTE:
The following CATH hierarchy description lines are typically found together 
(as found in the CathDomainList and CathNames files)

CATHCODE  2.40.10.10
CLASS     Mainly Beta
ARCH      Beta Barrel
TOPOL     Thrombin, subunit H
HOMOL     Trypsin-like serine proteases

The following domain sequence lines are typically found together 
(as found in the CathDomain Fasta Sequence File)

DLENGTH   87
DSEQH     >pdb|9lprA01
DSEQS     IVGGIEYSINNASLCSVGFSVTRGATKGFVTAGHCGTVNATARIGGAVVGTFAARVFPGNDRAWVSLTSA
DSEQS     QTLLLQPILSQYGLSLV

Segment sequence lines are always initiated with a 'SEGMENT' tag
and terminated with an 'ENDSEG' tag. The number of segments in the domain
always precedes the first segment using the 'NSEGMENTS' tag.

The following segment sequence lines are typically found together
(as found in the CathSegments Fasta Sequence File)

SEGMENT   9lprA01:1:2
SRANGE    START=16  STOP=115
SLENGTH   74
SSEQH     >pdb|9lprA01:1:2
SSEQS     IVGGIEYSINNASLCSVGFSVTRGATKGFVTAGHCGTVNATARIGGAVVGTFAARVFPGNDRAWVSLTSA
SSEQS     QTLL
ENDSEG


CATH Domain and Segment Naming Conventions
------------------------------------------
CATH Domain Names
-----------------
The domain names have seven characters (e.g. 1oaiA00).

CHARACTERS 1-4: PDB Code
The first 4 characters determine the PDB code e.g. 1oai

CHARACTER 5: Chain Character
This determines which PDB chain is represented.
Chain characters of zero ('0') indicate that the PDB file has no chain field.

CHARACTER 6-7: Domain Number
The domain number is a 2-figure, zero-padded number (e.g. '01', '02' ... '10', '11', '12'). Where the domain number is a double ZERO ('00') this indicates that the domain is a whole PDB chain with no domain chopping. 

CATH Segment Names
------------------
CATH segments (continuous regions of sequence within a domain) are described
adding colon separated numbers to the end of the domain name.
The first number is the sequential number of the segment.
The second number is the total number of segments in this domain.

1abcA01:1:2
xxxxxxxooooo

x = standard CATH six character domain name
o = segment information :ThisSegment:TotalSegments


----------------------------------------------------------------------------------------------------

CATH Domall File (CDF) Format 2.0
---------------------------------
The CATH Domall file describes domain boundaries for entries in the CATH database. 
All PDB chains in CATH that contain 1 or more domains have a CathDomall entry. Whole
chain domains can be identified where the number of domains is 1 and the
number of fragments is 0.

Comment lines start with a '#' character.

Segments are continuous sequence regions of domains.
Fragments are small regions of the protein chain that are excluded from the domain definition.

Column 1: Chain name (5 characters)
Column 2: Number of domains (formatted 'D%02d')
Column 3: Number of fragments (formatted 'F%02d')

The formatting of a CATH Domall file is best explained using examples.

Example CathDomall Entries
--------------------------

KEY:
N  = Number of segments
C  = Chain character
I  = Insert character/code ('-' indicates no insert character)
S  = Start PDB number
E  = End PDB number
NR = number of residues (fragment information only)

1chmA  D02 F00  1  A    2 - A  156 -  1  A  157 - A  402 -
                N |C    S I C    E I| N |C    S I C    E I|
               |<----Domain One---->|<-----Domain Two---->|
                  |<--Segment One-->|   |<--Segment One-->|

This translates to:
1chmA01 = Chain A; 2-156
1chmA02 = Chain A; 157-402

1cnsA  D02 F00  2  A    1 - A   87 -  A  146 - A  243 -  1  A   88 - A  145 -
                N |C    S I C    E I| C    S I C    E I| N |C    S I C    E I|
               |<--------------Domain One------------->|<-----Domain Two---->|
                  |<--Segment One-->|<---Segment Two-->|   |<--Segment One-->|

This translates to:
1cnsA01 = Chain A; 1-87, 146-243
1cnsA02 = Chain A; 88-145

Fragment Information
--------------------

Fragments are small regions of the protein chain that are not included
in the domain definition. These residue ranges are tagged on the end of the
segment information. The format is different from the segment range information.

1amg0  D02 F01  1  0    1 - 0  360 -  1  0  362 - 0  417 -  0  361 - 0  361 - (1)
                N |C    S I C    E I| N |C    S I C    E I| C    S I C    E I  NR|
               |<----Domain One---->|<-----Domain Two---->|<---Fragment One----->|
                  |<--Segment One-->|   |<--Segment One-->|

This translates to:
1amg001 = No chain character; 1-360
1amg002 = No chain character; 362-417
Fragment = 361

1bcmA0 D02 F02  1  A  257 - A  487 -  1  A  492 - A  559 -  A  488 - A  491 - (4)  A  560 - A  560 - (1)
                N |C    S I C    E I| N |C    S I C    E I| C    S I C    E I  NR| C    S I C    E I  NR|
               |<----Domain One---->|<-----Domain Two---->|<---Fragment One----->|<---Fragment Two----->|
                  |<--Segment One-->|   |<--Segment One-->|

This translates to:
1bcmA01 = Chain A; 257-487
1bcmA02 = Chain A; 492-559
Fragments = 488-491, 560


Cath Chain Names
----------------
The chain names have five characters (e.g. 1oaiA).

CHARACTERS 1-4: PDB Code
The first 4 characters determine the PDB code e.g. 1oai

CHARACTER 5: Chain Character
This determines which PDB chain is represented.
Chain characters of zero ('0') indicate that the PDB file has no chain field.

----------------------------------------------------------------------------------------------------

New to CLF Format 2.0
=====================
Domain ids are now 7 characters long to accommodate chains with more than 9 domains.

There is a new sequence family level at 60% sequence identity (column 7)

A count for each domain in a S100 sequence family is included (column 10) so that it is possible to represent domains uniquely using the CATH code.

C.A.T.H.S.O.L.I.D

C - Class
A - Architecture
T - Topology
H - Homologous Superfamily
S - Sequence Family (S35)
O - Orthogous Seqeuce Family (S60)
L - 'Like' Sequence Family (S95)
I - Identical (S100)
D - Domain (S100 count)

Clustering for the CLF uses a directed multi-linkage clustering algorithm and is order by increasing resolution, domain length and domain_id.


New to CNF Format 2.0
=====================
Zero padding has been removed from class and architecture node numbers.

Domain names are now 7 characters long rather than 6.

All nodes down to the Homologous Superfamily level are represented. Nodes
(homologous superfamily level only) that do not have a name are left blank.


New to CDDF Format 2.0
======================
Comments line is now optional

Domain names and segment names are based on the new CATH 7 character domain name format.

Where the HOMOL entry, naming the homologous superfamily level, is unnamed, the entry inherits the description from the parent topology level.


New to CDF Format 2.0
=====================
Now all chains that are in CATH are represented in the CATH Domall file whether
the chain has been chopped into domains or is a whole chain domain.

Whole chain domains can be identied where the number of domains is 1 and the
number of fragments is 0. 

Chain names are now represented by 5 characters instead of 6


New to CathChainList file
=========================
The CATH chain list now contains ALL chains that are in CATH rather than
just chopped chains. These chains are clustered in to sequence families
using a directed multi-linkage algorithm. Chain names use the full domain
name format of 7 characters where the domain number is '00' to indicate
'unchopped'.

----------------------------------------------------------------------------------------------------
END OF README
----------------------------------------------------------------------------------------------------
