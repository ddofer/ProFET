#!/usr/local/bin/bash
# ------------------------------------------------------------------------
# download_sequences.sh 
#
# Usage: download_keywords.sh
#
# To be executed from protodev machine (ssh protodev prior to execution)
# Written by: Nadav Rappoport
# Date: 23/02/2014
# ------------------------------------------------------------------------

# Run SQL query to download protein sequences
# -------------------------------------------
psql -d dev -h /tmp << SQL

\t\a

\f ','

\o SP_sequences.csv

select protidinsource, protseq
from proteins_release p, proteinsequences_release s 
where p.releaseid = 6 AND s.releaseid = 6 AND p.protid = s.protid AND EXISTS (select 1 from protein_actualsource_release pa WHERE pa.releaseid = 6 AND p.protid = pa.protid AND actual_sourceid = 1)
order by p.protid;

\o
\q
SQL

if [ $? -ne 0 ]
then
   echo "Could not download sequences"
   exit -1
fi

exit 0
