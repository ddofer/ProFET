#!/usr/local/bin/bash
# ------------------------------------------------------------------------
# download_domains.sh 
#
# Usage: download_domains.sh
#
# To be executed from protodev machine (ssh protodev prior to execution)
# Fields are
# 
# Written by: Nadav Rappoport
# Date: 23/02/2014
# ------------------------------------------------------------------------

# Run SQL query to download proteins' PFAM domains
# -------------------------------------------
psql -d dev -h /tmp << SQL

\t\a

-- Uncomment next line to get it in comma seperated format instead of tab delimited
--\f ','

\o SP_domains.tab

SELECT (SELECT protidinsource FROM proteins_release WHERE releaseid = 6 AND protid=p.protid) as protidinsource, p.motid, fromaa, toaa
FROM proteinmotif_release p, motifs_release m
WHERE p.releaseid = 6 and p.sysid = 12 AND m.releaseid = 6 AND m.sysid = p.sysid AND p.motifsourceid = 20 AND m.base_sysid = 12 AND m.motid = p.motid

-- Get only SwissProt (comment next line by '--' to remove this constrain)
--------------------------------------------------------------------------
 AND EXISTS (SELECT 1 FROM protein_actualsource_release pa WHERE pa.releaseid = 6 AND p.protid = pa.protid AND actual_sourceid = 1)
ORDER BY p.protid, fromaa, toaa, p.motid;

\o
\q
SQL

if [ $? -ne 0 ]
then
   echo "Could not download domains"
   exit -1
fi

exit 0
