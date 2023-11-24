#! /usr/bin/python3
# downloads CollecTRI data form OmniPath database
# Args:
# 1 taxon short name, human|mouse|rat
# 2 path to the data download directory (with a trailing '/')

import sys
import decoupler as dc
import omnipath as op

txn = sys.argv[1]
datdir = sys.argv[2] # with a trailing '/'
subdir = 'ctri/'
ext = '.csv'

#dcctri = dc.get_collectri(organism=txn, split_complexes=True)
#wpth = datdir + subdir + 'dc-ctri' + ext
#dcctri.to_csv(wpth, index=False) # 4 fields

opctri = op.interactions.CollecTRI.get(genesymbols=True, organism=txn, loops=True)
wpth = datdir + subdir + txn + ext
opctri.to_csv(wpth, index=False) # 17 fields
