#! /usr/bin/python3
# downloads CollecTRI data form OmniPath database
# Args:
# 1 taxon short name, human|mouse|rat TODO change to taxon IDs
# 2 path to the data download directory (with a trailing '/')

import sys
import decoupler as dc
import omnipath as op

txid = sys.argv[1]
datdir = sys.argv[2] # with a trailing '/'
taxa = {
	"9606": "human",
	"10090": "mouse",
	"10116": "rat"
}
txlbl = taxa[txid]
subdir = 'coltri/'
ext = '.csv'

#dcctri = dc.get_collectri(organism=txlbl, split_complexes=True)
#wpth = datdir + subdir + 'dc-ctri' + ext
#dcctri.to_csv(wpth, index=False) # 4 fields

opctri = op.interactions.CollecTRI.get(genesymbols=True, organism=txlbl, loops=True)
wpth = datdir + subdir + txid + ext
opctri.to_csv(wpth, index=False) # 17 fields
