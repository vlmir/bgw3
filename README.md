# BioGateway RDF generator
This is a code base implemented in Golang (making use of a couple of ancillary Python scrypts) for converting arbitrary datasets into RDF.
It will download data from a specified set of data sources and generate RDF N-Triples files using a unified data model.

# Installation

Clone the repository with:
git clone git@github.com:vlmir/bgw3.git
cd bgw3/src/adapters
From within each subdirectory run:
go build
Move all the generated binaries to an external directory, e.g. ~/bin

# Usage

mkdir myProject DAT OUT
Note: myProject DAT OUT can be ANY directories outside bgw3
Note: DAT & OUT will hold downloaded files & generated RDFs respectively
cd myProject
cp bgw3/misc/runAdaptors.sh ./
Note: myProject will contain all STDOUTs and STDERRs
Note: you may need to make temporary changes to runAdaptors.sh

./runAdaptors.sh <path to binaries> <path for downloading files> <path for RDFs> <path to a file mapping taxa and proteomes> <path to scripts> <OMIM version>

_Example_
./runAdaptors.sh ~/bin/ DAT/ OUT/ ~/repos/bgw3/misc/prm_txn.txt ~/repos/bgw3/misc/ 25 &

### Output
DAT: a directory will be created for each data source containing the downloaded files
OUT: a directory will be created for ewach RDF graph containing N-Triples files
The files are ready to be loaded into any implementation of RDF/triple store
The file bgw3/misc/loadBGW.sh can be used for uploading into a Virtuoso Open Source database
