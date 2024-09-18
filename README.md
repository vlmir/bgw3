# BioGateway RDF generator
This is a code base implemented in Golang (making use of a couple of ancillary Python scrypts) for converting arbitrary datasets into RDF.
It will download data from a specified set of data sources and generate RDF N-Triples files using a unified data model.

# Installation

Clone the repository with:
git clone git@github.com:vlmir/bgw3.git

cd bgw3/src/adapters

From within each subdirectory run:
go build

Move all the generated binaries e.g. to ~/bin

# Usage

mkdir myProject # can be any directory outside bgw3
cd myProject
mkdir DAT OUT # for holding downloaded files and generated RDFs, any names will do
cp bgw3/misc/runAdaptors.sh ./



./runAdaptors.sh <path to binaries> <path for downloading files> <path for RDFs> <path to a file mapping taxa and proteomes> <path to scripts> <OMIM version>

_Example_
./runAdaptors.sh ~/bin/ DAT/ OUT/ ~/repos/bgw3/misc/prm_txn.txt ~/repos/bgw3/misc/ 25 &
