# BioGateway RDF generator
This is a code base implemented in Golang (making use of a couple of ancillary Python scrypts) for converting arbitrary datasets into RDF.
It will download data from a specified set of data sources and generate RDF N-Triples files using a unified data model.

# Installation


Clone the repository with:
git clone git@github.com:vlmir/bgw3.git

cd bgw3/src
go test ./... # testing all source files, MUST pass

cd bgw3/src/adapters/

From within each subdirectory run:
go build

Move all the generated binaries to an external directory, e.g. ~/bin

# Usage


mkdir myProject DAT OUT
**Note: myProject DAT OUT may be ANY directories OUTSIDE bgw3
**Note: DAT & OUT will hold downloaded files & generated RDFs respectively

cd myProject
**Note: myProject will contain all STDOUTs and STDERRs

cp bgw3/misc/runAdaptors.sh ./
**Note: if you need to make changes to runAdaptors.sh it is more practical to do it here

Run the following command in the current directory to build a new set of RDF files:
./runAdaptors.sh path_to_binaries path_for_downloaded_files path_for_RDFs path_to_a_file_mapping_taxa_and_proteomes path_to_scripts OMIM_version

_Example_

./runAdaptors.sh ~/bin/ DAT/ OUT/ ~/repos/bgw3/misc/prm_txn.tsv ~/repos/bgw3/misc/ 25 &

### Output

subdirectories will be created in DAT for each data source filled with downloaded files

subdirectories will be created in OUT for each RDF graph filled with N-Triples files

The files are ready to be loaded into any implementation of RDF/triple store

The file bgw3/misc/loadBGW.sh can be used for uploading RDFs into a Virtuoso Open Source database
