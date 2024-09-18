#!/bin/bash
bindir=$1 # directory containing the binaries, with a trailing slash
datdir=$2 # directory for downloading files, with a trailing slash
bgwdir=$3 # directory for RDF files, with a trailing slash
prmtxn=$4 # path to a file providing mappings between taxa and proteomes
scrdir=$5 # directory holding scripts, with a trailing slash
version=$6 # only for downloading OMIM ontology from BioPortal, currently 25

## Downloading

#$1uniprot -d $datdir $bgwdir $prmtxn $scrdir > up.dout 2> up.derr &
#$1goa -d $datdir $bgwdir $prmtxn > goa.dout 2> goa.derr &
#$1intact -d $datdir $bgwdir $prmtxn > ia.dout 2> ia.derr &
#$1onto -d $datdir $bgwdir $version > onto.dout 2> onto.derr &
#$1signor -d $datdir $bgwdir > sig.dout 2> sig.derr &
#$1tfac2gene -d $datdir $bgwdir $scrdir > tf.dout 2> tf.derr &
wait # essential

## Exporting

#$1uniprot -e $datdir $bgwdir $prmtxn $scrdir > up.xout 2> up.xerr &
wait # essential
#$1goa -e $datdir $bgwdir $prmtxn > goa.xout 2> goa.xerr &
#$1intact -e $datdir $bgwdir $prmtxn > ia.xout 2> ia.xerr &
#$1onto -e $datdir $bgwdir $version > onto.xout 2> onto.xerr &
#$1signor -e $datdir $bgwdir > sig.xout 2> sig.xerr &
#$1tfac2gene -e $datdir $bgwdir $scrdir > tf.xout 2> tf.xerr &
wait # essential
