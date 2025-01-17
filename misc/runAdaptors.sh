#!/bin/bash
bindir=$1 # directory containing the binaries, with a trailing slash
datdir=$2 # directory for downloading files, with a trailing slash
bgwdir=$3 # directory for RDF files, with a trailing slash
prmtxn=$4 # path to a file providing mappings between taxa and proteomes
scrdir=$5 # directory holding scripts, with a trailing slash
version=$6 # only for downloading OMIM ontology from BioPortal

## Downloading

$1uniprot -d $datdir $bgwdir $prmtxn $scrdir > up.outd 2> up.errd &
$1goa -d $datdir $bgwdir $prmtxn > goa.outd 2> goa.errd &
$1intact -d $datdir $bgwdir $prmtxn > ia.outd 2> ia.errd &
$1onto -d $datdir $bgwdir $version > onto.outd 2> onto.errd &
$1signor -d $datdir $bgwdir > sig.outd 2> sig.errd &
$1tfac2gene -d $datdir $bgwdir $scrdir > tf.outd 2> tf.errd &

## special case
cp $5atregnet.tgz $2
tar -xzf $2atregnet.tgz
rm $2atregnet.tgz

wait # essential

## Exporting

$1uniprot -e $datdir $bgwdir $prmtxn $scrdir > up.outx 2> up.errx # sic, foreground!
$1goa -e $datdir $bgwdir $prmtxn > goa.outx 2> goa.errx &
$1intact -e $datdir $bgwdir $prmtxn > ia.outx 2> ia.errx &
$1onto -e $datdir $bgwdir $version > onto.outx 2> onto.errx &
$1signor -e $datdir $bgwdir > sig.outx 2> sig.errx &
$1tfac2gene -e $datdir $bgwdir $scrdir > tf.outx 2> tf.errx &

wait # essential

## special case
cp $5biolink-model.nt.gz $3onto
