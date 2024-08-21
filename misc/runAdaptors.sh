#!/bin/bash
bindir=$1 # with a trailing slash
datdir=$2 # with a trailing slash
bgwdir=$3 # with a trailing slash
prmtxn=$4
scrdir=$5 # with a trailing slash
year=$6 # only for downloadin omim

#$1goa -d $datdir $bgwdir $prmtxn > goa.dout 2> goa.derr &
#$1intact -d $datdir $bgwdir $prmtxn > ia.dout 2> ia.derr &
#$1onto -d $datdir $bgwdir $year > onto.dout 2> onto.derr &
#$1signor -d $datdir $bgwdir > sig.dout 2> sig.derr &
#$1tfac2gene -d $datdir $bgwdir $scrdir > tf.dout 2> tf.derr &
#$1uniprot -d $datdir $bgwdir $prmtxn $scrdir > up.dout 2> up.derr &
#wait

#$1uniprot -e $datdir $bgwdir $prmtxn $scrdir > up.xout 2> up.xerr &
wait
#$1goa -e $datdir $bgwdir $prmtxn > goa.xout 2> goa.xerr &
#$1intact -e $datdir $bgwdir $prmtxn > ia.xout 2> ia.xerr &
#$1onto -e $datdir $bgwdir $year > onto.xout 2> onto.xerr &
#$1signor -e $datdir $bgwdir $prmtxn > sig.xout 2> sig.xerr &
#$1tfac2gene -e $datdir $bgwdir $scrdir > tf.xout 2> tf.xerr &
wait
