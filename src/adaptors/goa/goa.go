package main

import (
	"flag"
	"github.com/vlmir/bgw3/src/dat4bgw"
	"github.com/vlmir/bgw3/src/rdf4bgw"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"time"
)

func main() {
	aP := flag.Bool("a", false, "[a]ll steps")
	dP := flag.Bool("d", false, "[d]ownload")
	eP := flag.Bool("e", false, "[e]xport")
	flag.Parse()
	if !flag.Parsed() {
		log.Fatalln("goa: Failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 3 {
		log.Fatalln("goa: Expecting more arguments than:", cnt)
	}
	log.Println("goa: Started with args:", args)
	datdir := args[0]                              // path to data directory (with a trailing '/')
	bgwdir := args[1]                              // path to rdf directory (with a trailing '/')
	rpthT := args[2]                               // path to a list of selected taxa and proteomes
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "_") // txnID -> proteomeID
	if err != nil {
		log.Fatalln("goa: Failed to make map:", rpthT, err)
	}
	n := len(txn2prm)
	if n == 0 {
		log.Fatalln("goa: Empty map:", rpthT)
	}

	if *aP || *dP {
		start := time.Now()
		if err := dat4bgw.SaveAllGaf(datdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("Done with goa in", time.Since(start))
	}

	if *aP || *eP {
		/*
		 */
		start := time.Now()
		ext := ".gaf"
		_, err = rdf4bgw.Prot2go(datdir, bgwdir, txn2prm, ext)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with goa in", time.Since(start))
		}
	}
}
