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
		log.Fatalln("tfac2gene: failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 4 {
		log.Fatalln("tfac2gene: Expecting more arguments than:", cnt)
	}
	log.Println("tfac2gene: Started with args:", args)
	datdir := args[0]                              // path to data directory (with a trailing '/')
	bgwdir := args[1]                              // path to rdf directory (with a trailing '/')
	rpthT := args[2]                               // path to a list of selected taxa and proteomes
	scrdir := args[3]                              // path to scripts directory
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "_") // txnID -> proteomeID
	if err != nil {
		log.Fatalln("tfac2gene: Failed to make map:", rpthT, err)
	}
	n := len(txn2prm)
	if n == 0 {
		log.Fatalln("tfac2gene: Empty map:", rpthT)
	}

	if *aP || *dP {
		start := time.Now()
		if err := dat4bgw.SaveAllTflink(datdir); err != nil {
			panic(err)
		}
		log.Println("tfac2gene: Downloaded tflink in", time.Since(start))
		start = time.Now()
		if err := dat4bgw.SaveAllColtri(datdir, scrdir); err != nil {
			panic(err)
		}
		log.Println("tfac2gene: Downloaded collectri in", time.Since(start))
	}

	if *aP || *eP {
		start := time.Now()
		if _, err := rdf4bgw.Tfac2gene(datdir, bgwdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("tfac2gene: Exported in", time.Since(start))
	}
}
