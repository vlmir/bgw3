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
		log.Fatalln("intact: Failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 3 {
		log.Fatalln("intact: Expecting more arguments than:", cnt)
	}
	log.Println("intact: Started with args:", args)
	datdir := args[0]                              // path to data directory (with a trailing '/')
	bgwdir := args[1]                              // path to rdf directory (with a trailing '/')
	rpthT := args[2]                               // path to a list of selected taxa and proteomes
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "_") // txnID -> proteomeID
	if err != nil {
		log.Fatalln("intact: Failed to make map:", rpthT, err)
	}
	n := len(txn2prm)
	if n == 0 {
		log.Fatalln("intact: Empty map:", rpthT)
	}

	if *aP || *dP {
		start := time.Now()
		if err := dat4bgw.SaveAllIntact(datdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("Done with intact in", time.Since(start))
	}

	if *aP || *eP {
		/*
		 */
		start := time.Now()
		rdf4bgw.Prot2prot(datdir, bgwdir, txn2prm)
		log.Println("Done with intact in", time.Since(start))
	}
}
