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
	datdir := args[0]                               // path to data directory (with a trailing '/')
	bgwdir := args[1]                               // path to rdf directory (with a trailing '/')
	rpthT := args[2]                                // path to a list of selected taxa and proteomes
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "\t") // Set2D: map txnID -> proteomeID
	if err != nil {
		panic(err) // sic
	}

	if *aP || *dP {
		start := time.Now()
		if err := dat4bgw.SaveAllIntact(datdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("intact: Downloaded in", time.Since(start))
	}

	if *aP || *eP {
		start := time.Now()
		if _, err := rdf4bgw.Prot2prot(datdir, bgwdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("intact: Exported in", time.Since(start))
	}
}
