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
		log.Fatalln("uniprot: Failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 4 {
		log.Fatalln("uniprot: Expecting more arguments than:", cnt)
	}
	log.Println("uniprot: Started with args:", args)
	datdir := args[0]                              // path to data directory (with a trailing '/')
	bgwdir := args[1]                              // path to rdf directory (with a trailing '/')
	rpthT := args[2]                               // path to a list of selected taxa and proteomes
	scrdir := args[3]                              // path to scripts directory
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "_") // txnID -> proteomeID
	if err != nil {
		log.Fatalln("uniprot: Failed to make map:", rpthT, err)
	}
	n := len(txn2prm)
	if n == 0 {
		log.Fatalln("uniprot: Empty map:", rpthT)
	}

	if *aP || *dP {
		start := time.Now()
		if err := dat4bgw.SaveAllIdmap(datdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("uniprot: Downloaded idmapping in", time.Since(start))
		start = time.Now() // including 'humsavar.txt'
		if err := dat4bgw.SaveAllUniprot(datdir, txn2prm, scrdir); err != nil {
			panic(err)
		}
		log.Println("uniprot: Downloaded uniprot in", time.Since(start))
	}

	if *aP || *eP {
		start := time.Now()
		// rdf4bgw.Geneprot(datdir, bgwdir, txn2prm) // MUST be run before the others !!!
		if err := rdf4bgw.Geneprot(datdir, bgwdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("uniprot: Exported gene, prot, xmap in", time.Since(start))
		start = time.Now()
		if _, err := rdf4bgw.Gene2phen(datdir, bgwdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("uniprot: Exported gene2phen in", time.Since(start))
		start = time.Now()
		if _, err := rdf4bgw.Ortho(datdir, bgwdir, txn2prm); err != nil {
			panic(err)
		}
		log.Println("uniprot: Exported ortho in", time.Since(start))
	}
}
