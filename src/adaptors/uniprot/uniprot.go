package main

import (
	"flag"
	"github.com/vlmir/bgw3/src/dat4bgw"
	"github.com/vlmir/bgw3/src/rdf4bgw"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"os"
	"path/filepath"
	"time"
)

func main() {
	aP := flag.Bool("a", false, "[a]ll steps")
	dP := flag.Bool("d", false, "[d]ownload")
	eP := flag.Bool("e", false, "[e]xport")
	flag.Parse()
	if !flag.Parsed() {
		log.Fatalln("failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 4 {
		log.Fatalln("Expecting more arguments than ", cnt)
	}
	log.Println("Started with args:", args)
	datdir := args[0]                              // path to data directory (with a trailing '/')
	bgwdir := args[1]                              // path to rdf directory (with a trailing '/')
	rpthT := args[2]                               // path to a list of selected taxa and proteomes
	scrdir := args[3]                              // path to scripts directory
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "_") // txnID -> proteomeID
	if err != nil {
		log.Fatalln("main:", err)
	}
	n := len(txn2prm)
	if n == 0 {
		log.Fatalln("main:Empty map:", rpthT)
	}
	log.Println("txn2prm:", n)

	if *aP || *dP {
		start := time.Now()
		subdir := "idmapping/"
		if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
			panic(err)
		}
		dat4bgw.SaveAllIdmap(datdir, txn2prm)
		log.Println("Done with idmapping in", time.Since(start))

		subdir = "uniprot/"
		if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
			panic(err)
		}

		uri := "http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
		wpth := datdir + "uniprot/9606.var"
		if _, err := dat4bgw.HttpFile(uri, wpth); err != nil {
			panic(err)
		}

		dat4bgw.SaveAllUniprot(datdir, txn2prm, scrdir)
		log.Println("Done with uniprot in", time.Since(start))
	}

	if *aP || *eP {
		/*
		 */
		start := time.Now()
		rdf4bgw.Geneprot(datdir, bgwdir, txn2prm) // MUST be run before the others !!!
		log.Println("Done with Geneprot in", time.Since(start))
		start = time.Now()
		rdf4bgw.Gene2phen(datdir, bgwdir, txn2prm)
		log.Println("Done with Gene2phen in", time.Since(start))
	}
}
