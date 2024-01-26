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
	if cnt < 3 {
		log.Fatalln("Expecting more arguments than ", cnt)
	}
	log.Println("Started with args:", args)
	datdir := args[0]                              // path to data directory (with a trailing '/')
	bgwdir := args[1]                              // path to rdf directory (with a trailing '/')
	rpthT := args[2]                               // path to a list of selected taxa and proteomes
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
		subdir := "signor/"
		if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
			panic(err)
		}
		dat4bgw.SaveOneSignor("9606", datdir)
		log.Println("Done with signor in", time.Since(start))
	}

	if *aP || *eP {
		start := time.Now()
		rdf4bgw.Rgr2trg(datdir, bgwdir, txn2prm)
		log.Println("Done with rgr2trg in", time.Since(start))
		start = time.Now()
		rdf4bgw.Reg2pway(datdir, bgwdir, txn2prm)
		log.Println("Done with reg2pway in", time.Since(start))
	}
}
