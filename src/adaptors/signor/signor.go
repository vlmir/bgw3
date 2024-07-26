package main

import (
	"flag"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/dat4bgw"
	"github.com/vlmir/bgw3/src/rdf4bgw"
	"log"
	"time"
)

func main() {
	aP := flag.Bool("a", false, "[a]ll steps")
	dP := flag.Bool("d", false, "[d]ownload")
	eP := flag.Bool("e", false, "[e]xport")
	flag.Parse()
	if !flag.Parsed() {
		log.Fatalln("signor: Failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 2 {
		log.Fatalln("signor: Expecting more arguments than:", cnt)
	}
	log.Println("signor: Started with args:", args)
	datdir := args[0] // path to data directory (with a trailing '/')
	bgwdir := args[1] // path to rdf directory (with a trailing '/')

	if *aP || *dP {
		start := time.Now()
		if err := dat4bgw.SaveOneSignor("9606", datdir); err != nil {
			panic(err)
		}
		log.Println("signor: Downloaded in", time.Since(start))
	}

	if *aP || *eP {
		start := time.Now()
		if _, err := rdf4bgw.Reg2targ(datdir, bgwdir, bgw.Taxa4regs); err != nil {
			panic(err)
		}
		log.Println("signor: Exported reg2targ in", time.Since(start))
		start = time.Now()
		if _, err := rdf4bgw.Reg2pway(datdir, bgwdir, bgw.Taxa4regs); err != nil {
			panic(err)
		}
		log.Println("signor: Exported reg2pway in", time.Since(start))
	}
}
