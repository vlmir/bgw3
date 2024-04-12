package main

import (
	"flag"
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
		log.Fatalln("onto: Failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 3 {
		log.Fatalln("onto: Expecting more arguments than:", cnt)
	}
	log.Println("onto: Started with args:", args)
	datdir := args[0] // path to data directory (with a trailing '/')
	bgwdir := args[1] // path to rdf directory (with a trailing '/')
	year := args[2]   //

	if *aP || *dP {
		start := time.Now()
		if err := dat4bgw.SaveAllOnto(datdir, year); err != nil {
			panic(err)
		}
		log.Println("onto: Downloaded in", time.Since(start))
	}

	if *aP || *eP {
		start := time.Now()
		if err := rdf4bgw.Onto(datdir, bgwdir); err != nil {
			panic(err)
		}
		log.Println("onto: Exported in", time.Since(start))
	}
}
