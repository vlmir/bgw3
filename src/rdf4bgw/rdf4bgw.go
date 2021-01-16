package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/export"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/semweb"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"time"
)

// should stay here, needs to be passed to tfac2gene()
var uris4tftg = rdf.Uris4tftg

func geneprot(datdir, bgwdir string, txn2prm util.Set2D) (ntg, ntp int, err error) {
	ntg = 0
	ntp = 0
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tgeneprot for:", txid)
		ext := ".upt"
		subdir := "uniprot/"
		rpthu := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
		subdir = "idmapping/"
		prmids := txn2prm[txid].Keys()
		prmid := fmt.Sprintf("%s_%s", prmids[0], txid)
		ext = ".idmapping"
		rpthi := fmt.Sprintf("%s%s%s%s", datdir, subdir, prmid, ext) // read
		ext = ".nt"
		subdir = "prot/"
		wpthp := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write prot nt
		subdir = "gene/"
		wpthg := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write gene nt
		ext = ".json"
		subdir = "xmap/"
		wpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write BGW map json
		/////////////////////////////////////////////////////////////////////////////
		idmkeys := bgw.Upkeys
		upac2xrf, err := parse.UpIdMap(rpthi, idmkeys)
		util.CheckE(err)
		/////////////////////////////////////////////////////////////////////////////
		//rpthi = datdir + "intact.lst" // for filtering, TODO eliminate the hard coding
		dat4rdf, err := parse.UpTab(rpthu, upac2xrf, txn2prm)
		util.CheckE(err)
		/////////////////////////////////////////////////////////////////////////////
		// passing pointers, seems slightly faster, at most by 10%
		//dat4rdf.Upac = &upac2xrf
		nlg, nlp, err := export.GeneProt(dat4rdf, wpthg, wpthp, wpthx)
		util.CheckE(err)
		ntg += nlg
		ntp += nlp
	}
	return ntg, ntp, nil
} // geneprot()

func tfac2gene(datdir, bgwdir string, txn2prm util.Set2D, uris4tftg map[string]string) (int, error) {
	keys := []bgw.Column{
		{0, ":", 0, "--", 0, ""},
		{0, ":", 1, "--", 0, ""},
	}
	vals := []bgw.Column{
		{1, "|", 0, "|", 0, "uniprot"},
		{2, "|", 0, "|", 0, "ncbig"},
		{3, "|", 0, ";", 0, "pubmed"},
		{4, "|", 0, ";", 0, "confidence"},
		{5, "|", 0, ";", 0, "mode"},
	}
	graph := "tfac2gene"
	wdir := fmt.Sprintf("%s%s/", bgwdir, graph)
	nln := 0
	for txid := range txn2prm {
		if txid != "9606" {
			continue
		} // for now
		log.Println("\n\ttfac2gene for:", txid)
		subdir := "tftg/"
		rdir := fmt.Sprintf("%s%s%s/", datdir, subdir, txid)
		ext := ".nt"
		wpth := fmt.Sprintf("%s%s%s", wdir, txid, ext)
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		dat4one := make(map[string]util.Set3D)
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		util.CheckE(err)
		//upac2bgw := xmap.Upac
		//gene2bgw := xmap.Lblg
		//gene2bgw := xmap.Ncbig
		for src, uri := range uris4tftg {
			ext := ".f2g"
			rpth := fmt.Sprintf("%s%s%s", rdir, src, ext)
			duos, err := parse.GetSetFromTab(rpth, keys, vals)
			if err != nil {
				msg := fmt.Sprintf("rdf4bgw.go:main.tfac2gene():%s: %s", err, src)
				log.Println(msg)
				//continue // SIC!
			}
			duos.Add(src, "uri", uri)
			dat4one[src] = make(util.Set3D)
			dat4one[src] = duos
		} //src
		if len(dat4one) == 0 {
			msg := fmt.Sprintf("rdf4bgw.go:main.tfac2gene():%s: NoData", txid)
			log.Println(msg)
		} // TODO delete, makes no sense
		//n, err := export.Tfac2gene(dat4one, upac2bgw, gene2bgw, wpth)
		n, err := export.Tfac2gene(dat4one, xmap, wpth)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.tfac2gene():%s: %s", err, txid)
			log.Println(msg)
			continue
		}
		nln += n
	} //txid
	return nln, nil
} //tfac2gene

func gene2phen(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	nln := 0
	for txid := range txn2prm {
		if txid != "9606" {
			continue
		}
		log.Println("\n\tgene2phen for:", txid)
		subdir := "uniprot/"
		ext := ".var"
		rpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
		subdir = "gene2phen/"
		ext = ".nt"
		wpth := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext)
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		util.CheckE(err)
		gsym2bgw := xmap.Lblg
		/////////////////////////////////////////////////////////////////////////////
		//duos, err := parse.UpVar(rpth, gsym2bgw)
		duos, err := parse.UpVar(rpth)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.gene2phen():%s: %s", err, txid)
			log.Println(msg)
		}
		/////////////////////////////////////////////////////////////////////////////
		nts, err := export.Gene2phen(duos, gsym2bgw, wpth)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.gene2phen():%s: %s", err, txid)
			log.Println(msg)
			continue
		}
		nln += nts
	}
	return nln, nil
} // gene2phen()

func prot2go(datdir, bgwdir string, txn2prm util.Set2D, fx string) (int, error) {
	nln := 0
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tprot2go for:", txid)
		subdir := "goa/"
		rpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, fx) // read Goa data
		subdir = "xmap/"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json") // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		util.CheckE(err)
		upac2bgw := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		bps := make(util.Set3D)
		ccs := make(util.Set3D)
		mfs := make(util.Set3D)
		if fx == ".gpa" {
			//out := parse.Gpa(rpth, upac2bgw) // TODO
		} else if fx == ".gaf" {
			//bps, ccs, mfs, err = parse.Gaf(rpth, upac2bgw)
			bps, ccs, mfs, err = parse.Gaf(rpth)
		} else {
			msg := fmt.Sprintf("%s: UnrecognizedFileExtension: %s", txid, fx)
			panic(errors.New(msg))
		}
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.prot2go():%s: %s", err, txid)
			log.Println(msg)
		}
		/////////////////////////////////////////////////////////////////////////////
		xpth := ""
		ext := ".nt"
		subdir = "prot2bp/"
		xpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntp, err := export.Prot2go(bps, upac2bgw, xpth)
		util.CheckE(err)
		subdir = "prot2cc/"
		xpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntc, err := export.Prot2go(ccs, upac2bgw, xpth)
		util.CheckE(err)
		subdir = "prot2mf/"
		xpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntf, err := export.Prot2go(mfs, upac2bgw, xpth)
		util.CheckE(err)
		nts := ntp + ntc + ntf
		if nts == 0 {
			msg := fmt.Sprintf("%s: NoTriples", txid)
			panic(errors.New(msg))
		}
		nln += nts
	}
	return nln, nil
} // prot2go()

func prot2prot(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	nlt := 0
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tprot2prot for:", txid)
		subdir := "intact/"
		ext := ".mit"
		rpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext) // read IntAct tsv
		subdir = "prot2prot/"
		ext = ".nt"
		wpth := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write ppi nt
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		util.CheckE(err)
		upac2bgw := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		//duos, err := parse.MiTab(rpth, upac2bgw)
		duos, err := parse.MiTab(rpth)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.prot2prot():%s: %s", err, txid)
			log.Println(msg)
		} // NoData
		/////////////////////////////////////////////////////////////////////////////
		nts, err := export.Prot2prot(duos, upac2bgw, wpth)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.prot2prot():%s: %s", err, txid)
			log.Println(msg)
			continue
		}
		nlt += nts
	}
	return nlt, nil
} // prot2prot()

/////////////////////////////////////////////////////////////////////////////
func ortho(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	log.Println("\n\tortho for:", "all")
	idmkeys := bgw.Orthokeys
	nln := 0
	for _, txidL := range txn2prm.Keys() {
		for _, txidR := range txn2prm.Keys() {
			if txidL >= txidR {
				continue
			} // skipping symmetrical and digonal
			txids := [2]string{txidL, txidR}
			duos, err := parse.OrthoDuo(datdir, txidL, txidR, txn2prm, idmkeys)
			if err != nil {
				msg := fmt.Sprintf("rdf4bgw.go:main.ortho():%s:%v", err, txids)
				log.Println(msg)
			} // NoData
			/////////////////////////////////////////////////////////////////////////////
			up2bgw := make(util.Set3D)
			subdir := "xmap/"
			ext := ".json"
			// building up2bgw for one pair of taxa
			for _, txid := range txids {
				rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
				xmap := bgw.NewXmap()
				err := xmap.Unmarshal(rpthx)
				util.CheckE(err)
				upac2bgw := xmap.Upac
				for upac, all := range upac2bgw {
					for src, one := range all {
						for id, _ := range one {
							up2bgw.Add(upac, src, id)
						}
					}
				}
			}
			subdir = "ortho/"
			ext = ".nt"
			file := fmt.Sprintf("%s%s%s%s", txidL, "-", txidR, ext)
			xpth := fmt.Sprintf("%s%s%s", bgwdir, subdir, file)
			nts, err := export.Ortho(duos, up2bgw, xpth)
			if err != nil {
				msg := fmt.Sprintf("rdf4bgw.go:main.ortho():%s:%v", err, txids)
				log.Println(msg)
				continue
			}
			nln += nts
		} // txidR
	} // txidL
	return nln, nil
} // end of orhto

////////////////////////////////////////////////////////////////////////////
func main() {
	start := time.Now()
	aP := flag.Bool("a", false, "export [a]ll")
	eP := flag.Bool("e", false, "export gene and protein [e]ntities")
	iP := flag.Bool("i", false, "export molecular [i]nteractions")
	dP := flag.Bool("d", false, "export [d]isease associations")
	gP := flag.Bool("g", false, "export [g]ene ontology annotations")
	rP := flag.Bool("r", false, "export [r]egulatory associations")
	oP := flag.Bool("o", false, "export [o]rthology relateions")
	pP := flag.String("p", "./proteomes.pls", "[p]roteome list")
	tP := flag.String("t", "./taxa.tls", "selected [t]axa")
	var n int
	flag.Parse()
	if !flag.Parsed() {
		msg := fmt.Sprintf("flag.Parse(): Failed")
		panic(errors.New(msg))
	}
	args := flag.Args()
	n = len(args)
	m := 2
	if n < m {
		msg := fmt.Sprintf("Want at least: %d args have: %d", m, n)
		panic(errors.New(msg))
	}
	datdir := args[0] // path to data directory (with a trailing '/')
	bgwdir := args[1] // path to rdf directory (with a trailing '/')
	rpthP := *pP      // path to list of RefProts
	rpthT := *tP      // path to a list of selected taxa
	log.Println("Started rdf4bgw with args:", args)
	/////////////////////////////////////////////////////////////////////////////
	tx2pm, err := util.MakeMap(rpthP, 1, 0, "_")
	util.CheckE(err)
	log.Println("AllRefProteomes:", len(tx2pm))
	////////////////////////////////////////////////////////////////////////////
	taxa4bgw, err := util.MakeMap(rpthT, 0, 1, ".")
	util.CheckE(err)
	log.Println("TaxaForBgw:", len(taxa4bgw))
	/////////////////////////////////////////////////////////////////////////////
	txn2prm := make(util.Set2D)
	for txid, _ := range taxa4bgw {
		prms, ok := tx2pm[txid]
		if !ok {
			continue
		} // filtering by ReferenceProteomes
		if l := len(prms); l != 1 {
			msg := fmt.Sprintf("main:%s: %d ReferenceProteomes", txid, l)
			panic(errors.New(msg))
		} // 20200531 download: never happens
		txn2prm.Add(txid, prms.Keys()[0])
	}
	if len(txn2prm) == 0 {
		msg := fmt.Sprintf("%s: NoProteomes", rpthP)
		panic(errors.New(msg))
	}
	/////////////////////////////////////////////////////////////////////////////
	if *aP || *eP {
		mystart := time.Now()
		geneprot(datdir, bgwdir, txn2prm) // MUST be run before the others !!!
		log.Println("Done with geneprot in", time.Since(mystart))
	}
	if *aP || *dP {
		mystart := time.Now()
		gene2phen(datdir, bgwdir, txn2prm)
		log.Println("Done with gene2phen in", time.Since(mystart))
	}
	if *aP || *gP {
		mystart := time.Now()
		ext := ".gaf"
		_, err = prot2go(datdir, bgwdir, txn2prm, ext)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with prot2go in", time.Since(mystart))
		}
	}
	if *aP || *iP {
		mystart := time.Now()
		prot2prot(datdir, bgwdir, txn2prm)
		log.Println("Done with prot2prot in", time.Since(mystart))
	}
	if *aP || *rP {
		mystart := time.Now()
		tfac2gene(datdir, bgwdir, txn2prm, uris4tftg)
		log.Println("Done with tfac2gene in", time.Since(mystart))
	}
	if *aP || *oP {
		mystart := time.Now()
		_, err := ortho(datdir, bgwdir, txn2prm)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with ortho in", time.Since(mystart))
		}
	}
	log.Println("Done with rdf4bgw in", time.Since(start))
}
