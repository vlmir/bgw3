package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/semweb"
	"github.com/vlmir/bgw3/src/export"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"time"
)

// should stay here, needs to be passed to tfac2gene()
var uris4tftg = rdf.Uris4tftg
var msg = ""

func geneprot(datdir, bgwdir string, txn2pom util.Set2D) (ntg, ntp int, err error){
	ntg = 0
	ntp = 0
	for txid := range txn2pom {
		ext := ".upt"
		subdir := "uniprot/"
		rpthu := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
		subdir = "idmapping/"
		pomes := txn2pom[txid].Keys()
		pome := fmt.Sprintf("%s_%s", pomes[0], txid)
		ext = ".idmapping"
		rpthi := fmt.Sprintf("%s%s%s%s", datdir, subdir, pome, ext) // read
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
		upac2xrf, err := parse.Upidmap(rpthi, idmkeys)
		if err != nil {
			panic(err)
		}
		/////////////////////////////////////////////////////////////////////////////
		//rpthi = datdir + "intact.lst" // for filtering, TODO eliminate the hard coding
		dat4rdf, err := parse.Updat(rpthu, upac2xrf)
		if err != nil {
			panic(err)
		}
		updat := *dat4rdf.Udat
		txns := *dat4rdf.Txns
		if len(updat) == 0 || len(txns) == 0 {
			err := errors.New(fmt.Sprintf("%s: Not enough data", rpthu))
			panic(err)
		}
		if len(txns) > 1 {
			err := errors.New(fmt.Sprintf("%s: Multiple taxa: %v", rpthu, txns))
			panic(err)
		}
		/////////////////////////////////////////////////////////////////////////////
		// passing pointers, seems slightly faster, at most by 10%
		dat4rdf.Upac = &upac2xrf
		nlg, nlp, err := export.GeneProt(dat4rdf, wpthg, wpthp, wpthx)
		if err != nil {
			panic(err)
		}
		if nlg == 0 || nlp == 0 {
			err := errors.New(fmt.Sprintf("Taxon %s: No triples", txid))
			panic(err)
		}
		ntg += nlg
		ntp += nlp
	}
	return ntg, ntp, nil
}

func tfac2gene(datdir, bgwdir string, txn2pom util.Set2D, uris4tftg map[string]string) (int, error) {
	keys := []bgw.Column{
		{0, ":", 0, "--", 0, ""},
		{0, ":", 1, "--", 0, ""},
	}
	vals := []bgw.Column{
		{1, "|", 0, "|", 0, "uniprot"},
		{3, "|", 0, ";", 0, "pubmed"},
		{4, "|", 0, ";", 0, "confidence"},
		{5, "|", 0, ";", 0, "mode"},
	}
	graph := "tfac2gene"
	wdir := fmt.Sprintf("%s%s/", bgwdir, graph)
	nln := 0
	for txid := range txn2pom {
		if txid != "9606" {
			continue
		} // for now
		subdir := "tftg/"
		rdir := fmt.Sprintf("%s%s%s/", datdir, subdir, txid)
		ext := ".nt"
		wpth := fmt.Sprintf("%s%s%s", wdir, txid, ext)
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		dat4txn := make(map[string]util.Set3D)
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		if err != nil {
			panic(err)
		}
		upac2bgw := xmap.Upac
		gsym2bgw := xmap.Gsymb
		for src, uri := range uris4tftg {
			dat4txn[src] = make(util.Set3D)
			ext := ".f2g"
			rpth := fmt.Sprintf("%s%s%s", rdir, src, ext)
			duos, err := parse.Tab2set(rpth, keys, vals)
			if err != nil {
				panic(err)
			}
			duos.Add(src, "uri", uri)
			dat4txn[src] = duos
		} //src
		n, err := export.Tfac2gene(dat4txn, upac2bgw, gsym2bgw, wpth)
		if err != nil {
			panic(err)
		}
		if n == 0 {
			msg = fmt.Sprintf("%s: NoTriples", txid)
			panic(errors.New(msg))
		}
		nln += n
	} //txid
	return nln, nil
} //tfac2gene

func gene2phen(datdir, bgwdir string, txn2pom util.Set2D) {
	for txid := range txn2pom {
		if txid != "9606" {
			continue
		}
		subdir := "uniprot/"
		ext := ".var"
		rpthu := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
		subdir = "gene2phen/"
		ext = ".nt"
		wpth := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext)
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		if err != nil {
			panic(err)
		}
		gsym2bgw := xmap.Gsymb
		/////////////////////////////////////////////////////////////////////////////
		duos := parse.Upvar(rpthu, gsym2bgw)
		if len(duos) == 0 {
			msg = fmt.Sprintf("%s: NoData", rpthu)
			panic(errors.New(msg))
		}
		/////////////////////////////////////////////////////////////////////////////
		nts, err := export.Gene2phen(duos, gsym2bgw, wpth)
		if err != nil {
			panic(err)
		}
		if nts == 0 {
			msg = fmt.Sprintf("%s: NoTriples", txid)
			panic(errors.New(msg))
		}
	}
}

func prot2onto(datdir, bgwdir string, txn2pom util.Set2D, fx string) error {
	for txid := range txn2pom {
		subdir := "goa/"
		pth0 := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, fx) // read Goa data
		subdir = "xmap/"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json") // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		if err != nil {
			panic(err)
		}
		upac2bgw := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		bps := make(util.Set3D)
		ccs := make(util.Set3D)
		mfs := make(util.Set3D)
		if fx == ".gpa" {
			bps = parse.Gpa(pth0, upac2bgw)
		} else if fx == ".gaf" {
			bps, ccs, mfs = parse.Gaf(pth0, upac2bgw)
		} else {
			msg = fmt.Sprintf("%s: UnrecognizedFileExtension: %s", txid, fx)
			panic(errors.New(msg))
		}
		cnt := len(bps) + len(ccs) + len(mfs)
		if cnt == 0 {
			msg = fmt.Sprintf("%s: NoData", pth0)
			panic(errors.New(msg))
		}
		/////////////////////////////////////////////////////////////////////////////
		xpth := ""
		ext := ".nt"
		subdir = "prot2bp/"
		xpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntp, err := export.Prot2go(bps, upac2bgw, xpth)
		if err != nil {
			panic(err)
		}
		subdir = "prot2cc/"
		xpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntc, err := export.Prot2go(ccs, upac2bgw, xpth)
		if err != nil {
			panic(err)
		}
		subdir = "prot2mf/"
		xpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntf, err := export.Prot2go(mfs, upac2bgw, xpth)
		if err != nil {
			panic(err)
		}
		nts := ntp + ntc + ntf
		if nts == 0 {
			msg = fmt.Sprintf("%s: NoTriples", txid)
			panic(errors.New(msg))
		}
	}
	return nil
}

func prot2prot(datdir, bgwdir string, txn2pom util.Set2D) {
	for txid := range txn2pom {
		subdir := "intact/"
		ext := ".mit"
		pth0 := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext) // read IntAct tsv
		subdir = "prot2prot/"
		ext = ".nt"
		pth1 := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write ppi nt
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := bgw.NewXmap()
		err := xmap.Unmarshal(rpthx)
		if err != nil {
			panic(err)
		}
		upac2bgw := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		duos := parse.Mitab(pth0, upac2bgw)
		if len(duos) == 0 {
			msg = fmt.Sprintf("main.prot2prot():parse.Mitab():%s: NoData", txid)
			fmt.Printf("%s\n", msg)
			continue
		}
		/////////////////////////////////////////////////////////////////////////////
		nts, err := export.Prot2prot(duos, upac2bgw, pth1)
		if err != nil || nts == 0 {
			panic(err)
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
func ortho(datdir, bgwdir string, txn2pom util.Set2D) (int, error) {
	idmkeys := bgw.Orthokeys
	nln := 0
	for txidL := range txn2pom {
		for txidR := range txn2pom {
			if txidL >= txidR {
				continue
			} // skipping symmetrical and digonal
			duos, err := parse.Orthoduo(datdir, txidL, txidR, txn2pom, idmkeys)
			if err != nil {
				panic(err)
			}
			up2bgw := make(util.Set3D)
			subdir := "xmap/"
			ext := ".json"
			txids := [2]string{txidL, txidR}
			for _, txid := range txids {
				rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
				xmap := bgw.NewXmap()
				err := xmap.Unmarshal(rpthx)
				if err != nil {
					panic(err)
				}
				upac2bgw := xmap.Upac
				for acc, all := range upac2bgw {
					for src, one := range all {
						for id, _ := range one {
							up2bgw.Add(acc, src, id)
						}
					}
				}
			}
			subdir = "ortho/"
			ext = ".nt"
			file := fmt.Sprintf("%s%s%s%s", txidL, "-", txidR, ext)
			xpth := fmt.Sprintf("%s%s%s", bgwdir, subdir, file)
			nts, err := export.Ortho(duos, up2bgw, xpth)
			if nts == 0 {
				// taxon duos without matches already skiped in parse.Orthoduo
				msg = fmt.Sprintf("%v: NoContent", txids)
				//panic(errors.New(msg))
				fmt.Printf("main.ortho():%s\n", msg)
				continue
			}// TODO see why this happens
			nln += nts
		} // txidR
	} // txidL
	return nln, nil
}

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
		msg = fmt.Sprintf("flag.Parse(): Failed")
		panic(errors.New(msg))
	}
	args := flag.Args()
	n = len(args)
	m := 2
	if n < m {
		msg = fmt.Sprintf("Want at least: %d args have: %d", m, n)
		panic(errors.New(msg))
	}
	pth0 := *pP     // path to list of RefProts
	pth1 := args[0] // path to data directory (with a trailing '/')
	pth2 := args[1] // path to rdf directory (with a trailing '/')
	ptht := *tP     // path to a list of selected taxa
	log.Println("Started rdf4bgw with args:", args)
	/////////////////////////////////////////////////////////////////////////////
	/// used only by geneprot() for reading idmapping files
	tx2pm, err := util.Makemap(pth0, 1, 0, "_")
	if err != nil {
		panic(err)
	}
	n = len(tx2pm)
	if n == 0 {
		msg = fmt.Sprintf("%s: NoProteomes", pth0)
		panic(errors.New(msg))
	}
	log.Println("tx2pm:", n)
	////////////////////////////////////////////////////////////////////////////
	taxa4bgw, err := util.Makemap(ptht, 0, 1, ".")
	if err != nil {
		panic(err)
	}
	n = len(taxa4bgw)
	if n == 0 {
		msg = fmt.Sprintf("%s: NoData", ptht)
		panic(errors.New(msg))
	}
	log.Println("taxa4bgw:", n)
	/////////////////////////////////////////////////////////////////////////////
	txn2pom := make(util.Set2D)
	for txid, _ := range(taxa4bgw){
		poms, ok := tx2pm[txid]
		if !ok {continue}// filtering by ReferenceProteomes
		if l := len(poms); l > 1 {
			msg := fmt.Sprintf("main:%s: %d ReferenceProteomes", txid, l)
			fmt.Printf("%s,\n", msg)
		}// 20200531 download: never happens
		for pom, _ := range(poms){
			txn2pom.Add(txid, pom)
		}
	}
	if len(txn2pom) == 0 {
		msg = fmt.Sprintf("%s: NoProteomes", pth0)
		panic(errors.New(msg))
	}
	fmt.Println(txn2pom)
	/////////////////////////////////////////////////////////////////////////////
	if *aP || *eP {
		mystart := time.Now()
		geneprot(pth1, pth2, txn2pom) // MUST be run before the others !!!
		log.Println("Done with geneprot in", time.Since(mystart))
	}
	if *aP || *dP {
		mystart := time.Now()
		gene2phen(pth1, pth2, txn2pom)
		log.Println("Done with gene2phen in", time.Since(mystart))
	}
	if *aP || *gP {
		mystart := time.Now()
		ext := ".gaf"
		err = prot2onto(pth1, pth2, txn2pom, ext)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with prot2onto in", time.Since(mystart))
		}
	}
	if *aP || *iP {
		mystart := time.Now()
		prot2prot(pth1, pth2, txn2pom)
		log.Println("Done with prot2prot in", time.Since(mystart))
	}
	if *aP || *rP {
		mystart := time.Now()
		tfac2gene(pth1, pth2, txn2pom, uris4tftg)
		log.Println("Done with tfac2gene in", time.Since(mystart))
	}
	if *aP || *oP {
		mystart := time.Now()
		_, err := ortho(pth1, pth2, txn2pom)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with ortho in", time.Since(mystart))
		}
	}
	log.Println("Done with rdf4bgw in", time.Since(start))
}
