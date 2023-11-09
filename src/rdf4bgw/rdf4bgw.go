package main

import (
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/export"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/semweb"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"os"
	"time"
)

func geneprot(datdir, bgwdir string, txn2prm util.Set2D) (err error) {
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
		var xmap bgw.Xmap
		xmap.New()
		err := export.Gene(rpthu, rpthi, wpthg, &xmap)
		util.CheckE(err)
		err = export.Prot(rpthu, rpthi, wpthp, &xmap)
		util.CheckE(err)
		// xmap export
		wfhX, err := os.Create(wpthx)
		util.CheckE(err)
		j, err := json.MarshalIndent(&xmap, "", " ")
		util.CheckE(err)
		wfhX.Write(j)
	}
	return nil
} // geneprot()

func reg2pway(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	log.Println("\n\reg2pway for:", "all")
	cnts := make(util.Set2D)
	var pdcks = []string{
		"reg2ptrg",
		"reg2ntrg",
		"reg2utrg",
	}
	for srck, _ := range rdf.Uris4rgrtrg {
		// define keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		ext := ""
		rpth := ""
		if srck == "signor" {
			keys, vals = bgw.SigPwaysParseConf()
			ext = ".tsv"
		}

		for txid := range txn2prm {
			if txid != "9606" {
				continue
			} // for now
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()
			if srck == "signor" {
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", "pathways", ext)
			}
			log.Println("Rdf4bgw.reg2pway(): processing", rpth)
			err := parse.Tab2struct(rpth, keys, vals, &d4b)
			if err != nil {
				log.Printf("%s%s", "reg2pway:parse.Tab2struct: ", err)
				continue // sic!
			}
			// d4b is now loaded with data
			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			ext := ".json"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
			err = xmap.Unmarshal(rpthx)
			util.CheckE(err)

			if srck == "signor" {
				// generating map signor-id -> entitity-ids

				sigmap := make(util.Set3D)
				rdir := fmt.Sprintf("%s%s%s", datdir, srck, "/")
				smpths := []string{
					rdir + "complexes.map",
					rdir + "families.map",
				}
				err := parse.Sig2up(sigmap, smpths)
				if err != nil {
					log.Printf("%s%s", "parse.Sig2up: ", err)
					continue // sic!
				}
				xmap.Signor = sigmap
			}

			d4b.Src = srck
			d4b.Taxid = txid
			err = export.SigPways(&d4b, &xmap, bgwdir)
			if err != nil {
				//panic(err)
				log.Println(err)
				continue
			}
			for _, pdck := range pdcks {
				cnts.Add(pdck, srck)
				cnts[pdck][srck] = d4b.Cnts[pdck][srck]
			}
		}
	}
	return cnts, nil
} // reg2pway

func rgr2trg(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	log.Println("\n\trgr2trg for:", "all")
	cnts := make(util.Set2D)
	var pdcks = []string{
		"reg2ptrg",
		"reg2ntrg",
		"reg2utrg",
	}
	for srck, _ := range rdf.Uris4rgrtrg {
		// define keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		ext := ""
		rpth := ""
		if srck == "signor" {
			keys, vals = bgw.SignorParseConf()
			ext = ".mi28"
		}

		for txid := range txn2prm {
			if txid != "9606" {
				continue
			} // for now
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()
			if srck == "signor" {
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", txid, ext)
			}
			log.Println("Rdf4bgw.rgr2trg(): processing", rpth)
			err := parse.Tab2struct(rpth, keys, vals, &d4b)
			if err != nil {
				log.Printf("%s%s", "rgr2trg:parse.Tab2struct: ", err)
				continue // sic!
			}
			// d4b is now loaded with data
			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			ext := ".json"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
			err = xmap.Unmarshal(rpthx)
			util.CheckE(err)

			if srck == "signor" {
				// generating map signor-id -> entitity-ids

				sigmap := make(util.Set3D)
				rdir := fmt.Sprintf("%s%s%s", datdir, srck, "/")
				smpths := []string{
					rdir + "complexes.map",
					rdir + "families.map",
				}
				err := parse.Sig2up(sigmap, smpths)
				if err != nil {
					log.Printf("%s%s", "parse.Sig2up: ", err)
					continue // sic!
				}
				xmap.Signor = sigmap
			}

			d4b.Src = srck
			d4b.Taxid = txid
			err = export.Rgr2trg(&d4b, &xmap, bgwdir)
			if err != nil {
				//panic(err)
				log.Println(err)
				continue
			}
			for _, pdck := range pdcks {
				cnts.Add(pdck, srck)
				cnts[pdck][srck] = d4b.Cnts[pdck][srck]
			}
		}
	}
	return cnts, nil
} // rgr2trg

func tfac2gene(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	log.Println("\n\ttfac2gene for:", "all")
	cnts := make(util.Set2D)
	var pdcks = []string{
		"reg2ptrg",
		"reg2ntrg",
		"reg2utrg",
	}
	for srck, _ := range rdf.Uris4tftg {
		// define keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		rpth := ""

		for txid := range bgw.Tflink {
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()
			if srck == "tflink" {
				keys, vals = bgw.TflinkParseConf()
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", txid, ".tsv")
			} else {
				keys, vals = bgw.TftgParseConf()
				rpth = fmt.Sprintf("%s%s%s%s%s%s", datdir, "static/", srck, "/", txid, ".f2g")
			}
			// log.Println("Rdf4bgw.tfac2gene(): processing", rpth)
			err := parse.Tab2struct(rpth, keys, vals, &d4b)
			if err != nil { // normal
				// log.Printf("%s%s", "tfac2gene:parse.Tab2struct: ", err)
				continue // sic!
			}
			/// d4b is now loaded with data

			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json") // read BGW map json
			err = xmap.Unmarshal(rpthx)
			util.CheckE(err)

			d4b.Src = srck
			d4b.Taxid = txid
			err = export.Tfac2gene(&d4b, &xmap, bgwdir)
			if err != nil {
				//panic(err)
				log.Println(err)
				continue
			}
			for _, pdck := range pdcks {
				cnts.Add(pdck, srck)
				cnts[pdck][srck] = d4b.Cnts[pdck][srck]
			}
		} // txid
	} // srck
	return cnts, nil
} // tfac2gene

func prot2prot(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	log.Println("\n\tprot2prot for:", "all")
	cnts := make(util.Set2D)
	var pdcks = []string{
		"tlp2tlp",
	}
	srck := "intact" // for now
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tprot2prot for:", txid)
		subdir := "intact/"
		ext := ".mi25"
		rpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext) // read IntAct tsv
		subdir = "prot2prot/"
		ext = ".nt"
		// wpth := fmt.Sprintf("%s%sintact-%s%s", bgwdir, subdir, txid, ext) // write ppi nt
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		var xmap bgw.Xmap
		xmap.New()
		err := xmap.Unmarshal(rpthx)
		util.CheckE(err)
		/////////////////////////////////////////////////////////////////////////////
		//duos, err := parse.MiTab(rpth, upac2bgw)
		var d4b bgw.Dat4bridge
		d4b.New()
		keys, vals := bgw.IntactParseConf()
		err = parse.Tab2struct(rpth, keys, vals, &d4b)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.prot2prot():%s: %s", err, txid)
			log.Println(msg)
		} // NoData
		/////////////////////////////////////////////////////////////////////////////
		d4b.Src = srck
		d4b.Taxid = txid
		err = export.Prot2prot(&d4b, &xmap, bgwdir)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.prot2prot():%s: %s", err, txid)
			log.Println(msg)
			continue
		}
		for _, pdck := range pdcks {
			cnts.Add(pdck, srck)
			cnts[pdck][srck] = d4b.Cnts[pdck][srck]
		}
	} // taxid
	return cnts, nil
} // prot2prot()

func gene2phen(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	// TODO interface similar to tfac2gene etc.
	log.Println("\n\tgene2phen for:", "all") // is not printed TODO
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
		var xmap bgw.Xmap
		xmap.New()
		err := xmap.Unmarshal(rpthx)
		util.CheckE(err)
		gsym2bgw := xmap.Lblg
		/////////////////////////////////////////////////////////////////////////////
		//duos, err := parse.UpVar(rpth, gsym2bgw) // the second arg is optional
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
	// TODO interface similar to tfac2gene etc.
	log.Println("\n\tprot2go for:", "all") // is not printed TODO
	nln := 0
	for _, txid := range txn2prm.Keys() {
		if txid == "9031" {
			continue
		}
		log.Println("\n\tprot2go for:", txid)
		subdir := "goa/"
		rpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, fx) // read Goa data
		subdir = "xmap/"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json") // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		var xmap bgw.Xmap
		xmap.New()
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
			log.Println("rdf4bgw.prot2go():", msg)
			continue // TODO test this
			// panic(errors.New(msg))
		}
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.prot2go():%s: %s", err, txid)
			log.Println(msg)
		}
		/////////////////////////////////////////////////////////////////////////////
		wpth := ""
		ext := ".nt"
		subdir = "prot2bp/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntp, err := export.Prot2go(bps, upac2bgw, wpth)
		util.CheckE(err)
		subdir = "prot2cc/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntc, err := export.Prot2go(ccs, upac2bgw, wpth)
		util.CheckE(err)
		subdir = "prot2mf/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntf, err := export.Prot2go(mfs, upac2bgw, wpth)
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

// ///////////////////////////////////////////////////////////////////////////
func ortho(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	// TODO interface similar to tfac2gene etc.
	log.Println("\n\tortho for:", "all")
	nln := 0
	for _, txidL := range txn2prm.Keys() {
		for _, txidR := range txn2prm.Keys() {
			if txidL >= txidR {
				continue
			} // skipping symmetrical and digonal
			txids := [2]string{txidL, txidR}
			duos, err := parse.OrthoDuo(datdir, txidL, txidR, txn2prm)
			if err != nil {
				msg := fmt.Sprintf("rdf4bgw.go:main.ortho():%s:%v", err, txids)
				log.Println(msg)
			} // NoData
			/////////////////////////////////////////////////////////////////////////////
			subdir := "ortho/"
			ext := ".nt"
			file := fmt.Sprintf("%s%s%s%s", txidL, "-", txidR, ext)
			wpth := fmt.Sprintf("%s%s%s", bgwdir, subdir, file)
			nts, err := export.Ortho(duos, wpth)
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

// //////////////////////////////////////////////////////////////////////////
func main() {
	// TODO for all functions: add taxa and proteome lists as arguments !
	start := time.Now()
	// pointers:
	aP := flag.Bool("a", false, "export [a]ll")
	eP := flag.Bool("e", false, "export gene and protein [e]ntities")
	iP := flag.Bool("i", false, "export molecular [i]nteractions")
	dP := flag.Bool("d", false, "export [d]isease associations")
	gP := flag.Bool("g", false, "export [g]ene ontology annotations")
	rP := flag.Bool("r", false, "export [r]egulatory associations")
	oP := flag.Bool("o", false, "export [o]rthology relateions")
	// tP := flag.String("t", "./prm_txn.txt", "selected [t]axa")

	var n int
	flag.Parse()
	if !flag.Parsed() {
		msg := fmt.Sprintf("flag.Parse(): Failed")
		panic(errors.New(msg))
	}
	args := flag.Args()
	n = len(args)
	m := 3
	if n < m {
		msg := fmt.Sprintf("Want at least: %d args have: %d", m, n)
		panic(errors.New(msg))
	}
	datdir := args[0] // path to data directory (with a trailing '/')
	bgwdir := args[1] // path to rdf directory (with a trailing '/')
	rpthT := args[2]  // path to a list of selected taxa and proteomes
	log.Println("Started rdf4bgw with args:", args)
	/////////////////////////////////////////////////////////////////////////////
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "_") // txnID -> proteomeID

	if len(txn2prm) == 0 {
		msg := fmt.Sprintf("%s: NoProteomes", rpthT)
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
		tfac2gene(datdir, bgwdir, txn2prm)
		log.Println("Done with tfac2gene in", time.Since(mystart))
		mystart = time.Now()
		rgr2trg(datdir, bgwdir, txn2prm)
		log.Println("Done with rgr2trg in", time.Since(mystart))
		mystart = time.Now()
		reg2pway(datdir, bgwdir, txn2prm)
		log.Println("Done with reg2pway in", time.Since(mystart))
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
