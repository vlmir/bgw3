// package main
package rdf4bgw

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
	"os/exec"
	"path/filepath"
	"time"
)

func rdfpipe(strs ...string) error {
	// TODO tests
	// TODO capture stdin and write as in:
	// https://www.sohamkamani.com/golang/exec-shell-command/
	var cmd *exec.Cmd
	// Attn: all flags separately!
	rpth := strs[0]
	ifmt := strs[1]
	ofmt := strs[2]
	wpth := strs[3]
	cmd = exec.Command("rdfpipe", "-i", ifmt, "-o", ofmt, rpth)
	out, err := cmd.Output()
	if err != nil {
		return err
	}
	wfh, err := os.Create(wpth)
	if err != nil {
		return err
	}
	defer wfh.Close()
	wfh.Write([]byte(out))
	return nil
}

func Geneprot(datdir, bgwdir string, txn2prm util.Set2D) (err error) {
	subdir := "gene/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
	subdir = "prot/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
	subdir = "xmap/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tGeneprot for:", txid)
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
		err := export.Gene(rpthi, wpthg, &xmap)
		if err != nil {
			return err
		}
		if err := util.Gzip(wpthg); err != nil {
			return err
		}

		err = export.Prot(rpthu, rpthi, wpthp, &xmap)
		if err != nil {
			return err
		}
		if err := util.Gzip(wpthp); err != nil {
			return err
		}

		wfhX, err := os.Create(wpthx)
		if err != nil {
			return err
		}
		j, err := json.MarshalIndent(&xmap, "", " ")
		if err != nil {
			return err
		}
		wfhX.Write(j)
	} // txid
	return nil
} // Geneprot()

func Reg2pway(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	subdir := "reg2pway/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
	log.Println("\n\trdf4bgw.Reg2pway for:", "all")
	cnts := make(util.Set2D)
	for srck, _ := range rdf.Uris4rgrtrg {
		// define keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		ext := ""
		rpth := ""
		if srck == "signor" {
			// keys, vals = bgw.SigPwaysParseConf()
			keys = bgw.SigPwaysConf.Keys
			vals = bgw.SigPwaysConf.Vals
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
			log.Println("rdf4bgw.Reg2pway(): processing", rpth)
			err := parse.Tab2struct(rpth, keys, vals, &d4b, "\t")
			if err != nil {
				msg := fmt.Sprintf("%s%s", "rdf4bgw.Reg2pway:parse.Tab2struct: ", err)
				fmt.Printf("%s\n", msg)
				continue // sic! TODO double check
			}
			// d4b is now loaded with data
			wpths := map[string]string{
				"reg2ptrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "p", srck, txid),
				"reg2ntrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "n", srck, txid),
				"reg2utrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "u", srck, txid),
			}
			d4b.Out = wpths
			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			ext := ".json"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
			err = xmap.Unmarshal(rpthx)
			if err != nil {
				return cnts, err
			}

			if srck == "signor" {
				// generating map signor-id -> entitity-ids
				sigmap, _ := parse.Sigmap(datdir)
				xmap.Signor = sigmap
			}

			d4b.Src = srck
			d4b.Taxid = txid
			err = export.SigPways(&d4b, &xmap, bgwdir)
			if err != nil {
				return cnts, err
			}
			for pdck, wpth := range wpths {
				if err := util.Gzip(wpth); err != nil {
					return cnts, err
				}
				cnts.Add(pdck, srck)
				cnts[pdck][srck] = d4b.Cnts[pdck][srck]
			}
		} // txid
	} // srck
	return cnts, nil
} // Reg2pway

func Reg2targ(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	subdir := "reg2targ/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
	log.Println("\n\trdf4bgw.Reg2targ for:", "all")
	cnts := make(util.Set2D)
	for srck, _ := range rdf.Uris4rgrtrg {
		// define keys and vals for parsing
		ext := ""
		rpth := ""
		if srck == "signor" {
			ext = ".mi28"
		}

		for txid := range txn2prm {
			if txid != "9606" {
				continue
			} // for now
			wpths := map[string]string{
				"reg2ptrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "p", srck, txid),
				"reg2ntrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "n", srck, txid),
				"reg2utrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "u", srck, txid),
			}

			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()
			if srck == "signor" {
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", txid, ext)
			}
			log.Println("rdf4bgw.Reg2targ(): processing", rpth)
			err := parse.Tab2struct(rpth, bgw.SignorConf.Keys, bgw.SignorConf.Vals, &d4b, "\t")
			if err != nil {
				msg := fmt.Sprintf("%s%s", "rdf4bgw.Reg2targ:parse.Tab2struct: ", err)
				fmt.Printf("%s\n", msg)
				return cnts, err
			}
			// d4b is now loaded with data
			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			ext := ".json"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
			err = xmap.Unmarshal(rpthx)
			if err != nil {
				return cnts, err
			}

			if srck == "signor" {
				// generating map signor-id -> entitity-ids
				sigmap, _ := parse.Sigmap(datdir)
				xmap.Signor = sigmap
			}

			d4b.Src = srck
			d4b.Taxid = txid
			d4b.Out = wpths
			err = export.Reg2targ(&d4b, &xmap, bgwdir)
			if err != nil {
				return cnts, err
			}
			for pdck, wpth := range wpths {
				err := util.Gzip(wpth)
				if err != nil {
					return cnts, err
				}
				cnts.Add(pdck, srck)
				cnts[pdck][srck] = d4b.Cnts[pdck][srck]
			}
		} // txid
	} // srck
	return cnts, nil
} // Reg2targ

func Tfac2gene(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	subdir := "tfac2gene/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
	log.Println("\n\trdf4bgw.Tfac2gene for:", "all")
	cnts := make(util.Set2D)
	for srck, _ := range rdf.Uris4tftg {
		// define keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		rpth := ""
		dlm := "\t"

		// looping over all taxa present in Tflink
		for txid := range bgw.Tflink {
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()
			if srck == "tflink" {
				keys = bgw.TflinkConf.Keys
				vals = bgw.TflinkConf.Vals
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", txid, ".tsv")
			}
			if srck == "coltri" {
				keys = bgw.ColtriConf.Keys
				vals = bgw.ColtriConf.Vals
				dlm = "," // re-defining
				if txid != "9606" {
					continue
				}
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", txid, ".csv")
			}
			// log.Println("rdf4bgw.Tfac2gene(): processing", rpth)
			err := parse.Tab2struct(rpth, keys, vals, &d4b, dlm)
			if err != nil { // normal
				msg := fmt.Sprintf("%s%s", "Tfac2gene:parse.Tab2struct: ", err)
				fmt.Printf("%s\n", msg)
				// return cnts, err // tests fail !!
				continue // next taxon
			}
			/// d4b is now loaded with data

			wpths := map[string]string{
				"reg2ptrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "p", srck, txid),
				"reg2ntrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "n", srck, txid),
				"reg2utrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "u", srck, txid),
			}
			d4b.Out = wpths
			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json") // read BGW map json
			err = xmap.Unmarshal(rpthx)
			if err != nil {
				return cnts, err
			}

			d4b.Src = srck
			d4b.Taxid = txid
			err = export.Tfac2gene(&d4b, &xmap, bgwdir)
			if err != nil {
				return cnts, err
			}
			for pdck, wpth := range wpths {
				err := util.Gzip(wpth)
				if err != nil {
					return cnts, err
				}
				cnts.Add(pdck, srck)
				cnts[pdck][srck] = d4b.Cnts[pdck][srck]
			}
		} // txid
	} // srck
	return cnts, nil
} // Tfac2gene

func Prot2prot(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
	subdir := "prot2prot/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
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
		//		wpth := fmt.Sprintf("%s%sintact-%s%s", bgwdir, subdir, txid, ext) // path to ppi nt
		wpth := fmt.Sprintf("%s%s/%s-%s.nt", bgwdir, subdir, srck, txid) // path to .nt file
		subdir = "xmap/"
		ext = ".json"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		var xmap bgw.Xmap
		xmap.New()
		err := xmap.Unmarshal(rpthx)
		if err != nil {
			return cnts, err
		}
		/////////////////////////////////////////////////////////////////////////////
		//duos, err := parse.MiTab(rpth, upac2bgw)
		var d4b bgw.Dat4bridge
		d4b.New()
		err = parse.Tab2struct(rpth, bgw.IntactConf.Keys, bgw.IntactConf.Vals, &d4b, "\t")
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.Prot2prot():%s: %s", err, txid)
			log.Println(msg)
		} // NoData
		/////////////////////////////////////////////////////////////////////////////
		d4b.Src = srck
		d4b.Taxid = txid
		err = export.Prot2prot(&d4b, &xmap, bgwdir)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go:main.Prot2prot():%s: %s", err, txid)
			log.Println(msg)
			// continue
			return cnts, err
		}
		for _, pdck := range pdcks {
			cnts.Add(pdck, srck)
			cnts[pdck][srck] = d4b.Cnts[pdck][srck]
		}
		if err = util.Gzip(wpth); err != nil {
			log.Println("util.Gzip(): Failed to gzip:", wpth)
			return cnts, err
		}
	} // taxid
	return cnts, nil
} // Prot2prot()

func Gene2phen(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	subdir := "gene2phen/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("rdf4bgw.Gene2phen: Failed to make: %s %s", subdir, err)
		return 0, errors.New(msg)
	}
	// TODO interface similar to Tfac2gene etc.
	cnt := 0
	for txid := range txn2prm {
		if txid != "9606" {
			continue
		}
		log.Println("\n\tGene2phen for:", txid)
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
		if err != nil {
			return 0, err
		}
		gsym2bgw := xmap.Lblg
		/////////////////////////////////////////////////////////////////////////////
		duos, err := parse.UpVar(rpth) // returns error if len(duos) == 0
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.Gene2phen(): parse.UpVar(%s) %s: ", rpth, err)
			return 0, errors.New(msg)
		}
		/////////////////////////////////////////////////////////////////////////////
		nrel, err := export.Gene2phen(duos, gsym2bgw, wpth)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.Gene2phen(): export.Gene2phen(_, _, %s): %s ", wpth, err)
			return 0, errors.New(msg)
		}
		cnt += nrel
		if err = util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("rdf4bgw.Gene2phen(): util.Gzip(%s): %s", wpth, err)
			return 0, errors.New(msg)
		}
	} // txid
	return cnt, nil
} // Gene2phen

func Prot2go(datdir, bgwdir string, txn2prm util.Set2D, fx string) (int, error) {
	subdir := "prot2bp/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		return 0, err
	}
	subdir = "prot2cc/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		return 0, err
	}
	subdir = "prot2mf/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		return 0, err
	}
	// TODO interface similar to Tfac2gene etc.
	log.Println("\n\tProt2go for:", "all") // is not printed TODO
	nln := 0
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tProt2go for:", txid)
		subdir := "goa/"
		rpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, fx) // read Goa data
		subdir = "xmap/"
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json") // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		var xmap bgw.Xmap
		xmap.New()
		err := xmap.Unmarshal(rpthx)
		if err != nil {
			return 0, err
		}
		upac2bgw := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		bps := make(util.Set3D)
		ccs := make(util.Set3D)
		mfs := make(util.Set3D)
		if fx == ".gpa" {
			//out := parse.Gpa(rpth) // TODO
		}
		if fx == ".gaf" {
			bps, ccs, mfs, err = parse.Gaf(rpth)
			if err != nil {
				msg := fmt.Sprintf("rdf4bgw.go:main.Prot2go():%s: %s", err, txid)
				log.Println(msg)
				return 0, err
			}
		}
		/////////////////////////////////////////////////////////////////////////////
		wpth := ""
		ext := ".nt"
		subdir = "prot2bp/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntp, err := export.Prot2go(bps, upac2bgw, wpth)
		if err = util.Gzip(wpth); err != nil {
			log.Println("util.Gzip(): Failed to gzip:", wpth)
			return 0, err
		}
		subdir = "prot2cc/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntc, err := export.Prot2go(ccs, upac2bgw, wpth)
		if err = util.Gzip(wpth); err != nil {
			log.Println("util.Gzip(): Failed to gzip:", wpth)
			return 0, err
		}
		subdir = "prot2mf/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntf, err := export.Prot2go(mfs, upac2bgw, wpth)
		if err = util.Gzip(wpth); err != nil {
			log.Println("util.Gzip(): Failed to gzip:", wpth)
			return 0, err
		}

		nts := ntp + ntc + ntf
		if nts == 0 {
			msg := fmt.Sprintf("%s: NoTriples", txid)
			return 0, errors.New(msg)
		}
		nln += nts
	} // txid
	return nln, nil
} // Prot2go()

// ///////////////////////////////////////////////////////////////////////////
func Ortho(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	subdir := "ortho/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Println(err)
	}
	// TODO interface similar to Tfac2gene etc.
	log.Println("\n\tortho for:", "all")
	cnt := 0
	for _, txidL := range txn2prm.Keys() {
		for _, txidR := range txn2prm.Keys() {
			if txidL >= txidR {
				continue
			} // skipping symmetrical and digonal
			duos, err := parse.OrthoDuo(datdir, txidL, txidR, txn2prm)
			if err != nil {
				// OrthoDuo returns error if no orthology data for one of the taxa (occurs)
				// OR no orthologues found for a pair of taxa (never occured so far)
				msg := fmt.Sprintf("rdf4bgw.Ortho(): parse.OrthoDuo(_, %s, %s, _): %s: ", txidL, txidR, err)
				fmt.Printf("%s\n", msg)
				continue // works as intended
			} // NoData
			/////////////////////////////////////////////////////////////////////////////
			subdir := "ortho/"
			ext := ".nt"
			file := fmt.Sprintf("%s%s%s%s", txidL, "-", txidR, ext)
			wpth := fmt.Sprintf("%s%s%s", bgwdir, subdir, file)
			nrel, err := export.Ortho(duos, wpth)
			if err != nil {
				msg := fmt.Sprintf("rdf4bgw.Ortho(): export.Ortho(_, %s): %s", wpth, err)
				return 0, errors.New(msg) // assuming duos is not empty
			}
			if err = util.Gzip(wpth); err != nil {
				msg := fmt.Sprintf("rdf4bgw.Ortho(): util.Gzip(%s): %s", wpth, err)
				return 0, errors.New(msg)
			}
			cnt += nrel
		} // txidR
	} // txidL
	return cnt, nil
} // Orhto

func Onto(datdir, bgwdir string) error {
	var ontos = map[string]string{
		/*
			"biolink-model": ".ttl", // not available anymore
		*/
		"omim":      ".ttl",
		"bfo":       ".owl",
		"go-basic":  ".owl",
		"mi":        ".owl",
		"ncbitaxon": ".owl",
		"ro":        ".owl",
		"sio":       ".owl",
	}
	subdir := "onto/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		log.Fatal(err)
	}
	for onto, ext := range ontos {
		rpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, onto, ext)
		wpth := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, onto, ".nt")
		ifmt := ""
		if ext == ".ttl" {
			ifmt = "turtle"
		} else {
			ifmt = "application/rdf+xml"
		}
		if err := rdfpipe(rpth, ifmt, "nt", wpth); err != nil {
			log.Println("rdf4bgw.rdfpipe(): Failed to convert: ", rpth)
			return err
		}
		if err := util.Gzip(wpth); err != nil {
			log.Println("util.Gzip(): Failed to gzip:", wpth)
			return err
		}
	}
	return nil
}

// //////////////////////////////////////////////////////////////////////////
func main() {
	// TODO for all functions: add taxa and proteome lists as arguments !
	// pointers:
	aP := flag.Bool("a", false, "export [a]ll")
	eP := flag.Bool("e", false, "export gene and protein [e]ntities")
	iP := flag.Bool("i", false, "export molecular [i]nteractions")
	dP := flag.Bool("d", false, "export [d]isease associations")
	gP := flag.Bool("g", false, "export [g]ene ontology annotations")
	rP := flag.Bool("r", false, "export [r]egulatory associations")
	oP := flag.Bool("o", false, "export [o]rthology relations")
	// tP := flag.String("t", "./prm_txn.txt", "selected [t]axa")
	start := time.Now()

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
		start := time.Now()
		Geneprot(datdir, bgwdir, txn2prm) // MUST be run before the others !!!
		log.Println("Done with Geneprot in", time.Since(start))
	}
	if *aP || *dP {
		start := time.Now()
		Gene2phen(datdir, bgwdir, txn2prm)
		log.Println("Done with Gene2phen in", time.Since(start))
	}
	if *aP || *gP {
		start := time.Now()
		ext := ".gaf"
		_, err = Prot2go(datdir, bgwdir, txn2prm, ext)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with Prot2go in", time.Since(start))
		}
	}
	if *aP || *iP {
		start := time.Now()
		Prot2prot(datdir, bgwdir, txn2prm)
		log.Println("Done with Prot2prot in", time.Since(start))
	}
	if *aP || *rP {
		start := time.Now()
		Tfac2gene(datdir, bgwdir, txn2prm)
		log.Println("Done with rdf4bgw.Tfac2gene in", time.Since(start))
		start = time.Now()
		Reg2targ(datdir, bgwdir, txn2prm)
		log.Println("Done with rdf4bgw.Reg2targ in", time.Since(start))
		start = time.Now()
		Reg2pway(datdir, bgwdir, txn2prm)
		log.Println("Done with rdf4bgw.Reg2pway in", time.Since(start))
	} // rP
	if *aP || *oP {
		start := time.Now()
		_, err := Ortho(datdir, bgwdir, txn2prm)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with rdf4bgw.Ortho in", time.Since(start))
		}
	}
	log.Println("Done with rdf4bgw in", time.Since(start))
} // main()
