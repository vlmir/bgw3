package rdf4bgw

import (
	"encoding/json"
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/export"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"os"
	"os/exec"
	"path/filepath"
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
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	wfh, err := os.Create(wpth)
	if err != nil {
		msg := fmt.Sprintf("%s: os.Create: %s", util.FN(0), err)
		return errors.New(msg)
	}
	defer wfh.Close()
	wfh.Write([]byte(out))
	return nil
}

func Geneprot(datdir, bgwdir string, txn2prm util.Set2D) (err error) {
	subdir := "gene/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return errors.New(msg)
	}
	subdir = "prot/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return errors.New(msg)
	}
	subdir = "xmap/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return errors.New(msg)
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

//func Reg2pway(datdir, bgwdir string, txn2prm util.Set2D) (util.Set2D, error) {
func Reg2pway(datdir, bgwdir string, taxa map[string][]string) (util.Set2D, error) {
	subdir := "reg2pway/"
	cnts := make(util.Set2D)
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err) // no need for FN(1)
		return cnts, errors.New(msg)
	}
	for srck, txids := range taxa {
		// keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		rpth := ""
		dlm := ""
		ext := ""

		if srck == "signor" {
			keys = bgw.SigPwaysConf.Keys
			vals = bgw.SigPwaysConf.Vals
			ext = ".tsv"
			dlm = "\t"
		}

		for _, txid := range txids {
			log.Println("\n\tReg2pway for:", txid)
			wpths := map[string]string{
				"reg2ptrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "p", srck, txid),
				"reg2ntrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "n", srck, txid),
				"reg2utrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "u", srck, txid),
			}
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()

			/// parsing
			if srck == "signor" {
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", "pathways", ext)
			}
			err := parse.Tab2struct(rpth, keys, vals, &d4b, dlm)
			if err != nil {
				// either failed to open rpth or no interactions extracted
				msg := fmt.Sprintf("%s: parse.Tab2struct: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			// d4b is now loaded with data

			/// exporting
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

func Reg2targ(datdir, bgwdir string, taxa map[string][]string) (util.Set2D, error) {
	subdir := "reg2targ/"
	cnts := make(util.Set2D)
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return cnts, errors.New(msg)
	}
	for srck, txids := range taxa {
		// keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		rpth := ""
		dlm := ""
		ext := ""

		if srck == "signor" {
			keys = bgw.SignorConf.Keys
			vals = bgw.SignorConf.Vals
			dlm = "\t"
			ext = ".mi28"
		}

		for _, txid := range txids {
			log.Println("\n\tReg2targ for:", txid)
			wpths := map[string]string{
				"reg2ptrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "p", srck, txid),
				"reg2ntrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "n", srck, txid),
				"reg2utrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "u", srck, txid),
			}
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()

			/// parsing
			rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", txid, ext)
			err := parse.Tab2struct(rpth, keys, vals, &d4b, dlm)
			if err != nil {
				// either failed to open rpth or no interactions extracted
				msg := fmt.Sprintf("%s: parse.Tab2struct: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			// d4b is now loaded with data

			/// exporting
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

func Tfac2gene(datdir, bgwdir string, taxa map[string][]string) (util.Set2D, error) {
	subdir := "tfac2gene/"
	cnts := make(util.Set2D)
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return cnts, errors.New(msg)
	}
	for srck, txids := range taxa {
		// keys and vals for parsing
		var vals []bgw.Column
		var keys []bgw.Column
		rpth := ""
		dlm := ""
		ext := ""

		if srck == "tflink" {
			keys = bgw.TflinkConf.Keys
			vals = bgw.TflinkConf.Vals
			dlm = "\t"
			ext = ".tsv"
		} else if srck == "atregnet" {
			keys = bgw.AtregnetConf.Keys
			vals = bgw.AtregnetConf.Vals
			dlm = "\t"
			ext = ".tsv"
		} else if srck == "coltri" {
			keys = bgw.ColtriConf.Keys
			vals = bgw.ColtriConf.Vals
			dlm = ","
			ext = ".csv"
		}

		for _, txid := range txids {
			log.Println("\n\tTfac2gene for:", txid)
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()
			/// parsing
			rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", txid, ext)
			err := parse.Tab2struct(rpth, keys, vals, &d4b, dlm)
			if err != nil {
				// either failed to open rpth or no interactions extracted
				msg := fmt.Sprintf("%s: parse.Tab2struct: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			/// d4b is now loaded with data

			/// exporting
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
				// currently if the header does not have requred number of triples
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
	cnts := make(util.Set2D)
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return cnts, errors.New(msg)
	}
	var pdcks = []string{
		"tlp2tlp",
	}
	srck := "intact" // for now
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tProt2prot for:", txid)
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
			msg := fmt.Sprintf("%s", err) // err includes: Prot2prot Tab2struct rpth
			return cnts, errors.New(msg)  // TODO check
		} // NoData
		/////////////////////////////////////////////////////////////////////////////
		d4b.Src = srck
		d4b.Taxid = txid
		err = export.Prot2prot(&d4b, &xmap, bgwdir)
		if err != nil {
			msg := fmt.Sprintf("rdf4bgw.go: Prot2prot():%s: %s", err, txid)
			return cnts, errors.New(msg) // TODO check
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
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
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
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return 0, errors.New(msg)
	}
	subdir = "prot2cc/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return 0, errors.New(msg)
	}
	subdir = "prot2mf/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return 0, errors.New(msg)
	}
	// TODO interface similar to Tfac2gene etc.
	nln := 0
	for _, txid := range txn2prm.Keys() {
		log.Println("\n\tProt2go for:", txid) // is not printed TODO
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
				msg := fmt.Sprintf("rdf4bgw.go: Prot2go():%s: %s", err, txid)
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
	// Note: currently no mapping to OrthoDB for 284812, 44689
	subdir := "ortho/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return 0, errors.New(msg)
	}
	// TODO interface similar to Tfac2gene etc.
	log.Println("\n\tOrtho")
	cnt := 0 // used in testing
	for _, txidL := range txn2prm.Keys() {
		for _, txidR := range txn2prm.Keys() {
			if txidL >= txidR {
				continue
			} // skipping symmetrical and digonal
			duos, err := parse.OrthoDuo(datdir, txidL, txidR, txn2prm)
			if err != nil {
				// OrthoDuo returns error if no orthology data for one of the taxa (occurs)
				// err passed from parse.Idmap()
				// OR no orthologues found for a pair of taxa (never occured so far)
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				fmt.Printf("%s, skipping\n", msg)
				continue // works as intended
			} // NoData
			/////////////////////////////////////////////////////////////////////////////
			subdir := "ortho/"
			ext := ".nt"
			file := fmt.Sprintf("%s%s%s%s", txidL, "-", txidR, ext)
			wpth := fmt.Sprintf("%s%s%s", bgwdir, subdir, file)
			nrel, err := export.Ortho(duos, wpth)
			// export.Ortho() retuns ALWAYS nil !! TODO
			// NB: by this time ALL duos are non-empty
			if err != nil {
				msg := fmt.Sprintf("rdf4bgw.Ortho(): export.Ortho(_, %s): %s", wpth, err)
				return 0, errors.New(msg)
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
	subdir := "onto/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: os.Mkdir: %s", util.FN(0), err)
		return errors.New(msg)
		//log.Fatal(err)
	}
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
			msg := fmt.Sprintf("%s: rdf4bgw.rdfpipe(%s): %s: ", util.FN(0), rpth, err)
			return errors.New(msg)
		}
		if err := util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("%s: util.Gzip(%s): %s: ", util.FN(0), wpth, err)
			return errors.New(msg)
		}
	}
	return nil
}
