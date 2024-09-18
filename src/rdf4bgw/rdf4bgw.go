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
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	defer wfh.Close()
	wfh.Write([]byte(out))
	return nil
}

func Geneprot(datdir, bgwdir string, txn2prm util.Set2D) (err error) {
	subdir := "gene/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	subdir = "prot/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	subdir = "xmap/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
		if err := export.Gene(rpthi, wpthg, &xmap); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}
		if err := util.Gzip(wpthg); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}

		if err := export.Prot(rpthu, rpthi, wpthp, &xmap); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}
		if err := util.Gzip(wpthp); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}

		wfhX, err := os.Create(wpthx) // sic!
		if err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}
		j, err := json.MarshalIndent(&xmap, "", " ") // sic!
		if err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}
		wfhX.Write(j)
	} // txid
	return nil
} // Geneprot()

func Reg2pway(datdir, bgwdir string, taxa map[string][]string) (util.Set2D, error) {
	subdir := "reg2targ/"
	cnts := make(util.Set2D)
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
		// generating map signor-id -> entitity-ids
		sigmap, err := parse.Sigmap(datdir)
		if err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return cnts, errors.New(msg)
		}

		for _, txid := range txids {
			log.Println("\n\tReg2pway for:", txid)
			wpths := map[string]string{
				"reg2ptrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "ppw", srck, txid),
				"reg2ntrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "npw", srck, txid),
				"reg2utrg": fmt.Sprintf("%s%s%s-%s-%s.nt", bgwdir, subdir, "upw", srck, txid),
			}
			var d4b bgw.Dat4bridge // one source, one taxon
			d4b.New()

			/// parsing
			if srck == "signor" {
				rpth = fmt.Sprintf("%s%s%s%s%s", datdir, srck, "/", "pathways", ext)
			}
			if err := parse.Tab2struct(rpth, keys, vals, &d4b, dlm); err != nil {
				// either failed to open rpth or no interactions extracted
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			// d4b is now loaded with data

			/// exporting
			d4b.Out = wpths
			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			ext := ".json"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext)
			if err := xmap.Unmarshal(rpthx); err != nil {
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}

			if srck == "signor" {
				xmap.Signor = sigmap
			}

			d4b.Src = srck
			d4b.Taxid = txid
			if err := export.SigPways(&d4b, &xmap, bgwdir); err != nil {
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			for pdck, wpth := range wpths {
				if err := util.Gzip(wpth); err != nil {
					msg := fmt.Sprintf("%s: %s", util.FN(0), err)
					return cnts, errors.New(msg)
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
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
		// generating map signor-id -> entitity-ids
		sigmap, err := parse.Sigmap(datdir)
		if err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return cnts, errors.New(msg)
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
			if err := parse.Tab2struct(rpth, keys, vals, &d4b, dlm); err != nil {
				// either failed to open rpth or no interactions extracted
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			// d4b is now loaded with data

			/// exporting
			var xmap bgw.Xmap
			xmap.New()
			subdir := "xmap/"
			ext := ".json"
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext)
			if err := xmap.Unmarshal(rpthx); err != nil {
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}

			if srck == "signor" {
				xmap.Signor = sigmap
			}

			d4b.Src = srck
			d4b.Taxid = txid
			d4b.Out = wpths
			if err := export.Reg2targ(&d4b, &xmap, bgwdir); err != nil {
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			for pdck, wpth := range wpths {
				if err := util.Gzip(wpth); err != nil {
					msg := fmt.Sprintf("%s: %s", util.FN(0), err)
					return cnts, errors.New(msg)
				}
				cnts.Add(pdck, srck)
				cnts[pdck][srck] = d4b.Cnts[pdck][srck]
			}
		} // txid
	} // srck
	return cnts, nil
} // Reg2targ

func Tfac2gene(datdir, bgwdir string, taxa map[string][]string) (util.Set2D, error) {
	// Tflink contains only generic interactions (no +/-) for 10116 6239 7227 7955 ->
	// 8 empty files saved
	subdir := "tfac2gene/"
	cnts := make(util.Set2D)
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
			if err := parse.Tab2struct(rpth, keys, vals, &d4b, dlm); err != nil {
				// either failed to open rpth or no interactions extracted
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
			rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json")
			if err := xmap.Unmarshal(rpthx); err != nil {
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}

			d4b.Src = srck
			d4b.Taxid = txid
			if err := export.Tfac2gene(&d4b, &xmap, bgwdir); err != nil {
				// if the header does not have requred number of triples
				// OR no data triples
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return cnts, errors.New(msg)
			}
			for pdck, wpth := range wpths {
				if err := util.Gzip(wpth); err != nil {
					msg := fmt.Sprintf("%s: %s", util.FN(0), err)
					return cnts, errors.New(msg)
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
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
		wpth := fmt.Sprintf("%s%s/%s-%s.nt", bgwdir, subdir, srck, txid) // path to .nt file
		subdir = "xmap/"
		ext = ".json"
		var xmap bgw.Xmap
		xmap.New()
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext)
		if err := xmap.Unmarshal(rpthx); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return cnts, errors.New(msg)
		}
		var d4b bgw.Dat4bridge
		d4b.New()
		if err := parse.Tab2struct(rpth, bgw.IntactConf.Keys, bgw.IntactConf.Vals, &d4b, "\t"); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return cnts, errors.New(msg) // TODO check
		} // NoData
		d4b.Src = srck
		d4b.Taxid = txid
		if err := export.Prot2prot(&d4b, &xmap, bgwdir); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return cnts, errors.New(msg) // TODO check
		}
		for _, pdck := range pdcks {
			cnts.Add(pdck, srck)
			cnts[pdck][srck] = d4b.Cnts[pdck][srck]
		}
		if err := util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return cnts, errors.New(msg) // TODO check
		}
	} // taxid
	return cnts, nil
} // Prot2prot()

func Gene2phen(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	subdir := "gene2phen/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
		var xmap bgw.Xmap
		xmap.New()
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext)
		if err := xmap.Unmarshal(rpthx); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}
		duos, err := parse.UpVar(rpth) // sic!
		if err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}
		gsym2bgw := xmap.Lblg
		nrel, err := export.Gene2phen(duos, gsym2bgw, wpth) // sic!
		if err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}
		cnt += nrel
		if err = util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}
	} // txid
	return cnt, nil // TODO should be return nil
} // Gene2phen

func Prot2go(datdir, bgwdir string, txn2prm util.Set2D, fx string) (int, error) {
	subdir := "prot2bp/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return 0, errors.New(msg)
	}
	subdir = "prot2cc/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return 0, errors.New(msg)
	}
	subdir = "prot2mf/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
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
		var xmap bgw.Xmap
		xmap.New()
		rpthx := fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ".json")
		if err := xmap.Unmarshal(rpthx); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}
		upac2bgw := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		bps := make(util.Set3D)
		ccs := make(util.Set3D)
		mfs := make(util.Set3D)
		var err error
		if fx == ".gpa" {
			//out := parse.Gpa(rpth) // TODO
		}
		if fx == ".gaf" {
			if bps, ccs, mfs, err = parse.Gaf(rpth); err != nil {
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return 0, errors.New(msg)
			}
		}
		/////////////////////////////////////////////////////////////////////////////
		wpth := ""
		ext := ".nt"
		subdir = "prot2bp/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntp, err := export.Prot2go(bps, upac2bgw, wpth)
		if err = util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}
		subdir = "prot2cc/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntc, err := export.Prot2go(ccs, upac2bgw, wpth)
		if err = util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}
		subdir = "prot2mf/"
		wpth = fmt.Sprintf("%s%s%s%s", bgwdir, subdir, txid, ext) // write goa nt
		ntf, err := export.Prot2go(mfs, upac2bgw, wpth)
		if err = util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return 0, errors.New(msg)
		}

		nts := ntp + ntc + ntf
		if nts == 0 {
			msg := fmt.Sprintf("%s: NoTriplesForTaxon: %s", util.FN(0), txid)
			return 0, errors.New(msg)
		}
		nln += nts
	} // txid
	return nln, nil // nln is used for tests
} // Prot2go()

// ///////////////////////////////////////////////////////////////////////////
func Ortho(datdir, bgwdir string, txn2prm util.Set2D) (int, error) {
	// currently OrthoDB is the only source
	// no mapping to OrthoDB for 284812, 44689
	// NO empty files saved anymore !!
	subdir := "ortho/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return 0, errors.New(msg)
	}
	// TODO interface similar to Tfac2gene etc.
	cnt := 0 // used in testing
	for _, txidL := range txn2prm.Keys() {
		outL, err := parse.OrthoSolo(datdir, txidL, txn2prm) // Set3D
		if err != nil {
			continue // skipping taxa with no orthology data
		}
		for _, txidR := range txn2prm.Keys() {
			outR, err := parse.OrthoSolo(datdir, txidR, txn2prm) // Set3D
			if err != nil {
				continue // skipping taxa wothno orthology data
			}
			// now there is data for both taxa
			if txidL >= txidR {
				continue
			} // skipping symmetrical and digonal
			duos := make(util.Set3D)
			for srck, _ := range bgw.Orthokeys {
				clusters := util.Shared(outL[srck].Keys(), outR[srck].Keys())
				if len(clusters) == 0 {
					continue // sic - some pairs of taxa may have no orthologues, sofar none
				}
				for _, xid := range clusters {
					for _, upLac := range outL[srck][xid].Keys() {
						for _, upRac := range outR[srck][xid].Keys() {
							duoid := fmt.Sprintf("uniprot!%s--uniprot!%s", upLac, upRac)
							duos.Add(duoid, srck, xid)
						}
					}
				} // xid
			} // srck
			/////////////////////////////////////////////////////////////////////////////
			subdir := "ortho/"
			ext := ".nt"
			txpair := fmt.Sprintf("%s-%s", txidL, txidR)
			file := fmt.Sprintf("%s%s", txpair, ext)
			wpth := fmt.Sprintf("%s%s%s", bgwdir, subdir, file)
			log.Println("\n\tOrtho for:", txpair)
			nrel, err := export.Ortho(duos, wpth)
			// export.Ortho() retuns ALWAYS nil !! TODO
			// NB: by this time ALL duos are non-empty
			if err != nil {
				msg := fmt.Sprintf("%s: %s", util.FN(0), err)
				return 0, errors.New(msg)
			}
			if err = util.Gzip(wpth); err != nil {
				msg := fmt.Sprintf("rdf4bgw.Ortho(): util.Gzip(%s): %s", wpth, err)
				return 0, errors.New(msg)
			}
			cnt += nrel
		} // txidR
	} // txidL
	if cnt == 0 {
		msg := fmt.Sprintf("%s: NoOrthogs", util.FN(0))
		return cnt, errors.New(msg)
	}
	return cnt, nil
} // Orhto

func Onto(datdir, bgwdir string) error {
	subdir := "onto/"
	if err := os.MkdirAll(filepath.Join(bgwdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
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
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}
		if err := util.Gzip(wpth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}
	}
	return nil
}
