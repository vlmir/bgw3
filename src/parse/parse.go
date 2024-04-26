package parse

import (
	"bufio"
	"encoding/csv"
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/util"
	"io"
	"log"
	"os"
	"strings"
)

func primaryKey(cells []string, keys []bgw.Column) string {
	// TODO tests
	// cells: values after splitting on the primary delimiter
	pkbits := make([]string, 0) // sic, for appending values
	//var pkbits = []string{} // equvalent
	dlm := keys[0].Dlm2
	var pk string
	for _, bc := range keys {
		// bc: bgw.Column, bc.Key (the last field) may contain empty strings
		kcell := strings.TrimSpace(cells[bc.Ind1]) // any value may occur in the datb
		if kcell == "" {
			return ""
		}
		// special case, does occur
		if kcell == "-" {
			// skipping a component of the primary key, occurs in some IntAct files
			continue // otherwise tests fail for some sources other than  IntAct- why?? TODO
		}
		kbits := strings.Split(kcell, bc.Dlm1) // sub-fields
		if len(kbits) < bc.Ind2+1 {
			return ""
		}
		if bc.Key != "" && kbits[0] != bc.Key {
			return "" // filtering by the type of data specified by bc.Key
		}
		kval := strings.TrimSpace(kbits[bc.Ind2])
		if kval == "" {
			return ""
		}
		pkbits = append(pkbits, kval)
	} // the order of componenets in the key is deterministic, no variation from call to call
	if keys[0].Ind3 == -1 {
		// the clause is needed for Tab2set3D()
		// len(keys) == 1 for SigMaps  
		if len(pkbits) == 1 {
			pkbits = append(pkbits, pkbits[0])
		}
	}
	pk = strings.Join(pkbits, dlm)
	return pk // OK to return string
} // primaryKey

// addSubFields() splits a string (pval) on the delimiter
// specified in a bgw.Column (v) and adds the values to util.Set3D (subfields)
func addSubFields(pval, pk string, v bgw.Column, subfields util.Set3D) {
	// Dlm2: secondary delimiter
	svals := strings.Split(pval, v.Dlm2)
	for ind, sval := range svals {
		sval = strings.TrimSpace(svals[ind])
		if len(sval) == 0 {
			continue
		}
		if v.Ind2 >= 0 && ind != v.Ind2 {
			continue // only one subfield is used
		} // othgerwise all subfields are used
		// special case db:id
		// skipping dbs other than that specified in v.Key
		if v.Ind3 == -1 {
			src := strings.TrimSpace(svals[0]) // Attn: src != sval !!
			if src != v.Key {
				continue
			}
		}
		subfields.Add(pk, v.Key, sval) // the args are non-empty strings
	} // sval
	// no need to return any values
}

// Sig2up() parses mapping files provided by Signor and returns a map structure
func Sig2up(sigmap util.Set3D, pths []string) error {
	// not used anymore, kept just as an example of using encoding/csv
	for _, pth := range pths {
		// open the file
		csvfile, err := os.Open(pth)
		if err != nil {
			msg := fmt.Sprintf("parse.Sig2up(): os.Open(%s): %s", pth, err)
			return errors.New(msg)
		}

		// Parse the file
		r := csv.NewReader(csvfile)
		r.Comma = ';'
		r.Comment = '#'

		// Iterate through the records
		for {
			// Read each record from csv
			rec, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Fatal(err)
			}
			//rec[2] = strings.Replace(rec[2], ", ", "", -1) ??
			if len(rec[2]) == 0 {
				continue // no entities
			}
			// list of individual protein entities
			list := strings.Split(rec[2], ", ")
			// Signor ids and names
			sigmap.Add(rec[0], "lbl", rec[1])
			for _, item := range list {
				id := strings.TrimSpace(item)
				if len(id) == 0 {
					continue
				}
				sigmap.Add(rec[0], "ids", id)
			}
		}
	}
	return nil
}

// TODO generalize, improve error handling
func Sigmap(datdir string) (util.Set3D, error) {
	// replacement for Sig2up accomodating tab-delimited files
	// used only in rdf4bgw
	sigmap := make(util.Set3D)
	rdir := fmt.Sprintf("%s%s%s", datdir, "signor", "/")
	smpths := []string{
		rdir + "complexes.tsv",
		rdir + "families.tsv",
	}
	keys, vals := bgw.SigMapParseConf()
	// len(keys) == 1 len(vals) == 1
	for _, rpth := range smpths {
		out, err := Tab2set3D(rpth, keys, vals)
		if err != nil {
			return out, err
		}
		for k, v := range out {
			sigmap[k] = v
		}
	}
	return sigmap, nil
}

func Idmap(rpth string, idmkeys map[string]string, i1, i2, i3 int) (util.Set3D, error) {
	out := make(util.Set3D)
	fh, err := os.Open(rpth)
	if err != nil {
		msg := fmt.Sprintf("parse.Idmap(%s, %v, %d, %d, %d): os.Open(%s): %s", rpth, idmkeys, i1, i2, i3, rpth, err)
		return out, errors.New(msg)
	}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		if len(cells) != 3 {
			continue
		}
		_, ok := idmkeys[cells[1]] // filtering
		if !ok {
			continue
		} // filtering by idmkeys
		// ALL proteomes 2022-12-14: no '"' anymore, 28 occurences of "''"
		key1 := strings.Replace(cells[i1], "\"", "`", -1) // was present in 44689
		key2 := strings.Replace(cells[i2], "\"", "`", -1) // was present in 44689
		key3 := strings.Replace(cells[i3], "\"", "`", -1) // was present in 44689
		out.Add(key1, key2, key3)
	}
	if len(out) == 0 {
		msg := fmt.Sprintf("parse.Idmap(%s, %v, %d, %d, %d): EmptyMap", rpth, idmkeys, i1, i2, i3)
		return out, errors.New(msg)
	}

	return out, nil
}

// vals.Ind1 - column index (split on an arbitrary string)
// vals.Dlm1 - primary separator for multiple values, all values used
// vals.Ind2 - the index of the value to use after splitting on Dlm2, if < 0 all values
// the 3 fields below used only by addSubFields()
// vals.Dlm2 - secondary separator for multiple values
// vals.Key - the string to be used as the seondary key in the output maps
// vals.Ind3 - integer used for controling the output
// if Ind3 == -1 vals.Key is used for filteering the values by addSubFields()
func Tab2struct(rpth string, keys, vals []bgw.Column, p *bgw.Dat4bridge, dlm string) (err error) {
	// NO empty values added to Dat4bridge !
	// the value in vals.Ind1 is split by vals.Dlm1, ALL subfilds are processed by addSubFields()
	// addSubFields() is NOT concerned with Ind1 and Dlm1
	// if Ind2 < 0 all subfields are used, otherwise only the one in vals.Ind2
	// if Ind3 == -1 subfields are filtered by vals.Key
	d4b := *p
	maxind := 0
	// finding the  maximal index in keys+vlals
	for _, one := range keys {
		n := one.Ind1
		if n > maxind {
			maxind = n
		}
		n = one.Ind2
		if n > maxind {
			maxind = n
		}
	}
	for _, one := range vals {
		n := one.Ind1
		if n > maxind {
			maxind = n
		}
		n = one.Ind2
		if n > maxind {
			maxind = n
		}
	}

	duos := d4b.Duos // to be filled with data
	fh, err := os.Open(rpth)
	if err != nil {
		return err
	}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	ln := 0 // current line number
	for scanner.Scan() {
		// by default scans for '\n'
		line := scanner.Text()
		ln++
		if len(line) == 0 {
			continue
		}
		if string(line[0]) == "#" {
			continue
		}
		cells := strings.Split(line, dlm) // fields
		if len(cells) < maxind+1 {
			msg := fmt.Sprintf("%s:%d: TooFetFields: want %d have %d", rpth, ln, maxind+1, len(cells))
			log.Println(msg)
			continue
		}
		/// primary key
		pk := primaryKey(cells, keys)
		if pk == "" {
			continue
		}
		/// values
		for _, v := range vals {
			// v: bgw.Column; specifies fields to be extracted
			cell := strings.TrimSpace(cells[v.Ind1])
			if (cell == "") || (cell == "-") {
				continue
			}
			// Dlm1: primary delimiter for multiple values
			// there is at least one value, may contain secondary delimiters
			for _, pval := range strings.Split(cell, v.Dlm1) {
				// subfields
				addSubFields(pval, pk, v, duos) // pval trimmed and not empty
			}
		} // one field
	} // one line
	if len(duos) == 0 {
		msg := fmt.Sprintf("parse.Tab2struct():%s: NoData", rpth)
		return errors.New(msg)
	}
	d4b.Duos = duos
	*p = d4b
	return nil
} // Tab2struct

// Tab2set3D is a generic parser for tab-delimeted files with sub-fields (up to 2 levels)
// Tab2set3D converts tabular data into a map keyed on an arbitrary combination of fields
// Arguments:
// 1. path to the input file (tab separated)
// 2. data structure specifying primary key
// 3. data structure specifying values
// Returns:
// 1. map primary key -> secondary key -< values -> counts
// 2. error
func Tab2set3D(rpth string, keys, vals []bgw.Column) (out util.Set3D, err error) {
	// TODO generalize for accepting any primary delimiter as in Tab2struct()
	// used by export.Gene(), export.Prot(), parse.Sigmap()
	// NO empty values added to Set3D
	maxind := 0
	for _, one := range keys {
		n := one.Ind1
		if n > maxind {
			maxind = n
		}
		n = one.Ind2
		if n > maxind {
			maxind = n
		}
	}
	for _, one := range vals {
		n := one.Ind1
		if n > maxind {
			maxind = n
		}
		n = one.Ind2
		if n > maxind {
			maxind = n
		}
	}

	out = make(util.Set3D)
	fh, err := os.Open(rpth)
	if err != nil {
		return out, err
	}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	ln := 0 // current line number
	for scanner.Scan() {
		// by default scans for '\n'
		line := scanner.Text()
		ln++
		if len(line) == 0 {
			continue
		}
		if string(line[0]) == "#" {
			continue
		}
		cells := strings.Split(line, "\t") // fields
		if cells[0] == "Entry" {           // header line TODO generalize
			continue
		}
		if len(cells) < maxind+1 {
			msg := fmt.Sprintf("parse.Tab2set3D(%s, _, _): line: %d TooFewFields", rpth, ln)
			fmt.Printf("%s\n", msg)
			continue
		}
		/// primary key
		pk := primaryKey(cells, keys)
		if pk == "" {
			continue
		}
		/// values
		for _, v := range vals {
			// v: bgw.Column; specify fields to be extracted
			cell := strings.TrimSpace(cells[v.Ind1]) // value of the selected column
			if (cell == "") || (cell == "-") {
				continue
			}
			// Dlm1: primary delimiter for multiple values
			// there is at least one value, may contain secondary delimiters
			for _, pval := range strings.Split(cell, v.Dlm1) {
				addSubFields(pval, pk, v, out) // pval trimmed and not empty
			}
		} // one field
	} // one line
	if len(out) == 0 {
		msg := fmt.Sprintf("parse.Tab2set3D():%s: NoData", rpth)
		return out, errors.New(msg)
	}
	return out, nil
} // Tab2set3D()

func UpVar(rpth string) (duos util.Set3D, err error) {
	// NO empty values in Set3D
	/*
	 AARS	  P49588	 VAR_063527  p.Arg329His	Disease	   rs267606621 Charcot-Marie-Tooth disease 2N(CMT2N) [MIM:613287]
	  0  => 0 gene name
	 10 => 10 up acc // always present and single
	 21 => 21
	 33 => 33 aa change
	 48 => 48 type
	 62 => 57 '-' possible
	 77 => 72 description '-' possible
	*/
	duos = make(util.Set3D)
	fhR, err := os.Open(rpth)
	if err != nil {
		return duos, err
	}
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	nsL := "hgncsymb"
	nsR := "omim"
	for scanner.Scan() { // by default scans for '\n'
		line := scanner.Text()
		if len(line) < 74 {
			continue
		} // skipping lines with '-' in the last field SIC!
		if strings.TrimSpace(line[48:56]) != "Disease" {
			//	continue // Disease => LP/P; seems superfluous anyway
		}
		upca := strings.TrimSpace(line[10:20])
		symG := strings.TrimSpace(line[0:9])
		oriL := symG
		idL := fmt.Sprintf("%s%s%s", nsL, "!", oriL)
		dfn := strings.TrimSpace(line[72:])
		// only one definition per line
		if dfn == "" {
			continue
		} // seems redundant
		bits := strings.Split(dfn, "[MIM:") // sic, returns 1 value
		// Attn: multiple ':' may occur !!!
		// 20200531: only for "OXCT1     P55809":
		// Succinyl-CoA:3-oxoacid CoA transferase deficiency (SCOTD) [MIM:245050]
		if len(bits) < 2 {
			continue
		} // no MIM ID
		if len(bits) > 2 {
			msg := fmt.Sprintf("parse.UpVar():%s:%s: MultiMIMs, skipped", upca, symG)
			fmt.Printf("%s\n", msg)
			continue
		} // normally should never happen
		// removing the trailing ']':
		oriR := strings.TrimSuffix(bits[1], "]")
		idR := fmt.Sprintf("%s%s%s", nsR, "!", oriR)
		if idR == "" {
			continue
		}
		pairid := fmt.Sprintf("%s%s%s", idL, "--", idR) // idL and idR trimmed
		duos.Add(pairid, "dfn", strings.TrimSpace(bits[0]))
		duos.Add(pairid, "upca", upca) // upca trimmed
	}
	if len(duos) == 0 {
		msg := fmt.Sprintf("parse.UpVar(%s): NoData", rpth)
		return duos, errors.New(msg)
	}
	return duos, nil
} // UpVar

/*
fields:
4:values:
"" # very many !
NOT
NOT|colocalizes_with
NOT|contributes_to
colocalizes_with
contributes_to
*/
func Gaf(rpth string) (bp, cc, mf util.Set3D, err error) {
	// TODO replace with Tab2struct()
	ourppys := map[string]string{
		"C": "gp2cc",
		"F": "gp2mf",
		"P": "gp2bp",
	}
	srckL := "uniprot"
	srckR := "obo"
	bp = make(util.Set3D)
	cc = make(util.Set3D)
	mf = make(util.Set3D)
	fhR, err := os.Open(rpth)
	if err != nil {
		msg := fmt.Sprintf("parse.Gaf(): os.Open(%s): %s", rpth, err)
		return bp, cc, mf, errors.New(msg)
	}
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		upac := ""
		if cells[0] == "UniProtKB" {
			upac = cells[1]
		} else {
			continue
		}
		goid := strings.Replace(cells[4], ":", "_", 1)
		if cells[3] == "NOT" {
			continue
		}
		ppy, ok := ourppys[cells[8]] // aspect (C|F|P)
		if !ok {
			continue
		}
		idL := fmt.Sprintf("%s%s%s", srckL, "!", upac)
		idR := fmt.Sprintf("%s%s%s", srckR, "!", goid)
		pairid := fmt.Sprintf("%s%s%s", idL, "--", idR)
		refs := make([]string, 0)
		for _, ref := range strings.Split(cells[5], "|") {
			bits := strings.Split(ref, ":")
			if bits[0] == "PMID" {
				refid := fmt.Sprintf("%s%s%s", "pubmed", "!", bits[1])
				refs = append(refs, refid)
			}
		}
		goc := cells[6] // GO evidence code
		switch ppy {
		case "gp2bp":
			bp.Add(pairid, "ppy", ppy)
			bp.Add(pairid, "goc", goc)
			for _, ref := range refs {
				bp.Add(pairid, "ref", ref)
			}
		case "gp2mf":
			mf.Add(pairid, "ppy", ppy)
			mf.Add(pairid, "goc", goc)
			for _, ref := range refs {
				mf.Add(pairid, "ref", ref)
			}
		case "gp2cc":
			cc.Add(pairid, "ppy", ppy)
			cc.Add(pairid, "goc", goc)
			for _, ref := range refs {
				cc.Add(pairid, "ref", ref)
			}
		}
	}
	if len(bp)+len(cc)+len(mf) == 0 {
		msg := fmt.Sprintf("parse.Gaf():%s: NoData", rpth)
		err := errors.New(msg)
		return bp, cc, mf, err
	}
	return bp, cc, mf, nil
} // Gaf()

/*
!   name                         required? cardinality             GAF column #
1   DB                           yes       1                         1
2   DB_Object_ID                 yes       1                         2 / 17
3   Qualifier                    no        0 or greater              4
4   GO ID                        yes       1                         5
5   DB:Reference                 yes       1 or greater              6
6   ECO evidence code            yes       1                         7 + 6 (GO evidence code + reference)
7   With/From                    no        0 or greater              8
8   Interacting taxon ID         no        0 or 1                   13
9   Date                         yes       1                        14
10   Assigned_by                  yes       1                        15
11   Annotation Extension         no        0 or greater             16
12   Annotation Properties        no        0 or 1                   n/a
3:values:
NOT|colocalizes_with
NOT|contributes_to
NOT|enables
NOT|involved_in
NOT|part_of
acts_upstream_of
acts_upstream_of_negative_effect
acts_upstream_of_or_within
acts_upstream_of_or_within_negative_effect
acts_upstream_of_or_within_positive_effect
acts_upstream_of_positive_effect
colocalizes_with
contributes_to
enables
involved_in
is_active_in
part_of

	gpakeys := map[int][]string{
		0: []string{"db", "|"},
		1: []string{"acc", "|"},
		2: []string{"qlrs", "|"}, // Attn. entries like 'NOT|enables'
		3: []string{"goid", "|"},
		4: []string{"refs", "|"}, // e.g. PMID:11463391, GO_REF:0000024 not all have pubmed refs!
		5: []string{"eco", "|"},
	}
*/
func Gpa(rpth string) (duos util.Set3D) {
	// TODO replace with Tab2struct()
	ourppys := map[string]string{
		"part_of":     "gp2cc",
		"enables":     "gp2mf",
		"involved_in": "gp2bp",
	}
	srckL := "uniprot"
	srckR := "obo"
	duos = make(util.Set3D)
	fhR, err := os.Open(rpth)
	if err != nil {
		msg := fmt.Sprintf("parse.Gpa(): os.Open(%s): %s", rpth, err)
		panic(errors.New(msg))
	}
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		upac := ""
		if cells[0] == "UniProtKB" {
			upac = cells[1]
		} else {
			continue
		}
		goid := strings.Replace(cells[3], ":", "_", 1)
		qlrs := strings.Split(cells[2], "|")
		//if qlrs[0] == "NOT" { continue }
		if len(qlrs) != 1 {
			continue
		}
		ppy, ok := ourppys[qlrs[0]]
		if !ok {
			continue
		} // filtering py rels
		idL := fmt.Sprintf("%s%s%s", srckL, "!", upac)
		idR := fmt.Sprintf("%s%s%s", srckR, "!", goid)
		pairid := fmt.Sprintf("%s%s%s", idL, "--", idR)
		eco := strings.Replace(cells[5], ":", "_", 1)
		refs := make([]string, 0)
		for _, ref := range strings.Split(cells[4], "|") {
			bits := strings.Split(ref, ":")
			if bits[0] == "PMID" {
				refid := fmt.Sprintf("%s%s%s", "pubmed", "!", bits[1])
				refs = append(refs, refid)
			}
		}
		duos.Add(pairid, "ppy", ppy)
		duos.Add(pairid, "eco", eco)
		for _, ref := range refs {
			duos.Add(pairid, "ref", ref)
		}
	}
	return duos
} // Gpa()

// orthosolo() extracts orthology relations from UniProt idmappings for one taxon (sic!)
// orthosolo() is used by OrthoDuo()
// retrns a map for 1 taxon, primary key sourse dabase label
func orthosolo(datdir, txid string, txn2prm util.Set2D) (util.Set3D, error) {
	idmkeys := bgw.Orthokeys
	solos := make(util.Set3D)
	subdir := "idmapping/"
	ext := ".idmapping"
	prmids := txn2prm[txid].Keys() // normally 1, yet multiple may occur
	if len(prmids) > 1 {
		msg := fmt.Sprintf("parse.orthosolo(_, %s, _): MultipleProteomes: %s", txid, prmids)
		fmt.Printf("%s\n", msg)
	} // never occors
	for _, prmid := range prmids {
		prmid := fmt.Sprintf("%s%s%s", prmid, "_", txid)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, prmid, ext) // read
		dat, err := Idmap(pth, idmkeys, 1, 2, 0)
		if err != nil {
			msg := fmt.Sprintf("parse.orthosolo(_, %s, _):  %s", txid, err)
			return solos, errors.New(msg)
		}
		for idmk, all := range dat {
			for xid, one := range all {
				for upac, _ := range one {
					solos.Add(idmk, xid, upac)
				}
			} // xid: external cluster ID
		} // idmk: e.g. KO, OrthoDB
	}
	count := 0
	for idmk, _ := range solos {
		count += len(solos[idmk])
	}
	if count == 0 {
		msg := fmt.Sprintf("parse.orthosolo(): NoDataForTaxon: %s", txid)
		return solos, errors.New(msg)
	} // no orthology in any of the sources
	return solos, nil
} // orthosolo()

// OrthoDuo() extracts orthology relations from UniProt idmappings for a pair of taxa
// returns a map, keys:
// primary: relation label
// secondary: source database label
// tertiary: relation identifier in the source database
// func OrthoDuo(datdir, txidL, txidR string, txn2prm util.Set2D, idmkeys map[string]string) (util.Set3D, error) {
func OrthoDuo(datdir, txidL, txidR string, txn2prm util.Set2D) (util.Set3D, error) {
	idmkeys := bgw.Orthokeys
	duos := make(util.Set3D)
	outL, err := orthosolo(datdir, txidL, txn2prm)
	if err != nil {
		return duos, err
	} // NoDate for one taxaon
	outR, err := orthosolo(datdir, txidR, txn2prm)
	if err != nil {
		return duos, err
	} // NoDate for one taxaon
	for idmk, _ := range idmkeys {
		allL, ok := outL[idmk]
		if !ok {
			continue
		}
		allR, ok := outR[idmk]
		if !ok {
			continue
		} // now data from 'idmk' is present in both taxa
		for id, oneL := range allL {
			oneR, ok := allR[id]
			if !ok {
				continue
			} // only shared clusters
			for upLac, _ := range oneL {
				for upRac, _ := range oneR {
					if upLac >= upRac {
						continue
					} // skipping diagonal and symmetrical
					duoid := fmt.Sprintf("uniprot!%s--uniprot!%s", upLac, upRac)
					duos.Add(duoid, idmk, id)
				}
			}
		} // id: external cluster ID
	}
	if len(duos) == 0 {
		msg := fmt.Sprintf("parse.OrthoDuo(): NoDataForTaxa: %s--%s", txidL, txidR)
		return duos, errors.New(msg)
	} // 20200531: none Note: no filtering by BGW
	return duos, nil
} // OrthoDuo()
