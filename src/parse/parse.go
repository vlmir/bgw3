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

func checkID(id string, filter util.Set3D, counter util.Set1D) (out int) {
	out = 0
	_, ok := filter[id]
	if !ok {
		out++
		_, ok := counter[id]
		if !ok {
			out++
			counter[id]++
		}
	}
	return out
}

// addSubFields() splits a value in a tab-delimited file on the delimiter
// specified in bgw.Column and adds the values to util.Set3D
func addSubFields(pval, pk string, v bgw.Column, out util.Set3D) {
	// Dlm2: secondary delimiter
	svals := strings.Split(pval, v.Dlm2)
	ind2 := v.Ind2
	sk := v.Key // secondary key explicitely specified
	for i, sval := range svals {
		if ind2 >= 0 && i != ind2 {
			continue // only one subfield is used
		}
		sval = strings.TrimSpace(svals[i])
		if len(sval) == 0 {
			continue
		}
		// db:id
		// skipping dbs other that specified in v.Key
		if v.Ind3 == -1 {
			if strings.TrimSpace(svals[0]) != sk {
				continue
			}
		}
		out.Add(pk, sk, sval)
	}
}

// Sig2up() parses mapping files provided by Signor and returns a map structure
func Sig2up(sigmap util.Set3D, pths []string) error {
	for _, pth := range pths {

		// Open the file
		csvfile, err := os.Open(pth)
		if err != nil {
			log.Fatalln("Couldn't open the csv file", err)
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

// used only in orthosolo()
func Idmap(rpth string, srcs map[string]string, i1, i2, i3 int) (util.Set3D, error) {
	out := make(util.Set3D)
	fh, err := os.Open(rpth)
	//util.CheckE(err)
	if err != nil {
		msg := fmt.Sprintf("parse.Idmap(%s, %v, %d, %d, %d):", rpth, srcs, i1, i2, i3)
		log.Println(msg)
		return out, err
	}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		if len(cells) != 3 {
			continue
		}
		_, ok := srcs[cells[1]] // filtering
		if !ok {
			continue
		} // filtering by srcs
		// ALL proteomes 2022-12-14: no '"' anymore, 28 occurences of "''"
		key1 := strings.Replace(cells[i1], "\"", "`", -1) // was present in 44689
		key2 := strings.Replace(cells[i2], "\"", "`", -1) // was present in 44689
		key3 := strings.Replace(cells[i3], "\"", "`", -1) // was present in 44689
		out.Add(key1, key2, key3)
	}
	return out, nil
}

// val.Ind1 - column index
// val.Dlm1 - primary separator of multiple values
// val.Dlm2 - secondary separator of multiple values
// val.Ind2 - the index of the value to use after splitting on Dlm2, if < 0 all values
// val.Key - the string to be used as the seondary key in the output maps
// val.Ind3 - integer used for conntroleing the output
// if Ind3 == -1 val.Key is used for filteering the values
func Tab2struct(rpth string, keys, vals []bgw.Column, p *bgw.Dat4bridge) (err error) {
	// the value in val.Ind1 is split by val.Dlm1, ALL subfilds are processed by addSubFields()
	// if Ind2 < 0 all subfields are used, otherwise only the one in val.Ind2
	// if Ind3 == -1 subfields are filtered by val.Key
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
	// util.CheckE(err)
	if err != nil {
		return err
	}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	ln := 0 // current line number
	for scanner.Scan() {
		// by default scans for '\n'
		// TODO generalize, consider using csv.NewReader as in Sig2up()
		line := scanner.Text()
		ln++
		if len(line) == 0 {
			continue
		}
		if string(line[0]) == "#" {
			continue
		}
		cells := strings.Split(line, "\t") // fields
		if len(cells) < maxind+1 {
			msg := fmt.Sprintf("%s:%d: TooFetFields: want %d have %d", rpth, ln, maxind+1, len(cells))
			//panic(errors.New(msg))
			log.Println(msg)
			continue
		}
		/// primary key
		var pk string // the primary key to be used in the output map
		for i, k := range keys {
			// i: index, k: bgw.Column
			items := strings.Split(cells[k.Ind1], k.Dlm1) // sub-fields
			spk := strings.TrimSpace(items[k.Ind2])
			util.CheckStrings(spk)
			if i == 0 {
				pk = spk
			} else {
				// joining on Dlm2
				pk = fmt.Sprintf("%s%s%s", pk, k.Dlm2, spk)
			}
		} // the order of componenets in the key is deterministic, no variation from call to call
		util.CheckStrings(pk)
		/// values
		for _, v := range vals {
			// v: bgw.Column; specify fields to be extracted
			cell := strings.TrimSpace(cells[v.Ind1])
			if (cell == "") || (cell == "-") {
				continue
			}
			// Dlm1: primary delimiter for multiple values
			// there is at least one value, may contain secondary delimiters
			for _, pval := range strings.Split(cell, v.Dlm1) {
				// subfields
				addSubFields(pval, pk, v, duos)
			}
		} // one field
	} // one line
	if len(duos) == 0 {
		msg := fmt.Sprintf("parse.Tab2set3D():%s: NoData", rpth)
		return errors.New(msg)
	}
	d4b.Duos = duos
	*p = d4b
	return nil
} // Tab2struct()

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
	// util.CheckE(err)
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
		if len(cells) < maxind+1 {
			msg := fmt.Sprintf("%s:%d: TooFewFields", rpth, ln)
			panic(errors.New(msg))
		}
		/// primary key
		var pk string // the primary key to be used in the output map
		for i, k := range keys {
			// i: index, k: bgw.Column
			items := strings.Split(cells[k.Ind1], k.Dlm1) // sub-fields
			spk := strings.TrimSpace(items[k.Ind2])
			util.CheckStrings(spk)
			if i == 0 {
				pk = spk
			} else {
				// joining on Dlm2
				pk = fmt.Sprintf("%s%s%s", pk, k.Dlm2, spk)
			}
		} // the order of componenets in the key is deterministic, no variation from call to call
		util.CheckStrings(pk)
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
				addSubFields(pval, pk, v, out)
			}
		} // one field
	} // one line
	if len(out) == 0 {
		msg := fmt.Sprintf("parse.Tab2set3D():%s: NoData", rpth)
		return out, errors.New(msg)
	}
	return out, nil
} // Tab2set3D()

// parses the output of the system 'ls' utility
// any string can be used as a basename-extension separator
/*
func Basenames(rpth, dlm string) (set util.Set2D, err error) {
	//set := make(util.Set2D)
	keyind := 0
	valind := 1
	fh, err := os.Open(rpth)
	if err != nil {return set, err}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	for scanner.Scan() {// by default scans for '\n'
		cells := strings.Split(scanner.Text(), "/")
		file := cells[len(cells)-1]
		bits := strings.Split(file, dlm) // normally "."
		if len(bits) < 2 { continue }
		basename := bits[keyind]
		set.Add(basename, bits[valind])
	}
	return set, nil
}
*/

// UpIdMap filters a file of the form db1id\tdb2label\tdb2id by values of db2lablel
func UpIdMap(rpth string, idmkeys map[string]string) (upac2xrf util.Set3D, err error) {
	// used one only in rdf4bgw.go pass upac2xrf to parse.UpTab():
	// needed for retrieving all iso-form accessions in export.GeneProt()!
	upac2xrf = make(util.Set3D)
	rfh, err := os.Open(rpth)
	util.CheckE(err)
	defer rfh.Close()
	scanner := bufio.NewScanner(rfh)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		if len(cells) != 3 {
			continue
		} // just a sanity check
		key, ok := idmkeys[cells[1]]
		if !ok {
			continue
		} // filtering by idmkeys
		upac := cells[0]
		// ALL proteomes 2022-12-14: no '"' anymore, 28 occurences of "''"
		xrf := strings.Replace(cells[2], "\"", "`", -1) // was present in 44689
		upac2xrf.Add(upac, key, xrf)
		// mapping canonical accesssion and isoforms
		// TODO see how to get rid of this
		bits := strings.Split(upac, "-")
		if len(bits) == 2 {
			upca := bits[0]
			upac2xrf.Add(upca, "upac", upac)
		} // only for iso-forms
	}
	return upac2xrf, nil
}

/*
9606:95744 entries no proteome
Proteome examples other than empty (73928 entries):
UP000005640: Chromosome X
UP000005640: Chromosome Y
UP000005640: Mitochondrion
UP000005640: Unplaced // ro be filtered out
Multiple chromosems:
	mironov arch2 ~/g/p/d/uniprot> grep ""UP000005640: Chromosome"" 9606.tsv | cut -f 8 *.tsv | grep "",""
UP000005640: Chromosome 14, UP000005640: Chromosome 19, UP000005640: Chromosome 2
UP000005640: Chromosome 1, UP000005640: Chromosome 5
UP000005640: Chromosome 1, UP000005640: Chromosome 12, UP000005640: Chromosome 6
UP000005640: Chromosome 1, UP000005640: Chromosome 17
*/
// Entry   Entry name      Organism        Organism ID     Protein names   Proteomes       PubMed ID       Annotation
func UpTab(rpth string, upac2xrf util.Set3D, txn2prm util.Set2D) (out bgw.Dat4rdf, err error) {
	// Attn: the input file MUST hane no blank lines ??
	allPs := make(util.Set3D)
	allTs := make(util.Set3D)
	allGs := make(util.Set3D)
	fhR, err := os.Open(rpth)
	util.CheckE(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	notInRefProt := make(util.Set1D)
	noGeneName := make(util.Set1D)
	noChrom := make(util.Set1D)
	multiProtChrom := make(util.Set1D)
	multiGene := make(util.Set1D)
	ln := 0 // line counter
	for scanner.Scan() {
		/// by default scans for '\n'
		// one line per UniProt canonical accession
		cells := strings.Split(scanner.Text(), "\t")
		// the # of fields is currently 11 but only 10 in test files !!
		// the last field (GeneID) is not currently used
		upca := cells[0]
		if upca == "Entry" {
			continue
		} // skipping the header line
		if len(cells) < 10 {
			msg := fmt.Sprintf("parse.UpTab():%s:line:%d: TooFewFields", rpth, ln)
			panic(errors.New(msg))
		} // blanc lines skipped as well
		ln++
		_, ok := upac2xrf[upca] // RefProt ACs
		if !ok {
			notInRefProt[upca]++
			continue
		} // filtering by RefProt
		//v34 upid := cells[1]
		txid := cells[5]
		prmid := txn2prm[txid].Keys()[0] // a single value
		oneP := make(util.Set2D)
		// Pairs Proteome: Chromosome separated by ", "
		// A single Reference Proteome expected
		prm2chr := make(util.Set2D)
		/// Attn: multiple entries do occur due to diff in chromosomes !!
		if len(cells[7]) == 0 {
			panic("PrmChrEmpty") // never happens for Reference Proteoms
		} // at least one pair RefProtID: ChromoID
		prmchrs := strings.Split(cells[7], ", ") // proteome: chromosome pairs
		/// Processing all Proteome:Chromosome values
		for _, prmchr := range prmchrs {
			bits := strings.Split(prmchr, ": ") // proteome: chromosome
			if bits[0] != prmid {
				continue
			} // only the reference proteomes now
			if len(bits) < 2 {
				// happens in 367110 and 284812
				// "Linkage Group I" "Linkage Group II" etc in 367110
				// e.g.UP000001805: Chromosome 1, Linkage Group I, UP000001805: Unassembled WGS sequence
				// "gap-filling sequence" in 284812
				// e.g. UP000002485: Chromosome II, gap-filling sequence
				// at least one pair contains RefProtID
				// one more problem in 36329:
				// UP000001450: Chromosome: api 20200531: 10
				// UP000001450: Chromosome: mit 20200531: 2
				// alongside with
				// Apicoplast A/B 20200531: 30
				// UP000001450: Mitochondrion 20200531: 3
				continue // go to the next pair
			}
			chr := bits[1]
			chrbits := strings.Split(chr, " ")
			chrbit := chrbits[0]
			/*
				trgs := []string{"Chromosome", "Mitochondrion", "Chloroplast", "Apicoplast", "Plasmid"}
				cnt := 0
				for _, trg := range trgs {
					if trg == chrbit {
						cnt++
					}
				}
			*/
			chrid := ""
			switch chrbit {
			case "Chromosome": // "e.g. 'Chromosome 1'"
				shortbit := strings.ToLower(chrbit[0:3])
				if len(bits) == 2 {
					chrnbr := chrbits[1]
					chrid = strings.Join([]string{shortbit, chrnbr}, "-")
				} else if len(bits) == 3 {
					chrlbl := bits[2]
					chrid = fmt.Sprintf("%s-%s", shortbit, chrlbl)
				}
			case "Mitochondrion":
				chrid = strings.ToLower(chrbits[0][0:6])
			case "Chloroplast":
				chrid = strings.ToLower(chrbits[0][0:6])
			case "Apicoplast":
				// "Apicoplast A" "Apicoplat B" in 36329
				chrid = fmt.Sprintf("%s-%s", strings.ToLower(chrbits[0][0:5]), chrbits[1])
			case "Plasmid":
				// "Plasmid 2-micron" in 559292 S.cerevisiae
				lastbit := strings.Replace(chrbits[1], "-", "", -1)
				chrid = fmt.Sprintf("%s-%s", strings.ToLower(chrbits[0][0:4]), lastbit[0:3])
			}
			if chrid != "" {
				prm2chr.Add(prmid, chrid)
			}
		} // all Proteome:Chromosome values for the given UP entry
		if len(prm2chr) == 0 {
			// junk like 'Unassembled WGS sequence', 'Geneome assembly', 'Unplaced' etc
			noChrom[upca]++
			//v34 chrid := fmt.Sprintf("dummy") // unmapped entries allowed as of 2020-08-19
			//v34 prm2chr.Add(prmid, chrid)
		}
		// a SINGLE proteome from now on !!

		// proteome-chromosome pairs
		chrs := prm2chr[prmid].Keys()
		if len(chrs) > 1 {
			// happens in all the 25 taxa; 10230 in 367110
			multiProtChrom[upca]++
			// msg := fmt.Sprintf("parse.UpTab():%s:%s: MultiChrom: %d %v", txid, upca, len(chrs), chrs)
			// fmt.Printf("%s\n", msg)
		}
		// Primary Gene Names separated by "; ":
		// empty strings for many entries
		// numerous entries like: "14-3-3 protein", "20 amino acid ORF", "7 ER" etc
		gnmbag := strings.Split(cells[2], "; ")
		// Groups of " " separated GeneSynonyms, separated by "; "
		// empty string and something like"; ; ;" occur as values
		gsnmbag := strings.Split(cells[3], "; ")
		if len(gsnmbag) != len(gnmbag) {
			msg := fmt.Sprintf("parse.UpTab():%s:%s: Mismatch: %v:%v", txid, upca, gnmbag, gsnmbag)
			panic(errors.New(msg))
		} // never happens
		// removing empty strings, occur in 3702 for example
		var gnms, gsnms []string
		for i, gnm := range gnmbag {
			if gnm == "" {
				continue
			}
			gnms = append(gnms, gnm)
			gsnms = append(gsnms, gsnmbag[i]) // empty strings likely
		}
		if len(gnms) == 0 {
			noGeneName[upca]++
			msg := fmt.Sprintf("parse.UpTab():%s:%s: NoGene", txid, upca)
			fmt.Printf("%s\n", msg)
			//v334 gnms = append(gnms, upid) // allowing entries without Gene Name as of 200-07-19
		}

		if len(gnms) > 1 {
			/*
				for i, _ := range gnms {
					gnms[i] = cells[1]
				} // replacing GeneNames with UniProtID
			*/
			multiGene[upca]++
			// msg := fmt.Sprintf("parse.UpTab():%s:%s: MultiGene: %d %v", txid, upca, len(gnms), gnms)
			// fmt.Printf("%s\n", msg)
		}
		if len(chrs) > 1 && len(gnms) == 1 {
			// msg := fmt.Sprintf("parse.UpTab():%s:%s: OneGeneMultiChrom: %s %v", txid, upca, gnms[0], chrs)
			// fmt.Printf("%s\n", msg)
		}
		if len(chrs) > 1 && len(gnms) > 1 {
			//v34 chrs = []string{"multi"}
			// msg := fmt.Sprintf("parse.UpTab():%s:%s: MultiGeneMultiChrom: %d %v", txid, upca, len(gnms), chrs)
			// fmt.Printf("%s\n", msg)
		}
		for i, gnm := range gnms {
			allGs.Add(gnm, "upca", upca)
			if len(chrs) == 1 {
				allGs.Add(gnm, "chr", chrs[0])
			}

			oneP.Add("gnm", gnm)
			if len(gsnms) != len(gnms) {
				continue
			} // occurs if an artificial GeneName added. e.g. if UP provides no primary gene name
			for _, gsnm := range strings.Split(gsnms[i], " ") {
				if gsnm == "" {
					continue
				} // no synonyms for a particulare GeneName
				allGs.Add(gnm, "gsnm", gsnm)
				oneP.Add("gsnm", gsnm)
			}
		}
		oneP.Add("prmid", prmid)
		//		for chr := range prm2chr[prmid] {
		for _, chr := range chrs {
			chr = strings.Replace(chr, " ", "-", 1)
			oneP.Add("chr", chr)
			allTs.Add(txid, "chr", chr)
			if len(gnms) == 1 {
				allGs.Add(gnms[0], "chr", chr)
			}
		}
		oneP.Add("upid", cells[1])
		oneP.Add("spnm", cells[4])
		oneP.Add("txid", txid)
		oneP.Add("pdfn", strings.Replace(cells[6], "\"", "''", -1))
		//oneP.Add("pubmed", cells[8])
		for _, pmid := range strings.Split(cells[8], "; ") {
			if pmid == "" {
				continue
			}
			oneP.Add("pubmed", pmid)
		}
		oneP.Add("score", strings.Split(cells[9], " ")[0])
		allTs.Add(txid, "spnm", cells[4]) // NEEDED !!
		allTs.Add(txid, "prmid", prmid)
		allPs[upca] = oneP
	} // one UniProt canonical accession
	/////////////////////////////////////////////////////////////////////////////
	txids := allTs.Keys()
	msg := ""
	if len(txids) != 1 {
		msg = fmt.Sprintf("parse.UpTab(): MaldefinedTaxon: %v", txids)
		return out, errors.New(msg)
	}
	if len(allPs) == 0 {
		msg = fmt.Sprintf("parse.UpTab():%s: NoProts", txids)
		return out, errors.New(msg)
	}
	if len(allGs) == 0 {
		msg = fmt.Sprintf("parse.UpTab():%s: NoGenes", txids)
		return out, errors.New(msg)
	}
	msg = fmt.Sprintf("parse.UpTab():%v: ProcessedEntries: %d", txids, ln)
	log.Println(msg)
	// skipped:
	if len(notInRefProt) != 0 {
		msg = fmt.Sprintf("parse.UpTab():%v: Skipped: NotInRefProt: %d", txids, len(notInRefProt))
		log.Println(msg)
	}
	// not skipped:
	if len(noGeneName) != 0 {
		msg = fmt.Sprintf("parse.UpTab():%v: Warning: NoGeneName: %d", txids, len(noGeneName))
		log.Println(msg)
	}
	if len(noChrom) != 0 {
		msg = fmt.Sprintf("parse.UpTab():%v: Warning: NoChrom: %d", txids, len(noChrom))
		log.Println(msg)
	}
	if len(multiProtChrom) != 0 {
		// occurs in all the 25 taxa, normally limited number though goes upto 10230 in 367110
		// 15 canonical accessions in 9606
		msg = fmt.Sprintf("parse.UpTab():%v: Warning: MultiChrom: %d", txids, len(multiProtChrom))
		log.Println(msg)
	}
	if len(multiGene) != 0 {
		msg = fmt.Sprintf("parse.UpTab():%v: Warning: MultiGene: %d", txids, len(multiGene))
		log.Println(msg)
	}
	msg = fmt.Sprintf("parse.UpTab():%s: AcceptedEntries: %d", txids, len(allPs))
	log.Println(msg)
	out.Upac = &upac2xrf // checked in rdf4bgw
	// all checked above
	out.Txns = &allTs
	out.Udat = &allPs
	out.Gnm = &allGs
	return out, nil
} // UpTab()

func UpVar(rpth string, filters ...util.Set3D) (duos util.Set3D, err error) {
	// filters is used only for filtering
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
	nsL := "hgncsymb"
	nsR := "omim"
	fhR, err := os.Open(rpth)
	util.CheckE(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	duos = make(util.Set3D)
	notInBgw := make(util.Set1D)
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
		if len(filters) == 1 {
			if out := checkID(oriL, filters[0], notInBgw); out != 0 {
				continue
			}
		}
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
		pairid := fmt.Sprintf("%s%s%s", idL, "--", idR)
		duos.Add(pairid, "dfn", bits[0])
		duos.Add(pairid, "upca", upca)
	}
	if len(filters) == 1 {
		msg := fmt.Sprintf("parse.UpVar():%s: NotInBgw: %d", rpth, len(notInBgw))
		log.Println(msg)
		//fmt.Printf("%s\n", msg)
	}
	if len(duos) == 0 {
		msg := fmt.Sprintf("parse.UpVar():%s: NoData", rpth)
		return duos, errors.New(msg)
	}
	return duos, nil
} // UpVar()

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
func Gaf(rpth string, filters ...util.Set3D) (bp, cc, mf util.Set3D, err error) {
	ourppys := map[string]string{
		"C": "gp2cc",
		"F": "gp2mf",
		"P": "gp2bp",
	}
	srcL := "uniprot"
	srcR := "obo"
	bp = make(util.Set3D)
	cc = make(util.Set3D)
	mf = make(util.Set3D)
	fhR, err := os.Open(rpth)
	util.CheckE(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	notInBgw := make(util.Set1D)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		upac := ""
		if cells[0] == "UniProtKB" {
			upac = cells[1]
		} else {
			continue
		}
		// filtering by RefProt
		if len(filters) == 1 {
			if out := checkID(upac, filters[0], notInBgw); out != 0 {
				continue
			}
		}
		goid := strings.Replace(cells[4], ":", "_", 1)
		if cells[3] == "NOT" {
			continue
		}
		ppy, ok := ourppys[cells[8]] // aspect (C|F|P)
		if !ok {
			continue
		}
		idL := fmt.Sprintf("%s%s%s", srcL, "!", upac)
		idR := fmt.Sprintf("%s%s%s", srcR, "!", goid)
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
	if len(filters) == 1 {
		msg := fmt.Sprintf("parse.Gaf():%s: NotInBgw: %d", rpth, len(notInBgw))
		log.Println(msg)
		//fmt.Printf("%s\n", msg)
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
func Gpa(rpth string, filters ...util.Set3D) (duos util.Set3D) {
	ourppys := map[string]string{
		"part_of":     "gp2cc",
		"enables":     "gp2mf",
		"involved_in": "gp2bp",
	}
	srcL := "uniprot"
	srcR := "obo"
	duos = make(util.Set3D)
	fhR, err := os.Open(rpth)
	util.CheckE(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	notInBgw := make(util.Set1D)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		upac := ""
		if cells[0] == "UniProtKB" {
			upac = cells[1]
		} else {
			continue
		}
		// filtering by RefProt
		if len(filters) == 1 {
			if out := checkID(upac, filters[0], notInBgw); out != 0 {
				continue
			}
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
		idL := fmt.Sprintf("%s%s%s", srcL, "!", upac)
		idR := fmt.Sprintf("%s%s%s", srcR, "!", goid)
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
	if len(filters) == 1 {
		msg := fmt.Sprintf("parse.Gpa():%s: NotInBgw: %d", rpth, len(notInBgw))
		log.Println(msg)
		//fmt.Printf("%s\n", msg)
	}
	return duos
} // Gpa()

// TODO  use Tab2set3D()
func MiTab(rpth string, filters ...util.Set3D) (duos util.Set3D, err error) {
	tsvkeys := map[int][]string{
		// NB: 0:1 duos are NOT inique !
		//0: []string{"upacA", "|"}, // 9606: all single
		//1: []string{"upacB", "|"}, // 9606: all single
		6:  {"mtds", "|"}, // 9606: all single
		8:  {"pubmedIds", "|"},
		11: {"intactTypes", "|"}, // 9606: all single
		13: {"intactIds", "|"},   // intact:EBI-20559053|imex:IM-26397-1
		14: {"cnfvals", "|"},     // author score:D|intact-miscore:0.37
		//18: []string{"xprlAs", "|"}, // 9606: all single
		//19: []string{"xprlBs", "|"}, // 9606: all single
		//28: []string{"hosts", "|"}, // TODO expunge
	}
	duos = make(util.Set3D)
	fhR, err := os.Open(rpth)
	util.CheckE(err)
	defer fhR.Close()
	srcL := "uniprot"
	srcR := "uniprot"
	scanner := bufio.NewScanner(fhR)
	notInBgw := make(util.Set1D)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		/// processing IDs
		idAs := strings.Split(cells[0], "|")
		idBs := strings.Split(cells[1], "|")
		if len(idAs) == 0 {
			continue
		}
		if len(idBs) == 0 {
			continue
		}
		upacLs := util.X1type(idAs, "uniprotkb", ":")
		if len(upacLs) == 0 {
			continue
		}
		upacA := upacLs[0]
		upacRs := util.X1type(idBs, "uniprotkb", ":")
		upacB := ""
		if len(upacRs) == 0 {
			upacB = upacA // self-interaction
		} else {
			upacB = upacRs[0]
		}
		// filtering by RefProt
		if len(filters) == 1 {
			if out := checkID(upacA, filters[0], notInBgw); out != 0 {
				continue
			}
			if out := checkID(upacB, filters[0], notInBgw); out != 0 {
				continue
			}
		}
		// re-assigning
		cells[0] = upacA
		cells[1] = upacB
		idL := fmt.Sprintf("%s%s%s", srcL, "!", upacA)
		idR := fmt.Sprintf("%s%s%s", srcR, "!", upacB)
		pairid := fmt.Sprintf("%s%s%s", idL, "--", idR)
		/// filling up with the data
		for i, cell := range cells {
			slice, ok := tsvkeys[i]
			if !ok {
				continue
			} // filtering by tsvkeys
			key1 := slice[0] // the key to be used downstream
			dlm := slice[1]
			bits := strings.Split(cell, dlm)
			for _, bit := range bits {
				duos.Add(pairid, key1, bit)
			}
		}
	}
	if len(filters) == 1 {
		msg := fmt.Sprintf("parse.MiTab():%s: NotInBgw: %d", rpth, len(notInBgw))
		log.Println(msg)
		//fmt.Printf("%s\n", msg)
	}
	if len(duos) == 0 {
		msg := fmt.Sprintf("parse.MiTab():%s: NoData", rpth)
		return duos, errors.New(msg)
	}
	return duos, nil
} // MiTab()

func orthosolo(datdir, txid string, txn2prm util.Set2D, idmkeys map[string]string) (util.Set3D, error) {
	solos := make(util.Set3D)
	subdir := "idmapping/"
	ext := ".idmapping"
	prmids := txn2prm[txid].Keys() // normally 1, yet multiple may occur
	for _, prmid := range prmids {
		prmid := fmt.Sprintf("%s%s%s", prmid, "_", txid)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, prmid, ext) // read
		//dat, err := util.FilterByValues(pth, idmkeys, 1, 1, 0)
		dat, err := Idmap(pth, idmkeys, 1, 2, 0)
		util.CheckE(err)
		for src, all := range dat {
			for id, one := range all {
				for upac, _ := range one {
					solos.Add(src, id, upac)
				}
			}
		} // src; e.g. KO, OrthoDB
	}
	count := 0
	for src, _ := range solos {
		count += len(solos[src])
	}
	if count == 0 {
		msg := fmt.Sprintf("parse.orthosolo:%s: NoData", txid)
		return solos, errors.New(msg)
	} // no orthology in any of the sources
	return solos, nil
} // orthosolo()

func OrthoDuo(datdir, txidL, txidR string, txn2prm util.Set2D, idmkeys map[string]string) (util.Set3D, error) {
	duos := make(util.Set3D)
	outL, err := orthosolo(datdir, txidL, txn2prm, idmkeys)
	if err != nil {
		return duos, err
	} // NoDate for one taxaon
	outR, err := orthosolo(datdir, txidR, txn2prm, idmkeys)
	if err != nil {
		return duos, err
	} // NoDate for one taxaon
	for src, _ := range idmkeys {
		allL, ok := outL[src]
		if !ok {
			continue
		}
		allR, ok := outR[src]
		if !ok {
			continue
		} // now data from 'src' is present in both taxa
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
					duos.Add(duoid, src, id)
				}
			}
		}
	}
	if len(duos) == 0 {
		msg := fmt.Sprintf("parse.OrthoDuo():%s--%s: NoData", txidL, txidR)
		return duos, errors.New(msg)
	} // 20200531: none Note: no filtering by BGW
	return duos, nil
} // OrthoDuo()
