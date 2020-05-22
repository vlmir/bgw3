package parse

import (
	"bufio"
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"os"
	"strings"
)

func check(e error) {
	if e != nil {
		panic(e)
	}
}

// Tab2set converts tabular data into a map keyed on an arbitrary combination of fields
// Arguments:
// 1. path to the input file (tab separated)
// 2. data structure specifying primary key
// 3. data structure specifying values
// Returns:
// 1. map primary key -> secondary key -< values -> counts
// 2. error
func Tab2set(rpth string, keys, vals []bgw.Column) (out util.Set3D, err error) {
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
	check(err)
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	ln := 0
	for scanner.Scan() { // by default scans for '\n'
		line := scanner.Text()
		ln++
		if len(line) == 0 {
			continue
		}
		if string(line[0]) == "#" {
			continue
		}
		cells := strings.Split(line, "\t")
		var ks string
		if len(cells) < maxind+1 {
			msg := fmt.Sprintf("%s:%d: tooFewtFields", rpth, ln)
			panic(errors.New(msg))
		}
		for i, k := range keys {
			items := strings.Split(cells[k.Ind1], k.Dlm1)
			if i == 0 {
				ks = items[k.Ind2]
			} else {
				ks = fmt.Sprintf("%s%s%s", ks, k.Dlm2, items[k.Ind2])
			}
		}
		for _, v := range vals {
			items := strings.Split(cells[v.Ind1], v.Dlm1) // can be just one
			for _, item := range items {
				parts := strings.Split(item, v.Dlm2) // can be just one
				if v.Ind3 == -1 {                    // db:id
					key := v.Key
					if len(key) == 0 {
						key = parts[0]
					}
					out.Add(ks, key, parts[v.Ind2])
				} else {
					for _, part := range parts {
						if len(part) == 0 {
							continue
						}
						out.Add(ks, v.Key, part)
					}
				}
			}
		}
	}
	return out, nil
}

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

func Idmap(rpth string, srcs map[string]string, ind1, ind2, ind3 int) (util.Set3D, error) {
	out := make(util.Set3D)
	fh, err := os.Open(rpth)
	check(err)
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
		key1 := strings.Replace(cells[ind1], "\"", "''", -1) // present e.g. in 44689
		key2 := strings.Replace(cells[ind2], "\"", "''", -1) // present e.g. in 44689
		key3 := strings.Replace(cells[ind3], "\"", "''", -1) // present e.g. in 44689
		out.Add(key1, key2, key3)
	}
	return out, nil
}

func Upidmap(pthR string, idmkeys map[string]string) (upac2xrf, gnm2upca util.Set3D, err error) {
	upac2xrf = make(util.Set3D)
	gnm2upca = make(util.Set3D)
	fh1, err := os.Open(pthR)
	check(err)
	defer fh1.Close()
	scanner := bufio.NewScanner(fh1)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		if len(cells) != 3 {
			continue
		}
		key, ok := idmkeys[cells[1]]
		if !ok {
			continue
		} // filtering by idmkeys
		upac := cells[0]
		// TODO see if the replacement below is acceptable
		xrf := strings.Replace(cells[2], "\"", "''", -1) // present e.g. in 44689
		upac2xrf.Add(upac, key, xrf)
		upca := strings.Split(upac, "-")[0]
		upac2xrf.Add(upca, "upac", upac)
		upac2xrf.Add(upac, "upca", upca)
		if key == "gnm" {
			gnm2upca.Add(xrf, "upca", upca)
		}
	}
	return upac2xrf, gnm2upca, nil
}

/*
Fields: Dated!
Entry   Entry name      Gene names  (primary )  Gene names  (synonym )  Organism        Organism IDProtein names   Proteomes       Annotation      PubMed ID
Attn: gene Names space separated and some contain internal space - useless
TODO contact UniProt
*/
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
func Updat(pthR string, upacs util.Set3D) (updat, txns util.Set3D, err error) {
	updat = make(util.Set3D)
	txns = make(util.Set3D)
	fhR, err := os.Open(pthR)
	check(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	multiPome := 0
	noProteome := 0
	multiProtChrom := 0
	irregularProtChrom := 0
	notInRefProt := 0
	for scanner.Scan() {
		/// by default scans for '\n'
		// one line per UniProt canonical accession
		cells := strings.Split(scanner.Text(), "\t")
		if len(cells) < 8 {
			continue// does not happen
		}
		txid := cells[3]
		upca := cells[0]
		_, ok := upacs[upca] // RefProt ACs
		if !ok {
			notInRefProt++
			continue
		} // filtering by RefProt
		onedat := make(util.Set2D)
		ptm2chr := make(util.Set2D)
		/// Attn: multiple entries do occur due to diff in chromosomes !!
		if len(cells[5]) == 0 {
			continue// does not happen
		}
		ptmchrs := strings.Split(cells[5], ", ") // proteome: chromosome
		if len(ptmchrs) == 0 {
			continue// does not happen
		}
		if len(ptmchrs) > 1 {
			multiProtChrom++// happens in all the 25 taxa; 10230 in 367110
		}
		/// Processing all Proteome:Chromosome values
		for _, ptmchr := range ptmchrs {
			bits := strings.Split(ptmchr, ": ") // proteome: chromosome
			if len(bits) != 2 {
				//fmt.Printf("parse.Updat:%s: IrregularProtChrom:%v: skipped\n", txid, bits)
				// "Linkage Group I" "Linkage Group II" etc in 367110 
				// "gap-filling sequence" in 284812
				irregularProtChrom++
				continue// happens in 367110 and 284812
			}
			pomeid := bits[0]
			come := bits[1]
			chrmbits := strings.Split(come, " ")
			chrmbit := chrmbits[0]
			trgs := []string{"Chromosome", "Mitochondrion", "Chloroplast", "Apicoplast", "Plasmid"}
			cnt := 0
			for _, trg := range trgs {
				if trg == chrmbit {
					cnt++
				}
			}
			if cnt != 1 {
				//fmt.Printf("parse.Updat:%s: IrregularChrom:%v: skipped\n", txid, come)
			// junk like 'Unassembled WGS sequence', 'Geneome assembly', 'Unplaced' etc
				continue
			}
			comeid := ""
			switch chrmbit {
			case "Chromosome": // "e.g. 'Chromosome 1'"
				chrnbr := chrmbits[1]
				if len(chrnbr) == 1 {
					// chrnbr = strings.Join([]string{"0", chrnbr}, "")
				}
				comeid = strings.Join([]string{strings.ToLower(chrmbits[0][0:3]), chrnbr}, "-")
			case "Mitochondrion":
				comeid = strings.ToLower(chrmbits[0][0:6])
			case "Chloroplast":
				comeid = strings.ToLower(chrmbits[0][0:6])
			case "Apicoplast":
				// "Apicoplast A" "Apicoplat B" in 36329
				comeid = fmt.Sprintf("%s-%s", strings.ToLower(chrmbits[0][0:5]), chrmbits[1])
			case "Plasmid":
				// "Plasmid 2-micron" in 559292 S.cerevisiae
				lastbit := strings.Replace(chrmbits[1], "-", "", -1)
				comeid = fmt.Sprintf("%s-%s", strings.ToLower(chrmbits[0][0:4]), lastbit[0:3])
			}
			if comeid != "" {
				ptm2chr.Add(pomeid, comeid)
			}
		}// one Proteome:Chromosome value
		if len(ptm2chr) > 1 {
			multiPome++
		}
		if len(ptm2chr) == 0 {
			noProteome++
			continue
		}
		pome := ptm2chr.Keys()[0]
		onedat.Add("pome", pome)
		for come := range ptm2chr[pome] {
			come = strings.Replace(come, " ", "-", 1)
			onedat.Add("come", come)
			txns.Add(txid, "come", come)
		}
		onedat.Add("upid", cells[1])
		onedat.Add("spnm", cells[2])
		onedat.Add("txn", txid)
		onedat.Add("pdfn", strings.Replace(cells[4], "\"", ",,", 1))
		onedat.Add("pubmed", cells[6])
		bits := strings.Split(cells[7], " ")
		onedat.Add("score", bits[0])
		txns.Add(txid, "spnm", cells[2]) // NEEDED !!
		txns.Add(txid, "pome", pome)
		updat[upca] = onedat
	}// one UniProt canonical accession
	txids := txns.Keys()
	msg := ""
	if len(txids) != 1 {
		msg = fmt.Sprintf("MultipleTaxa: %v", txids)
		panic(errors.New(msg))
	}
	if noProteome != 0 {
		msg = fmt.Sprintf("parse.Updat():%v: noProteome, skipped: %d", txids, noProteome)
		fmt.Printf("%s\n", msg)
	}
	if irregularProtChrom != 0 {
		// 10231 in 367110, 5 in 284812
		msg = fmt.Sprintf("parse.Updat():%v: irregularProtChrom, skipped: %d", txids, irregularProtChrom)
		fmt.Printf("%s\n", msg)
	}
	if notInRefProt != 0 {
		msg = fmt.Sprintf("parse.Updat():%v: notInRefProt, skipped: %d", txids, notInRefProt)
		fmt.Printf("%s\n", msg)
	}
	if multiProtChrom != 0 {
		// occurs in all the 25 taxa, normally limited number though goes upto 10230 in 367110
		msg = fmt.Sprintf("parse.Updat():%v: Warning: multiProtChrom : %d", txids, multiProtChrom)
		fmt.Printf("%s\n", msg)
	}
	if multiPome != 0 {
		// only once in 4588
		msg = fmt.Sprintf("parse.Updat():%v: Warning: multiPome : %d", txids, multiPome)
		fmt.Printf("%s\n", msg)
	}
	return updat, txns, nil
}

func Upvar(pthR string, gsym2bgw util.Set3D) (duos util.Set3D) {
	// gsym2bgw is used only for filtering
	/*
	 AARS	  P49588	 VAR_063527  p.Arg329His	Disease	   rs267606621 Charcot-Marie-Tooth disease 2N(CMT2N) [MIM:613287]
	  0 gene name
	 10 up acc // always present and single
	 21
	 33 aa change
	 48 type
	 62 '-' possible
	 77 description '-' possible
	*/
	nsL := "hgncsymb"
	nsR := "omim"
	fhR, err := os.Open(pthR)
	check(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	duos = make(util.Set3D)
	notInBgw := 0
	for scanner.Scan() { // by default scans for '\n'
		line := scanner.Text()
		if len(line) < 79 {
			continue
		} // skipping lines with '-' in the last field
		typ := strings.TrimSpace(line[48:61])
		if typ != "Disease" {
			continue
		} // filtering, now all lines longer 78
		//oriL := strings.TrimSpace(line[10:20])
		oriL := strings.TrimSpace(line[0:9]) // Gene Symbol
		_, ok := gsym2bgw[oriL]
		if !ok {
			notInBgw++
			continue
		}
		idL := fmt.Sprintf("%s%s%s", nsL, "!", oriL)
		// TODO consider splitting on '[' and ']'
		upca := strings.TrimSpace(line[10:20])
		nmR := strings.TrimSpace(line[77:])
		if nmR == "" {
			continue
		}
		bits := strings.Split(nmR, ":") // sic, returns 1 value
		if len(bits) < 2 {
			continue
		} // no MIM ID
		oriR := strings.TrimSpace(strings.Split(bits[1], "]")[0])
		idR := fmt.Sprintf("%s%s%s", nsR, "!", oriR)
		if idR == "" {
			continue
		}
		pairid := fmt.Sprintf("%s%s%s", idL, "--", idR)
		duos.Add(pairid, "nmR", nmR)
		duos.Add(pairid, "upca", upca)
	}
	log.Println("parse.Upvar(): notInBgw:", notInBgw)
	return duos
}

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
func Gaf(pthR string, upac2bgw util.Set3D) (bp, cc, mf util.Set3D) {
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
	fhR, err := os.Open(pthR)
	check(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	notInBgw := 0
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		upac := ""
		if cells[0] == "UniProtKB" {
			upac = cells[1]
		} else {
			continue
		}
		// filtering by RefProt
		_, ok := upac2bgw[upac]
		if !ok {
			notInBgw++
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
	log.Println("parse.Gaf(): notInBgw:", notInBgw)
	return bp, cc, mf
}

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
func Gpa(pthR string, upac2bgw util.Set3D) (duos util.Set3D) {
	ourppys := map[string]string{
		"part_of":     "gp2cc",
		"enables":     "gp2mf",
		"involved_in": "gp2bp",
	}
	srcL := "uniprot"
	srcR := "obo"
	duos = make(util.Set3D)
	fhR, err := os.Open(pthR)
	check(err)
	defer fhR.Close()
	scanner := bufio.NewScanner(fhR)
	notInBgw := 0
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), "\t")
		upac := ""
		if cells[0] == "UniProtKB" {
			upac = cells[1]
		} else {
			continue
		}
		// filtering by RefProt
		_, ok := upac2bgw[upac]
		if !ok {
			notInBgw++
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
	log.Println("parse.Gpa(): notInBgw:", notInBgw)
	return duos
}

func Mitab(pthR string, upac2bgw util.Set3D) (duos util.Set3D) {
	tsvkeys := map[int][]string{
		// NB: 0:1 duos are NOT inique !
		//0: []string{"upacA", "|"}, // 9606: all single
		//1: []string{"upacB", "|"}, // 9606: all single
		6:  {"mtds", "|"}, // 9606: all single
		8:  {"pubids", "|"},
		11: {"inactps", "|"}, // 9606: all single
		13: {"inacids", "|"}, // intact:EBI-20559053|imex:IM-26397-1
		14: {"cnfvals", "|"}, // author score:D|intact-miscore:0.37
		//18: []string{"xprlAs", "|"}, // 9606: all single
		//19: []string{"xprlBs", "|"}, // 9606: all single
		//28: []string{"hosts", "|"}, // TODO expunge
	}
	duos = make(util.Set3D)
	fhR, err := os.Open(pthR)
	check(err)
	defer fhR.Close()
	srcL := "uniprot"
	srcR := "uniprot"
	scanner := bufio.NewScanner(fhR)
	notInBgw := 0
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
		_, ok := upac2bgw[upacA]
		if !ok {
			notInBgw++
			continue
		}
		_, ok = upac2bgw[upacB]
		if !ok {
			notInBgw++
			continue
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
	log.Println("parse.Mitab(): notInBgw:", notInBgw)
	return duos
}

// Parses a tab delimited file containing relations between human transcription factors
// and their known target genes extracted from multiple sources (Miguel)
// Arguments:
// 1. path to the data file
// 2. mapping from UniProt Accessions to BGW IDs
// 2. mapping from HGNC Symbols to BGW IDs
// Returns:
// 1. data structure containing core data
// 2. data structure containing metadata

func orthosolo(datdir, txid string, tx2pm util.Set2D, idmkeys map[string]string) (util.Set3D, error) {
	out := make(util.Set3D)
	subdir := "idmapping/"
	ext := ".idmapping"
	pomes := tx2pm[txid].Keys() // normally 1, yet multiple may occur
	for _, pome := range pomes {
		pome := fmt.Sprintf("%s%s%s", pome, "_", txid)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, pome, ext) // read
		dat, err := Idmap(pth, idmkeys, 1, 2, 0)
		check(err)
		for src, all := range dat {
			for id, one := range all {
				for upac, _ := range one {
					out.Add(src, id, upac)
				}
			}
		}
	}
	return out, nil
}
func Orthoduo(datdir, txidL, txidR string, tx2pm util.Set2D, idmkeys map[string]string) (util.Set3D, error) {
	out := make(util.Set3D)
	outL, err := orthosolo(datdir, txidL, tx2pm, idmkeys)
	check(err)
	outR, err := orthosolo(datdir, txidR, tx2pm, idmkeys)
	check(err)
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
					out.Add(duoid, src, id)
				}
			}
		}
	}
	return out, nil
}
