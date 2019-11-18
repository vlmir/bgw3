package parse

import (
	"aux"
	"bgw/util"
	"bufio"
	"fmt"
	"os"
	"strings"
)

// parses the output of the system 'ls' utility
// any string can be used as a basename-extension separator
/*
func Basenames(pth, dlm string) (set aux.Set2D, err error) {
	//set := make(aux.Set2D)
	keyind := 0
	valind := 1
	fh, err := os.Open(pth)
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

func Idmap(pth string, srcs map[string]string, ind1, ind2, ind3 int) (aux.Set3D, error) {
	set := make(aux.Set3D)
	fh, err := os.Open(pth)
	if err != nil {
		return set, err
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
		key1 := strings.Replace(cells[ind1], "\"", "''", -1) // present e.g. in 44689
		key2 := strings.Replace(cells[ind2], "\"", "''", -1) // present e.g. in 44689
		key3 := strings.Replace(cells[ind3], "\"", "''", -1) // present e.g. in 44689
		set.Add(key1, key2, key3)
	}
	return set, nil
}

func Upidmap(pth1 string, idmkeys map[string]string) (upcas aux.Set2D, upacs, gnms aux.Set3D, err error) {
	upcas = make(aux.Set2D)
	upacs = make(aux.Set3D)
	gnms = make(aux.Set3D)
	fh1, err := os.Open(pth1)
	if err != nil {
		return upcas, upacs, gnms, err
	}
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
		xrf := strings.Replace(cells[2], "\"", "", -1) // present e.g. in 44689
		upacs.Add(upac, key, xrf)
		upca := strings.Split(upac, "-")[0]
		//aux.Add2keys(upcas, upca, upac)
		upcas.Add(upca, upac)
		upacs.Add(upca, "upac", upac)
		if key == "gnm" {
			gnms.Add(xrf, "upac", upac)
		}
	}
	return upcas, upacs, gnms, nil
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
func Updat(pth0 string, upcas aux.Set2D) (updat, txns aux.Set3D, err error) {
	updat = make(aux.Set3D)
	txns = make(aux.Set3D)
	fh0, err := os.Open(pth0)
	if err != nil {
		err = fmt.Errorf("%s%s", "parse.Updat:os.Open:", err)
		return updat, txns, err
	}
	defer fh0.Close()
	scanner := bufio.NewScanner(fh0)
	Unplaced := 0
	Genome := 0
	WGS := 0
	multiChrom := 0
	multiPome := 0
	notInBgw := 0
	/// by default scans for '\n'
	for scanner.Scan() {
		cells := strings.Split(scanner.Text(), "\t")
		if len(cells) < 8 {
			notInBgw++
			continue
		}
		txid := cells[3]
		upca := cells[0]
		_, ok := upcas[upca] // RefProt ACs
		if !ok {
			notInBgw++
			continue
		} // filtering by RefProt
		onedat := make(aux.Set2D)
		pomes := make(aux.Set2D)
		/// Attn: multiple entries do occur due to diff in chromosomes !!
		items := strings.Split(cells[5], ", ") // proteome: chromosome
		for _, item := range items {
			bits := strings.Split(item, ": ") // proteome: chromosome
			if len(bits) != 2 {
				continue
			}
			pomeid := bits[0]
			come := bits[1]
			// other vals: 'Unassembled WGS sequence', 'Geneome assembly'
			chrmbits := strings.Split(come, " ")
			if len(chrmbits) > 2 {
				WGS++
				continue
			} // skipping 'Unassembled WGS sequence'
			// a simple minded way:
			//come = strings.Replace(come, " ", "-", -1)
			//pomes.Add(pomeid, come)
			// alternatovely:
			chrmbit := chrmbits[0]
			if chrmbit == "Unplaced" {
				Unplaced++
				continue
			}
			if chrmbit == "Genome" {
				Genome++
				continue
			}
			trgs := []string{"Chromosome", "Mitochondrion", "Chloroplast"}
			cnt := 0
			for _, trg := range trgs {
				if trg == chrmbit {
					cnt++
				}
			}
			if cnt > 1 {
				multiChrom++
			}
			if cnt != 1 {
				continue
			}
			comeid := ""
			if len(chrmbits) == 2 { // "e.g. 'Chromosome 1'"
				chrnbr := chrmbits[1]
				if len(chrnbr) == 1 {
					chrnbr = strings.Join([]string{"0", chrnbr}, "")
				}
				comeid = strings.Join([]string{strings.ToLower(chrmbits[0][0:3]), chrnbr}, "-")
			} else if len(chrmbits) == 1 {
				comeid = strings.ToLower(chrmbits[0][0:6]) // e.g. Mitochondrion
			} else {
				continue
			}
			if comeid == "" {
				continue
			}
			pomes.Add(pomeid, comeid)
		}
		if len(pomes) > 1 {
			multiPome++
		}
		if len(pomes) != 1 {
			notInBgw++
			continue
		}
		pome := pomes.Keys()[0]
		onedat.Add("pome", pome)
		for come := range pomes[pome] {
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
	}
	fmt.Println("parse.Updat():notInBgw:", notInBgw)
	fmt.Println("parse.Updat():Unplaced:", Unplaced)
	fmt.Println("parse.Updat():Genome:", Genome)
	fmt.Println("parse.Updat():WGS:", WGS)
	fmt.Println("parse.Updat():multiChrom:", multiChrom) // all: 0
	fmt.Println("parse.Updat():multiPome:", multiPome)   // all: 0, 1, 65
	return updat, txns, nil
}

func Upvar(pth0 string, gsmap aux.Set3D) (pairs aux.Set3D) {
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
	srcL := "uniprot"
	srcR := "omim"
	fh0, err := os.Open(pth0)
	if err != nil {
		return
	}
	defer fh0.Close()
	scanner := bufio.NewScanner(fh0)
	pairs = make(aux.Set3D)
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
		_, ok := gsmap[oriL]
		if !ok {
			notInBgw++
			continue
		}
		idL := fmt.Sprintf("%s%s%s", srcL, "!", oriL)
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
		idR := fmt.Sprintf("%s%s%s", srcR, "!", oriR)
		if idR == "" {
			continue
		}
		pairid := fmt.Sprintf("%s%s%s", idL, "--", idR)
		pairs.Add(pairid, "nmR", nmR)
		pairs.Add(pairid, "upca", upca)
	}
	fmt.Println("parse.Upvar():notInBgw:", notInBgw)
	return pairs
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
func Gaf(pth0 string, upac2bgw aux.Set3D) (bp, cc, mf aux.Set3D) {
	ourppys := map[string]string{
		"C": "gp2cc",
		"F": "gp2mf",
		"P": "gp2bp",
	}
	srcL := "uniprot"
	srcR := "obo"
	bp = make(aux.Set3D)
	cc = make(aux.Set3D)
	mf = make(aux.Set3D)
	fh0, err := os.Open(pth0)
	if err != nil {
		return
	}
	defer fh0.Close()
	scanner := bufio.NewScanner(fh0)
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
	fmt.Println("parse.Gaf():notInBgw:", notInBgw)
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
func Gpa(pth0 string, upac2bgw aux.Set3D) (pairs aux.Set3D) {
	ourppys := map[string]string{
		"part_of":     "gp2cc",
		"enables":     "gp2mf",
		"involved_in": "gp2bp",
	}
	srcL := "uniprot"
	srcR := "obo"
	pairs = make(aux.Set3D)
	fh0, err := os.Open(pth0)
	if err != nil {
		return
	}
	defer fh0.Close()
	scanner := bufio.NewScanner(fh0)
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
		pairs.Add(pairid, "ppy", ppy)
		pairs.Add(pairid, "eco", eco)
		for _, ref := range refs {
			pairs.Add(pairid, "ref", ref)
		}
	}
	fmt.Println("parse.Gpa():notInBgw:", notInBgw)
	return pairs
}

func Mitab(pth0 string, upac2bgw aux.Set3D) (pairs aux.Set3D) {
	tsvkeys := map[int][]string{
		// NB: 0:1 pairs are NOT inique !
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
	pairs = make(aux.Set3D)
	fh0, err := os.Open(pth0)
	if err != nil {
		return
	}
	defer fh0.Close()
	srcL := "uniprot"
	srcR := "uniprot"
	scanner := bufio.NewScanner(fh0)
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
		upacLs := aux.X1type(idAs, "uniprotkb", ":")
		if len(upacLs) == 0 {
			continue
		}
		upacA := upacLs[0]
		upacRs := aux.X1type(idBs, "uniprotkb", ":")
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
				pairs.Add(pairid, key1, bit)
			}
		}
	}
	fmt.Println("parse.Mitab():notInBgw:", notInBgw)
	return pairs
}

func Tftg(pth0 string, upac2bgw aux.Set3D, gsmap aux.Set3D) (aux.Set3D, util.Meta) {
	pairs := make(aux.Set3D)
	meta := util.NewMeta()
	fh0, err := os.Open(pth0)
	if err != nil {
		return pairs, meta
	} // needed in case of unnamed returned vars
	defer fh0.Close()
	scanner := bufio.NewScanner(fh0)
	pairs = make(aux.Set3D)
	notInBgw := 0
	for scanner.Scan() { // by default scans for '\n'
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}
		cells := strings.Split(line, "\t")
		if len(cells) < 28 {
			fmt.Println("Truncated data line for:", cells[0])
			continue
		}
		bits := strings.Split(cells[0], ":") // TF:TG
		tfnm := bits[0]                      // including AP1, NFKB, etc !
		tgnm := bits[1]
		_, ok := gsmap[tgnm]
		if !ok {
			notInBgw++
			continue
		}
		pairid := fmt.Sprintf("%s%s%s", tfnm, "--", tgnm)
		bits = strings.Split(cells[1], "|") // UP ACs, currently for AP1 and NFKB
		if len(bits) == 0 {
			continue
		}
		for _, val := range bits {
			_, ok := upac2bgw[val]
			if !ok {
				notInBgw++
				continue
			} // essentially filtering by RefProt
			pairs.Add(pairid, "upca", val)
		}
		/*
			Metadata
		*/
		src := "extri"
		bits = strings.Split(cells[4], ";")
		for _, val := range bits {
			if val != "" {
				meta.Refs.Add(pairid, src, val)
			}
		}
		bits = strings.Split(cells[3], ";")
		for _, val := range bits {
			if val != "" {
				meta.Cnfs.Add(pairid, src, val)
			}
		}
		src = "htri"
		bits = strings.Split(cells[8], ";")
		for _, val := range bits {
			if val != "" {
				meta.Refs.Add(pairid, src, val)
			}
		}
		bits = strings.Split(cells[9], ";")
		for _, val := range bits {
			if val != "" {
				meta.Cnfs.Add(pairid, src, val)
			}
		}
		src = "trrust"
		bits = strings.Split(cells[12], ";")
		for _, val := range bits {
			if val != "" {
				meta.Refs.Add(pairid, src, val)
			}
		}
		bits = strings.Split(cells[11], ";")
		for _, val := range bits {
			if val != "" {
				meta.Signs.Add(pairid, src, val)
			}
		}
		src = "signor"
		bits = strings.Split(cells[27], ";")
		for _, val := range bits {
			if val != "" {
				meta.Refs.Add(pairid, src, val)
			}
		}
		bits = strings.Split(cells[26], ";")
		for _, val := range bits {
			if val != "" {
				meta.Signs.Add(pairid, src, val)
			}
		}
		src = "tfacts"
		bits = strings.Split(cells[17], ";")
		for _, val := range bits {
			if val != "" {
				meta.Refs.Add(pairid, src, val)
			}
		}
		bits = strings.Split(cells[14], ";")
		for _, val := range bits {
			if val != "" {
				meta.Signs.Add(pairid, src, val)
			}
		}
		bits = strings.Split(cells[18], ";")
		for _, val := range bits {
			if val != "" {
				meta.Cnfs.Add(pairid, src, val)
			}
		}
		src = "goa"
		bits = strings.Split(cells[20], ";")
		for _, val := range bits {
			if val != "" {
				meta.Signs.Add(pairid, src, val)
			}
		}
		src = "intact"
		bits = strings.Split(cells[22], ";")
		for _, val := range bits {
			if val != "" {
				meta.Refs.Add(pairid, src, val)
			}
		}

	}
	fmt.Println("parse.Tftg():notInBgw:", notInBgw)
	return pairs, meta
}
