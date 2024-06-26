package export

// TODO GLOBAL SIO_000253 -> evidenceOrigin

import (
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/semweb"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"os"
	"strings"
)

func newFH(wpth string) *os.File {
	fh, err := os.Create(wpth)
	if err != nil {
		msg := fmt.Sprintf("%s: %s: os.Create(): %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	return fh
}

func counter(s []string, c util.Set2D, a, d, v string) (l int) {
	// used only in Gene2phen and Prot2go
	l = len(s)
	if l == 0 {
		c.Add(d, v)
	} else {
		c.Add(a, v)
	}
	return l
}

func gene2bgwg(genes []string, xmap bgw.Xmap) []string {
	// used only in Tfac2gene
	set := make(util.Set2D)
	for _, gene := range genes {
		for _, bgwg := range xmap.Ncbig[gene]["bgwg"].Keys() {
			set.Add(bgwg, gene)
		} // by NCBI gene id
		if len(set) == 0 {
			for _, bgwg := range xmap.Lblg[gene]["bgwg"].Keys() {
				set.Add(bgwg, gene)
			} // by gene symbol
		}
		if len(set) == 0 {
			for _, bgwg := range xmap.Syng[gene]["bgwg"].Keys() {
				set.Add(bgwg, gene)
			} // by gene synonym
		}
		if len(set) == 0 {
			for _, bgwg := range xmap.Ensg[gene]["bgwg"].Keys() {
				set.Add(bgwg, gene)
			} // by Ensembl gene ID
		}
	}
	return set.Keys()
}

func gene2bgwp(genes []string, xmap bgw.Xmap) []string {
	// used only in Tfac2gene
	set := make(util.Set2D)
	for _, gene := range genes {
		for _, bgwp := range xmap.Ncbig[gene]["bgwp"].Keys() {
			set.Add(bgwp, gene)
		} // by NCBI gene id
		if len(set) == 0 {
			for _, bgwp := range xmap.Ensg[gene]["bgwp"].Keys() {
				set.Add(bgwp, gene)
			} // by Ensembl gene ID
		}
	}
	return set.Keys()
}

func upac2bgwp(upacs []string, xmap bgw.Xmap) []string {
	// used in Reg2targ, SigPways, Tfac2gene
	set := make(util.Set2D)
	for _, upac := range upacs {
		upca := strings.Split(upac, "-")[0]
		for _, bgwp := range xmap.Upac[upca]["bgwp"].Keys() {
			set.Add(bgwp, upac)
		}
	}
	return set.Keys()
}

func Prot2prot(d *bgw.Dat4bridge, x *bgw.Xmap, wdir string) error {
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"sub2cls",
		"stm2sbj",
		"stm2obj",
		"stm2pdc",
		"tlp2tlp",
		"ins2cls",
		"sth2src",
		"sth2evd",
		"sth2mtd",
		"sub2set",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"evd2lvl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	header, nln := rdf.Capita(keys4b)
	if nln < 28 { // 23? // special case
		msg := fmt.Sprintf("export.Prot2prot(): rdf.Capita(%v): MalformedHeader", keys4b)
		panic(errors.New(msg))
	}
	flags := make(util.Set1D) // for printing the header only once per file
	var srcs = map[string]string{
		"uniprot": "http://uniprot.org/uniprot",
		"intact":  "http://identifiers.org/intact",
	}
	nss := rdf.Nss // BGW URI name spaces
	d4b := *d
	srck := d4b.Src
	txid := d4b.Taxid
	sdir := "prot2prot"
	wpth := fmt.Sprintf("%s%s/%s-%s.nt", wdir, sdir, srck, txid)
	wfh := newFH(wpth)
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	srcU := rdf.FormU(srcs[srck])

	duoNS := "http://rdf.biogateway.eu/prot-prot/"
	pdck := "tlp2tlp"
	xmap := *x
	cnts := d4b.Cnts // Set2D

	duos := d4b.Duos
	for _, duoid := range duos.Keys() {
		duo := duos[duoid]
		idA := ""
		idB := ""
		ids := duo["uniprotkb"].Keys()
		if len(ids) > 2 {
			continue
		} // sanity check
		if len(ids) == 1 {
			idA = strings.Split(ids[0], "-")[0] // UP canonical accession
			idB = idA
		} else {
			idA = strings.Split(ids[0], "-")[0] // UP canonical accession
			idB = strings.Split(ids[1], "-")[0] // UP canonical accession
		} // ids sorted
		bgwAs := xmap.Upac[idA]["bgwp"].Keys() // needed for filtering by RefProt
		bgwBs := xmap.Upac[idB]["bgwp"].Keys() // needed for filtering by RefProt
		if len(bgwAs) == 0 || len(bgwBs) == 0 {
			continue
		} // filterong by Reference Proteome
		if len(bgwAs) > 1 {
			msg := fmt.Sprintf("export.prot2prot(): Multiple BGW ids: idA: %s bgwAs: %v", idA, bgwAs)
			log.Println(msg)
		} // normally should never happen
		if len(bgwBs) > 1 {
			msg := fmt.Sprintf("export.prot2prot(): Multiple BGW ids: idB: %s bgwBs: %v", idB, bgwBs)
			log.Println(msg)
		} // normally should never happen
		cnts.Add(pdck, srck)
		bgwA := bgwAs[0]
		bgwB := bgwBs[0]
		duoid = fmt.Sprintf("uniprot!%s--uniprot!%s", bgwA, bgwB) // redefining
		duoU := rdf.CompU(duoNS, duoid)
		uriA := rdf.CompU(nss["uniprot"], bgwA)
		uriB := rdf.CompU(nss["uniprot"], bgwB)
		sb.WriteString(rdf.FormT(duoU, ourUs["ins2cls"], rdf.CompU(nss["owl"], "Class")))
		sb.WriteString(rdf.FormT(duoU, ourUs["sub2cls"], ourUs["stm"]))
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2lbl"], rdf.FormL(duoid)))
		duodfn := fmt.Sprintf("%s %s %s", bgwA, rdf.Opys[pdck][2], bgwB)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2dfn"], rdf.FormL(duodfn)))
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(nss["rdf"], "predicate"), ourUs[pdck]))
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(nss["rdf"], "subject"), uriA))
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(nss["rdf"], "object"), uriB))
		sb.WriteString(rdf.FormT(uriA, ourUs[pdck], uriB))
		sb.WriteString(rdf.FormT(uriB, ourUs[pdck], uriA))

		/// INSTANCES
		insid := fmt.Sprintf("%s%s%s", duoid, "#", "intact")
		insU := rdf.CompU(duoNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], duoU))
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		/// metadata, aggregated for all subjects and objects
		/// intaracton IDs
		for _, key := range duo["intact"].Keys() {
			myU := rdf.CompU(nss["intact"], key)
			sb.WriteString(rdf.FormT(insU, ourUs["sub2set"], myU)) // Attn: change prop
		}
		/// publications
		for _, key := range duo["pubmed"].Keys() {
			myU := rdf.CompU(nss["pubmed"], key)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], myU))
		}
		/// confidence values
		for _, key := range duo["intact-miscore"].Keys() {
			sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(key)))
		}
		/// interaction types
		for _, key := range duo["typeABid"].Keys() {
			myU := rdf.CompU(nss["obo"], strings.Replace(key, ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], myU))
		}
		/// detection methods
		for _, key := range duo["mtd"].Keys() {
			myU := rdf.CompU(nss["obo"], strings.Replace(key, ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], myU))
		}

		if cnts[pdck][srck] > 0 {
			flags.Add(pdck)
			if flags[pdck] == 1 {
				wfh.Write([]byte(header))
			}
			wfh.Write([]byte(sb.String()))
			flags[pdck]++
		}
		sb.Reset()
	} // duoid
	return nil
} // Prot2prot

// TODO re-implement for a single property (pdck)
func Reg2targ(d *bgw.Dat4bridge, x *bgw.Xmap, wdir string) error {
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"sub2cls",
		"stm2sbj",
		"stm2obj",
		"stm2pdc",
		"reg2ptrg",
		"reg2ntrg",
		"reg2utrg",
		"ins2cls",
		"sth2src",
		"sth2evd",
		"sth2mtd",
		"mi2bp",
		"sth2rlm", // 'interaction type' TODO find a better property - rdf:type ?
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"evd2lvl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	header, nln := rdf.Capita(keys4b)
	if nln < 34 {
		msg := fmt.Sprintf("export.Reg2targ(): rdf.Capita(%v): MalformedHeader", keys4b)
		panic(errors.New(msg))
	}
	d4b := *d
	xmap := *x
	srck := d4b.Src
	txid := d4b.Taxid
	wpths := d4b.Out
	cnts := d4b.Cnts // Set2D
	// just sortcuts
	p2t := "reg2ptrg"
	n2t := "reg2ntrg"
	u2t := "reg2utrg"
	modes := make(util.Set3D)
	modes.Add("signor", "up-regulates", p2t)
	modes.Add("signor", "down-regulates", n2t)

	fhs := make(map[string]*os.File)
	fhs[n2t] = newFH(wpths[n2t])
	fhs[p2t] = newFH(wpths[p2t])
	fhs[u2t] = newFH(wpths[u2t])
	for _, wfh := range fhs {
		defer wfh.Close()
	}
	ourUs := rdf.FmtURIs(keys4b)
	nss := rdf.Nss // BGW URI name spaces
	//srcns := nss[srck] // fully qualified namespace
	srcU := rdf.FormU("http://signor.uniroma2.it")
	rdfns := nss["rdf"]
	// graphns := fmt.Sprintf("%s%s", nss["bgw"], "graph/")

	duos := d4b.Duos
	// future extensions:
	// psi-mi:"MI:0328"(small molecule)
	// psi-mi:"MI:2258"(xenobiotic) # chemical?
	mytypes := []string{"MI:0314", "MI:1304", "MI:0326", "MI:0328"} // for filetering
	myRegTypes := []string{}                                        // for filetering
	entAns := nss["uniprot"]                                        // sic, never changes
	entBns := nss["uniprot"]                                        // sic, never changes
	bgwNS := "http://rdf.biogateway.eu/"
	flags := make(util.Set1D) // for printing the header only once per file
	for _, duokey := range duos.Keys() {
		duo := duos[duokey]
		oriAid := duo["Aid"].Keys()[0] // justly assuming a single value
		oriBid := duo["Bid"].Keys()[0] // justly assuming a single value
		Abits := strings.Split(oriAid, "\"")
		Bbits := strings.Split(oriBid, "\"")
		// chemicals are in the 'else' clause
		if len(Abits) == 1 {
			bits := strings.Split(oriAid, ":")
			oriAid = bits[1]
		} else {
			oriAid = strings.Replace(Abits[1], ":", "_", 1)
		}
		if len(Bbits) == 1 {
			bits := strings.Split(oriBid, ":")
			oriBid = bits[1]
			if bits[0] == "uniprotkb" {
				oriBid = strings.Split(bits[1], "-")[0]
			}
		} else {
			oriBid = strings.Replace(Bbits[1], ":", "_", 1)
		}
		var oriAids []string // UP ACs for one duokey
		var oriBids []string // UP ACs for one duokey
		// filtering by the level of regulation
		if len(myRegTypes) != 0 {
			if len(util.Shared(myRegTypes, duo["reglevelid"].Keys())) == 0 {
				continue
			}
		}
		typeAs := duo["typeAid"].Keys()
		typeBs := duo["typeBid"].Keys()
		// filtering by the type of entities
		if len(util.Shared(mytypes, typeAs)) == 0 {
			continue
		}
		if len(util.Shared(mytypes, typeBs)) == 0 {
			continue
		}
		// skipping pairs with ambiguous entyty types
		if len(typeAs) > 1 {
			msg := fmt.Sprintf("export.Reg2targ(): multiple entity A types: %s:%s, skipping: ", duokey, typeAs)
			fmt.Printf("%s\n", msg)
			continue // 2024-04: 87
		}
		if len(typeBs) > 1 {
			msg := fmt.Sprintf("export.Reg2targ(): multiple entity B types: %s:%s, skipping: ", duokey, typeBs)
			fmt.Printf("%s\n", msg)
			continue // 2024-04: 0
		}
		typeA := typeAs[0] // assuming a single vwlue
		typeB := typeBs[0] // assuming a single value
		/// expanding families and complexes
		if typeA == "MI:0314" || typeA == "MI:1304" {
			oriAids = xmap.Signor[oriAid]["ids"].Keys() // all UP accessions
			if len(oriAids) == 0 {
				continue
			}
		}
		if typeB == "MI:0314" || typeB == "MI:1304" {
			oriBids = xmap.Signor[oriBid]["ids"].Keys()
			if len(oriBids) == 0 {
				continue
			}
		}
		// for chemicals? TODO confirm
		if len(oriAids) == 0 {
			oriAids = append(oriAids, oriAid)
		}
		if len(oriBids) == 0 {
			oriBids = append(oriBids, oriBid)
		}

		// converting oriAids to Bgwids
		bgwAids := upac2bgwp(oriAids, xmap)
		// converting oriBids to Bgwids
		bgwBids := upac2bgwp(oriBids, xmap)
		if len(bgwAids) == 0 || len(bgwBids) == 0 {
			continue
		}

		//  assigning predicate depending of  the direction of regulation
		pdcks := make(util.Set1D) // all predicate keys for one duokey, u2t first
		pdcks.Add(u2t)
		for _, onemod := range duo["modelbl"].Keys() {
			onemod = util.StripParQuots(onemod)
			onemod = strings.Split(onemod, " ")[0]
			// keys for predicates
			for _, pdck := range modes[srck][onemod].Keys() {
				pdcks.Add(pdck)
			}
		} // onemod

		for _, pdck := range pdcks.Keys() {
			clsns := fmt.Sprintf("%s%s/", bgwNS, pdck)
			var sb strings.Builder
			for _, entAid := range bgwAids {
				for _, entBid := range bgwBids {
					cnts.Add(pdck, srck)
					clsid := fmt.Sprintf("uniprot!%s--uniprot!%s", entAid, entBid)
					clsU := rdf.CompU(clsns, clsid)
					entAU := rdf.CompU(entAns, entAid)
					entBU := rdf.CompU(entBns, entBid)
					sb.WriteString(rdf.FormT(clsU, ourUs["ins2cls"], rdf.CompU(nss["owl"], "Class")))
					sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
					sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					clsdfn := fmt.Sprintf("%s %s %s", entAid, rdf.Opys[pdck][2], entBid)
					sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "predicate"), ourUs[pdck]))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "subject"), entAU))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "object"), entBU))
					sb.WriteString(rdf.FormT(entAU, ourUs[pdck], entBU))

					/// INSTANCES
					insid := fmt.Sprintf("%s#%s", clsid, srck)
					insU := rdf.CompU(clsns, insid)
					sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
					for _, id := range duo["reglevelid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["mi2bp"], rdf.CompU(nss["obo"], id)))
					}
					for _, id := range duo["typeABid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["sth2rlm"], rdf.CompU(nss["obo"], id))) // TODO find a better prpperty
					}
					for _, key := range duo["pubmed"].Keys() {
						pubmedU := rdf.CompU(nss["pubmed"], key)
						sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], pubmedU))
					}
					for _, key := range duo["score"].Keys() {
						sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(key)))
					}
					/// detection methods
					for _, key := range duo["mtd"].Keys() {
						myU := rdf.CompU(nss["obo"], strings.Replace(key, ":", "_", 1))
						sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], myU))
					}

				} // entBid
			} // entAid
			if cnts[pdck][srck] > 0 {
				flags.Add(pdck)
				fh := fhs[pdck]
				if flags[pdck] == 1 {
					fh.Write([]byte(header))
				}
				fh.Write([]byte(sb.String()))
				flags[pdck]++
			}
		} // pdck
	} // duokey
	log.Println("export.Reg2targ():", srck, txid, cnts)
	return nil
} // Reg2targ

// TODO re-implement for a single property (pdck)
func Tfac2gene(d *bgw.Dat4bridge, x *bgw.Xmap, wdir string) error {
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"sub2cls",
		"stm2sbj",
		"stm2obj",
		"stm2pdc",
		"reg2ptrg",
		"reg2ntrg",
		"reg2utrg",
		"ins2cls",
		"sth2src",
		"sth2evd",
		"sth2mtd",
		"mi2bp",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"evd2lvl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	header, nln := rdf.Capita(keys4b)
	if nln < 32 {
		msg := fmt.Sprintf("export.Tfac2gene(): rdf.Capita(%v): MalformedHeader", keys4b)
		panic(errors.New(msg))
	}
	d4b := *d
	xmap := *x
	srck := d4b.Src
	txid := d4b.Taxid
	wpths := d4b.Out
	cnts := d4b.Cnts // Set2D
	// just sortcuts
	p2t := "reg2ptrg"
	n2t := "reg2ntrg"
	u2t := "reg2utrg"
	modes := make(util.Set3D) // used only for compiling sets of predicates for each duo
	modes.Add("tflink", "activator", p2t)
	modes.Add("tflink", "repressor", n2t)
	modes.Add("atregnet", "Activation", p2t)
	modes.Add("atregnet", "Repression", n2t)
	modes.Add("atregnet", "Represson", n2t)
	modes.Add("coltri", "pos", p2t)
	modes.Add("coltri", "neg", n2t)
	fhs := make(map[string]*os.File)
	// newFH(wpth) panics if fails to open wpth
	fhs[n2t] = newFH(wpths[n2t])
	fhs[p2t] = newFH(wpths[p2t])
	fhs[u2t] = newFH(wpths[u2t])
	for _, wfh := range fhs {
		defer wfh.Close()
	}
	ourUs := rdf.FmtURIs(keys4b)
	nss := rdf.Nss // BGW URI name spaces
	srcU := rdf.FormU(rdf.Uris4tftg[srck])
	rdfns := nss["rdf"]
	entAns := nss["uniprot"] // sic, never changes
	entBns := nss["gene"]    // sic, never changes
	bgwNS := "http://rdf.biogateway.eu/"
	flags := make(util.Set1D) // for printing the header only once per file

	duos := d4b.Duos
	for _, duokey := range duos.Keys() {
		duo := duos[duokey]
		var oriAids []string
		var oriBids []string
		if srck == "coltri" {
			oriAids = duo["Aupca"].Keys()
			oriBids = duo["Bglbl"].Keys()
		} else if srck == "atregnet" {
			oriAids = duo["Agloc"].Keys()
			oriBids = duo["Bgloc"].Keys()
		} else {
			oriAids = duo["uniprot"].Keys() // UP canonical ACs for one duokey
			oriBids = duo["ncbig"].Keys()   // NCBI GeneID, single
		}

		if len(oriAids) == 0 || len(oriBids) == 0 {
			continue
		}

		// converting oriAids to Bgwids
		bgwAids := upac2bgwp(oriAids, xmap)
		if srck == "atregnet" {
			bgwAids = gene2bgwp(oriAids, xmap)
		}
		// converting oriBids to Bgwids; BGW genes
		bgwBids := gene2bgwg(oriBids, xmap)
		if len(bgwAids) == 0 || len(bgwBids) == 0 {
			continue
		}

		// special treatment for Collectri data:
		// 2 columns:is_stimulation & is_inhibition, True||Fals, single values
		for _, mode := range modes["coltri"].Keys() {
			if len(duo[mode].Keys()) != 0 {
				val := duo[mode].Keys()[0]
				if val == "True" {
					duo.Add("mode", mode) // pos||neg
				}
			}
		}
		//  assigning predicate depending of  the direction of regulation
		pdcks := make(util.Set1D) // all predicate keys for one duokey, u2t first
		pdcks.Add(u2t)
		for _, onemod := range duo["mode"].Keys() {
			// keys for predicates
			for _, pdck := range modes[srck][onemod].Keys() {
				pdcks.Add(pdck)
			}
		} // onemod

		for _, pdck := range pdcks.Keys() {
			clsns := fmt.Sprintf("%s%s/", bgwNS, pdck)
			var sb strings.Builder
			for _, entAid := range bgwAids {
				entAU := rdf.CompU(entAns, entAid)
				for _, entBid := range bgwBids {
					cnts.Add(pdck, srck)
					clsid := fmt.Sprintf("uniprot!%s--gene!%s", entAid, entBid)
					clsU := rdf.CompU(clsns, clsid)
					entBU := rdf.CompU(entBns, entBid)
					sb.WriteString(rdf.FormT(clsU, ourUs["ins2cls"], rdf.CompU(nss["owl"], "Class")))
					sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
					sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					clsdfn := fmt.Sprintf("%s %s %s", entAid, rdf.Opys[pdck][2], entBid)
					sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "predicate"), ourUs[pdck]))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "subject"), entAU))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "object"), entBU))
					sb.WriteString(rdf.FormT(entAU, ourUs[pdck], entBU))
					id := "MI_2247" // transcriptional regulation
					sb.WriteString(rdf.FormT(clsU, ourUs["mi2bp"], rdf.CompU(nss["obo"], id)))

					/// INSTANCES
					insid := fmt.Sprintf("%s#%s", clsid, srck)
					insU := rdf.CompU(clsns, insid)
					sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
					for _, id := range duo["mtdid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], rdf.CompU(nss["obo"], id)))
					} // only for tflink
					// clean up of the mess in the data
					for _, key := range duo["pubmed"].Keys() {
						b := util.IsDigital(key)
						if !b {
							continue
						}
						pubmedU := rdf.CompU(nss["pubmed"], key)
						sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], pubmedU))
					}
					for _, key := range duo["score"].Keys() {
						if srck == "tflink" {
							if key == "Small-scale.evidence:Yes" {
								key = "High"
							} else {
								key = "Low"
							}
						} else if srck == "atregnet" {
							if key == "Confirmed" || key == "confirmed" {
								key = "High"
							} else {
								key = "Low"
							}
						}
						sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(key)))
					} // scores
				} // entBid
			} // entAid
			if cnts[pdck][srck] > 0 {
				flags.Add(pdck)
				fh := fhs[pdck]
				if flags[pdck] == 1 {
					fh.Write([]byte(header))
				} // header is written once per file and only if there are duos
				fh.Write([]byte(sb.String())) // some files are empty
				flags[pdck]++
			}
		} // pdck
	} // duokey
	log.Println("export.Tfac2gene():", srck, txid, cnts)
	return nil
} // Tfac2gene

// TODO re-implement for a single property (pdck)
func SigPways(d *bgw.Dat4bridge, x *bgw.Xmap, wdir string) error {
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"sub2cls",
		"stm2sbj",
		"stm2obj",
		"stm2pdc",
		"reg2ptrg",
		"reg2ntrg",
		"reg2utrg",
		"reg2dtrg", // pways specific
		"reg2itrg", // pways specific
		"ins2cls",
		"sth2src",
		"sth2evd",
		"mi2loc",    // 'occurs in'; pways specific
		"step2pway", // 'part of'; pways specific
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"evd2lvl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	header, nln := rdf.Capita(keys4b)
	if nln < 36 {
		msg := fmt.Sprintf("export.SigPways(): rdf.Capita(%v): MalformedHeader", keys4b)
		panic(errors.New(msg))
	}
	d4b := *d
	xmap := *x
	srck := d4b.Src
	txid := d4b.Taxid
	wpths := d4b.Out
	cnts := d4b.Cnts // Set2D
	p2t := "reg2ptrg"
	n2t := "reg2ntrg"
	u2t := "reg2utrg"
	modes := make(util.Set3D)
	modes.Add("signor", "up-regulates", p2t)
	modes.Add("signor", "down-regulates", n2t)

	fhs := make(map[string]*os.File)
	fhs[n2t] = newFH(wpths[n2t])
	fhs[p2t] = newFH(wpths[p2t])
	fhs[u2t] = newFH(wpths[u2t])
	for _, wfh := range fhs {
		defer wfh.Close()
	}
	ourUs := rdf.FmtURIs(keys4b)
	nss := rdf.Nss // BGW URI name spaces
	srcU := rdf.FormU("http://signor.uniroma2.it")
	rdfns := nss["rdf"]
	// graphns := fmt.Sprintf("%s%s", nss["bgw"], "graph/")

	duos := d4b.Duos
	mytypes := []string{"complex", "protein", "proteinfamily"} // for filetering
	entAns := nss["uniprot"]                                   // sic, never changes
	entBns := nss["uniprot"]                                   // sic, never changes
	bgwNS := "http://rdf.biogateway.eu/"
	flags := make(util.Set1D) // for printing the header only once per file
	for _, duokey := range duos.Keys() {
		duo := duos[duokey]
		// if duo["txid"].Keys()[0] != txid {continue} // filtering by host
		oriAid := duo["Aid"].Keys()[0] // assuming a single value
		oriBid := duo["Bid"].Keys()[0] // assuming a single value

		var oriAids []string // UP ACs for one duokey
		var oriBids []string // UP ACs for one duokey
		typeAs := duo["typeAlbl"].Keys()
		typeBs := duo["typeBlbl"].Keys()
		// currently 3 types for both A and B entities
		// filtering by the type of entities
		if len(util.Shared(mytypes, typeAs)) == 0 {
			continue
		}
		if len(util.Shared(mytypes, typeBs)) == 0 {
			continue
		}
		// skipping pairs with ambiguous entyty types - never happens
		if len(typeAs) > 1 {
			msg := fmt.Sprintf("export.SigPways(): multiple entity A types: %s:%s, skipping: ", duokey, typeAs)
			fmt.Printf("%s\n", msg)
			continue // should never happen TODO return error
		}
		if len(typeBs) > 1 {
			msg := fmt.Sprintf("export.SigPways(): multiple entity B types: %s:%s, skipping: ", duokey, typeBs)
			fmt.Printf("%s\n", msg)
			continue // should never happen TODO return error
		}
		typeA := typeAs[0]
		typeB := typeBs[0]
		/// expanding families and complexes
		if typeA == "complex" || typeA == "proteinfamily" {
			oriAids = xmap.Signor[oriAid]["ids"].Keys() // all UP accessions
			if len(oriAids) == 0 {
				continue
			}
		}
		if typeB == "complex" || typeB == "proteinfamily" {
			oriBids = xmap.Signor[oriBid]["ids"].Keys() // all UP accessions
			if len(oriBids) == 0 {
				continue
			}
		}
		if typeA == "protein" {
			oriAids = append(oriAids, oriAid)
		}
		if typeB == "protein" {
			oriBids = append(oriBids, oriBid)
		}
		// now oriAids oriBids are all UP AC

		// converting oriAids to Bgwids
		bgwAids := upac2bgwp(oriAids, xmap)
		// converting oriBids to Bgwids
		bgwBids := upac2bgwp(oriBids, xmap)
		if len(bgwAids) == 0 || len(bgwBids) == 0 {
			continue
		}

		//  assigning predicate depending of  the direction of regulation
		pdcks := make(util.Set1D) // all predicate keys for one duokey, u2t first
		pdcks.Add(u2t)
		for _, onemod := range duo["modelbl"].Keys() {
			onemod = strings.Split(onemod, " ")[0]
			// keys for predicates
			for _, pdck := range modes[srck][onemod].Keys() {
				pdcks.Add(pdck)
			}
		} // onemod

		for _, pdck := range pdcks.Keys() {
			clsns := fmt.Sprintf("%s%s/", bgwNS, pdck)
			var sb strings.Builder
			for _, entAid := range bgwAids {
				for _, entBid := range bgwBids {
					cnts.Add(pdck, srck)
					clsid := fmt.Sprintf("uniprot!%s--uniprot!%s", entAid, entBid)
					clsU := rdf.CompU(clsns, clsid)
					entAU := rdf.CompU(entAns, entAid)
					entBU := rdf.CompU(entBns, entBid)
					sb.WriteString(rdf.FormT(clsU, ourUs["ins2cls"], rdf.CompU(nss["owl"], "Class")))
					sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
					sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					clsdfn := fmt.Sprintf("%s %s %s", entAid, rdf.Opys[pdck][2], entBid)
					sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "predicate"), ourUs[pdck]))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "subject"), entAU))
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfns, "object"), entBU))
					sb.WriteString(rdf.FormT(entAU, ourUs[pdck], entBU))

					/// INSTANCES
					insid := fmt.Sprintf("%s#%s", clsid, srck)
					insU := rdf.CompU(clsns, insid)
					sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
					for _, id := range duo["cellid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["mi2loc"], rdf.CompU(nss["obo"], id)))
					}
					for _, id := range duo["tissueid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["mi2loc"], rdf.CompU(nss["obo"], id)))
					}
					for _, key := range duo["pubmed"].Keys() {
						oneU := rdf.CompU(nss["pubmed"], key)
						sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], oneU))
					}
					for _, key := range duo["pwayid"].Keys() {
						oneU := rdf.CompU(nss["sigpway"], key)
						sb.WriteString(rdf.FormT(insU, ourUs["step2pway"], oneU))
					}
					for _, key := range duo["score"].Keys() {
						sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(key)))
					}
					for _, key := range duo["isdirect"].Keys() {
						if key == "t" {
							sb.WriteString(rdf.FormT(entAU, ourUs["reg2dtrg"], entBU))
						} else if key == "f" {
							sb.WriteString(rdf.FormT(entAU, ourUs["reg2itrg"], entBU))
						}
					}
				} // entBid
			} // entAid
			if cnts[pdck][srck] > 0 {
				flags.Add(pdck)
				fh := fhs[pdck]
				if flags[pdck] == 1 {
					fh.Write([]byte(header))
				}
				fh.Write([]byte(sb.String()))
				flags[pdck]++
			}
		} // pdck
	} // duokey
	log.Println("export.SigPways():", srck, txid, cnts)
	return nil
} // SigPways

// Note: no isoforms in this graph
func Gene2phen(duos, gsym2bgw util.Set3D, wpth string) (int, error) {
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"sub2cls",
		"stm2sbj",
		"stm2obj",
		"stm2pdc",
		"gn2phn",
		"ins2cls",
		"sth2src",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	header, nln := rdf.Capita(keys4b)
	if nln < 20 {
		msg := fmt.Sprintf("export.Gene2phen(): rdf.Capita(%v): MalformedHeader", keys4b)
		panic(errors.New(msg))
	}
	nss := rdf.Nss // BGW URI name spaces
	srck := "uniprot"
	var srcs = map[string]string{
		"uniprot": "http://uniprot.org/uniprot",
	}
	srcU := rdf.FormU(srcs[srck])
	clsU := rdf.CompU(nss["owl"], "Class")
	// gene2phen graph ini
	wfh := newFH(wpth)
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	sb.WriteString(header)
	nln = 0
	stmNS := "http://rdf.biogateway.eu/gene-phen/"
	rdfNS := nss["rdf"]
	cnt := make(util.Set2D) // genes absent in BGW
	cntD := 0               // accepted duos
	for _, duoid := range duos.Keys() {
		duo := duos[duoid]
		bits := strings.Split(duoid, "--")
		idL := bits[0]
		idR := bits[1]
		oriL := strings.Split(idL, "!")[1] // Gene Name
		bgwLs := gsym2bgw[oriL]["bgwg"].Keys()
		cntLs := len(bgwLs)
		if cntLs != 1 {
			msg := fmt.Sprintf("export.Gene2phen():%s: bgwLs: %v", oriL, bgwLs)
			fmt.Printf("%s\n", msg)
		} // 2303: 1 missing BGW gene, likely due to RefProt filtering
		// filtering, superfluous, done anyway by looping over bgwLs
		l := 0
		if l = counter(bgwLs, cnt, "addG", "dropG", oriL); l == 0 {
			continue
		}
		oriR := strings.Split(idR, "!")[1] // OMIM ID
		duoU := rdf.CompU(stmNS, duoid)
		sb.WriteString(rdf.FormT(duoU, ourUs["ins2cls"], clsU))
		sb.WriteString(rdf.FormT(duoU, ourUs["sub2cls"], ourUs["stm"]))
		clslbl := fmt.Sprintf("%s--%s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
		upcas := duo["upca"].Keys()
		if len(upcas) > 1 {
			// msg := fmt.Sprintf("export.Gene2phen():%s: upcas: %v", duoid, upcas)
			// fmt.Printf("%s\n", msg)
		} // 20200531: 2 using symG in the key, all the 4 accs with a single dfn
		// 20200531: 69 with > 1 dfns using symG in the key; the same MIM ID indeed
		dfns := duo["dfn"].Keys()
		if len(dfns) != 1 {
			msg := fmt.Sprintf("export.Gene2phen():%s:%v: %v dfns: %v", duoid, upcas, len(dfns), dfns)
			fmt.Printf("%s\n", msg)
		} // 230303: 74
		dfn := strings.Join(dfns, "; ")
		clsdfn := fmt.Sprintf("Association between gene %s and disease %v", oriL, dfn)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))

		// an attempt to get desease names in the App, to no avail
		uriR := rdf.CompU(nss["omim"], oriR)
		sb.WriteString(rdf.FormT(uriR, ourUs["sth2lbl"], rdf.FormL(dfn)))

		pdc := "gn2phn"
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		// multiple subjects (never happens)
		for _, bgwL := range bgwLs { // sorted above
			uriL := rdf.CompU(nss["gene"], bgwL)
			sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "subject"), uriL))
			sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
		}
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "object"), uriR))

		/// INSTANCES
		insid := fmt.Sprintf("%s%s%s", duoid, "#", srck)
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], duoU))
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		bytes := []byte(sb.String())
		if len(bytes) != 0 {
			wfh.Write(bytes)
			sb.Reset()
			cntD++
		}
	} // duoid
	msg := ""
	if cntD == 0 {
		msg = fmt.Sprintf("export.Prot2phen(): NoDuos") // sic!
		panic(errors.New(msg))
	}
	msg = fmt.Sprintf("export.Gene2phen(): Pairs: added: %d dropped: %d", cntD, len(duos)-cntD)
	log.Println(msg)
	msg = fmt.Sprintf("export.Gene2phen(): Genes: added: %d dropped: %d", len(cnt["addG"]), len(cnt["dropG"]))
	log.Println(msg)
	return cntD, nil
} // Gene2phen

// arg1: output of parse.Gaf or parse.Gpa
// arg2: mapping fromn UniProt accession to BGW IDs generated by GeneProt()
// arg3: path for exporting the RDF file
func Prot2go(duos, upac2bgw util.Set3D, wpth string) (int, error) {
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"sub2cls",
		"stm2sbj",
		"stm2obj",
		"stm2pdc",
		"gp2bp",
		"gp2cc",
		"gp2mf",
		"ins2cls",
		"sth2src",
		"sth2evd",
		"sth2mtd",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	header, nln := rdf.Capita(keys4b)
	if nln < 28 {
		msg := fmt.Sprintf("export.Prot2go(): rdf.Capita(%v): MalformedHeader", keys4b)
		panic(errors.New(msg))
	}
	nss := rdf.Nss // BGW URI name spaces
	/*
		Attn: no isoforms in GAF files 1
	*/
	gosubs := map[string]string{
		"gp2cc": "cellular component",
		"gp2mf": "molecular function",
		"gp2bp": "biological process",
	}
	srcU := rdf.FormU(nss["goa"])
	clsU := rdf.CompU(nss["owl"], "Class")
	// prot2bp prot2cc prot2mf graph ini
	wfh := newFH(wpth)
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	sb.WriteString(header)
	nln = 0

	stmNS := "http://rdf.biogateway.eu/prot-onto/"
	rdfNS := nss["rdf"]
	count := make(util.Set3D)

	cnt := make(util.Set2D) // count proteins absent in BGW
	cntD := 0
	for _, duoid := range duos.Keys() {
		duo := duos[duoid]
		//for duoid, duo := range duos {
		ppys := duo["ppy"].Keys()
		if l := len(ppys); l != 1 { // unnecessary, may help debugging
			msg := fmt.Sprintf("export.Prot2go():%s: Want 1 property have: %d: %v", duoid, l, ppys)
			panic(errors.New(msg))
		}
		refs := duo["ref"].Keys()
		/// Class level
		duoU := rdf.CompU(stmNS, duoid)

		bits := strings.Split(duoid, "--")
		idL := bits[0]
		idR := bits[1]
		oriL := strings.Split(idL, "!")[1] // UP AC
		oriR := strings.Split(idR, "!")[1] // GO ID
		bgwLs := upac2bgw[oriL]["bgwp"].Keys()
		if l := counter(bgwLs, cnt, "addP", "dropP", oriL); l == 0 {
			continue
		}

		sb.WriteString(rdf.FormT(duoU, ourUs["ins2cls"], clsU))
		nln++
		sb.WriteString(rdf.FormT(duoU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		oboid := strings.Replace(oriR, "_", ":", 1)
		clslbl := fmt.Sprintf("%s--%s", oriL, oboid)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
		nln++
		pdc := ppys[0]
		clsdfn := fmt.Sprintf("Association between protein %s and %s %s", oriL, gosubs[pdc], oboid)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		nln++
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++

		uriR := rdf.CompU(nss["obo"], oriR)
		// multiple subjects
		for _, bgwL := range bgwLs { // sorted above
			if len(bgwLs) > 1 {
				count.Add("oriL", oriL, bgwL)
			}
			uriL := rdf.CompU(nss["uniprot"], bgwL)
			sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
			nln++
		}
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "object"), uriR))
		nln++

		/// INSTANCES
		insid := fmt.Sprintf("%s%s%s", duoid, "#", "goa")
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], duoU))
		nln++
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		nln++
		for _, ref := range util.X1type(refs, "pubmed", "!") {
			myU := rdf.CompU(nss["pubmed"], ref)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], myU))
			nln++
		}
		for _, eco := range duo["eco"].Keys() { // for GPA files
			myU := rdf.CompU(nss["obo"], eco)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], myU))
			nln++
		}
		for _, goc := range duo["goc"].Keys() { // for GAF files
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], rdf.FormL(goc)))
			nln++
		}
		cntD++
		wfh.Write([]byte(sb.String()))
		sb.Reset()
	} // duoid
	msg := ""
	if cntD == 0 {
		msg = fmt.Sprintf("export.Prot2go(): NoDuos")
		panic(errors.New(msg))
	}
	msg = fmt.Sprintf("export.Prot2go(): Pairs: added: %d dropped: %d", cntD, len(duos)-cntD)
	log.Println(msg)
	msg = fmt.Sprintf("export.Prot2go(): Prots: added: %d dropped: %d", len(cnt["addP"]), len(cnt["dropP"]))
	log.Println(msg)
	return nln, nil
} // Prot2go

// Arg1: orthology data for one pair of taxa, output of parse.OrthoDuo(), non empty
// Arg2: path for writing RDF file
func Ortho(duos util.Set3D, wpth string) (int, error) {
	// TODO defer writing the header until the end of the main loop (duoid)
	// duos: output of parse.OrthoDuos()
	// Note: no UP isoforms in this graph; only RefProt canonical accessions
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"sub2cls",
		"stm2sbj",
		"stm2obj",
		"stm2pdc",
		"orl2orl",
		"ins2cls",
		"sth2src",
		"sub2set",
		"mbr2lst",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	header, nln := rdf.Capita(keys4b)
	if nln < 24 {
		msg := fmt.Sprintf("export.Ortho(): rdf.Capita(%v): MalformedHeader", keys4b)
		panic(errors.New(msg))
	}
	nss := rdf.Nss // BGW URI name spaces
	var srcs = map[string]string{
		"uniprot":   "http://uniprot.org/uniprot",
		"keggortho": "http://identifiers.org/kegg.orthology",
		"orthodb":   "https://www.orthodb.org",
	}
	clsU := rdf.CompU(nss["owl"], "Class")
	// ortho graph ini
	wfh := newFH(wpth)
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	sb.WriteString(header)
	///////////////////////////////////////////////////////////////////////////////
	stmNS := "http://rdf.biogateway.eu/ortho/"
	rdfNS := nss["rdf"]
	idmkeys := bgw.Orthokeys // currently only "OrthoDB": "orthodb", TODO move here?
	cntD := 0                // number of orthology relations for a pair of taxa
	nln = 0                  // number of lines written for a pair of taxa
	for _, duoid := range duos.Keys() {
		duo := duos[duoid]
		duoU := rdf.CompU(stmNS, duoid)
		bits := strings.Split(duoid, "--")
		oriL := strings.Split(bits[0], "!")[1] // UniProt Canonical Accession
		oriR := strings.Split(bits[1], "!")[1] // UniProt Canonical Accession
		sb.WriteString(rdf.FormT(duoU, ourUs["ins2cls"], clsU))
		sb.WriteString(rdf.FormT(duoU, ourUs["sub2cls"], ourUs["stm"]))
		clslbl := fmt.Sprintf("%s--%s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
		clsdfn := fmt.Sprintf("Pair of orthologous proteins %s and %s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		pdc := "orl2orl"
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		uriL := rdf.CompU(nss["uniprot"], oriL)
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "subject"), uriL))
		uriR := rdf.CompU(nss["uniprot"], oriR)
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "object"), uriR))
		sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
		sb.WriteString(rdf.FormT(uriR, ourUs[pdc], uriL))

		/// INSTANCES
		for _, idmk := range duo.Keys() {
			srck, ok := idmkeys[idmk]
			if !ok {
				continue
			} // needed! ?
			insid := fmt.Sprintf("%s%s%s", duoid, "#", srck)
			insU := rdf.CompU(stmNS, insid)
			sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], duoU))
			sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
			srcU := rdf.FormU(srcs[srck])
			sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
			// looping over orthology clusters:
			for _, setid := range duo[idmk].Keys() {
				setU := rdf.FormU(fmt.Sprintf("%s%s", nss[srck], setid))
				sb.WriteString(rdf.FormT(insU, ourUs["sub2set"], setU)) // part_of
				sb.WriteString(rdf.FormT(uriL, ourUs["mbr2lst"], setU))
				sb.WriteString(rdf.FormT(uriR, ourUs["mbr2lst"], setU))
			}
		}
		bytes := []byte(sb.String())
		if len(bytes) != 0 {
			wfh.Write(bytes)
			sb.Reset()
			cntD++
		}
	} // duoid
	return cntD, nil
} // Ortho

func Gene(rpthI, wpth string, p *bgw.Xmap) error {
	xmap := *p
	// Reference proteome IDs:
	// Idmap returns errors if fails to open the file OR the output map is empty
	xrf2upac, err := parse.Idmap(rpthI, bgw.Upkeys, 1, 2, 0) // RefProts only
	if err != nil {
		msg := fmt.Sprintf("export.Gene():  xrf2upac: %s", err)
		return errors.New(msg)
	}
	upac2xrf, err := parse.Idmap(rpthI, bgw.Upkeys, 0, 1, 2) // RefProt only
	if err != nil {
		msg := fmt.Sprintf("export.Gene(): upac2xrf: %s", err)
		return errors.New(msg)
	}
	// the 4 maps below are RefProt limited
	// building the maps:
	upca2upac := make(util.Set2D) // all accessions including iso-forms
	gnm2upca := make(util.Set2D)  // gene names from UP idmapping
	gnm2upid := make(util.Set2D)  // gene names from UP idmapping
	gnm2gsnm := make(util.Set2D)  // gene synonyms from UP idmapping
	for upac, xrfs := range upac2xrf {
		bits := strings.Split(upac, "-")
		upca := bits[0]
		upca2upac.Add(upca, upac)
		upid := ""
		// upids are NOT attached to usoforms
		if len(xrfs["UniProtKB-ID"]) == 1 {
			upid = xrfs["UniProtKB-ID"].Keys()[0]
		}
		// gene names & synonyms are NOT attached to usoforms
		for _, nmG := range xrfs["Gene_Name"].Keys() {
			gnm2upca.Add(nmG, upca)
			gnm2upid.Add(nmG, upid)
			for _, synG := range xrfs["Gene_Synonym"].Keys() {
				gnm2gsnm.Add(nmG, synG)
			}
		}
	}

	// txids are NOT attached to usoforms
	txid := xrf2upac["NCBI_TaxID"].Keys()[0]
	nss := rdf.Nss                         // BGW URI name spaces
	txnU := rdf.CompU(nss["ncbitx"], txid) // taxon URI
	xkeys := []string{
		"Ensembl",
		"EnsemblGenome",
		"GeneID",
	} // for selecting xrefs
	keys4g := make(util.SliceSet)
	// keys of object properties ('gene' graph)
	keys4g["Opys"] = []string{
		"gn2gp",
		"be2txn",
		"ins2cls",
		// "mbr2lst",
		"sth2clm",
		"sth2src",
		"sub2cls",
	}
	// keys of annotation properties ('gene' graph)
	keys4g["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"sth2syn",
	}
	// keys of parental classes ('gene' graph)
	keys4g["Prns"] = []string{
		"gn",
	}
	clsU := rdf.CompU(nss["owl"], "Class")
	/////////////////////////////////////////////////////////////////////////////
	wfh := newFH(wpth)
	defer wfh.Close()
	// 'gene' graph ini
	var sbG strings.Builder
	gnUs := rdf.FmtURIs(keys4g) // URIs for 'gene' graph
	header, nln := rdf.Capita(keys4g)
	if nln < 20 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}
	// main loop
	for _, nmG := range gnm2upca.Keys() {
		// Note: gene names and synonyms are associated only with canonical accessions
		upcas := gnm2upca[nmG].Keys() // for one gene, Refprot only, multiple possible
		if len(upcas) == 0 {
			continue // there must be at least one upca for each gene name
		}
		upacSet := make(util.Set1D) // all types of accessions, UNIQUE
		for _, upca := range upcas {
			for _, upac := range upca2upac[upca].Keys() {
				upacSet.Add(upac)
			}
		} // all accessions from RefProt
		upacs := upacSet.Keys()
		if len(upacs) == 0 {
			continue // nust contain at least one upca
		}
		// collecting selected xrefs for one gene:
		xrfs4oneg := make(util.Set2D) // needed, used for synonyms & instances
		for _, upac := range upacs {
			xrfs := upac2xrf[upac]
			// xkey as used in UniProt idmapping files
			for _, xkey := range xkeys { // slice of strings, local
				for _, xrf := range xrfs[xkey].Keys() {
					xrfs4oneg.Add(xkey, xrf)
				}
			}
		}
		upids := gnm2upid[nmG].Keys() // used only in gene defs
		if len(upids) == 0 {
			continue
		}
		lblG := strings.Replace(nmG, " ", "_", -1)
		clsG := strings.Join([]string{txid, lblG}, "/")
		clsGU := rdf.CompU(nss["gene"], clsG)
		sbG.WriteString(rdf.FormT(clsGU, gnUs["ins2cls"], clsU))
		sbG.WriteString(rdf.FormT(clsGU, gnUs["sub2cls"], gnUs["gn"]))
		sbG.WriteString(rdf.FormT(clsGU, gnUs["be2txn"], txnU))
		xmap.Lblg.Add(lblG, "bgwg", clsG)
		xmap.Bgwg.Add(clsG, "lblg", lblG)
		for _, upca := range upcas {
			clsP := upca
			clsPU := rdf.CompU(nss["uniprot"], clsP)
			clsGU := rdf.CompU(nss["gene"], clsG)
			sbG.WriteString(rdf.FormT(clsGU, gnUs["gn2gp"], clsPU))
			xmap.Bgwp.Add(clsP, "bgwg", clsG)
			xmap.Bgwg.Add(clsG, "bgwp", clsP)
		}
		dfnG := fmt.Sprintf("gene %s/%s encoding %s", txid, lblG, upids)
		sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2dfn"], rdf.FormL(dfnG)))
		sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2lbl"], rdf.FormL(nmG)))
		for _, synG := range gnm2gsnm[nmG].Keys() {
			sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2syn"], rdf.FormL(synG)))
			xmap.Syng.Add(synG, "bgwg", clsG)
			xmap.Bgwg.Add(clsG, "syng", synG)
		}
		// synonyms & instances
		for _, xkey := range xrfs4oneg.Keys() {
			var xrfU string
			for _, xrf := range xrfs4oneg[xkey].Keys() {
				if xkey == "EnsemblGenome" {
					xrfU = rdf.CompU(nss[bgw.Ensomes[txid]], xrf)
				} else {
					xrfU = rdf.CompU(nss[bgw.Upkeys[xkey]], xrf)
				}
				sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2clm"], xrfU))
				insG := fmt.Sprintf("%s%s%s", clsG, "#", xrf)
				insGU := rdf.CompU(nss["gene"], insG)
				sbG.WriteString(rdf.FormT(insGU, gnUs["ins2cls"], clsGU))
				lbl := fmt.Sprintf("%s#%s", xkey, xrf)
				sbG.WriteString(rdf.FormT(insGU, gnUs["sth2lbl"], rdf.FormL(lbl)))
				// TODO generalise
				for _, upac := range upacs {
					// including iso-forms
					upca := strings.Split(upac, "-")[0]
					if xkey == "GeneID" {
						xmap.Ncbig.Add(xrf, "bgwp", upca)
						xmap.Ncbig.Add(xrf, "bgwg", clsG)
						xmap.Bgwg.Add(clsG, "ncbig", xrf)
					}
					if xkey == "Ensembl" {
						xmap.Ensg.Add(xrf, "bgwp", upca)
						xmap.Ensg.Add(xrf, "bgwg", clsG)
						xmap.Bgwg.Add(clsG, "ensg", xrf)
					}
					if xkey == "EnsemblGenome" {
						xmap.Ensg.Add(xrf, "bgwp", upca)
						xmap.Ensg.Add(xrf, "bgwg", clsG)
						xmap.Bgwg.Add(clsG, "ensg", xrf)
					}
				}
			}
		}
	} // lblG
	outG := sbG.String()
	if len(outG) == 0 {
		msg := fmt.Sprintf("%s: CalledBy: %s: NoDataForTaxon: %s", util.FN(0), util.FN(1), txid)
		return errors.New(msg)
	}
	wfh.Write([]byte(header))
	wfh.Write([]byte(outG))
	sbG.Reset()
	return nil
} // Gene()

func Prot(rpthUP, rpthI, wpth string, p *bgw.Xmap) error {
	xmap := *p
	// Reference proteome IDs:
	// Idmap returns errors if fails to open the file OR the output map is empty
	xrf2upac, err := parse.Idmap(rpthI, bgw.Upkeys, 1, 2, 0) // Set3D
	if err != nil {
		msg := fmt.Sprintf("export.Prot():  xrf2upac: %s", err)
		return errors.New(msg)
	}
	upac2xrf, err := parse.Idmap(rpthI, bgw.Upkeys, 0, 1, 2) // Set3D
	if err != nil {
		msg := fmt.Sprintf("export.Prot():  upac2xrf: %s", err)
		return errors.New(msg)
	}
	upca2upac := make(util.Set2D)
	// keys: canonical accessions
	// vals: all accessions including canonical and iso-forms
	for _, upac := range upac2xrf.Keys() {
		bits := strings.Split(upac, "-")
		upca := bits[0]
		upca2upac.Add(upca, upac)
	}

	nss := rdf.Nss // BGW URI name spaces
	txid := xrf2upac["NCBI_TaxID"].Keys()[0]

	keys4p := make(util.SliceSet)
	// keys of object properties ('prot' graph)
	keys4p["Opys"] = []string{
		"be2txn",
		"ins2cls",
		"sth2clm",
		"sth2evd",
		// "sth2src", // special case
		"sub2cls",
	}
	// keys of annotation properties ('prot' graph)
	keys4p["Apys"] = []string{
		"evd2lvl",
		"sth2dfn",
		"sth2lbl",
		"sth2syn",
	}
	// keys of parental classes ('prot' graph)
	keys4p["Prns"] = []string{
		"tlp",
	}
	clsU := rdf.CompU(nss["owl"], "Class")
	/////////////////////////////////////////////////////////////////////////////
	wfh := newFH(wpth)
	defer wfh.Close()
	// prot graph ini
	var sbP strings.Builder
	gpUs := rdf.FmtURIs(keys4p) // URIs for 'prot' graph
	header, nln := rdf.Capita(keys4p)
	if nln < 20 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}
	/////////////////////////////////////////////////////////////////////////////

	allUPs, err := parse.Tab2set3D(rpthUP, bgw.UpdatConf.Keys, bgw.UpdatConf.Vals)
	if err != nil {
		msg := fmt.Sprintf("export.Prot(): parse.Tab2set3D(%s, _, _): allUPs: EmptySet", rpthUP)
		return errors.New(msg)
	}
	txnU := rdf.CompU(nss["ncbitx"], txid) // taxon URI
	xkeys := []string{
		"Ensembl_PRO",
		"RefSeq",
	}
	// main loop
	for _, upca := range upca2upac.Keys() {
		dfnP := ""                      // descriptor
		upacs := upca2upac[upca].Keys() // including upcas
		// filtering xrefs by xkeys:
		xrfs4onep := make(util.Set2D)
		for _, upac := range upacs {
			xrfs := upac2xrf[upac] // Set2D
			for _, xkey := range xkeys {
				for _, xrf := range xrfs[xkey].Keys() {
					xrfs4onep.Add(xkey, xrf)
				}
			}
		}
		uparc2upac := make(util.Set2D) // used for instances
		for _, upac := range upacs {
			uparc2upac.Add(upac2xrf[upac]["UniParc"].Keys()[0], upac)
		}

		clsP := upca
		/// Some sanity checks
		oneP, ok := allUPs[upca] // UniProt entry
		if !ok {
			continue
		}
		if len(oneP) < 5 {
			msg := fmt.Sprintf("export.Prot():%s:%s: TooFewDataFields: %d", txid, upca, len(oneP))
			fmt.Printf("%s\n", msg)
			continue
		} // all filds present and have at least one entry
		oneXs, ok := upac2xrf[upca] // xrefs for the canonical accession only
		if !ok {
			// never happens
			continue
		}
		if len(oneXs) < 3 {
			msg := fmt.Sprintf("export.Prot():%s:%s: TooFewDataFields: %d", txid, upca, len(oneXs))
			fmt.Printf("%s\n", msg)
			continue
		} // all filds present and have at least one entry
		// not all have pubmed refs, thus 5 iso 6
		upids := oneP["upid"].Keys()
		if len(upids) != 1 {
			continue
		}
		lblP := upids[0]
		pdfns := oneP["pdfns"].Keys()
		if len(pdfns) == 1 {
			dfnP = strings.TrimSpace(pdfns[0])
		} else {
			//msg := fmt.Sprintf("export.Prot():%s:%s: NoDefinition", txid, upca)
			//fmt.Printf("%s\n", msg)
			dfnP = "Unspecified" // 90 in 25 species
		} // 20200531: 124
		oriU := rdf.CompU(nss["uniprot"], upca)
		sbP.WriteString(rdf.FormT(oriU, gpUs["ins2cls"], clsU))
		sbP.WriteString(rdf.FormT(oriU, gpUs["sub2cls"], gpUs["tlp"]))
		sbP.WriteString(rdf.FormT(oriU, gpUs["be2txn"], txnU))
		sbP.WriteString(rdf.FormT(oriU, gpUs["sth2dfn"], rdf.FormL(dfnP)))
		sbP.WriteString(rdf.FormT(oriU, gpUs["sth2lbl"], rdf.FormL(lblP)))
		xmap.Lblp.Add(lblP, "bgwp", clsP)
		xmap.Bgwp.Add(clsP, "lblp", lblP)

		// synonyms
		synPs := make(util.Set2D)
		for _, val := range oneXs["Gene_Name"].Keys() {
			synPs.Add(val, upca)
		}
		for _, val := range oneXs["Gene_Synonym"].Keys() {
			synPs.Add(val, upca)
		}
		synPs.Add(upca, upca)
		for _, synP := range synPs.Keys() {
			sbP.WriteString(rdf.FormT(oriU, gpUs["sth2syn"], rdf.FormL(synP)))
			xmap.Synp.Add(synP, "bgwp", clsP)
			xmap.Bgwp.Add(clsP, "synp", synP)
		}

		score := oneP["score"].Keys()[0]
		sbP.WriteString(rdf.FormT(oriU, gpUs["evd2lvl"], rdf.FormL(string(score))))
		pubmeds := oneP["pubmed"].Keys()
		for _, pubmed := range pubmeds { // sorted above
			pubmedU := rdf.CompU(nss["pubmed"], pubmed)
			sbP.WriteString(rdf.FormT(oriU, gpUs["sth2evd"], pubmedU))
		}
		xmap.Upac.Add(upca, "bgwp", clsP) // for backward compatibility

		// Instances
		for _, uparc := range uparc2upac.Keys() {
			insU := rdf.CompU(nss["uniparc"], uparc)
			sbP.WriteString(rdf.FormT(insU, gpUs["ins2cls"], oriU))
			lbl := fmt.Sprintf("instance of %s", lblP)
			sbP.WriteString(rdf.FormT(insU, gpUs["sth2lbl"], rdf.FormL(lbl)))
			xmap.Bgwp.Add(clsP, "uparc", uparc) // full set
			for _, upac := range uparc2upac[uparc].Keys() {
				sbP.WriteString(rdf.FormT(insU, gpUs["sth2syn"], rdf.FormL(upac)))
				xmap.Upac.Add(upca, "upac", upac) // full set
				xmap.Upac.Add(upac, "uparc", uparc)
			}
		}
		for _, xkey := range xrfs4onep.Keys() {
			for _, xrf := range xrfs4onep[xkey].Keys() {
				xrfU := rdf.CompU(nss[bgw.Upkeys[xkey]], xrf)
				sbP.WriteString(rdf.FormT(oriU, gpUs["sth2clm"], xrfU))
				// TODO generalise
				if xkey == "RefSeq" {
					xmap.Rfsq.Add(xrf, "bgwp", clsP)
					xmap.Bgwp.Add(clsP, "rfsq", xrf)
				}
				if xkey == "Ensembl_PRO" {
					xmap.Ensp.Add(xrf, "bgwp", clsP)
					xmap.Bgwp.Add(clsP, "ensp", xrf)
				}
			}
		}
	} // upca
	/////////////////////////////////////////////////////////////////////////////
	outP := sbP.String()
	if len(outP) == 0 {
		msg := fmt.Sprintf("%s: CalledBy: %s: NoDataForTaxon: %s", util.FN(0), util.FN(1), txid)
		return errors.New(msg)
	}
	wfh.Write([]byte(header))
	wfh.Write([]byte(outP))
	sbP.Reset()
	return nil
} // Prot()
