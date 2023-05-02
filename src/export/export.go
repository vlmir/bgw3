package export

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

func counter(s []string, c util.Set2D, a, d, v string) (l int) {
	l = len(s)
	if l == 0 {
		c.Add(d, v)
	} else {
		c.Add(a, v)
	}
	return l
}

func gene2bgwg(genes []string, xmap bgw.Xmap) []string {
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
	}
	return set.Keys()
}

func upac2bgwp(upacs []string, xmap bgw.Xmap) []string {
	set := make(util.Set2D)
	for _, upac := range upacs {
		upca := strings.Split(upac, "-")[0]
		for _, bgwp := range xmap.Upac[upca]["bgwp"].Keys() {
			set.Add(bgwp, upac)
		}
	}
	return set.Keys()
}

func Tfac2gene(d *bgw.Dat4bridge, x *bgw.Xmap, wdir string) error {
	d4b := *d
	xmap := *x
	src := d4b.Src
	txid := d4b.Taxid
	cnts := d4b.Cnts // Set2D
	p2t := "reg2ptrg"
	n2t := "reg2ntrg"
	u2t := "reg2utrg"
	modes := make(util.Set3D)
	modes.Add("signor", "UP", p2t)
	modes.Add("signor", "DOWN", n2t)
	modes.Add("cytreg", "Activation", p2t)
	modes.Add("cytreg", "Repression", n2t)
	modes.Add("trrust", "Activation", p2t)
	modes.Add("trrust", "Repression", n2t)
	modes.Add("geredb", "positive", p2t)
	modes.Add("geredb", "negative", n2t)
	modes.Add("tfacts", "UP", p2t)
	modes.Add("tfacts", "DOWN", n2t)
	modes.Add("ntnu", "+", p2t)
	modes.Add("ntnu", "-", n2t)
	modes.Add("tflink", "activator", p2t)
	modes.Add("tflink", "repressor", n2t)
	sdir := "tfac2gene"
	wpths := map[string]string{
		p2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "p", src, txid),
		n2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "n", src, txid),
		u2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "u", src, txid),
	}

	fh4p, err := os.Create(wpths[p2t])
	util.CheckE(err)
	defer fh4p.Close()
	fh4n, err := os.Create(wpths[n2t])
	util.CheckE(err)
	defer fh4n.Close()
	fh4u, err := os.Create(wpths[u2t])
	util.CheckE(err)
	defer fh4u.Close()
	fhs := make(map[string]*os.File)
	fhs[p2t] = fh4p
	fhs[n2t] = fh4n
	fhs[u2t] = fh4u
	keys4b := bgw.Keys4rgrs()
	ourUs := rdf.FmtURIs(keys4b)
	nss := rdf.Nss // BGW URI name spaces
	//srcns := nss[src] // fully qualified namespace
	srcU := rdf.FormU(rdf.Uris4tftg[src])
	rdfns := nss["rdf"]
	// graphns := fmt.Sprintf("%s%s", nss["bgw"], "graph/")
	header, nln := rdf.Capita(keys4b)
	if nln < 24 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}

	duos := d4b.Duos
	entAns := nss["uniprot"] // sic, never changes
	entBns := nss["gene"]    // sic, never changes
	flags := make(util.Set1D) // for printing the header only once per file
	for _, duokey := range duos.Keys() {
		duo := duos[duokey]
		oriAB := strings.Split(duokey, "--")
		oriBid := oriAB[1] // gene symbol

		var allAids []string // UP canonical ACs for one duokey
		var allBids []string // gene symbols for one duokey
		allAids = duo["uniprot"].Keys()
		allBids = duo["ncbig"].Keys() // NCBI GeneID, single
		allBids = append(allBids, oriBid)
		// now allAids are all UP AC

		// converting allAids to Bgwids
		bgwAids := upac2bgwp(allAids, xmap)
		// converting allBids to Bgwids; BGW genes
		bgwBids := gene2bgwg(allBids, xmap)
		if len(bgwAids) == 0 || len(bgwBids) == 0 {
			continue
		}

		//  assigning predicate depending of  the direction of regulation
		pdcks := []string{u2t} // all predicate keys for one duokey, u2t first
		for _, onemod := range duo["mode"].Keys() {
			// keys for predicates
			for _, pdck := range modes[src][onemod].Keys() {
				pdcks = append(pdcks, pdck)
			}
		} // onemod
		for _, pdck := range pdcks {
			clsns := fmt.Sprintf("%s%s/", nss["bgw"], pdck)
			var sb strings.Builder
			//sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
			for _, entAid := range bgwAids {
				for _, entBid := range bgwBids {
					cnts.Add(pdck, src)
					clsid := fmt.Sprintf("uniprot!%s--gene!%s", entAid, entBid)
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
					insU := fmt.Sprintf("%s%s%s", clsU, "#", src)
					sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
					id := "MI:2247" // transcriptional regulation
					id = strings.Replace(id, ":", "_", 1)
					sb.WriteString(rdf.FormT(insU, ourUs["gp2bp"], rdf.CompU(nss["obo"], id)))
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
						bits := strings.Split(key, ":")
						if len(bits) == 2 {
							if bits[1] == "Yes" {
								key = "High"
							} else {
								key = "Low"
							}
						}
						sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(key)))
					}
				} // entBid
			} // entAid
			if cnts[pdck][src] > 0 {
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
	log.Println("export.Tfac2gene():", src, txid, cnts)
	return nil
} // Tfac2gene

func SigPways(d *bgw.Dat4bridge, x *bgw.Xmap, wdir string) error {
	d4b := *d
	xmap := *x
	src := d4b.Src
	txid := d4b.Taxid
	cnts := d4b.Cnts // Set2D
	p2t := "reg2ptrg"
	n2t := "reg2ntrg"
	u2t := "reg2utrg"
	modes := make(util.Set3D)
	modes.Add("signor", "up-regulates", p2t)
	modes.Add("signor", "down-regulates", n2t)
	sdir := "sigpways"
	wpths := map[string]string{
		p2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "p", src, txid),
		n2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "n", src, txid),
		u2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "u", src, txid),
	}

	fh4p, err := os.Create(wpths[p2t])
	util.CheckE(err)
	defer fh4p.Close()
	fh4n, err := os.Create(wpths[n2t])
	util.CheckE(err)
	defer fh4n.Close()
	fh4u, err := os.Create(wpths[u2t])
	util.CheckE(err)
	defer fh4u.Close()
	fhs := make(map[string]*os.File)
	fhs[p2t] = fh4p
	fhs[n2t] = fh4n
	fhs[u2t] = fh4u
	keys4b := bgw.Keys4rgrs()
	ourUs := rdf.FmtURIs(keys4b)
	nss := rdf.Nss // BGW URI name spaces
	//srcns := nss[src] // fully qualified namespace
	srcU := rdf.FormU(rdf.Uris4tftg[src])
	rdfns := nss["rdf"]
	// graphns := fmt.Sprintf("%s%s", nss["bgw"], "graph/")
	header, nln := rdf.Capita(keys4b)
	if nln < 24 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}

	duos := d4b.Duos
	mytypes := []string{"complex", "protein", "proteinfamily"} // for filetering
	entAns := nss["uniprot"]                                   // sic, never changes
	entBns := nss["uniprot"]                                   // sic, never changes
	flags := make(util.Set1D) // for printing the header only once per file
	for _, duokey := range duos.Keys() {
		duo := duos[duokey]
		// if duo["txid"].Keys()[0] != txid {continue} // filtering by host
		oriAid := duo["Aid"].Keys()[0] // assuming a single value
		oriBid := duo["Bid"].Keys()[0] // assuming a single value

		var allAids []string // UP ACs for one duokey
		var allBids []string // UP ACs for one duokey
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
			fmt.Printf("%s%s %s %s", "export.SigPways(): multiple entity A types, skipping: ", src, duokey, typeAs)
			continue
		}
		if len(typeBs) > 1 {
			fmt.Printf("%s%s %s %s", "export.SigPways(): multiple entity B types, skipping: ", src, duokey, typeBs)
			continue
		}
		typeA := typeAs[0]
		typeB := typeBs[0]
		/// expanding families and complexes
		if typeA == "complex" || typeA == "proteinfamily" {
			allAids = xmap.Signor[oriAid]["ids"].Keys() // all UP accessions
			if len(allAids) == 0 {
				continue
			}
		}
		if typeB == "complex" || typeB == "proteinfamily" {
			allBids = xmap.Signor[oriBid]["ids"].Keys() // all UP accessions
			if len(allBids) == 0 {
				continue
			}
		}
		if typeA == "protein" {
			allAids = append(allAids, oriAid)
		}
		if typeB == "protein" {
			allBids = append(allBids, oriBid)
		}
		// now allAids allBids are all UP AC

		// converting allAids to Bgwids
		bgwAids := upac2bgwp(allAids, xmap)
		// converting allBids to Bgwids
		bgwBids := upac2bgwp(allBids, xmap)
		if len(bgwAids) == 0 || len(bgwBids) == 0 {
			continue
		}

		//  assigning predicate depending of  the direction of regulation
		pdcks := []string{u2t} // all predicate keys for one duokey, u2t first
		for _, onemod := range duo["modelbl"].Keys() {
			onemod = strings.Split(onemod, " ")[0]
			// keys for predicates
			for _, pdck := range modes[src][onemod].Keys() {
				pdcks = append(pdcks, pdck)
			}
		} // onemod

		for _, pdck := range pdcks {
			clsns := fmt.Sprintf("%s%s/", nss["bgw"], pdck)
			var sb strings.Builder
			//sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
			for _, entAid := range bgwAids {
				for _, entBid := range bgwBids {
					cnts.Add(pdck, src)
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
					insU := fmt.Sprintf("%s%s%s", clsU, "#", src)
					sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
					for _, id := range duo["cellid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["pcs2loc"], rdf.CompU(nss["obo"], id)))
					}
					for _, id := range duo["tissueid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["pcs2loc"], rdf.CompU(nss["obo"], id)))
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
			if cnts[pdck][src] > 0 {
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
	log.Println("export.SigPways():", src, txid, cnts)
	return nil
} // SigPways

func Rgr2trg(d *bgw.Dat4bridge, x *bgw.Xmap, wdir string) error {
	d4b := *d
	xmap := *x
	src := d4b.Src
	txid := d4b.Taxid
	cnts := d4b.Cnts // Set2D
	p2t := "reg2ptrg"
	n2t := "reg2ntrg"
	u2t := "reg2utrg"
	modes := make(util.Set3D)
	modes.Add("signor", "up-regulates", p2t)
	modes.Add("signor", "down-regulates", n2t)
	sdir := "reg2targ"
	wpths := map[string]string{
		p2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "p", src, txid),
		n2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "n", src, txid),
		u2t: fmt.Sprintf("%s%s/%s-%s-%s.nt", wdir, sdir, "u", src, txid),
	}

	fh4p, err := os.Create(wpths[p2t])
	util.CheckE(err)
	defer fh4p.Close()
	fh4n, err := os.Create(wpths[n2t])
	util.CheckE(err)
	defer fh4n.Close()
	fh4u, err := os.Create(wpths[u2t])
	util.CheckE(err)
	defer fh4u.Close()
	fhs := make(map[string]*os.File)
	fhs[p2t] = fh4p
	fhs[n2t] = fh4n
	fhs[u2t] = fh4u
	keys4b := bgw.Keys4rgrs()
	ourUs := rdf.FmtURIs(keys4b)
	nss := rdf.Nss // BGW URI name spaces
	//srcns := nss[src] // fully qualified namespace
	srcU := rdf.FormU(rdf.Uris4tftg[src])
	rdfns := nss["rdf"]
	// graphns := fmt.Sprintf("%s%s", nss["bgw"], "graph/")
	header, nln := rdf.Capita(keys4b)
	if nln < 24 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}

	duos := d4b.Duos
	// future extensions:
	// psi-mi:"MI:0328"(small molecule)
	// psi-mi:"MI:2258"(xenobiotic) # chemical?
	mytypes := []string{"MI:0314", "MI:1304", "MI:0326", "MI:0328"} // for filetering
	myRegTypes := []string{}                                        // for filetering
	entAns := nss["uniprot"]                                        // sic, never changes
	entBns := nss["uniprot"]                                        // sic, never changes
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
		var allAids []string // UP ACs for one duokey
		var allBids []string // UP ACs for one duokey
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
		// skipping pairs with ambiguous entyty types - never happens
		if len(typeAs) > 1 {
			fmt.Printf("%s%s %s %s", "export.Rgr2trg(): multiple entity A types, skipping: ", src, duokey, typeAs)
			continue
		}
		if len(typeBs) > 1 {
			fmt.Printf("%s%s %s %s", "export.Rgr2trg(): multiple entity B types, skipping: ", src, duokey, typeBs)
			continue
		}
		typeA := typeAs[0] // assuming a single vwlue
		typeB := typeBs[0] // assuming a single value
		/// expanding families and complexes
		if typeA == "MI:0314" || typeA == "MI:1304" {
			allAids = xmap.Signor[oriAid]["ids"].Keys() // all UP accessions
			if len(allAids) == 0 {
				continue
			}
		}
		if typeB == "MI:0314" || typeB == "MI:1304" {
			allBids = xmap.Signor[oriBid]["ids"].Keys()
			if len(allBids) == 0 {
				continue
			}
		}
		// for chemicals? TODO confirm
		if len(allAids) == 0 {
			allAids = append(allAids, oriAid)
		}
		if len(allBids) == 0 {
			allBids = append(allBids, oriBid)
		}

		// converting allAids to Bgwids
		bgwAids := upac2bgwp(allAids, xmap)
		// converting allBids to Bgwids
		bgwBids := upac2bgwp(allBids, xmap)
		if len(bgwAids) == 0 || len(bgwBids) == 0 {
			continue
		}

		//  assigning predicate depending of  the direction of regulation
		pdcks := []string{u2t} // all predicate keys for one duokey, u2t first
		for _, onemod := range duo["modelbl"].Keys() {
			onemod = util.StripParQuots(onemod)
			onemod = strings.Split(onemod, " ")[0]
			// keys for predicates
			for _, pdck := range modes[src][onemod].Keys() {
				pdcks = append(pdcks, pdck)
			}
		} // onemod

		for _, pdck := range pdcks {
			clsns := fmt.Sprintf("%s%s/", nss["bgw"], pdck)
			var sb strings.Builder
			//sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
			for _, entAid := range bgwAids {
				for _, entBid := range bgwBids {
					cnts.Add(pdck, src)
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
					insU := fmt.Sprintf("%s%s%s", clsU, "#", src)
					// inslbl := fmt.Sprintf("%s%s%s", clsid, " from ", src) // changed
					sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clsid)))
					sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
					for _, id := range duo["reglevelid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["gp2bp"], rdf.CompU(nss["obo"], id)))
					}
					for _, id := range duo["typeABid"].Keys() {
						id = strings.Replace(id, ":", "_", 1)
						sb.WriteString(rdf.FormT(insU, ourUs["sth2rlm"], rdf.CompU(nss["obo"], id)))
					}
					for _, key := range duo["pubmed"].Keys() {
						pubmedU := rdf.CompU(nss["pubmed"], key)
						sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], pubmedU))
					}
					for _, key := range duo["score"].Keys() {
						sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(key)))
					}
				} // entBid
			} // entAid
			if cnts[pdck][src] > 0 {
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
	log.Println("export.Rgr2trg():", src, txid, cnts)
	return nil
} // Rgr2trg

// func Gene(rpthI string, wpthG string, p *bgw.Xmap) error {
func Gene(rpthU, rpthI, wpthG string, p *bgw.Xmap) error {
	xmap := *p
	// Reference proteome IDs:
	xrf2upac, _ := parse.Idmap(rpthI, bgw.Upkeys, 1, 2, 0)
	upac2xrf, _ := parse.Idmap(rpthI, bgw.Upkeys, 0, 1, 2)
	upca2upac := make(util.Set2D) // TODO double-check the use of
	for _, upac := range upac2xrf.Keys() {
		bits := strings.Split(upac, "-")
		if len(bits) == 2 {
			upca := bits[0]
			upca2upac.Add(upca, upac)
		} else {
			upca2upac.Add(upac, upac)
		}
	}
	keys, vals := bgw.UpdatParseConf()
	allUPs, err := parse.Tab2set3D(rpthU, keys, vals) // keys UP canonical accessions, not all RefProt
	allUPs.Check()                                    // 8 linker keys:
	gnm2upca := make(util.Set2D)
	for _, upca := range allUPs.Keys() {
		for _, gnm := range allUPs[upca]["gnms"].Keys() {
			gnm2upca.Add(gnm, upca)
		}
	}

	txid := xrf2upac["NCBI_TaxID"].Keys()[0]
	nss := rdf.Nss                         // BGW URI name spaces
	txnU := rdf.CompU(nss["ncbitx"], txid) // taxon URI
	xkeys := []string{
		"Ensembl",
		"EnsemblGenome",
		"GeneID",
	}
	keys4g := make(util.SliceSet)
	// keys of object properties ('gene' graph)
	keys4g["Opys"] = []string{
		"gn2gp",
		"gn2txn",
		"ins2cls",
		// "mbr2lst",
		"sth2clm",
		"sth2ori",
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
	srcU := rdf.FormU(nss["uniprot"])
	graphUg := "<http://rdf.biogateway.eu/graph/gene>"
	/////////////////////////////////////////////////////////////////////////////
	wfhG, err := os.Create(wpthG)
	util.CheckE(err)
	defer wfhG.Close()
	// 'gene' graph ini
	var sbG strings.Builder
	gnUs := rdf.FmtURIs(keys4g) // URIs for 'gene' graph
	header, _ := rdf.Capita(keys4g)
	sbG.WriteString(header)
	sbG.WriteString(rdf.FormT(graphUg, gnUs["sth2src"], srcU))
	nsG := rdf.FormU(nss["gene"])
	sbG.WriteString(rdf.FormT(nsG, gnUs["ins2cls"], clsU))
	dfn := "Set of genes in Biogateway"
	sbG.WriteString(rdf.FormT(nsG, gnUs["sth2dfn"], rdf.FormL(dfn)))
	lbl := "gene"
	sbG.WriteString(rdf.FormT(nsG, gnUs["sth2lbl"], rdf.FormL(lbl)))
	// cnt := make(util.Set2D)
	for _, gnm := range gnm2upca.Keys() {
		upcas := gnm2upca[gnm].Keys() // Attn: not limited to RefProt !
		upacs4oneg := make(util.Set3D)
		for _, upca := range upcas {
			upacs := upca2upac[upca].Keys() // including upcas, RefProt limited
			if len(upacs) == 0 {
				continue
			}
			for _, upac := range upacs {
				upacs4oneg.Add("upacs", upac, upca)
				upacs4oneg.Add("upcas", upca, upac)
			}
		}
		upcas = upacs4oneg["upcas"].Keys()
		upacs := upacs4oneg["upacs"].Keys()
		if len(upcas) == 0 || len(upacs) == 0 {
			continue
		} // all accessions now from RefProt
		// canonical accessions for a gene, filtered by Reference Proteome
		// Note: gene names and synonyms are associated only with canonical accessions
		xrfs4oneg := make(util.Set2D)
		upids4oneg := make(util.Set2D)
		gsnms4oneg := make(util.Set2D)
		for _, upca := range upcas {
			oneUP, ok := allUPs[upca]
			if !ok {
				continue
			}
			for _, upid := range oneUP["upid"].Keys() {
				upids4oneg.Add(upid, upca)
			}
			if len(upids4oneg.Keys()) == 0 {
				continue
			}
			for _, gsnm := range allUPs[upca]["gsnms"].Keys() {
				if len(gsnm) < 3 {
					continue
				}
				if gsnm != gnm {
					gsnms4oneg.Add(gsnm, upca)
				}
			}
		} // upids4oneg and gsnms4oneg constructed
		if len(upids4oneg) == 0 {
			continue
		}
		upids := upids4oneg.Keys()
		lblG := strings.Replace(gnm, " ", "_", -1)
		clsG := strings.Join([]string{txid, lblG}, "/")
		clsGU := rdf.CompU(nss["gene"], clsG)
		sbG.WriteString(rdf.FormT(clsGU, gnUs["ins2cls"], clsU))
		sbG.WriteString(rdf.FormT(clsGU, gnUs["sub2cls"], gnUs["gn"]))
		sbG.WriteString(rdf.FormT(clsGU, gnUs["gn2txn"], txnU))
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
		sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2lbl"], rdf.FormL(gnm)))
		for _, synG := range gsnms4oneg.Keys() {
			sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2syn"], rdf.FormL(synG)))
			xmap.Syng.Add(synG, "bgwg", clsG)
			xmap.Bgwg.Add(clsG, "syng", synG)
		}
		for _, upac := range upacs {
			xrfs := upac2xrf[upac]
			for _, xkey := range xkeys {
				for _, xrf := range xrfs[xkey].Keys() {
					xrfs4oneg.Add(xkey, xrf)
				}
			}
		}
		// xkey as used in UniProt idmapping files
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
				if xkey == "GeneID" {
					xmap.Ncbig.Add(xrf, "bgwg", clsG)
					xmap.Bgwg.Add(clsG, "ncbig", xrf)
				}
				if xkey == "Ensembl" {
					xmap.Ensg.Add(xrf, "bgwg", clsG)
					xmap.Bgwg.Add(clsG, "ensg", xrf)
				}
				if xkey == "EnsemblGenome" {
					xmap.Ensg.Add(xrf, "bgwg", clsG)
					xmap.Bgwg.Add(clsG, "ensg", xrf)
				}
			}
		}
	} // lblG
	wfhG.Write([]byte(sbG.String()))
	sbG.Reset()
	return nil
} // Gene()

func Prot(rpthU, rpthI, wpthP string, p *bgw.Xmap) error {
	xmap := *p
	// Reference proteome IDs:
	xrf2upac, _ := parse.Idmap(rpthI, bgw.Upkeys, 1, 2, 0)
	upac2xrf, _ := parse.Idmap(rpthI, bgw.Upkeys, 0, 1, 2)
	upca2upac := make(util.Set2D) // TODO double-check the use of
	for _, upac := range upac2xrf.Keys() {
		bits := strings.Split(upac, "-")
		if len(bits) == 2 {
			upca := bits[0]
			upca2upac.Add(upca, upac)
		} else {
			upca2upac.Add(upac, upac)
		}
	}

	nss := rdf.Nss // BGW URI name spaces
	txid := xrf2upac["NCBI_TaxID"].Keys()[0]

	keys4p := make(util.SliceSet)
	// keys of object properties ('prot' graph)
	keys4p["Opys"] = []string{
		"gp2txn",
		"ins2cls",
		"sth2eqv",
		// "mbr2lst",
		"sth2clm",
		"sth2evd",
		"sth2ori",
		"sth2src",
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
	srcU := rdf.FormU(nss["uniprot"])
	graphUp := "<http://rdf.biogateway.eu/graph/prot>"
	/////////////////////////////////////////////////////////////////////////////
	wfhP, err := os.Create(wpthP)
	util.CheckE(err)
	defer wfhP.Close()
	// prot graph ini
	var sbP strings.Builder
	gpUs := rdf.FmtURIs(keys4p) // URIs for 'prot' graph
	header, _ := rdf.Capita(keys4p)
	sbP.WriteString(header)
	sbP.WriteString(rdf.FormT(graphUp, gpUs["sth2src"], srcU))
	//srcUa = "<https://www.uniprot.org/uniprot/?query=organism:" + txid + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab>"
	/////////////////////////////////////////////////////////////////////////////

	keys, vals := bgw.UpdatParseConf()
	allUPs, err := parse.Tab2set3D(rpthU, keys, vals)
	allUPs.Check()                         // 8 linker keys:
	txnU := rdf.CompU(nss["ncbitx"], txid) // taxon URI
	xkeys := []string{
		"Ensembl_PRO",
		"RefSeq",
	}
	for _, upca := range upca2upac.Keys() {
		dfnP := ""                      // descriptor
		upacs := upca2upac[upca].Keys() // including upcas
		xrfs4onep := make(util.Set2D)
		for _, upac := range upacs {
			xrfs := upac2xrf[upac]
			for _, xkey := range xkeys {
				for _, xrf := range xrfs[xkey].Keys() {
					xrfs4onep.Add(xkey, xrf)
				}
			}
		}
		upi2upac := make(util.Set2D) // used for instances
		for _, upac := range upacs {
			upi2upac.Add(upac2xrf[upac]["UniParc"].Keys()[0], upac)
		}

		clsP := upca
		/// Some sanity checks
		oneP, ok := allUPs[upca]
		if !ok {
			continue
		}
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
		if len(oneP) < 5 {
			msg := fmt.Sprintf("export.Prot():%s:%s: TooFewDataFields: %d", txid, upca, len(oneP))
			fmt.Printf("%s\n", msg)
			continue
		} // all filds present and have at least one entry
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
		sbP.WriteString(rdf.FormT(oriU, gpUs["gp2txn"], txnU))
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
		for _, upi := range upi2upac.Keys() {
			insU := rdf.CompU(nss["uniparc"], upi)
			sbP.WriteString(rdf.FormT(insU, gpUs["ins2cls"], oriU))
			lbl := fmt.Sprintf("instance of %s", lblP)
			sbP.WriteString(rdf.FormT(insU, gpUs["sth2lbl"], rdf.FormL(lbl)))
			xmap.Bgwp.Add(clsP, "upi", upi) // full set
			for _, upac := range upi2upac[upi].Keys() {
				sbP.WriteString(rdf.FormT(insU, gpUs["sth2syn"], rdf.FormL(upac)))
				xmap.Upac.Add(upca, "upac", upac) // full set
				xmap.Upac.Add(upac, "upi", upi)
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
	wfhP.Write([]byte(sbP.String()))
	sbP.Reset()
	return nil
} // Prot()

// GeneProt wrtes RDF files for genes (graph 'gene') and proteins (graph 'prot')
// Arguments:
// 1. data structure aggregating data from UniProt for a single taxon
// 2; path for writing 'gene' file
// e. path for writing 'prot' file
// 4. path for writing cross-references
// Retturns:
// number of lines written to 'gene' file
// number of lines written to 'prot' file
// error
//func GeneProt(dat4rdf bgw.Dat4rdf, wpthG, wpthP, wpthX string) (int, int, error) {
//	/*
//		Fields in allUPs:
//		Add("pome",pome)
//		Add("chr",chr)
//		Add("upid",cells[1])
//		Add("spnm",cells[2])
//		Add("txid",txid)
//		Add("dfnP",strings.Replace(cells[4],"\"",",,",1))
//		Add("pubmed",cells[6])
//		Add("score",bits[0])
//	*/
//	// generated by parse.UpTab()
//	allUPs := *dat4rdf.Udat // util.Set3D, keys: Uniprot canonical accessions
//	allUPs.Check()          // 8 linker keys:
//	// generated by parse.UpTab()
//	allTs := *dat4rdf.Txns // util.Set3D, taxonomic data
//	if len(allTs) != 1 {
//		msg := fmt.Sprintf("TaxonUndefined")
//		panic(errors.New(msg))
//	}
//	txid := allTs.Keys()[0]
//	allGs := *dat4rdf.Gnm // util.Set3D, keys: UniProt gene names
//	allGs.Check()         // 2 linkers: "upac", "gsnm"
//	// generated by parse.UpIdMap() in rdf4bgw
//	upac2xrf := *dat4rdf.Upac // util.Set3D, keys: UniProt accessions
//	upac2xrf.Check()
//
//	nss := rdf.Nss // BGW URI name spaces
//	keys4g := make(util.SliceSet)
//	// keys of object properties ('gene' graph)
//	keys4g["Opys"] = []string{
//		"gn2gp",
//		"gn2txn",
//		"ins2cls",
//		// "mbr2lst",
//		"sth2clm",
//		"sth2ori",
//		"sth2src",
//		"sub2cls",
//	}
//	// keys of annotation properties ('gene' graph)
//	keys4g["Apys"] = []string{
//		"sth2dfn",
//		"sth2lbl",
//		"sth2syn",
//	}
//	// keys of parental classes ('gene' graph)
//	keys4g["Prns"] = []string{
//		"gn",
//	}
//	keys4p := make(util.SliceSet)
//	// keys of object properties ('prot' graph)
//	keys4p["Opys"] = []string{
//		"gp2txn",
//		"ins2cls",
//		"sth2eqv",
//		// "mbr2lst",
//		"sth2clm",
//		"sth2evd",
//		"sth2ori",
//		"sth2src",
//		"sub2cls",
//	}
//	// keys of annotation properties ('prot' graph)
//	keys4p["Apys"] = []string{
//		"evd2lvl",
//		"sth2dfn",
//		"sth2lbl",
//		"sth2syn",
//	}
//	// keys of parental classes ('prot' graph)
//	keys4p["Prns"] = []string{
//		"tlp",
//	}
//	clsU := rdf.CompU(nss["owl"], "Class")
//	srcU := rdf.FormU(nss["uniprot"])
//	graphUg := "<http://rdf.biogateway.eu/graph/gene>"
//	graphUp := "<http://rdf.biogateway.eu/graph/prot>"
//	dfn := "" // descriptor
//	lbl := "" // short label
//	nlng := 0 // number of lines 'gene' graph
//	nlnp := 0 // number of lines 'prot' graph
//	/////////////////////////////////////////////////////////////////////////////
//	wfhG, err := os.Create(wpthG)
//	util.CheckE(err)
//	defer wfhG.Close()
//	wfhP, err := os.Create(wpthP)
//	util.CheckE(err)
//	defer wfhP.Close()
//	// 'gene' graph ini
//	var sbG strings.Builder
//	gnUs := rdf.FmtURIs(keys4g) // URIs for 'gene' graph
//	header, ng := rdf.Capita(keys4g)
//	sbG.WriteString(header)
//	nlng += ng
//	sbG.WriteString(rdf.FormT(graphUg, gnUs["sth2src"], srcU))
//	nlng++
//	//srcUa = "<https://www.uniprot.org/uniprot/?query=organism:" + txid + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab>"
//	//sbG.WriteString(rdf.FormT("<http://rdf.biogateway.eu/graph/gene>", gnUs["sth2ori"], srcUa)); nlng++
//	nsG := rdf.FormU(nss["gene"])
//	sbG.WriteString(rdf.FormT(nsG, gnUs["ins2cls"], clsU))
//	nlng++
//	dfn = "Set of genes in Biogateway"
//	sbG.WriteString(rdf.FormT(nsG, gnUs["sth2dfn"], rdf.FormL(dfn)))
//	nlng++
//	sbG.WriteString(rdf.FormT(nsG, gnUs["sth2lbl"], rdf.FormL("gene")))
//	nlng++
//	// prot graph ini
//	var sbP strings.Builder
//	gpUs := rdf.FmtURIs(keys4p) // URIs for 'prot' graph
//	header, np := rdf.Capita(keys4p)
//	sbP.WriteString(header)
//	nlnp += np
//	sbP.WriteString(rdf.FormT(graphUp, gpUs["sth2src"], srcU))
//	nlnp++
//	//srcUa = "<https://www.uniprot.org/uniprot/?query=organism:" + txid + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab>"
//	//sbP.WriteString(rdf.FormT("<http://rdf.biogateway.eu/graph/prot>", gpUs["sth2ori"], srcUa)); nlnp++
//	nsP := rdf.FormU(nss["uniprot"])
//	sbP.WriteString(rdf.FormT(nsP, gpUs["ins2cls"], clsU))
//	nlnp++
//	dfn = "Set of translation products in Biogateway"
//	sbP.WriteString(rdf.FormT(nsP, gpUs["sth2dfn"], rdf.FormL(dfn)))
//	nlnp++
//	sbP.WriteString(rdf.FormT(nsP, gpUs["sth2lbl"], rdf.FormL("prot")))
//	nlnp++
//	/////////////////////////////////////////////////////////////////////////////
//	upca2bgwg := make(util.Set3D)         // keys: UniProt canonical accessions
//	bgwg2bgwp := make(util.Set3D)         // keys: BGW gene IDs
//	spnm := allTs[txid]["spnm"].Keys()[0] // biological species name
//	// removing alternative names:
//	bits := strings.Split(spnm, "(")
//	spnm = strings.TrimSpace(bits[0])
//	txnU := rdf.CompU(nss["ncbitx"], txid) // taxon URI
//	txnGU := rdf.CompU(nss["gene"], txid)  // BGW gene set URI
//	sbG.WriteString(rdf.FormT(txnGU, gnUs["ins2cls"], clsU))
//	nlng++
//	sbG.WriteString(rdf.FormT(txnGU, gnUs["sub2cls"], nsG))
//	nlng++
//	// sbG.WriteString(rdf.FormT(txnGU, gnUs["gn2txn"], txnU))
//	// nlng++
//	dfn = fmt.Sprintf("Set of %s genes in Biogateway", spnm)
//	sbG.WriteString(rdf.FormT(txnGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
//	nlng++
//	lbl = strings.Join([]string{"gene", txid}, "/")
//	sbG.WriteString(rdf.FormT(txnGU, gnUs["sth2lbl"], rdf.FormL(lbl)))
//	nlng++
//	txnPU := rdf.CompU(nss["uniprot"], txid)
//	sbP.WriteString(rdf.FormT(txnPU, gpUs["ins2cls"], clsU))
//	nlnp++
//	sbP.WriteString(rdf.FormT(txnPU, gpUs["sub2cls"], nsP))
//	nlnp++
//	dfn = fmt.Sprintf("Set of %s translation products in Biogateway", spnm)
//	sbP.WriteString(rdf.FormT(txnPU, gpUs["sth2dfn"], rdf.FormL(dfn)))
//	nlnp++
//	lbl = strings.Join([]string{"prot", txid}, "/")
//	sbP.WriteString(rdf.FormT(txnPU, gpUs["sth2lbl"], rdf.FormL(lbl)))
//	nlnp++
//	/* TODO
//	if nlng < 35 || nlnp < 37 {
//		msg := fmt.Sprintf("MalformedHeader")
//		panic(errors.New(msg))
//	}
//	*/
//	nlng = 0
//	nlnp = 0
//	/////////////////////////////////////////////////////////////////////////////
//	for _, chr := range allTs[txid]["chr"].Keys() {
//		chrid := strings.Join([]string{txid, chr}, "/")
//		chrGU := rdf.CompU(nss["gene"], chrid)
//		chrPU := rdf.CompU(nss["uniprot"], chrid)
//		sbG.WriteString(rdf.FormT(chrGU, gnUs["ins2cls"], clsU))
//		nlng++
//		sbG.WriteString(rdf.FormT(chrGU, gnUs["sub2cls"], txnGU))
//		nlng++
//		dfn := fmt.Sprintf("Set of genes residing in %s %s", spnm, chr)
//		sbG.WriteString(rdf.FormT(chrGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
//		nlng++
//		lbl = strings.Join([]string{"gene", chrid}, "/")
//		sbG.WriteString(rdf.FormT(chrGU, gnUs["sth2lbl"], rdf.FormL(lbl)))
//		nlng++
//		sbP.WriteString(rdf.FormT(chrPU, gpUs["ins2cls"], clsU))
//		nlnp++
//		sbP.WriteString(rdf.FormT(chrPU, gpUs["sub2cls"], txnPU))
//		nlnp++
//		dfn = fmt.Sprintf("Set of tranlation products encoded by %s %s", spnm, chr)
//		sbP.WriteString(rdf.FormT(chrPU, gpUs["sth2dfn"], rdf.FormL(dfn)))
//		nlnp++
//		lbl = strings.Join([]string{"prot", chrid}, "/")
//		sbP.WriteString(rdf.FormT(chrPU, gpUs["sth2lbl"], rdf.FormL(lbl)))
//		nlnp++
//	}
//	/////////////////////////////////////////////////////////////////////////////
//	wfhX, err := os.Create(wpthX)
//	util.CheckE(err)
//	defer wfhX.Close()
//	var xmap bgw.Xmap
//	xmap.New()
//
//	/////////////////////////////////////////////////////////////////////////////
//	/// Genes ///
//	// allGs checked
//	// allGs: ALL GeneNams in the Reference Protome!
//	cnt := make(util.Set2D)
//	for _, lblG := range allGs.Keys() {
//		oneG, ok := allGs[lblG]
//		if !ok {
//			panic("Impossible!")
//		}
//		oneG.Check() // not empty
//		clsG := strings.Join([]string{txid, lblG}, "/")
//		clsGU := rdf.CompU(nss["gene"], clsG)
//		// canonical accessions for a gene, filtered by Reference Proteome in parse.UpTab:
//		// Note: gene names and synonyms are associated only with canonical accessions
//		upcas := oneG["upca"].Keys()
//		if len(upcas) == 0 {
//			msg := fmt.Sprintf("export.GeneProt():%s:%s: NoUpcas", txid, lblG)
//			panic(errors.New(msg))
//		}
//		for _, upca := range upcas {
//			xrfs, ok := upac2xrf[upca]
//			if !ok {
//				// never happens
//				msg := fmt.Sprintf("export.GeneProt():%s:%s: NoXrefs", txid, upca)
//				panic(errors.New(msg))
//			}
//			m := 1 // at least "uniparc" - SIC!
//			if l := len(xrfs); l < m {
//				// never happens
//				msg := fmt.Sprintf("export.GeneProt():%s:%s: NoXrefs", txid, upca)
//				panic(errors.New(msg))
//			} // 20200531: 4278 entries having only UnniParc refs
//		}
//		sbG.WriteString(rdf.FormT(clsGU, gnUs["ins2cls"], clsU))
//		nlng++
//		for _, chr := range oneG["chr"].Keys() { // SIC! already filtered
//			chrid := strings.Join([]string{txid, chr}, "/")
//			chrGU := rdf.CompU(nss["gene"], chrid)
//			sbG.WriteString(rdf.FormT(clsGU, gnUs["sub2cls"], chrGU))
//			nlng++
//		}
//		sbG.WriteString(rdf.FormT(clsGU, gnUs["sub2cls"], gnUs["gn"]))
//		nlng++
//		sbG.WriteString(rdf.FormT(clsGU, gnUs["gn2txn"], txnU))
//		nlng++
//		dfn := fmt.Sprintf("%s gene %s", spnm, lblG)
//		sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
//		nlng++
//		sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2lbl"], rdf.FormL(lblG)))
//		nlng++
//		/// gene synonyms ///
//		// includes original GeneNames for entries with multipls GeneNames for now
//		for _, synG := range oneG["gsnm"].Keys() {
//			if synG == lblG {
//				continue
//			}
//			sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2syn"], rdf.FormL(synG)))
//			nlng++
//			xmap.Syng.Add(synG, "bgwg", clsG)
//			xmap.Syng.Add(clsG, "syng", synG)
//		}
//		xmap.Lblg.Add(lblG, "bgwg", clsG)
//		xmap.Lblg.Add(clsG, "gnm", lblG)
//		xrfs4oneg := make(util.Set2D)
//		for _, upca := range upcas {
//			clsP := fmt.Sprintf("%s%s%s", txid, "/", upca)
//			clsPU := rdf.CompU(nss["uniprot"], clsP)
//			clsGU := rdf.CompU(nss["gene"], clsG)
//			sbG.WriteString(rdf.FormT(clsGU, gnUs["gn2gp"], clsPU))
//			xmap.Bgwp.Add(clsP, "bgwg", clsG)
//			xmap.Bgwg.Add(clsG, "bgwp", clsP)
//			nlng++ // NB: => gene graph
//			xrfs, ok := upac2xrf[upca]
//			ensgs := xrfs["ensgene"].Keys()
//			ensoms := xrfs["ensom"].Keys()
//			ncbigs := xrfs["ncbigene"].Keys()
//			if len(ensgs)+len(ensoms)+len(ncbigs) == 0 {
//				cnt.Add("noGx", lblG)
//			} // 20200531: 44691
//			oneP, ok := allUPs[upca]
//			// skipping accessions filtered out in parse.UpTab():
//			if !ok {
//				panic("Redundant?")
//			} // probably redundant
//			upca2bgwg.Add(upca, "bgwg", clsG)
//			oriU := rdf.CompU(nss["uniprot"], upca)
//			sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2ori"], oriU))
//			nlng++
//			/// id mapping ///
//			// multiple BGW GeneIDs per GeneName are possible, e.g.:
//			//<http://rdf.biogateway.eu/gene/9606/multi/CALM1>
//			// <http://rdf.biogateway.eu/gene/9606/chr-14/CALM1>
//			if len(oneP["gnm"]) == 1 {
//				for _, ensg := range ensgs {
//					xrfs4oneg.Add("ensgene", ensg)
//					xmap.Ensg.Add(ensg, "bgwg", clsG)
//					xmap.Bgwg.Add(clsG, "ensgene", ensg)
//				}
//				for _, ensom := range ensoms {
//					xrfs4oneg.Add(bgw.Ensomes[txid], ensom)
//					xmap.Ensg.Add(ensom, "bgwg", clsG)
//					xmap.Bgwg.Add(clsG, "ensom", ensom)
//				}
//				for _, ncbig := range ncbigs {
//					xrfs4oneg.Add("ncbigene", ncbig)
//					xmap.Ncbig.Add(ncbig, "bgwg", clsG)
//					xmap.Bgwg.Add(clsG, "ncbigene", ncbig)
//				}
//			}
//		} // upca
//		for _, src := range xrfs4oneg.Keys() {
//			for _, xrf := range xrfs4oneg[src].Keys() {
//				xrfU := rdf.CompU(nss[src], xrf)
//				sbG.WriteString(rdf.FormT(clsGU, gnUs["sth2clm"], xrfU))
//				nlng++
//				insG := fmt.Sprintf("%s%s%s", clsG, "#", xrf)
//				insGU := rdf.CompU(nss["gene"], insG)
//				sbG.WriteString(rdf.FormT(insGU, gpUs["ins2cls"], clsGU))
//				nlng++
//				lbl := fmt.Sprintf("%s#%s", src, xrf)
//				sbG.WriteString(rdf.FormT(insGU, gpUs["sth2lbl"], rdf.FormL(lbl)))
//				nlng++
//				sbG.WriteString(rdf.FormT(insGU, gnUs["sth2clm"], xrfU))
//				nlng++
//			}
//		}
//	} // gene name
//	/////////////////////////////////////////////////////////////////////////////
//	/// Proteins ///
//	for _, upca := range allUPs.Keys() {
//
//		/// Some sanity checks
//		oneP, ok := allUPs[upca]
//		if !ok {
//			panic("Impossible!")
//		}
//		//for upca, oneP := range allUPs {
//		xrfsC, ok := upac2xrf[upca]
//		if !ok {
//			// never happens
//			msg := fmt.Sprintf("export.GeneProt():%s:%s: NoXrfs", txid, upca)
//			panic(errors.New(msg))
//		}
//		// not all have synonyms, pubmed refs, thus 8 iso 10
//		if len(oneP) < 6 {
//			msg := fmt.Sprintf("export.GeneProt():%s:%s: TooFewDataFields", txid, upca)
//			fmt.Printf("%s\n", msg)
//			panic(errors.New(msg)) // TODO:3.4: see why tis fails
//		} // all filds present and have at least one entry
//		lblGs := oneP["gnm"].Keys()
//		l := len(lblGs)
//		if l == 0 {
//			msg := fmt.Sprintf("export.GeneProt():%s:%s: NoGeneName", txid, upca)
//			fmt.Printf("%s\n", msg)
//		}
//		upid := oneP["upid"].Keys()[0]
//		/*
//			if l > 1 {
//				// 267 accessions in 25 taxa;
//				// max 24 GeneNames in 7227 for P02283 (His2B), P02299 (His3), P84040 (His4)
//				// P04745  AMY1_HUMAN      AMY1A; AMY1B; AMY1C: single Chr  UP000005640: Chromosome 1
//				msg := fmt.Sprintf("export.GeneProt():%s:%s: %d GeneNames: %v", txid, upca, l, lblGs)
//				fmt.Printf("%s\n", msg)
//			}
//		*/
//		/// protein synonyms
//		// not really correct TODO
//		/*
//			synPs := oneP["symG"].Keys()
//			for _, gnm := range lblGs {
//				synPs = append(synPs, gnm)
//			}
//		*/
//
//		// removing alternative definitions:
//		bits := strings.Split(oneP["dfnP"].Keys()[0], "(")
//		dfnP := strings.TrimSpace(bits[0])
//		if len(dfnP) == 0 {
//			//msg := fmt.Sprintf("export.GeneProt():%s:%s: NoDefinition", txid, upca)
//			//fmt.Printf("%s\n", msg)
//			dfnP = "Unspecified" // 90 in 25 species
//		} // 20200531: 124
//		score := oneP["score"].Keys()[0]
//		pubmeds := oneP["pubmed"].Keys()
//		// upca2bgwg constructed above, no isoforms herein
//		clsGs := upca2bgwg[upca]["bgwg"].Keys() // sorted
//		if l := len(clsGs); l > 1 {
//			// 20200531: 267
//			// 20200531:  none if multiple genes aggregated
//			// msg := fmt.Sprintf("export.GeneProt():%s:%s: %d BgwGeneIDs: %v", txid, upca, l, clsGs)
//			// fmt.Printf("%s\n", msg)
//		}
//		oriU := rdf.CompU(nss["uniprot"], upca)
//		clsP := fmt.Sprintf("%s%s%s", txid, "/", upca) // NEEDED!
//		clsPU := rdf.CompU(nss["uniprot"], clsP)          // NEEDED!
//		chrs := oneP["chr"].Keys()                     //v34 ALL chromosomes
//		if len(chrs) > 1 {
//			/*
//				msg := fmt.Sprintf("export.GeneProt():%s:%s: chrs: %d ENSG: %v", txid, upca, len(chrs), ensgs)
//				fmt.Printf("%s\n", msg)
//				msg = fmt.Sprintf("export.GeneProt():%s:%s: chrs: %d NCBIG: %v", txid, upca, len(chrs), ncbigs)
//				fmt.Printf("%s\n", msg)
//			*/
//		}
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["ins2cls"], clsU))
//		nlnp++
//		for _, chr := range chrs {
//			chrid := strings.Join([]string{txid, chr}, "/")
//			chrPU := rdf.CompU(nss["uniprot"], chrid)
//			sbP.WriteString(rdf.FormT(clsPU, gpUs["sub2cls"], chrPU))
//			nlnp++
//		}
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["sub2cls"], gpUs["tlp"]))
//		nlnp++
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["gp2txn"], txnU))
//		nlnp++
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2eqv"], oriU))
//		nlnp++
//		// dfn := fmt.Sprintf("Protein class %s from %s", upid, spnm) //v34
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2dfn"], rdf.FormL(dfnP)))
//		nlnp++
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2lbl"], rdf.FormL(upid)))
//		nlnp++
//		// synonyms
//		var synPs []string // per upac
//		synPs = append(synPs, upca)
//		for _, lblG := range lblGs {
//			synPs = append(synPs, lblG)
//		}
//		for _, clsG := range clsGs { // sorted above
//			bgwg2bgwp.Add(clsG, "bgwp", clsP)
//		} // clsGs
//		for _, synP := range synPs {
//			sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2syn"], rdf.FormL(synP)))
//			nlnp++
//		}
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2ori"], oriU))
//		nlnp++
//		sbP.WriteString(rdf.FormT(clsPU, gpUs["evd2lvl"], rdf.FormL(string(score))))
//		nlnp++                           // conversion needed?
//		for _, pubmed := range pubmeds { // sorted above
//			pubmedU := rdf.CompU(nss["pubmed"], pubmed)
//			sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2evd"], pubmedU))
//			nlnp++
//		}
//		xmap.Upac.Add(upca, "bgwp", clsP)
//		upiCs := upac2xrf[upca]["uniparc"].Keys()
//		if l := len(upiCs); l != 1 {
//			msg := fmt.Sprintf("export.GeneProt():%s:%s: %d UpArcs: %v", txid, upca, l, upiCs)
//			panic(errors.New(msg))
//		}
//		upiC := upiCs[0]
//
//		// iso-forms
//		// xrfsC - all xrefs for the current canonical accession
//		// walking over iso-forms, NO canonical accesssion
//		upi2upac := make(util.Set2D)
//		for _, upac := range xrfsC["upac"].Keys() {
//			if upac == "" {
//				panic("upacEmpty")
//			}
//			upiIs := upac2xrf[upac]["uniparc"].Keys()
//			if l := len(upiIs); l != 1 {
//				msg := fmt.Sprintf("export.GeneProt():%s:%s: %d UpArcs: %v", txid, upac, l, upiIs)
//				panic(errors.New(msg))
//			}
//			upi2upac.Add(upiIs[0], upac)
//		}
//		_, ok = upi2upac[upiC]
//		if !ok {
//			upac := fmt.Sprintf("%s-0", upca)
//			upi2upac.Add(upiC, upac)
//		}
//
//		for _, upi := range upi2upac.Keys() {
//			upacs := upi2upac[upi].Keys()
//			if l := len(upacs); l != 1 {
//				msg := fmt.Sprintf("export.GeneProt():%s:%s: %d MultiUpac: %v", txid, upca, l, upacs)
//				panic(errors.New(msg))
//			}
//			upac := upacs[0]
//			lblP := fmt.Sprintf("uniprot#%s", upac)
//			insP := fmt.Sprintf("%s%s%s", clsP, "#", upi)
//			insPU := rdf.CompU(nss["uniprot"], insP)
//			sbP.WriteString(rdf.FormT(insPU, gpUs["ins2cls"], clsPU))
//			nlnp++
//			sbP.WriteString(rdf.FormT(insPU, gpUs["sth2lbl"], rdf.FormL(lblP)))
//			nlnp++
//			// sbP.WriteString(rdf.FormT(insPU, gpUs["sth2ori"], oriU))
//			//nlnp++
//			/// id mapping ///
//			// NB: multiple ENSPs and even NPs per UPAC are common
//			xmap.Upac.Add(upac, "bgwp", insP)
//			xmap.Bgwp.Add(insP, "upac", upac)
//			xrf := ""
//			//			if upac == fmt.Sprintf("%s-0", upca) {
//			//				// no iso-forms or the set of iso-forms excludes the canonical sequence
//			//				upac = upca
//			//				sbP.WriteString(rdf.FormT(insPU, gpUs["sth2clm"], oriU))
//			//				nlnp++
//			//			} else {
//			//				xrf = fmt.Sprintf("%s%s%s", upca, "#", upac)
//			//				xrfU := rdf.CompU(nss["uniprot"], xrf)
//			//				sbP.WriteString(rdf.FormT(insPU, gpUs["sth2clm"], xrfU))
//			//				nlnp++
//			//			}
//			xrf = upi
//			xrfU := rdf.CompU(nss["uniparc"], xrf)
//			sbP.WriteString(rdf.FormT(insPU, gpUs["sth2clm"], xrfU))
//			nlnp++
//			for _, xrf = range upac2xrf[upac]["enspro"].Keys() {
//				enspU := rdf.CompU(nss["enspro"], xrf)
//				sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2clm"], enspU))
//				nlnp++
//				sbP.WriteString(rdf.FormT(insPU, gpUs["sth2clm"], enspU))
//				nlnp++
//				xmap.Ensp.Add(xrf, "bgwp", insP)
//				xmap.Bgwp.Add(insP, "enspro", xrf)
//			}
//			for _, xrf = range upac2xrf[upac]["refseq"].Keys() {
//				rfsqU := rdf.CompU(nss["refseq"], xrf)
//				sbP.WriteString(rdf.FormT(clsPU, gpUs["sth2clm"], rfsqU))
//				nlnp++
//				sbP.WriteString(rdf.FormT(insPU, gpUs["sth2clm"], rfsqU))
//				nlnp++
//				xmap.Rfsq.Add(xrf, "bgwp", insP)
//				xmap.Bgwp.Add(insP, "refseq", xrf)
//			}
//		} // upac2xrf
//	} // prot
//	for bgwP, xrefs := range xmap.Bgwp { // no need to sort
//		// 20200531: canonical accessions have the same UPI as one of the iso-forms
//		upacs := xrefs["upac"].Keys()
//		if l := len(upacs); l > 2 {
//			msg := fmt.Sprintf("export.GeneProt():%s: %d UniProtAccessions: %v", bgwP, l, upacs)
//			fmt.Printf("%s\n", msg)
//		}
//		// 20200531: 159
//		// mostly worm and fly entries, typically 1 TrEMBL, 2 SP: canonical + 1 isoform
//		// the rest (27) human, up to 8 canonical accessions
//		// 20200531: 813320 with a single accession
//
//	}
//	/////////////////////////////////////////////////////////////////////////////
//	// gene -> prot relations
//	/*
//		for _, clsG := range bgwg2bgwp.Keys() {
//			idPmap := bgwg2bgwp[clsG]
//			//for clsG, idPmap := range bgwg2bgwp { // keep this
//			clsGU := rdf.CompU(nss["gene"], clsG)
//			for _, bgwP := range idPmap["bgwp"].Keys() {
//				bgwPU := rdf.CompU(nss["uniprot"], bgwP)
//				sbG.WriteString(rdf.FormT(clsGU, gnUs["gn2gp"], bgwPU))
//				nlng++ // NB: => gene graph
//			}
//		}
//	*/
//	/////////////////////////////////////////////////////////////////////////////
//	if nlng == 0 {
//		msg := fmt.Sprintf("export.GeneProt(): Genes: NoContent")
//		return nlng, nlnp, errors.New(msg)
//	}
//	if nlnp == 0 {
//		msg := fmt.Sprintf("export.GeneProt(): Prots: NoContent")
//		return nlng, nlnp, errors.New(msg)
//	}
//	wfhG.Write([]byte(sbG.String()))
//	sbG.Reset()
//	wfhP.Write([]byte(sbP.String()))
//	sbP.Reset()
//	// xmap export
//	j, err := json.MarshalIndent(&xmap, "", " ")
//	util.CheckE(err)
//	wfhX.Write(j)
//	/////////////////////////////////////////////////////////////////////////////
//	if l := len(cnt["noGx"]); l > 0 {
//		log.Println("export.GeneProt(): Warning: NoGeneXref:", l)
//	}
//	return nlng, nlnp, nil
//} // GeneProt

// Note: no isoforms in this graph
func Gene2phen(duos, gsym2bgw util.Set3D, wpth string) (int, error) {
	nss := rdf.Nss // BGW URI name spaces
	nln := 0
	src := "uniprot"
	srcU := rdf.FormU(nss[src])
	if len(srcU) == 0 {
		msg := fmt.Sprintf("export.Gene2phen():UnknownNamespace: %s", src)
		panic(errors.New(msg))
	}
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"gn2phn",
		"ins2cls",
		"sth2ori",
		"sth2src",
		"sub2cls",
	}
	keys4b["Apys"] = []string{
		"sth2vrs",
		"sth2dfn",
		"sth2lbl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	clsU := rdf.CompU(nss["owl"], "Class")
	// gene2phen graph ini
	wfh, err := os.Create(wpth)
	if err != nil {
		panic(err)
	}
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	graphU := "<http://rdf.biogateway.eu/graph/gene2phen>"
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	myori := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
	oriU := rdf.FormU(myori)
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2ori"], oriU))
	nln++
	if nln < 18 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}
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
		nln++
		sb.WriteString(rdf.FormT(duoU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		clslbl := fmt.Sprintf("%s--%s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
		nln++
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
		nln++

		// an attempt to get desease names in the App, to no avail
		uriR := rdf.CompU(nss["omim"], oriR)
		sb.WriteString(rdf.FormT(uriR, ourUs["sth2lbl"], rdf.FormL(dfn)))
		nln++

		pdc := "gn2phn"
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		// multiple subjects (never happens)
		for _, bgwL := range bgwLs { // sorted above
			uriL := rdf.CompU(nss["gene"], bgwL)
			sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
			nln++
		}
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "object"), uriR))
		nln++

		/// INSTANCES
		insid := fmt.Sprintf("%s%s%s", duoid, "#", src)
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], duoU))
		nln++
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		nln++
		// metadata
		for _, upca := range upcas { // sorted above
			myU := rdf.CompU(nss["uniprot"], upca)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
			nln++
		}
		cntD++
		wfh.Write([]byte(sb.String()))
		sb.Reset()
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
	return nln, nil
} // Gene2phen

// arg1: output of parse.Gaf or parse.Gpa
// arg2: mapping fromn UniProt accession to BGW IDs generated by GeneProt()
// arg3: path for exporting the RDF file
func Prot2go(duos, upac2bgw util.Set3D, wpth string) (int, error) {
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
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"gp2bp",
		"gp2cc",
		"gp2mf",
		"ins2cls",
		"sth2evd",
		"sth2mtd",
		"sth2ori",
		"sth2src",
		"sub2cls",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	clsU := rdf.CompU(nss["owl"], "Class")
	// prot2bp prot2cc prot2mf graph ini
	nln := 0
	wfh, err := os.Create(wpth)
	util.CheckE(err)
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	bits := strings.Split(wpth, "/")
	graph := bits[len(bits)-2]
	graphU := fmt.Sprintf("<http://rdf.biogateway.eu/graph/%s>", graph)
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	/*
		myori := "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?taxonId=" + "9606" // TODO implement
		oriU := rdf.FormU(myori)
		sb.WriteString(rdf.FormT(graphU, ourUs["sth2ori"], oriU)); nln++
	*/
	if nln < 25 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}
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
		// metadata
		myU := rdf.CompU(nss["goa"], oriL)
		sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
		nln++
		// refs sorted above. X1type() preserves the order
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

func Prot2prot(duos, upac2bgw util.Set3D, wpth string) (int, error) {
	nss := rdf.Nss // BGW URI name spaces
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"ins2cls",
		"sth2evd",
		"sth2mtd",
		"sth2src",
		"tlp2tlp",
		"sub2cls",
		"sub2set",
	}
	keys4b["Apys"] = []string{
		"evd2lvl",
		"sth2dfn",
		"sth2lbl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	clsU := rdf.CompU(nss["owl"], "Class")
	srcU := rdf.FormU(nss["intact"])
	// prot2prot graph ini
	nln := 0
	wfh, err := os.Create(wpth)
	util.CheckE(err)
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	graphU := "<http://rdf.biogateway.eu/graph/prot2prot>"
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	/*
		myori := "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:" + "9606" // TODO implement
		oriU := rdf.FormU(myori)
		sb.WriteString(rdf.FormT(graphU, ourUs["sth2ori"], oriU)); nln++
	*/
	if nln < 23 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}
	nln = 0

	stmNS := "http://rdf.biogateway.eu/prot-prot/"
	rdfNS := nss["rdf"]
	count := make(util.Set3D)
	cnt := make(util.Set2D) // count proteins absent in BGW
	cntD := 0
	for _, duoid := range duos.Keys() {
		duo := duos[duoid]
		//for duoid, duo := range duos {
		refs := duo["pubmedIds"].Keys()
		refs = util.X1type(refs, "pubmed", ":")

		bits := strings.Split(duoid, "--")
		idL := bits[0]
		idR := bits[1]
		oriL := strings.Split(idL, "!")[1]                        // UP accession
		oriL = strings.Split(oriL, "-")[0]                        // UP canonical accession
		oriR := strings.Split(idR, "!")[1]                        // UP accession
		oriR = strings.Split(oriR, "-")[0]                        // UP canonical accession
		duoid = fmt.Sprintf("uniprot!%s--uniprot!%s", oriL, oriR) // redefining
		duoU := rdf.CompU(stmNS, duoid)
		bgwLs := upac2bgw[oriL]["bgwp"].Keys() // needed for filtering by RefProt
		bgwRs := upac2bgw[oriR]["bgwp"].Keys() // needed for filtering by RefProt
		if l := counter(bgwLs, cnt, "addP", "dropP", oriL); l == 0 {
			continue
		}
		if l := counter(bgwRs, cnt, "addP", "dropP", oriR); l == 0 {
			continue
		}
		sb.WriteString(rdf.FormT(duoU, ourUs["ins2cls"], clsU))
		nln++
		sb.WriteString(rdf.FormT(duoU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		clslbl := fmt.Sprintf("%s--%s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
		nln++
		clsdfn := fmt.Sprintf("Pair of molecular interactors %s and %s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		nln++
		pdc := "tlp2tlp"
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		for _, bgwL := range bgwLs { // sorted above
			if len(bgwLs) > 1 {
				count.Add("oriL", oriL, bgwL)
			}
			uriL := rdf.CompU(nss["uniprot"], bgwL)
			sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			for _, bgwR := range bgwRs { // sorted above
				if len(bgwRs) > 1 {
					count.Add("oriR", oriR, bgwR)
				}
				uriR := rdf.CompU(nss["uniprot"], bgwR)
				sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "object"), uriR))
				nln++
				sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
				nln++
				sb.WriteString(rdf.FormT(uriR, ourUs[pdc], uriL))
				nln++
			}
		}

		/// INSTANCES
		insid := fmt.Sprintf("%s%s%s", duoid, "#", "intact")
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], duoU))
		nln++
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		nln++
		/// metadata, aggregated for all subjects and objects
		/// intaracton IDs
		for _, key := range duo["intactIds"].Keys() {
			bits := strings.Split(key, ":")
			if bits[0] != "intact" {
				continue
			}
			myU := rdf.CompU(nss["intact"], bits[1])
			//sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
			sb.WriteString(rdf.FormT(insU, ourUs["sub2set"], myU)) // Attn: change prop
			nln++
		}
		/// publications
		for _, item := range refs {
			myU := rdf.CompU(nss["pubmed"], item)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], myU))
			nln++
		}
		/// confidence values
		for _, key := range duo["cnfvals"].Keys() {
			bits := strings.Split(key, ":")
			if bits[0] != "intact-miscore" {
				continue
			}
			sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(bits[1])))
			nln++
		}
		/// interaction types
		for _, key := range duo["intactTypes"].Keys() {
			bits := strings.Split(key, "\"")
			if bits[0] != "psi-mi:" {
				continue
			}
			myU := rdf.CompU(nss["obo"], strings.Replace(bits[1], ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], myU))
			nln++
		}
		/// detection methods
		for _, key := range duo["mtds"].Keys() {
			bits := strings.Split(key, "\"")
			if bits[0] != "psi-mi:" {
				continue
			}
			myU := rdf.CompU(nss["obo"], strings.Replace(bits[1], ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], myU))
			nln++
		}
		cntD++
		wfh.Write([]byte(sb.String()))
		sb.Reset()
	} // duoid
	msg := ""
	if cntD == 0 {
		msg = fmt.Sprintf("export.Prot2prot(): NoDuos")
		panic(errors.New(msg))
	}
	msg = fmt.Sprintf("export.Prot2prot(): Pairs: added: %d dropped: %d", cntD, len(duos)-cntD)
	log.Println(msg)
	msg = fmt.Sprintf("export.Prot2prot(): Prots: added: %d dropped: %d", len(cnt["addP"]), len(cnt["dropP"]))
	log.Println(msg)
	if nln == 0 {
		msg = fmt.Sprintf("export.Prot2prot():%s: NoContent", wpth)
		return nln, errors.New(msg)
	}
	return nln, nil
} // Prot2prot

func Ortho(duos util.Set3D, wpth string) (int, error) {
	// Note: no UP isoforms in this graph; only RefProt canonical accessions
	nss := rdf.Nss // BGW URI name spaces
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"orl2orl",
		"mbr2lst",
		"ins2cls",
		"sth2src",
		"sub2cls",
		"sub2set",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	clsU := rdf.CompU(nss["owl"], "Class")
	// ortho graph ini
	nln := 0
	wfh, err := os.Create(wpth)
	util.CheckE(err)
	defer wfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	graphU := "<http://rdf.biogateway.eu/graph/ortho>"
	srcU := rdf.FormU(nss["uniprot"])
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	if nln < 17 {
		msg := fmt.Sprintf("MalformedHeader")
		panic(errors.New(msg))
	}
	nln = 0
	///////////////////////////////////////////////////////////////////////////////
	stmNS := "http://rdf.biogateway.eu/ortho/"
	rdfNS := nss["rdf"]
	idmkeys := bgw.Orthokeys
	cntD := 0
	for _, duoid := range duos.Keys() {
		duo := duos[duoid]
		duoU := rdf.CompU(stmNS, duoid)
		bits := strings.Split(duoid, "--")
		oriL := strings.Split(bits[0], "!")[1] // UniProt Canonical Accession
		oriR := strings.Split(bits[1], "!")[1] // UniProt Canonical Accession
		sb.WriteString(rdf.FormT(duoU, ourUs["ins2cls"], clsU))
		nln++
		sb.WriteString(rdf.FormT(duoU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		clslbl := fmt.Sprintf("%s--%s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
		nln++
		clsdfn := fmt.Sprintf("Pair of orthologous proteins %s and %s", oriL, oriR)
		sb.WriteString(rdf.FormT(duoU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		nln++
		pdc := "orl2orl"
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		uriL := rdf.CompU(nss["uniprot"], oriL)
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "subject"), uriL))
		nln++
		uriR := rdf.CompU(nss["uniprot"], oriR)
		sb.WriteString(rdf.FormT(duoU, rdf.CompU(rdfNS, "object"), uriR))
		nln++
		sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
		nln++
		sb.WriteString(rdf.FormT(uriR, ourUs[pdc], uriL))
		nln++

		/// INSTANCES
		for _, srck := range duo.Keys() {
			src, ok := idmkeys[srck]
			if !ok {
				continue
			} // needed!
			insid := fmt.Sprintf("%s%s%s", duoid, "#", src)
			insU := rdf.CompU(stmNS, insid)
			sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], duoU))
			nln++
			sb.WriteString(rdf.FormT(insU, ourUs["sth2lbl"], rdf.FormL(clslbl)))
			srcU := rdf.FormU(nss[src])
			sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
			nln++
			// looping over orthology clusters:
			for _, setid := range duo[srck].Keys() {
				prefix := ""
				if srck == "OrthoDB" {
					prefix = "?query="
				}
				lastbit := fmt.Sprintf("%s%s", prefix, setid)
				setU := rdf.FormU(fmt.Sprintf("%s%s", nss[src], lastbit))
				sb.WriteString(rdf.FormT(insU, ourUs["sub2set"], setU)) // part_of
				nln++
				sb.WriteString(rdf.FormT(uriL, ourUs["mbr2lst"], setU))
				nln++
				sb.WriteString(rdf.FormT(uriR, ourUs["mbr2lst"], setU))
				nln++
			}
		}
		cntD++
		wfh.Write([]byte(sb.String()))
		sb.Reset()
	} // duoid
	msg := ""
	msg = fmt.Sprintf("export.Ortho(): Pairs: added: %d dropped: %d", cntD, len(duos)-cntD)
	log.Println(msg)
	/*
		msg = fmt.Sprintf("export.Ortho(): Prots: added: %d dropped: %d", len(cnt["addP"]), len(cnt["dropP"]))
		log.Println(msg)
	*/
	if nln == 0 {
		msg := fmt.Sprintf("export.Ortho(): NoContent")
		return nln, errors.New(msg)
	}
	return nln, nil
} // Ortho
