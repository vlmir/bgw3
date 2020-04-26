package export

import (
	"encoding/json"
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/semweb"
	"github.com/vlmir/bgw3/src/util" // pkg 'util'
	"os"
	"strings"
)

func GeneProt(dat4rdf bgw.Dat4rdf, xpthP string, xpthG string, wpthX string) (int, int, error) {
	nss := rdf.NameSpaces()
	keys4g := make(util.SliceSet)
	keys4g["Opys"] = []string{
		"gn2gp",
		"gn2txn",
		"mbr2lst",
		"sth2clm",
		"sth2ori",
		"sth2src",
		"sub2cls",
		"sub2set",
	}
	keys4g["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"sth2syn",
	}
	keys4g["Prns"] = []string{
		"bag",
		"gn",
	}
	keys4p := make(util.SliceSet)
	keys4p["Opys"] = []string{
		"gp2txn",
		"mbr2lst",
		"sth2clm",
		"sth2evd",
		"sth2ori",
		"sth2src",
		"sub2cls",
		"sub2set",
	}
	keys4p["Apys"] = []string{
		"evd2lvl",
		"sth2dfn",
		"sth2lbl",
		"sth2syn",
	}
	keys4p["Prns"] = []string{
		"bag",
		"tlp",
	}
	nlng := 0
	nlnp := 0
	txns := *dat4rdf.Txns
	if len(txns) == 0 {
		panic(errors.New(fmt.Sprintf("%s", "No taxon ID!"))) // happens!
	}
	txn := txns.Keys()[0]
	updat := *dat4rdf.Udat
	if len(updat) == 0 {
		panic(errors.New(fmt.Sprintf("%s%s", "No data for taxon: ", txn)))
	}
	upacs := *dat4rdf.Upac
	upcas := *dat4rdf.Upca
	gn2acs := *dat4rdf.Gnm
	cnts := make(util.Set2D)
	cnts["txns"] = make(util.Set1D)
	cnts["chrs"] = make(util.Set1D)
	cnts["bgwgs"] = make(util.Set1D)
	cnts["bgwps"] = make(util.Set1D)
	srcU := rdf.FormU(nss["uniprot"])
	graphUg := "<http://biogateway.eu/graph/gene>"
	graphUp := "<http://biogateway.eu/graph/prot>"
	//srcUa := "" // actual source NB: oriU is used downstream
	dfn := ""
	lbl := ""
	/////////////////////////////////////////////////////////////////////////////
	xfhG, err := os.Create(xpthG)
	if err != nil {
		panic(err)
	}
	defer xfhG.Close()
	xfhP, err := os.Create(xpthP)
	if err != nil {
		panic(err)
	}
	defer xfhP.Close()
	// gene graph ini
	var sbG strings.Builder
	gnUs := rdf.FmtURIs(keys4g)// URIs for the graph 'gene'
	header, ng := rdf.Capita(keys4g)
	sbG.WriteString(header)
	nlng += ng
	sbG.WriteString(rdf.FormT(graphUg, gnUs["sth2src"], srcU))
	nlng++
	//srcUa = "<https://www.uniprot.org/uniprot/?query=organism:" + txn + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab>"
	//sbG.WriteString(rdf.FormT("<http://biogateway.eu/graph/gene>", gnUs["sth2ori"], srcUa)); nlng++
	GU := rdf.FormU(nss["gene"])
	sbG.WriteString(rdf.FormT(GU, gnUs["sub2cls"], gnUs["bag"]))
	nlng++
	dfn = strings.Join([]string{"The set of genes in Biogateway", "."}, " ")
	sbG.WriteString(rdf.FormT(GU, gnUs["sth2dfn"], rdf.FormL(dfn)))
	nlng++
	sbG.WriteString(rdf.FormT(GU, gnUs["sth2lbl"], rdf.FormL("gene")))
	nlng++
	xfhG.Write([]byte(sbG.String()))
	sbG.Reset()
	// prot graph ini
	var sbP strings.Builder
	gpUs := rdf.FmtURIs(keys4p)
	header, np := rdf.Capita(keys4p)
	sbP.WriteString(header)
	nlnp += np
	sbP.WriteString(rdf.FormT(graphUp, gpUs["sth2src"], srcU))
	nlnp++
	//srcUa = "<https://www.uniprot.org/uniprot/?query=organism:" + txn + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab>"
	//sbP.WriteString(rdf.FormT("<http://biogateway.eu/graph/prot>", gpUs["sth2ori"], srcUa)); nlnp++
	PU := rdf.FormU(nss["prot"])
	sbP.WriteString(rdf.FormT(PU, gpUs["sub2cls"], gpUs["bag"]))
	nlnp++
	dfn = strings.Join([]string{"The set of translation products in Biogateway", "."}, " ")
	sbP.WriteString(rdf.FormT(PU, gpUs["sth2dfn"], rdf.FormL(dfn)))
	nlnp++
	sbP.WriteString(rdf.FormT(PU, gpUs["sth2lbl"], rdf.FormL("prot")))
	nlnp++
	xfhP.Write([]byte(sbP.String()))
	sbP.Reset()
	/////////////////////////////////////////////////////////////////////////////
	up2bgw := make(util.Set3D)
	gnm2bgw := make(util.Set3D)
	gene2prot := make(util.Set3D)
	spnm := txns[txn]["spnm"].Keys()[0]
	txnU := rdf.CompU(nss["ncbitx"], txn)
	txnGU := rdf.CompU(nss["gene"], txn)
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sub2cls"], gnUs["bag"]))
	nlng++
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sub2set"], GU))
	nlng++
	sbG.WriteString(rdf.FormT(txnGU, gnUs["gn2txn"], txnU))
	nlng++
	dfn = strings.Join([]string{"The set of", spnm, "genes in Biogateway", "."}, " ")
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
	nlng++
	lbl = strings.Join([]string{"gene", txn}, "/")
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sth2lbl"], rdf.FormL(lbl)))
	nlng++
	txnPU := rdf.CompU(nss["prot"], txn)
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sub2cls"], gpUs["bag"]))
	nlnp++
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sub2set"], PU))
	nlnp++
	sbP.WriteString(rdf.FormT(txnPU, gpUs["gp2txn"], txnU))
	nlnp++
	dfn = strings.Join([]string{"The set of", spnm, "translation products in Biogateway", "."}, " ")
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sth2dfn"], rdf.FormL(dfn)))
	nlnp++
	lbl = strings.Join([]string{"prot", txn}, "/")
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sth2lbl"], rdf.FormL(lbl)))
	nlnp++
	/////////////////////////////////////////////////////////////////////////////
	for chr := range txns[txn]["come"] {
		chrid := strings.Join([]string{txn, chr}, "/")
		chrGU := rdf.CompU(nss["gene"], chrid)
		chrPU := rdf.CompU(nss["prot"], chrid)
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sub2cls"], gnUs["bag"]))
		nlng++
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sub2set"], txnGU))
		nlng++
		dfn := strings.Join([]string{"The set of genes residing in", spnm, chr, "."}, " ")
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
		nlng++
		lbl = strings.Join([]string{"gene", chrid}, "/")
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sth2lbl"], rdf.FormL(lbl)))
		nlng++
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sub2cls"], gpUs["bag"]))
		nlnp++
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sub2set"], txnPU))
		nlnp++
		dfn = strings.Join([]string{"The set of tranlation products encoded by", spnm, chr, "."}, " ")
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sth2dfn"], rdf.FormL(dfn)))
		nlnp++
		lbl = strings.Join([]string{"prot", chrid}, "/")
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sth2lbl"], rdf.FormL(lbl)))
		nlnp++
	}

	/////////////////////////////////////////////////////////////////////////////
	wfhX, err := os.Create(wpthX)
	if err != nil {
		panic(err)
	}
	defer wfhX.Close()
	xmap := bgw.NewXmap()

	/////////////////////////////////////////////////////////////////////////////
	for lblG, mapG := range gn2acs {
		for _, upca := range mapG["upac"].Keys() {
			oriU := rdf.CompU(nss["uniprot"], upca)
			for chr := range updat[upca]["come"] {
				chrid := strings.Join([]string{txn, chr}, "/")
				chrGU := rdf.CompU(nss["gene"], chrid)
				idG := fmt.Sprintf("%s%s%s%s%s", txn, "/", chr, "/", lblG)
				up2bgw.Add(upca, "bgwg", idG)
				gnm2bgw.Add(lblG, "bgwg", idG)
				ourGU := rdf.CompU(nss["gene"], idG)
				snms := make(map[string]int)
				for gac := range gn2acs[lblG]["upac"] {
					for snm := range upacs[gac]["gsnm"] {
						snms[snm]++ // only gene syns at this point
					}
				}
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sub2cls"], gnUs["gn"]))
				nlng++
				sbG.WriteString(rdf.FormT(ourGU, gnUs["gn2txn"], txnU))
				nlng++
				sbG.WriteString(rdf.FormT(ourGU, gnUs["mbr2lst"], chrGU))
				nlng++
				dfn := strings.Join([]string{"Gene", lblG, "from", spnm, chr, "."}, " ")
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
				nlng++
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2lbl"], rdf.FormL(lblG)))
				nlng++
				for snm := range snms {
					if snm == lblG {
						continue
					}
					sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2syn"], rdf.FormL(snm)))
					nlng++
				}
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2ori"], oriU))
				nlng++ // UP CA
				/// id mapping ///
				xmap.Gsymb.Add(lblG, "bgwg", idG)
				ensgs := make(map[string]int)
				for gac := range gn2acs[lblG]["upac"] {
					for ensg := range upacs[gac]["ensg"] {
						ensgs[ensg]++
					}
				}
				for ensg := range ensgs {
					ensgU := rdf.CompU(nss["ensg"], ensg)
					sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2clm"], ensgU))
					nlng++
					xmap.Ensg.Add(ensg, "bgwg", idG)
					xmap.Bgwg.Add(idG, "ensg", ensg)
				}
				ncbigs := make(map[string]int)
				for gac := range gn2acs[lblG]["upac"] {
					for ncbig := range upacs[gac]["ncbig"] {
						ncbigs[ncbig]++
					}
				}
				for ncbig := range ncbigs {
					ncbigU := rdf.CompU(nss["ncbig"], ncbig)
					sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2clm"], ncbigU))
					nlng++
					xmap.Ncbig.Add(ncbig, "bgwg", idG)
					xmap.Bgwg.Add(idG, "ncbig", ncbig)
				}
			}
		}
	} // end of gene
	/////////////////////////////////////////////////////////////////////////////
	for upca, onemap := range updat {
		oriU := rdf.CompU(nss["uniprot"], upca)
		upid := onemap["upid"].Keys()[0]
		bits := strings.Split(upid, "_")
		if len(bits) != 2 {
			continue
		}
		lblP := bits[0]
		pdfn := onemap["pdfn"].Keys()[0]
		score := onemap["score"].Keys()[0]
		idGs := up2bgw[upca]["bgwg"]
		refmap := onemap["pubmed"]
		pubmeds := strings.Split(refmap.Keys()[0], "; ")
		for idG := range idGs {
			ourBU := rdf.CompU(nss["prot"], idG)
			bits := strings.Split(idG, "/")
			chr := bits[1]
			lblG := bits[2]
			chrid := strings.Join([]string{txn, chr}, "/")
			chrPU := rdf.CompU(nss["prot"], chrid)
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sub2cls"], gpUs["bag"]))
			nlnp++
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sub2set"], chrPU))
			nlnp++
			dfn = strings.Join([]string{"The set of tranlation products encoded by gene", lblG, "residing in", spnm, chr, "."}, " ")
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sth2dfn"], rdf.FormL(dfn)))
			nlnp++
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sth2lbl"], rdf.FormL(idG)))
			nlnp++
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sth2ori"], oriU))
			nlnp++

			gsnms := upacs[upca]["gsnm"].Keys()
			psnms := append(gsnms, lblG)
			for upac := range upcas[upca] {
				asnms := append(psnms, upac)
				uparcs := upacs[upac]["uparc"].Keys() // normally strictly one per UPAC
				if len(uparcs) != 1 {
					continue
				}
				idP := fmt.Sprintf("%s%s%s", idG, "/", uparcs[0])
				gene2prot.Add(idG, "bgwp", idP)
				ourPU := rdf.CompU(nss["prot"], idP)
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sub2cls"], gpUs["tlp"]))
				nlnp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["mbr2lst"], ourBU))
				nlnp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2dfn"], rdf.FormL(pdfn)))
				nlnp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2lbl"], rdf.FormL(lblP)))
				nlnp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["gp2txn"], txnU))
				nlnp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2ori"], oriU))
				nlnp++
				for _, snm := range asnms {
					if snm == lblP {
						continue
					}
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2syn"], rdf.FormL(snm)))
					nlnp++
				}
				sbP.WriteString(rdf.FormT(ourPU, gpUs["evd2lvl"], rdf.FormL(string(score))))
				nlnp++ // conversion needed?
				for _, pubmed := range pubmeds {
					if pubmed == "" {
						continue
					} // TODO see why this occurs
					pubmedU := rdf.CompU(nss["pubmed"], pubmed)
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2evd"], pubmedU))
					nlnp++
				}
				/// id mapping ///
				// NB: multiple ENSPs and even NPs per UPAC are common
				xmap.Upac.Add(upac, "bgwp", idP)
				xmap.Bgwp.Add(idP, "upac", upac)
				ensps := upacs[upac]["ensp"]
				for ensp := range ensps {
					enspU := rdf.CompU(nss["ensp"], ensp)
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2clm"], enspU))
					nlnp++
					xmap.Ensp.Add(ensp, "bgwp", idP)
					xmap.Bgwp.Add(idP, "ensp", ensp)
				}
				rfsqs := upacs[upac]["rfsq"]
				for rfsq := range rfsqs {
					rfsqU := rdf.CompU(nss["rfsq"], rfsq)
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2clm"], rfsqU))
					nlnp++
					xmap.Rfsq.Add(rfsq, "bgwp", idP)
					xmap.Bgwp.Add(idP, "rfsq", rfsq)
				}
			} // end of upacs
		} // end of idGs
	} // end of prot
	/////////////////////////////////////////////////////////////////////////////
	for idG, idPmap := range gene2prot {
		ourGU := rdf.CompU(nss["gene"], idG)
		for idP := range idPmap["bgwp"] {
			ourPU := rdf.CompU(nss["prot"], idP)
			sbG.WriteString(rdf.FormT(ourGU, gnUs["gn2gp"], ourPU))
			nlng++ // NB: => gene graph
		}
	}
	/////////////////////////////////////////////////////////////////////////////
	xfhG.Write([]byte(sbG.String()))
	sbG.Reset()
	xfhP.Write([]byte(sbP.String()))
	sbP.Reset()
	j, err := json.MarshalIndent(&xmap, "", " ")
	if err != nil {
		//log.Fatalln("export.GeneProt:json.MarshalIndent:", err)
		panic(err)
	}
	wfhX.Write(j)
	return nlng, nlnp, nil
}

func Upvar(duos util.Set3D, upac2bgw util.Set3D, gsmap util.Set3D, xpth string) (int, error) {
	nss := rdf.NameSpaces()
	nln := 0
	mysrc := "uniprot"
	srcU := rdf.FormU(nss[mysrc])
	if len(srcU) == 0 {
		panic(errors.New(fmt.Sprintf("%s%s", "Unknown namespase :", mysrc)))
	}
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"gp2phn",
		"ins2cls",
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
	// gene2phen graph ini
	xfh, err := os.Create(xpth)
	if err != nil {
		panic(err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	graphU := "<http://biogateway.eu/graph/gene2phen>"
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	myori := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
	oriU := rdf.FormU(myori)
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2ori"], oriU))
	nln++
	xfh.Write([]byte(sb.String()))
	sb.Reset()

	stmNS := "http://rdf.biogateway.eu/gene-phen/"
	rdfNS := nss["rdf"]
	count := make(util.Set3D)

	for clsid, onemap := range duos {
		bits := strings.Split(clsid, "--")
		idL := bits[0]
		idR := bits[1]
		oriL := strings.Split(idL, "!")[1] // Gene Name
		ourLs := gsmap[oriL]["bgwg"].Keys()
		if len(ourLs) == 0 {
			fmt.Println("export.Upvar:Warning:", oriL, "ourLs", ourLs) // 9606:0
			continue
		}
		oriR := strings.Split(idR, "!")[1]                                    // MIM ID
		clsid = fmt.Sprintf("%s%s%s%s%s", "gene!", oriL, "--", "omim!", oriR) // redifining
		clsU := rdf.CompU(stmNS, clsid)
		sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
		nln++

		onedfn := fmt.Sprintf("%s%s%s%s", "Association between gene ", oriL, " and disease MIM:", oriR)
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(onedfn)))
		nln++
		upcas := onemap["upca"].Keys()
		upca := upcas[0]
		nmRs := onemap["nmR"].Keys()
		// multiple labels// indeed multiple names in the source
		for _, nmR := range nmRs {
			if len(nmRs) > 1 {
				count.Add("oriR", oriR, nmR)
			} // 9606:68
			//onedfn := fmt.Sprintf("%s%s%s%s", "Association between gene ", oriL, " and disease ", nmR)
			//sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(onedfn))); nln++
		}
		pdc := "gp2phn"
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		uriR := rdf.CompU(nss["omim"], oriR)
		// multiple subjects
		for _, ourL := range ourLs {
			if len(ourLs) > 1 {
				count.Add("oriL", oriL, ourL)
			} // 9606:11
			uriL := rdf.CompU(nss["gene"], ourL)
			sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
			nln++
		}
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "object"), uriR))
		nln++
		insid := fmt.Sprintf("%s%s%s", clsid, "#", mysrc)
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
		nln++
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		nln++
		// metadata
		myU := rdf.CompU(nss["uniprot"], upca)
		sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
		nln++
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	return nln, nil
}

// arg1: output of parse.Gaf or parse.Gpa
// arg2: mapping fromn UniProt accession to BGW IDs generated by GeneProt()
// arg3: path for exporting the RDF file
func Goa(duos util.Set3D, upac2bgw util.Set3D, xpth string) (int, error) {
	nss := rdf.NameSpaces()
	/*
		Attn: isoforms are present in gpa files but lost in rdf files TODO !!
		mironov@manjaro ~/g/g/s/b/prot2onto (master ⚡ → =)> cut -f 2 p53.gpa | grep P04637 | sort -u | wc -l
		8
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
	// prot2onto graph ini
	nln := 0
	xfh, err := os.Create(xpth)
	if err != nil {
		panic(err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	graphU := "<http://biogateway.eu/graph/prot2onto>"
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	/*
		myori := "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?taxonId=" + "9606" // TODO implement
		oriU := rdf.FormU(myori)
		sb.WriteString(rdf.FormT(graphU, ourUs["sth2ori"], oriU)); nln++
	*/
	xfh.Write([]byte(sb.String()))
	sb.Reset()

	stmNS := "http://rdf.biogateway.eu/prot-onto/"
	rdfNS := nss["rdf"]
	count := make(util.Set3D)

	for clsid, onemap := range duos {
		ppys := onemap["ppy"].Keys()
		if len(ppys) != 1 { // unnecessary, may help debugging
			continue
			fmt.Println("export.Goa:Warning:", clsid, "ppys", ppys) // 9606: 0
		}
		refs := onemap["ref"].Keys()
		/// Class level
		clsU := rdf.CompU(stmNS, clsid)
		sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
		nln++

		bits := strings.Split(clsid, "--")
		idL := bits[0]
		idR := bits[1]
		oriL := strings.Split(idL, "!")[1] // UP AC
		oriR := strings.Split(idR, "!")[1] // GO ID
		ourLs := upac2bgw[oriL]["bgwp"].Keys()

		pdc := ppys[0]
		oboid := strings.Replace(oriR, "_", ":", 1)
		clsdfn := fmt.Sprintf("%s%s%s%s%s%s", "Association between protein ", oriL, " and ", gosubs[pdc], " ", oboid)
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		nln++
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++

		uriR := rdf.CompU(nss["obo"], oriR)
		// multiple subjects
		for _, ourL := range ourLs {
			if len(ourLs) > 1 {
				count.Add("oriL", oriL, ourL)
			}
			uriL := rdf.CompU(nss["prot"], ourL)
			sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
			nln++
		}
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "object"), uriR))
		nln++
		/// Instance level
		insid := fmt.Sprintf("%s%s%s", clsid, "#", "goa")
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
		nln++
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		nln++
		// metadata
		myU := rdf.CompU(nss["goa"], oriL)
		sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
		nln++
		for _, ref := range util.X1type(refs, "pubmed", "!") {
			myU := rdf.CompU(nss["pubmed"], ref)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], myU))
			nln++
		}
		for _, eco := range onemap["eco"].Keys() { // for GPA files
			myU := rdf.CompU(nss["obo"], eco)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], myU))
			nln++
		}
		for _, goc := range onemap["goc"].Keys() { // for GAF files
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], rdf.FormL(goc)))
			nln++
		}
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	return nln, nil
}

func Mitab(duos, upac2bgw util.Set3D, xpth string) (int, error) {
	nss := rdf.NameSpaces()
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
	srcU := rdf.FormU(nss["intact"])
	// prot2prot graph ini
	nln := 0
	xfh, err := os.Create(xpth)
	if err != nil {
		panic(err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	graphU := "<http://biogateway.eu/graph/prot2prot>"
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	/*
		myori := "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:" + "9606" // TODO implement
		oriU := rdf.FormU(myori)
		sb.WriteString(rdf.FormT(graphU, ourUs["sth2ori"], oriU)); nln++
	*/
	xfh.Write([]byte(sb.String()))
	sb.Reset()

	stmNS := "http://rdf.biogateway.eu/prot-prot/"
	rdfNS := nss["rdf"]
	count := make(util.Set3D)
	for clsid, onemap := range duos {
		refs := onemap["pubids"].Keys()
		refs = util.X1type(refs, "pubmed", ":")
		//if len(refs) == 0 {continue}

		clsU := rdf.CompU(stmNS, clsid)
		sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
		nln++
		bits := strings.Split(clsid, "--")
		idL := bits[0]
		idR := bits[1]
		oriL := strings.Split(idL, "!")[1] // UP AC
		oriR := strings.Split(idR, "!")[1] // UP AC
		clsdfn := fmt.Sprintf("%s%s%s%s", "A pair of molecular interactors ", oriL, " and ", oriR)
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		nln++
		pdc := "tlp2tlp"
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		ourLs := upac2bgw[oriL]["bgwp"].Keys()
		ourRs := upac2bgw[oriR]["bgwp"].Keys()
		// multiple subjects and objects
		// 53 human UP ACs with multiple BGW IDs
		for _, ourL := range ourLs {
			if len(ourLs) > 1 {
				count.Add("oriL", oriL, ourL)
			}
			uriL := rdf.CompU(nss["prot"], ourL)
			sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			for _, ourR := range ourRs {
				if len(ourRs) > 1 {
					count.Add("oriR", oriR, ourR)
				}
				uriR := rdf.CompU(nss["prot"], ourR)
				sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "object"), uriR))
				nln++
				sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
				nln++
				sb.WriteString(rdf.FormT(uriR, ourUs[pdc], uriL))
				nln++
			}
		}
		insid := fmt.Sprintf("%s%s%s", clsid, "#", "intact")
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
		nln++
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
		nln++
		/// metadata, aggregated for all subjects and objects
		/// intaracton IDs
		for key := range onemap["inacids"] {
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
		// TODO use this patter for all?
		for _, item := range refs {
			myU := rdf.CompU(nss["pubmed"], item)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], myU))
			nln++
		}
		/// confidence values
		for key := range onemap["cnfvals"] {
			bits := strings.Split(key, ":")
			if bits[0] != "intact-miscore" {
				continue
			}
			sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(bits[1])))
			nln++
		}
		/// interaction types
		for key := range onemap["inactps"] {
			bits := strings.Split(key, "\"")
			if bits[0] != "psi-mi:" {
				continue
			}
			myU := rdf.CompU(nss["obo"], strings.Replace(bits[1], ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], myU))
			nln++
		}
		/// detection methods
		for key := range onemap["mtds"] {
			bits := strings.Split(key, "\"")
			if bits[0] != "psi-mi:" {
				continue
			}
			myU := rdf.CompU(nss["obo"], strings.Replace(bits[1], ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], myU))
			nln++
		}
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	return nln, nil
}

func Tftg(duos util.Set3D, meta bgw.Meta, upac2bgw, gsmap util.Set3D, xpth string) (int, error) {
	nss := rdf.NameSpaces()
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"ins2cls",
		"sth2src",
		"sub2cls",
		"rgr2trg",
		"sth2evd",
	}
	keys4b["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"evd2lvl",
		"sth2val",
	}
	keys4b["Prns"] = []string{
		"stm",
	}
	srcs := []string{
		"extri",
		"htri",
		"trrust",
		"tfacts",
		"signor",
		"intact",
		"goa",
	}
	// tfac2gene graph ini
	nln := 0
	xfh, err := os.Create(xpth)
	if err != nil {
		panic(err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, n := rdf.Capita(keys4b)
	sb.WriteString(header)
	nln += n
	graphU := "<http://biogateway.eu/graph/tfac2gene>"
	for _, src := range srcs {
		srcU := rdf.FormU(nss[src])
		if len(srcU) == 0 {
			panic(errors.New(fmt.Sprintf("%s%s", "Unknown namespase :", src)))
		}
		sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
		nln++
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()

	stmNS := "http://rdf.biogateway.eu/tfac-gene/"
	rdfNS := nss["rdf"]
	for clsid, onemap := range duos {
		clsU := rdf.CompU(stmNS, clsid)
		sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
		nln++
		bits := strings.Split(clsid, "--")
		oriLs := onemap["upca"].Keys()
		tfnm := bits[0]
		/*
			only NFKB and AP1
			if len(oriLs) != 1 {fmt.Println("export.Tftg:Warning:", tfnm, "oriLs", oriLs)}
		*/
		oriR := bits[1]
		onedfn := fmt.Sprintf("%s%s%s%s", "Regulation of gene ", oriR, " by transcription factor ", tfnm)
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(onedfn)))
		nln++
		pdc := "rgr2trg"
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		ourRs := gsmap[oriR]["bgwg"].Keys()
		/*
			9606:1110
			if len(ourRs) != 1 {fmt.Println("export.Tftg:Warning:", oriR, "ourRs", ourRs)}
		*/
		for _, oriL := range oriLs {
			ourLs := upac2bgw[oriL]["bgwp"].Keys()
			// 9606:14, 13 unique; only 3 UP ACCs: P62805, P69905, Q16385
			//if len(ourLs) = 0 {fmt.Println("export.Tftg:notInBgw:", oriL); continue} // 9606:14
			if len(ourLs) != 1 {
				fmt.Println("export.Tftg:Warning:", len(ourLs), " BGW IDs for:", oriL)
			} // 9606:14
			for _, ourL := range ourLs {
				uriL := rdf.CompU(nss["prot"], ourL)
				sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
				nln++
				for _, ourR := range ourRs {
					uriR := rdf.CompU(nss["gene"], ourR)
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "object"), uriR))
					nln++
					sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
					nln++
					for _, src := range srcs {
						insid := fmt.Sprintf("%s%s%s", clsid, "#", src)
						insU := rdf.CompU(stmNS, insid)
						srcU := rdf.FormU(nss[src])
						keys := meta.Refs[clsid][src].Keys() // PubMed  IDs
						if len(keys) == 0 {
							continue
						}
						sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
						nln++
						sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
						nln++
						for _, key := range keys {
							pubmedU := rdf.CompU(nss["pubmed"], key)
							sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], pubmedU))
							nln++
						}
						keys = meta.Cnfs[clsid][src].Keys() // confidence levels
						for _, key := range keys {
							sb.WriteString(rdf.FormT(insU, ourUs["evd2lvl"], rdf.FormL(key)))
							nln++
						}
						keys = meta.Signs[clsid][src].Keys() // refulation mode
						for _, key := range keys {
							switch key {
							case "UP":
								sb.WriteString(rdf.FormT(insU, ourUs["sth2val"], rdf.FormL("positive")))
								nln++
							case "DOWN":
								sb.WriteString(rdf.FormT(insU, ourUs["sth2val"], rdf.FormL("negative")))
								nln++
							}
						}
					}
				}
			}
		}
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	return nln, nil
}

func Ortho(duos, upac2bgw util.Set3D, xpth string) (int, error) {
	nss := rdf.NameSpaces()
	keys4b := make(util.SliceSet)
	keys4b["Opys"] = []string{
		"orl2orl",
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
	// ortho graph ini
	nln := 0
	xfh, err := os.Create(xpth)
	if err != nil {
		panic(err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := rdf.FmtURIs(keys4b)
	header, _ := rdf.Capita(keys4b)
	//header, _ := rdf.Header(ourUs, keys4b)
	sb.WriteString(header)
	//nln += n
	graphU := "<http://biogateway.eu/graph/ortho>"
	srcU := rdf.FormU(nss["uniprot"])
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	//nln++
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	///////////////////////////////////////////////////////////////////////////////
	stmNS := "http://rdf.biogateway.eu/ortho/"
	rdfNS := nss["rdf"]
	idmkeys := bgw.Orthokeys
	for clsid, duo := range duos {
		clsU := rdf.CompU(stmNS, clsid)
		bits := strings.Split(clsid, "--")
		oriL := bits[0]
		oriR := bits[1]
		ourLs := upac2bgw[oriL]["bgwp"].Keys()
		ourRs := upac2bgw[oriR]["bgwp"].Keys()
		if len(ourLs) == 0 {continue}
		if len(ourRs) == 0 {continue}
		sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
		nln++
		clsdfn := fmt.Sprintf("%s%s%s%s", "Pair of orthologous proteins ", oriL, " and ", oriR)
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		nln++
		pdc := "orl2orl"
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		// multiple subjects and objects
		// 53 human UP ACs with multiple BGW IDs
		for _, ourL := range ourLs {
			uriL := rdf.CompU(nss["prot"], ourL)
			sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			for _, ourR := range ourRs {
				uriR := rdf.CompU(nss["prot"], ourR)
				sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "object"), uriR))
				nln++
				sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
				nln++
				sb.WriteString(rdf.FormT(uriR, ourUs[pdc], uriL))
				nln++
			}
		}
		for key, sets := range duo {
			src := idmkeys[key]
			insid := fmt.Sprintf("%s%s%s", clsid, "#", src)
			insU := rdf.CompU(stmNS, insid)
			sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
			nln++
			srcU := rdf.FormU(nss[src])
			sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
			nln++
			for setid, _ := range sets {
				prefix := ""
				if key == "OrthoDB" {
					prefix = "?query="
				}
				lastbit := fmt.Sprintf("%s%s", prefix, setid)
				setU := rdf.FormU(fmt.Sprintf("%s%s", nss[src], lastbit))
				sb.WriteString(rdf.FormT(insU, ourUs["sub2set"], setU))
				nln++
			}
		}
		xfh.Write([]byte(sb.String()))
		sb.Reset()
	}
	return nln, nil
}
