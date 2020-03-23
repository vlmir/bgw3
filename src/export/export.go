package export

import (
	"github.com/vlmir/bgw3/src/utils" // pkg 'aux'
	"github.com/vlmir/bgw3/src/semweb"
	"github.com/vlmir/bgw3/src/ancil"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"strings"
)

func GeneProt(dat4rdf util.Dat4rdf, xpthP string, xpthG string, wpthX string, zeno rdf.Zeno) (int, int, error) {
	keys4g := make(aux.SliceSet)
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
	keys4p := make(aux.SliceSet)
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
	nlg := 0
	nlp := 0
	txns := *dat4rdf.Txns
	if len(txns) == 0 {
		return nlg, nlp, fmt.Errorf("%s", "export.GeneProt: No taxon ID!") // happens!
	}
	txn := txns.Keys()[0]
	updat := *dat4rdf.Udat
	if len(updat) == 0 {
		return nlg, nlp, fmt.Errorf("%s%s", "export.GeneProt: No data for taxon: ", txn)
	}
	upacs := *dat4rdf.Upac
	upcas := *dat4rdf.Upca
	gn2acs := *dat4rdf.Gnm
	cnts := make(aux.Set2D)
	cnts["txns"] = make(aux.Set1D)
	cnts["chrs"] = make(aux.Set1D)
	cnts["bgwgs"] = make(aux.Set1D)
	cnts["bgwps"] = make(aux.Set1D)
	srcU := rdf.FormU(zeno.Uris["uniprot"])
	//srcUa := "" // actual source NB: oriU is used downstream
	dfn := ""
	lbl := ""
	/////////////////////////////////////////////////////////////////////////////
	// gene graph ini
	xfhG, err := os.Create(xpthG)
	if err != nil {
		return nlg, nlp, fmt.Errorf("%s%s%s", "export.GeneProt:os.Create:", xpthG, err)
	}
	defer xfhG.Close()
	var sbG strings.Builder
	gnUs := make(map[string]string)
	header, ng := rdf.Header(gnUs, keys4g, zeno)
	sbG.WriteString(header)
	nlg += ng
	sbG.WriteString(rdf.FormT("<http://biogateway.eu/graph/gene>", gnUs["sth2src"], srcU))
	nlg++
	//srcUa = "<https://www.uniprot.org/uniprot/?query=organism:" + txn + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab>"
	//sbG.WriteString(rdf.FormT("<http://biogateway.eu/graph/gene>", gnUs["sth2ori"], srcUa)); nlg++
	GU := rdf.FormU(zeno.Uris["gene"])
	sbG.WriteString(rdf.FormT(GU, gnUs["sub2cls"], gnUs["bag"]))
	nlg++
	dfn = strings.Join([]string{"The set of genes in Biogateway", "."}, " ")
	sbG.WriteString(rdf.FormT(GU, gnUs["sth2dfn"], rdf.FormL(dfn)))
	nlg++
	sbG.WriteString(rdf.FormT(GU, gnUs["sth2lbl"], rdf.FormL("gene")))
	nlg++
	xfhG.Write([]byte(sbG.String()))
	sbG.Reset()
	// prot graph ini
	xfhP, err := os.Create(xpthP)
	if err != nil {
		return nlg, nlp, fmt.Errorf("%s%s%s", "export.GeneProt:os.Create:", xpthP, err)
	}
	defer xfhP.Close()
	var sbP strings.Builder
	gpUs := make(map[string]string)
	header, np := rdf.Header(gpUs, keys4p, zeno)
	sbP.WriteString(header)
	nlp += np
	sbP.WriteString(rdf.FormT("<http://biogateway.eu/graph/prot>", gpUs["sth2src"], srcU))
	nlp++
	//srcUa = "<https://www.uniprot.org/uniprot/?query=organism:" + txn + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab>"
	//sbP.WriteString(rdf.FormT("<http://biogateway.eu/graph/prot>", gpUs["sth2ori"], srcUa)); nlp++
	PU := rdf.FormU(zeno.Uris["prot"])
	sbP.WriteString(rdf.FormT(PU, gpUs["sub2cls"], gpUs["bag"]))
	nlp++
	dfn = strings.Join([]string{"The set of translation products in Biogateway", "."}, " ")
	sbP.WriteString(rdf.FormT(PU, gpUs["sth2dfn"], rdf.FormL(dfn)))
	nlp++
	sbP.WriteString(rdf.FormT(PU, gpUs["sth2lbl"], rdf.FormL("prot")))
	nlp++
	xfhP.Write([]byte(sbP.String()))
	sbP.Reset()
	/////////////////////////////////////////////////////////////////////////////
	up2bgw := make(aux.Set3D)
	gnm2bgw := make(aux.Set3D)
	gene2prot := make(aux.Set3D)
	spnm := txns[txn]["spnm"].Keys()[0]
	txnU := rdf.CompU(zeno.Uris["ncbitx"], txn)
	txnGU := rdf.CompU(zeno.Uris["gene"], txn)
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sub2cls"], gnUs["bag"]))
	nlg++
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sub2set"], GU))
	nlg++
	sbG.WriteString(rdf.FormT(txnGU, gnUs["gn2txn"], txnU))
	nlg++
	dfn = strings.Join([]string{"The set of", spnm, "genes in Biogateway", "."}, " ")
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
	nlg++
	lbl = strings.Join([]string{"gene", txn}, "/")
	sbG.WriteString(rdf.FormT(txnGU, gnUs["sth2lbl"], rdf.FormL(lbl)))
	nlg++
	txnPU := rdf.CompU(zeno.Uris["prot"], txn)
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sub2cls"], gpUs["bag"]))
	nlp++
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sub2set"], PU))
	nlp++
	sbP.WriteString(rdf.FormT(txnPU, gpUs["gp2txn"], txnU))
	nlp++
	dfn = strings.Join([]string{"The set of", spnm, "translation products in Biogateway", "."}, " ")
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sth2dfn"], rdf.FormL(dfn)))
	nlp++
	lbl = strings.Join([]string{"prot", txn}, "/")
	sbP.WriteString(rdf.FormT(txnPU, gpUs["sth2lbl"], rdf.FormL(lbl)))
	nlp++
	/////////////////////////////////////////////////////////////////////////////
	for chr := range txns[txn]["come"] {
		chrid := strings.Join([]string{txn, chr}, "/")
		chrGU := rdf.CompU(zeno.Uris["gene"], chrid)
		chrPU := rdf.CompU(zeno.Uris["prot"], chrid)
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sub2cls"], gnUs["bag"]))
		nlg++
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sub2set"], txnGU))
		nlg++
		dfn := strings.Join([]string{"The set of genes residing in", spnm, chr, "."}, " ")
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
		nlg++
		lbl = strings.Join([]string{"gene", chrid}, "/")
		sbG.WriteString(rdf.FormT(chrGU, gnUs["sth2lbl"], rdf.FormL(lbl)))
		nlg++
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sub2cls"], gpUs["bag"]))
		nlp++
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sub2set"], txnPU))
		nlp++
		dfn = strings.Join([]string{"The set of tranlation products encoded by", spnm, chr, "."}, " ")
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sth2dfn"], rdf.FormL(dfn)))
		nlp++
		lbl = strings.Join([]string{"prot", chrid}, "/")
		sbP.WriteString(rdf.FormT(chrPU, gpUs["sth2lbl"], rdf.FormL(lbl)))
		nlp++
	}

	/////////////////////////////////////////////////////////////////////////////
	wfhX, err := os.Create(wpthX)
	if err != nil {
		return nlg, nlp, fmt.Errorf("%s%s%s", "export.GeneProt:os.Create:", wpthX, err)
	}
	defer wfhX.Close()
	xmap := util.NewXmap()

	/////////////////////////////////////////////////////////////////////////////
	for lblG, mapG := range gn2acs {
		for _, upca := range mapG["upac"].Keys() {
			oriU := rdf.CompU(zeno.Uris["uniprot"], upca)
			for chr := range updat[upca]["come"] {
				chrid := strings.Join([]string{txn, chr}, "/")
				chrGU := rdf.CompU(zeno.Uris["gene"], chrid)
				idG := fmt.Sprintf("%s%s%s%s%s", txn, "/", chr, "/", lblG)
				up2bgw.Add(upca, "bgwg", idG)
				gnm2bgw.Add(lblG, "bgwg", idG)
				ourGU := rdf.CompU(zeno.Uris["gene"], idG)
				snms := make(map[string]int)
				for gac := range gn2acs[lblG]["upac"] {
					for snm := range upacs[gac]["gsnm"] {
						snms[snm]++ // only gene syns at this point
					}
				}
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sub2cls"], gnUs["gn"]))
				nlg++
				sbG.WriteString(rdf.FormT(ourGU, gnUs["gn2txn"], txnU))
				nlg++
				sbG.WriteString(rdf.FormT(ourGU, gnUs["mbr2lst"], chrGU))
				nlg++
				dfn := strings.Join([]string{"Gene", lblG, "from", spnm, chr, "."}, " ")
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2dfn"], rdf.FormL(dfn)))
				nlg++
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2lbl"], rdf.FormL(lblG)))
				nlg++
				for snm := range snms {
					if snm == lblG {
						continue
					}
					sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2syn"], rdf.FormL(snm)))
					nlg++
				}
				sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2ori"], oriU))
				nlg++ // UP CA
				/// id mapping ///
				xmap.Gsymb.Add(lblG, "bgwg", idG)
				ensgs := make(map[string]int)
				for gac := range gn2acs[lblG]["upac"] {
					for ensg := range upacs[gac]["ensg"] {
						ensgs[ensg]++
					}
				}
				for ensg := range ensgs {
					ensgU := rdf.CompU(zeno.Uris["ensg"], ensg)
					sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2clm"], ensgU))
					nlg++
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
					ncbigU := rdf.CompU(zeno.Uris["ncbig"], ncbig)
					sbG.WriteString(rdf.FormT(ourGU, gnUs["sth2clm"], ncbigU))
					nlg++
					xmap.Ncbig.Add(ncbig, "bgwg", idG)
					xmap.Bgwg.Add(idG, "ncbig", ncbig)
				}
			}
		}
	} // end of gene
	/////////////////////////////////////////////////////////////////////////////
	for upca, onemap := range updat {
		oriU := rdf.CompU(zeno.Uris["uniprot"], upca)
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
			ourBU := rdf.CompU(zeno.Uris["prot"], idG)
			bits := strings.Split(idG, "/")
			chr := bits[1]
			lblG := bits[2]
			chrid := strings.Join([]string{txn, chr}, "/")
			chrPU := rdf.CompU(zeno.Uris["prot"], chrid)
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sub2cls"], gpUs["bag"]))
			nlp++
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sub2set"], chrPU))
			nlp++
			dfn = strings.Join([]string{"The set of tranlation products encoded by gene", lblG, "residing in", spnm, chr, "."}, " ")
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sth2dfn"], rdf.FormL(dfn)))
			nlp++
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sth2lbl"], rdf.FormL(idG)))
			nlp++
			sbP.WriteString(rdf.FormT(ourBU, gpUs["sth2ori"], oriU))
			nlp++

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
				ourPU := rdf.CompU(zeno.Uris["prot"], idP)
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sub2cls"], gpUs["tlp"]))
				nlp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["mbr2lst"], ourBU))
				nlp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2dfn"], rdf.FormL(pdfn)))
				nlp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2lbl"], rdf.FormL(lblP)))
				nlp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["gp2txn"], txnU))
				nlp++
				sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2ori"], oriU))
				nlp++
				for _, snm := range asnms {
					if snm == lblP {
						continue
					}
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2syn"], rdf.FormL(snm)))
					nlp++
				}
				sbP.WriteString(rdf.FormT(ourPU, gpUs["evd2lvl"], rdf.FormL(string(score))))
				nlp++ // conversion needed?
				for _, pubmed := range pubmeds {
					if pubmed == "" {
						continue
					} // TODO see why this occurs
					pubmedU := rdf.CompU(zeno.Uris["pubmed"], pubmed)
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2evd"], pubmedU))
					nlp++
				}
				/// id mapping ///
				// NB: multiple ENSPs and even NPs per UPAC are common
				xmap.Upac.Add(upac, "bgwp", idP)
				xmap.Bgwp.Add(idP, "upac", upac)
				ensps := upacs[upac]["ensp"]
				for ensp := range ensps {
					enspU := rdf.CompU(zeno.Uris["ensp"], ensp)
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2clm"], enspU))
					nlp++
					xmap.Ensp.Add(ensp, "bgwp", idP)
					xmap.Bgwp.Add(idP, "ensp", ensp)
				}
				rfsqs := upacs[upac]["rfsq"]
				for rfsq := range rfsqs {
					rfsqU := rdf.CompU(zeno.Uris["rfsq"], rfsq)
					sbP.WriteString(rdf.FormT(ourPU, gpUs["sth2clm"], rfsqU))
					nlp++
					xmap.Rfsq.Add(rfsq, "bgwp", idP)
					xmap.Bgwp.Add(idP, "rfsq", rfsq)
				}
			} // end of upacs
		} // end of idGs
	} // end of prot
	/////////////////////////////////////////////////////////////////////////////
	for idG, idPmap := range gene2prot {
		ourGU := rdf.CompU(zeno.Uris["gene"], idG)
		for idP := range idPmap["bgwp"] {
			ourPU := rdf.CompU(zeno.Uris["prot"], idP)
			sbG.WriteString(rdf.FormT(ourGU, gnUs["gn2gp"], ourPU))
			nlg++ // NB: => gene graph
		}
	}
	/////////////////////////////////////////////////////////////////////////////
	xfhG.Write([]byte(sbG.String()))
	sbG.Reset()
	xfhP.Write([]byte(sbP.String()))
	sbP.Reset()
	j, err := json.MarshalIndent(&xmap, "", " ")
	if err != nil {
		log.Fatalln("export.GeneProt:json.MarshalIndent:", err)
	}
	wfhX.Write(j)
	return nlg, nlp, nil
}

func Upvar(duos aux.Set3D, upac2bgw aux.Set3D, gsmap aux.Set3D, xpth string, zeno rdf.Zeno) (int, error) {
	nln := 0
	mysrc := "uniprot"
	srcU := rdf.FormU(zeno.Uris[mysrc])
	if len(srcU) == 0 {
		return nln, fmt.Errorf("%s%s", "export.Upvar:Unknown namespase :", mysrc)
	}
	keys4b := make(aux.SliceSet)
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
		return nln, fmt.Errorf("%s%s%s", "export.Upvar:os.Create:", xpth, err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := make(map[string]string)
	header, n := rdf.Header(ourUs, keys4b, zeno)
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
	rdfNS := zeno.Uris["rdf"]
	count := make(aux.Set3D)

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
		uriR := rdf.CompU(zeno.Uris["omim"], oriR)
		// multiple subjects
		for _, ourL := range ourLs {
			if len(ourLs) > 1 {
				count.Add("oriL", oriL, ourL)
			} // 9606:11
			uriL := rdf.CompU(zeno.Uris["gene"], ourL)
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
		myU := rdf.CompU(zeno.Uris["uniprot"], upca)
		sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
		nln++
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	return nln, nil
}

/*
func Upvar(duos aux.Set3D, upac2bgw aux.Set3D,  xpth string, zeno rdf.Zeno) error {
	mysrc := "uniprot"
	srcU := rdf.FormU(zeno.Uris[mysrc])
	if len(srcU) == 0 {
		return fmt.Errorf("%s%s", "export.Upvar:Unknown namespase :", mysrc)
	}
	keys4b := make(aux.SliceSet)
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
	// prot2phen graph ini
	xfh, err := os.Create(xpth)
	if err != nil {
		return fmt.Errorf("%s%s%s", "export.Upvar:os.Create:", xpth, err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := make(map[string]string)
	sb.WriteString(rdf.Header(ourUs, keys4b, zeno)); nln++
	graphU := "<http://biogateway.eu/graph/prot2phen>"
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU)); nln++
	xfh.Write([]byte(sb.String()))
	sb.Reset()

	stmNS := "http://rdf.biogateway.eu/prot-phen/"
	rdfNS := zeno.Uris["rdf"]
	count := make(aux.Set3D)

	for clsid, onemap := range duos {
		clsU := rdf.CompU(stmNS, clsid)
		sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"])); nln++
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid))); nln++

		bits := strings.Split(clsid, "--")
		idL := bits[0]
		idR := bits[1]
		oriL := strings.Split(idL, "!")[1] // UP AC
		oriR := strings.Split(idR, "!")[1] // MIM ID
		nmRs := onemap["nmR"].Keys()
		/*
		if len(nmRs) != 1 {// indeed multiple names in the source
			fmt.Println("export.Upvar:Warning:", idR, "nmRs", nmRs) // 9606:68
		}
		// multiple labels
		for _, nmR := range nmRs {
			if len(nmRs) > 1 { count.Add("oriR", oriR, nmR) }
			onedfn := fmt.Sprintf("%s%s%s%s", "Association between protein ", oriL, " and disease ", nmR)
			sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(onedfn))); nln++
		}
		pdc := "gp2phn"
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc])); nln++
		uriR := rdf.CompU(zeno.Uris["omim"], oriR)
		ourLs := upac2bgw[oriL]["bgwp"].Keys()
		// multiple subjects
		for _, ourL := range(ourLs) {
			if len(ourLs) > 1 { count.Add("oriL", oriL, ourL) }
			uriL := rdf.CompU(zeno.Uris["prot"], ourL)
			sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL)); nln++
			sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR)); nln++
		}
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "object"), uriR)); nln++
		insid := fmt.Sprintf("%s%s%s", clsid, "#", "uniprot")
		insU := rdf.CompU(stmNS, insid)
		sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU)); nln++
		sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU)); nln++
		// metadata
		myU := rdf.CompU(zeno.Uris["uniprot"], oriL)
		sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU)); nln++
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	return nil
}
*/

func Goa(duos aux.Set3D, upac2bgw aux.Set3D, xpth string, zeno rdf.Zeno) (int, error) {
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
	srcU := rdf.FormU(zeno.Uris["goa"])
	keys4b := make(aux.SliceSet)
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
		return nln, fmt.Errorf("%s%s%s", "export.Goa:os.Create: ", xpth, err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := make(map[string]string)
	header, n := rdf.Header(ourUs, keys4b, zeno)
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
	rdfNS := zeno.Uris["rdf"]
	count := make(aux.Set3D)

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

		uriR := rdf.CompU(zeno.Uris["obo"], oriR)
		// multiple subjects
		for _, ourL := range ourLs {
			if len(ourLs) > 1 {
				count.Add("oriL", oriL, ourL)
			}
			uriL := rdf.CompU(zeno.Uris["prot"], ourL)
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
		myU := rdf.CompU(zeno.Uris["goa"], oriL)
		sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
		nln++
		for _, ref := range aux.X1type(refs, "pubmed", "!") {
			myU := rdf.CompU(zeno.Uris["pubmed"], ref)
			sb.WriteString(rdf.FormT(insU, ourUs["sth2evd"], myU))
			nln++
		}
		for _, eco := range onemap["eco"].Keys() { // for GPA files
			myU := rdf.CompU(zeno.Uris["obo"], eco)
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

func Mitab(duos aux.Set3D, upac2bgw aux.Set3D, xpth string, zeno rdf.Zeno) (int, error) {
	keys4b := make(aux.SliceSet)
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
	srcU := rdf.FormU(zeno.Uris["intact"])
	// prot2prot graph ini
	nln := 0
	xfh, err := os.Create(xpth)
	if err != nil {
		return nln, fmt.Errorf("%s%s%s", "export.Mitab():os.Create:", xpth, err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := make(map[string]string)
	header, n := rdf.Header(ourUs, keys4b, zeno)
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
	rdfNS := zeno.Uris["rdf"]
	count := make(aux.Set3D)
	for clsid, onemap := range duos {
		refs := onemap["pubids"].Keys()
		refs = aux.X1type(refs, "pubmed", ":")
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
			uriL := rdf.CompU(zeno.Uris["prot"], ourL)
			sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			for _, ourR := range ourRs {
				if len(ourRs) > 1 {
					count.Add("oriR", oriR, ourR)
				}
				uriR := rdf.CompU(zeno.Uris["prot"], ourR)
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
			myU := rdf.CompU(zeno.Uris["intact"], bits[1])
			//sb.WriteString(rdf.FormT(insU, ourUs["sth2ori"], myU))
			sb.WriteString(rdf.FormT(insU, ourUs["sub2set"], myU)) // Attn: change prop
			nln++
		}
		/// publications
		// TODO use this patter for all?
		for _, item := range refs {
			myU := rdf.CompU(zeno.Uris["pubmed"], item)
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
			myU := rdf.CompU(zeno.Uris["obo"], strings.Replace(bits[1], ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], myU))
			nln++
		}
		/// detection methods
		for key := range onemap["mtds"] {
			bits := strings.Split(key, "\"")
			if bits[0] != "psi-mi:" {
				continue
			}
			myU := rdf.CompU(zeno.Uris["obo"], strings.Replace(bits[1], ":", "_", 1))
			sb.WriteString(rdf.FormT(insU, ourUs["sth2mtd"], myU))
			nln++
		}
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()
	return nln, nil
}

func Tftg(duos aux.Set3D, meta util.Meta, upac2bgw, gsmap aux.Set3D, xpth string, zeno rdf.Zeno) (int, error) {
	keys4b := make(aux.SliceSet)
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
		return nln, fmt.Errorf("%s%s%s", "export.Tftg:os.Create:", xpth, err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := make(map[string]string)
	header, n := rdf.Header(ourUs, keys4b, zeno)
	sb.WriteString(header)
	nln += n
	graphU := "<http://biogateway.eu/graph/tfac2gene>"
	for _, src := range srcs {
		srcU := rdf.FormU(zeno.Uris[src])
		if len(srcU) == 0 {
			return nln, fmt.Errorf("%s%s", "export.Tftg:Unknown namespase :", src)
		}
		sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
		nln++
	}
	xfh.Write([]byte(sb.String()))
	sb.Reset()

	stmNS := "http://rdf.biogateway.eu/tfac-gene/"
	rdfNS := zeno.Uris["rdf"]
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
				uriL := rdf.CompU(zeno.Uris["prot"], ourL)
				sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
				nln++
				for _, ourR := range ourRs {
					uriR := rdf.CompU(zeno.Uris["gene"], ourR)
					sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "object"), uriR))
					nln++
					sb.WriteString(rdf.FormT(uriL, ourUs[pdc], uriR))
					nln++
					for _, src := range srcs {
						insid := fmt.Sprintf("%s%s%s", clsid, "#", src)
						insU := rdf.CompU(stmNS, insid)
						srcU := rdf.FormU(zeno.Uris[src])
						keys := meta.Refs[clsid][src].Keys() // PubMed  IDs
						if len(keys) == 0 {
							continue
						}
						sb.WriteString(rdf.FormT(insU, ourUs["ins2cls"], clsU))
						nln++
						sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
						nln++
						for _, key := range keys {
							pubmedU := rdf.CompU(zeno.Uris["pubmed"], key)
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

func Ortho(duos, upac2bgw aux.Set3D, xpth string, zeno rdf.Zeno) (int, error) {
	keys4b := make(aux.SliceSet)
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
		return nln, fmt.Errorf("%s%s%s", "export.Ortho():os.Create:", xpth, err)
	}
	defer xfh.Close()
	var sb strings.Builder
	ourUs := make(map[string]string)
	header, n := rdf.Header(ourUs, keys4b, zeno)
	sb.WriteString(header)
	nln += n
	graphU := "<http://biogateway.eu/graph/ortho>"
	srcU := rdf.FormU(zeno.Uris["uniprot"])
	sb.WriteString(rdf.FormT(graphU, ourUs["sth2src"], srcU))
	nln++
	xfh.Write([]byte(sb.String()))

	stmNS := "http://rdf.biogateway.eu/ortho/"
	rdfNS := zeno.Uris["rdf"]
	var idmkeys = map[string]string{
		"KO":      "keggortho",
		"OrthoDB": "orthodb",
	}
	for clsid, duo := range duos {
		clsU := rdf.CompU(stmNS, clsid)
		sb.WriteString(rdf.FormT(clsU, ourUs["sub2cls"], ourUs["stm"]))
		nln++
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2lbl"], rdf.FormL(clsid)))
		nln++
		bits := strings.Split(clsid, "--")
		oriL := bits[0]
		oriR := bits[1]
		clsdfn := fmt.Sprintf("%s%s%s%s", "A pair of orthologous proteins ", oriL, " and ", oriR)
		sb.WriteString(rdf.FormT(clsU, ourUs["sth2dfn"], rdf.FormL(clsdfn)))
		nln++
		pdc := "orl2orl"
		sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "predicate"), ourUs[pdc]))
		nln++
		ourLs := upac2bgw[oriL]["bgwp"].Keys()
		ourRs := upac2bgw[oriR]["bgwp"].Keys()
		// multiple subjects and objects
		// 53 human UP ACs with multiple BGW IDs
		for _, ourL := range ourLs {
			uriL := rdf.CompU(zeno.Uris["prot"], ourL)
			sb.WriteString(rdf.FormT(clsU, rdf.CompU(rdfNS, "subject"), uriL))
			nln++
			for _, ourR := range ourRs {
				uriR := rdf.CompU(zeno.Uris["prot"], ourR)
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
			srcU := rdf.FormU(zeno.Uris[src])
			sb.WriteString(rdf.FormT(insU, ourUs["sth2src"], srcU))
			nln++
			for setid, _ := range sets {
				prefix := ""
				if key == "OrthoDB" { prefix = "?query=" }
				lastbit := fmt.Sprintf("%s%s", prefix, setid )
				setU := rdf.FormU(fmt.Sprintf("%s%s", zeno.Uris[src], lastbit))
				sb.WriteString(rdf.FormT(insU, ourUs["sub2set"], setU))
				nln++
			}
		}
		fmt.Println(sb.String())
	}
	return nln, nil
}
