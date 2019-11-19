package main

import (
	"github.com/vlmir/bgw3/pkg/utils" // pkg 'aux'
	"github.com/vlmir/bgw3/pkg/semweb"
	"github.com/vlmir/bgw3/pkg/ancil"
	"github.com/vlmir/bgw3/pkg/export"
	"github.com/vlmir/bgw3/pkg/parse"
	"flag"
	"fmt"
	"log"
	"time"
)

/// tx2pm needed for opening idmapping files
/// mitmap needed for filtering
func geneprot(datdir, rdfdir string, mitmap aux.Set2D, zeno rdf.Zeno, tx2pm aux.Set2D) {
	skippedTxn := 0
	for txid := range mitmap {
		_, ok := tx2pm[txid]
		if !ok {
			skippedTxn++
			continue
		}
		fmt.Println("main,geneprot:inRefProt:txid:", txid)
		ext := ".upt"
		subdir := "uniprot/"
		pth0 := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext) // read IntAct tsv
		subdir = "idmapping/"
		pomes := tx2pm[txid].Keys()
		// now in the parser if len(pomes) != 1 {fmt.Println("main.geneprot:Warning:", txid, "pomes", pomes)}
		pome := fmt.Sprintf("%s%s%s", pomes[0], "_", txid)
		ext = ".idmapping"
		pth1 := fmt.Sprintf("%s%s%s%s", datdir, subdir, pome, ext) // read
		ext = ".nt"
		subdir = "prot/"
		pth2 := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // write prot nt
		subdir = "gene/"
		pth3 := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // write gene nt
		ext = ".json"
		subdir = "xmap/"
		pthx := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // read BGW map json
		var cnt int
		/////////////////////////////////////////////////////////////////////////////
		var idmkeys = map[string]string{
			"Gene_Name":    "gnm",
			"Gene_Synonym": "gsnm",
			"Ensembl_PRO":  "ensp",
			"Ensembl":      "ensg",
			"GeneID":       "ncbig",
			"RefSeq":       "rfsq",
			"UniParc":      "uparc",
		}
		upcas, upacs, gnms, err := parse.Upidmap(pth1, idmkeys)
		if err != nil {
			fmt.Println("main.geneprot:parse.Upidmap:", err)
			continue
		}
		cnt = len(gnms)
		if cnt == 0 {
			fmt.Println("main.geneprot: No data in:", pth1)
			continue
		}
		cnt = len(upcas)
		if cnt == 0 {
			fmt.Println("main.geneprot: No data in:", pth1)
			continue
		}
		cnt = len(upacs)
		if cnt == 0 {
			fmt.Println("main.geneprot: No data in:", pth1)
			continue
		}
		/////////////////////////////////////////////////////////////////////////////
		//pth1 = datdir + "intact.lst" // for filtering, TODO eliminate the hard coding
		updat, txns, err := parse.Updat(pth0, upcas)
		if err != nil {
			fmt.Println("main.geneprot:parse.Updat:", err)
			continue
		}
		keys := txns.Keys()
		if len(txns) == 0 {
			continue
		} // TODO proper error handling
		if len(txns) > 1 {
			fmt.Println("main.geneprot:Multiple taxa, skipping:", "txns", keys)
			continue
		}
		cnt = len(updat)
		if cnt == 0 {
			fmt.Println("main.geneprot: No data in:", pth0)
			continue
		}
		/////////////////////////////////////////////////////////////////////////////
		// passing pointers, seems slightly faster, at most by 10%
		var dat4rdf util.Dat4rdf
		dat4rdf.Udat = &updat
		dat4rdf.Txns = &txns
		dat4rdf.Upac = &upacs
		dat4rdf.Upca = &upcas
		dat4rdf.Gnm = &gnms
		ntg, ntp, err := export.GeneProt(dat4rdf, pth2, pth3, pthx, zeno)
		if err != nil {
			fmt.Println("main.geneprot(): ", err)
			continue
		}
		nts := ntg + ntp
		if nts == 0 {
			fmt.Println("main.geneprot(): No triples for:", txid)
			continue
		}
	}
	fmt.Println("main.geneprot:skippedTxn:", skippedTxn)
}

func tfac2gene(datdir, rdfdir string, txmap aux.Set2D, zeno rdf.Zeno) {
	for txid := range txmap {
		if txid != "9606" {
			continue
		}
		subdir := "tftg/"
		ext := ".f2g"
		pth0 := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext) // read IntAct tsv
		subdir = "tfac2gene/"
		ext = ".nt"
		pth1 := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // write ppi nt
		subdir = "xmap/"
		ext = ".json"
		pthx := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := util.NewXmap()
		err := xmap.Unmarshal(pthx)
		if err != nil {
			continue
		}
		upmap := xmap.Upac
		gsmap := xmap.Gsymb
		/////////////////////////////////////////////////////////////////////////////
		pairs, meta := parse.Tftg(pth0, upmap, gsmap)
		var cnt int
		cnt = len(pairs)
		if cnt == 0 {
			fmt.Println("main.tfac2gene: No data in:", pth0)
			continue
		}
		/////////////////////////////////////////////////////////////////////////////
		nts, err := export.Tftg(pairs, meta, upmap, gsmap, pth1, zeno)
		if err != nil {
			fmt.Println("main.tfac2gene(): ", err)
			continue
		}
		if nts == 0 {
			fmt.Println("main.tfac2gene(): No triples for:", txid)
			continue
		}
	}
}

func gene2phen(datdir, rdfdir string, txmap aux.Set2D, zeno rdf.Zeno) {
	for txid := range txmap {
		if txid != "9606" {
			continue
		}
		subdir := "uniprot/"
		ext := ".var"
		pth0 := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext) // read IntAct tsv
		subdir = "gene2phen/"
		ext = ".nt"
		pth1 := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // write ppi nt
		subdir = "xmap/"
		ext = ".json"
		pthx := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := util.NewXmap()
		err := xmap.Unmarshal(pthx)
		if err != nil {
			continue
		}
		upmap := xmap.Upac
		gsmap := xmap.Gsymb
		/////////////////////////////////////////////////////////////////////////////
		pairs := parse.Upvar(pth0, gsmap)
		var cnt int
		cnt = len(pairs)
		if cnt == 0 {
			fmt.Println("main.gene2phen: No data in:", pth0)
			continue
		}
		/////////////////////////////////////////////////////////////////////////////
		nts, err := export.Upvar(pairs, upmap, gsmap, pth1, zeno)
		if err != nil {
			fmt.Println("main.gene2phen(): ", err)
			continue
		}
		if nts == 0 {
			fmt.Println("main.gene2phen(): No triples for:", txid)
			continue
		}
	}
}

func prot2onto(datdir, rdfdir string, txmap aux.Set2D, zeno rdf.Zeno, fx string) error {
	for txid := range txmap {
		subdir := "goa/"
		pth0 := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, fx) // read Goa data
		subdir = "xmap/"
		pthx := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ".json") // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := util.NewXmap()
		err := xmap.Unmarshal(pthx)
		if err != nil {
			continue
		}
		upmap := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		bps := make(aux.Set3D)
		ccs := make(aux.Set3D)
		mfs := make(aux.Set3D)
		if fx == ".gpa" {
			bps = parse.Gpa(pth0, upmap)
		} else if fx == ".gaf" {
			bps, ccs, mfs = parse.Gaf(pth0, upmap)
		} else {
			err := fmt.Errorf("%s%s", "main.prot2onto(): unrecognized file extension:", fx)
			return err
		}
		cnt := len(bps) + len(ccs) + len(mfs)
		if cnt == 0 {
			fmt.Println("main.prot2onto(): No data in:", pth0)
			continue
		}
		/////////////////////////////////////////////////////////////////////////////
		xpth := ""
		xpth = fmt.Sprintf("%s%s%s%s", rdfdir, "prot2bp/", txid, ".nt") // write goa nt
		ntp, err := export.Goa(bps, upmap, xpth, zeno)
		if err != nil {
			fmt.Println("main.prot2onto(): ", err)
			continue
		}
		xpth = fmt.Sprintf("%s%s%s%s", rdfdir, "prot2cc/", txid, ".nt") // write goa nt
		ntc, err := export.Goa(ccs, upmap, xpth, zeno)
		if err != nil {
			fmt.Println("main.prot2onto(): ", err)
			continue
		}
		xpth = fmt.Sprintf("%s%s%s%s", rdfdir, "prot2mf/", txid, ".nt") // write goa nt
		ntf, err := export.Goa(mfs, upmap, xpth, zeno)
		if err != nil {
			fmt.Println("main.prot2onto(): ", err)
			continue
		}
		nts := ntp + ntc + ntf
		if nts == 0 {
			fmt.Println("main.prot2onto(): No triples for:", txid)
			continue
		}
	}
	return nil
}

func prot2prot(datdir, rdfdir string, txmap aux.Set2D, zeno rdf.Zeno) {
	for txid := range txmap {
		subdir := "intact/"
		ext := ".mit"
		pth0 := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext) // read IntAct tsv
		subdir = "prot2prot/"
		ext = ".nt"
		pth1 := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // write ppi nt
		subdir = "xmap/"
		ext = ".json"
		pthx := fmt.Sprintf("%s%s%s%s", rdfdir, subdir, txid, ext) // read BGW map json
		/////////////////////////////////////////////////////////////////////////////
		xmap := util.NewXmap()
		err := xmap.Unmarshal(pthx)
		if err != nil {
			continue
		}
		upmap := xmap.Upac
		/////////////////////////////////////////////////////////////////////////////
		pairs := parse.Mitab(pth0, upmap)
		var cnt int
		cnt = len(pairs)
		if cnt == 0 {
			fmt.Println("main.prot2prot: No data in:", pth0)
			continue
		}
		/////////////////////////////////////////////////////////////////////////////
		nts, err := export.Mitab(pairs, upmap, pth1, zeno)
		if err != nil {
			fmt.Println("main.prot2prot(): ", err)
			continue
		}
		if nts == 0 {
			fmt.Println("main.prot2prot(): No triples for:", txid)
			continue
		}
	}
}

func main() {
	start := time.Now()
	aP := flag.Bool("a", false, "export [a]ll")
	eP := flag.Bool("e", false, "export gene and protein [e]ntities")
	iP := flag.Bool("i", false, "export molecular [i]nteractions")
	dP := flag.Bool("d", false, "export [d]isease associations")
	gP := flag.Bool("g", false, "export [g]ene ontology annotations")
	rP := flag.Bool("r", false, "export [r]egulatory associations")
	lP := flag.String("l", "./proteomes.lst", "proteome [l]ist")
	uP := flag.String("u", "./zeno.json", "[u]ri map")
	sP := flag.String("s", "./taxa.ls", "[s]elected taxa")
	var n int
	flag.Parse()
	if !flag.Parsed() {
		log.Fatalln("main:flag.Parsed: No flags")
	}
	args := flag.Args()
	n = len(args)
	if n < 2 {
		log.Fatalln("main:Expected more arguments than ", n)
	}
	pth0 := *lP     // path to list of RefProts
	pth1 := args[0] // path to data directory (with a trailing '/')
	pth2 := args[1] // path to rdf directory (with a trailing '/')
	pthz := *uP     // read map json
	pthx := *sP     // path to a list of selected taxa
	log.Println("Started rdf4bgw with args:", args)
	/////////////////////////////////////////////////////////////////////////////
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pthz)
	/////////////////////////////////////////////////////////////////////////////
	/// used only by geneprot() for reading idmapping files
	tx2pm, err := aux.Makemap(pth0, 1, 0, "_")
	if err != nil {
		log.Fatalln("main:", err)
	}
	n = len(tx2pm)
	if n == 0 {
		log.Fatalln("main:Empty map:", pth0)
	}
	log.Println("tx2pm:", n)
	/////////////////////////////////////////////////////////////////////////////
	mitmap, err := aux.Makemap(pthx, 0, 1, ".")
	if err != nil {
		log.Fatalln("main:", err)
	}
	n = len(mitmap)
	if n == 0 {
		log.Fatalln("main:Empty map:", pthx)
	}
	log.Println("mitmap:", n)
	/////////////////////////////////////////////////////////////////////////////
	if *aP || *eP {
		mystart := time.Now()
		geneprot(pth1, pth2, mitmap, zeno, tx2pm) // MUST be run before the others !!!
		log.Println("Done with geneprot in", time.Since(mystart))
	}
	if *aP || *dP {
		mystart := time.Now()
		gene2phen(pth1, pth2, mitmap, zeno)
		log.Println("Done with gene2phen in", time.Since(mystart))
	}
	if *aP || *gP {
		mystart := time.Now()
		ext := ".gaf"
		err = prot2onto(pth1, pth2, mitmap, zeno, ext)
		if err != nil {
			log.Println(err)
		} else {
			log.Println("Done with prot2onto in", time.Since(mystart))
		}
	}
	if *aP || *iP {
		mystart := time.Now()
		prot2prot(pth1, pth2, mitmap, zeno)
		log.Println("Done with prot2prot in", time.Since(mystart))
	}
	if *aP || *rP {
		mystart := time.Now()
		tfac2gene(pth1, pth2, mitmap, zeno)
		log.Println("Done with tfac2gene in", time.Since(mystart))
	}
	log.Println("Done with rdf4bgw in", time.Since(start))
}
