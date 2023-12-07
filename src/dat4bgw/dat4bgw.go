package main

import (
	"bytes"
	"errors"
	"flag"
	"fmt"
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/util"
	"io"
	"log"
	"net/http"
	"os"
	"os/exec"
	"strings"
	"time"
)

// function Rwget() recursively identifies target files and downloads in a slecified location
// irreproducible set of downloaded files. TODO
// arg1: base URI of the source
// arg2: pattern for filtering
// download location
func Rwget(ns, mask, ddir string) error {
	var cmd *exec.Cmd
	// Attn: all flags separately!
	cmd = exec.Command("wget", "-q", "-r", "-nd", "-A", mask, "-P", ddir, ns)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

/// Single File Download ///

// the 3 functions below have comparable performance
// HttpFile() seems prone to failures in case of very large files or poor connections

// from: https://siongui.github.io/2018/03/04/go-run-wget-via-shell-command/
// NB: preserves the time stamp !!
func WgetFile(strs ...string) (err error) {
	n := len(strs)
	url := ""
	var cmd *exec.Cmd
	if n == 2 {
		url = strs[0]
		pth := strs[1]
		cmd = exec.Command("wget", "-q", url, "-O", pth) // Attn: all flags separately!
		// sending to background with "-b" or "&" does not work here !
		// set pth to "-" for sending the content to STDOUT
	} else {
		panic(errors.New(fmt.Sprintf("Expexting at least 2 args, got: %d", n)))
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

// HttpFile will download a url to a local file. It's efficient because it will
// write as it downloads and not load the whole file into memory.
// from: https://golangcode.com/download-a-file-from-a-url/
func HttpFile(url, wpth string) (*bytes.Buffer, error) {
	// 2023-10-19: used only for ontologies
	// Get the data; resp is a ponter to a struct; resp.Body: io.ReadCloser
	buf := new(bytes.Buffer)
	resp, err := http.Get(url)
	if err != nil {
		return buf, err
	}
	defer resp.Body.Close()
	if wpth == "" {
		buf.ReadFrom(resp.Body)
	} else {
		// Create the file
		out, err := os.Create(wpth)
		if err != nil {
			return buf, err
		}
		defer out.Close()
		// Write the body to file
		_, err = io.Copy(out, resp.Body)
	}
	return buf, err
}

func GetFile(uri, key0, val0, wpth string) (*bytes.Buffer, error) {
	// 2023-10-19: used only for OMIM
	// Get the data; resp is a ponter to a struct; resp.Body: io.ReadCloser
	// TODO generalize
	// the field 'Header' is needed only for gpa files
	// with key0="Accept" and val0="text/gpad" works with saveOne*(), why ??
	buf := new(bytes.Buffer)
	req, err := http.NewRequest("GET", uri, nil)
	if err != nil {
		return buf, err
	}
	req.Header.Set(key0, val0)
	client := &http.Client{}
	resp, err := client.Do(req)
	if err != nil {
		return buf, err
	}
	defer resp.Body.Close()
	if wpth == "" {
		buf.ReadFrom(resp.Body)
	} else {
		out, err := os.Create(wpth)
		if err != nil {
			return buf, err
		}
		defer out.Close()
		_, err = io.Copy(out, resp.Body)
	}
	return buf, err
}

/// Single Taxon Download ///

func saveOneIdmap(txid, pome, datdir string) error {
	// Never use the 'ftp:' protocol hear !!!
	ns := "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
	subdir := "idmapping/"
	ext := ".idmapping.gz"
	file := fmt.Sprintf("%s_%s%s", pome, txid, ext)
	wpth := fmt.Sprintf("%s%s%s", datdir, subdir, file)
	uri := fmt.Sprintf("%s%s/%s", ns, pome, file)
	// HttpFile() cannot be used here!
	if err := WgetFile(uri, wpth); err != nil {
		log.Println("saveOneIdmap(): Warning: Failed to download data for:", txid, err)
		return err
	}
	return nil
}

func saveOneUniprot(txid, datdir, script string) error {
	var cmd *exec.Cmd
	cmd = exec.Command("/usr/bin/python3", script, txid, datdir, "&")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func saveOneCtri(txlbl, datdir, script string) error {
	var cmd *exec.Cmd
	cmd = exec.Command(script, txlbl, datdir, "&")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func saveOneSignor(txid string, datdir string) error {
// TODO see if all the files are limited to human data
	subdir := "signor/"
	ext := ".mi28"
	uri := "https://signor.uniroma2.it/getData.php?type=causalTab"
	wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := WgetFile(uri, wpth); err != nil {
		log.Println("saveOneSignor(): Warning: Failed to download data for:", txid, err)
		return err
	}
	ext = ".tsv"
	uri = "https://signor.uniroma2.it/API/getComplexData.php"
	wpth = fmt.Sprintf("%s%s%s%s", datdir, subdir, "complexes", ext)
	if err := WgetFile(uri, wpth); err != nil {
		log.Println("saveOneSignor(): Warning: Failed to download protein complexes: ", err)
		return err
	}
	uri = "https://signor.uniroma2.it/API/getProteinFamilyData.php"
	wpth = fmt.Sprintf("%s%s%s%s", datdir, subdir, "families", ext)
	if err := WgetFile(uri, wpth); err != nil {
		log.Println("saveOneSignor(): Warning: Failed to download protein families: ", err)
		return err
	}
	uri = "https://signor.uniroma2.it/getPathwayData.php?relations"
	wpth = fmt.Sprintf("%s%s%s%s%s", datdir, subdir, "pathways", txid, ext)
	if err := WgetFile(uri, wpth); err != nil {
		log.Println("saveOneSignor(): Warning: Failed to download pathway relations for:", txid, err)
		return err
	}
	return nil
}

func saveOneIntact(txid string, datdir string) error {
	uri := "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:" + txid
	subdir := "intact/"
	ext := ".mi25"
	wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	// NB: HttpFile() tends to fail here with large files !
	if err := WgetFile(uri, wpth); err != nil {
		//if _, err := GetFile(uri, "Accept", "text", wpth); err != nil {
		log.Println("saveOneIntact(): Warning: Failed to download data for:", txid, err)
		return err
	}
	return nil
}

func saveOneGaf(txid string, datdir string, gafpome string) error {
	// NB: this site is NOT recognised by https !!
	// Note: txid is used only in saved file names
	uri := "http://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/" + gafpome
	subdir := "goa/"
	ext := ".gaf"
	wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	// the 3 functions generate identical files in about the same time
	//if _, err := HttpFile(uri, wpth); err != nil {
	if err := WgetFile(uri, wpth); err != nil {
		//if _, err := GetFile(uri, "Accept", "text", wpth); err != nil {
		log.Println("txid:", txid)
		panic(err)
	}
	return nil
}

func saveOneGpa(txid string, datdir string) error {
	// Attn: the result is limited to 10000 annotations !!
	uri := "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?taxonId=" + txid
	subdir := "goa/"
	ext := ".gpa"
	wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	// NB: GetFile() MUST be used here !!
	// "text/gaf" for GAF files
	// Header is mandatory
	if _, err := GetFile(uri, "Accept", "text/gpad", wpth); err != nil {
		log.Println("saveOneGpa(): Warning: Failed to download data for:", txid, err)
		return err
	}
	return nil
}

func saveOneTflink(uri, txid, datdir string) error {
	subdir := "tflink/"
	ext := ".tsv.gz"
	wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := WgetFile(uri, wpth); err != nil {
		log.Println("saveOneTflink(): Warning: Failed to download data for:", txid, err)
		return err
	}
	return nil
}

/// Multiple Taxa Download ///

func saveAllCtdb(datdir string) error {
	// no filtering by taxon - data for all species in single files
	var uris = map[string]string{
		"chem2gene": "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz",
		"chem2dise": "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz",
		"gene2dise": "https://ctdbase.org/reports/CTD_genes_diseases.tsv.gz",
		"gene2path": "https://ctdbase.org/reports/CTD_genes_pathways.tsv.gz",
		"dise2path": "https://ctdbase.org/reports/CTD_diseases_pathways.tsv.gz",
		"chem2phen": "https://ctdbase.org/reports/CTD_pheno_term_ixns.tsv.gz",
	}
	subdir := "ctdb/"
	ext := ".tsv.gz"
	for lbl := range uris {
		uri := uris[lbl]
		wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, lbl, ext)
		if err := WgetFile(uri, wpth); err != nil {
			log.Println("saveOneCtdb(): Warning: Failed to download data for:", lbl, err)
			return err
		}
	} // for lbl
	return nil
}

func saveAllIdmap(datdir string, txn2prm util.Set2D) {
	for txid := range txn2prm {
		pome := txn2prm[txid].Keys()[0]
		if err := saveOneIdmap(txid, pome, datdir); err != nil {
			log.Println("Warning: Failed to download data for:", txid, err)
		}
	}
}

func saveAllCtri(datdir, scrdir string) {
	var taxa = map[string]string{
		"9606": "human",
	}
	script := scrdir + "download1ctri.py"
	for txid := range taxa {
		if err := saveOneCtri(taxa[txid], datdir, script); err != nil {
			log.Println("Warning: Failed to download data for:", txid, err)
		}
	}
}

func saveAllUniprot(datdir string, txn2prm util.Set2D, scrdir string) {
	script := scrdir + "download1up.py"
	for txid := range txn2prm {
		if err := saveOneUniprot(txid, datdir, script); err != nil {
			log.Println("Warning: Failed to download data for:", txid, err)
		}
	}
}

func saveAllIntact(datdir string, txn2prm util.Set2D) {
	for txid := range txn2prm {
		saveOneIntact(txid, datdir)
	}
}

func saveAllGaf(datdir string, txn2prm, gafmap util.Set2D) {
	for txid := range txn2prm {
		gpomes := gafmap[txid].Keys()
		if len(gpomes) != 1 {
			continue
		}
		saveOneGaf(txid, datdir, gpomes[0])
	}
}

func saveAllGpa(datdir string, txn2prm util.Set2D) {
	for txid := range txn2prm {
		saveOneGpa(txid, datdir)
	}
}

func saveAllTflink(datdir string) {
	ns := "https://cdn.netbiol.org/tflink/download_files/TFLink_"
	for txid := range bgw.Tflink {
		uri := ns + bgw.Tflink[txid]
		saveOneTflink(uri, txid, datdir)
	}
}

/// Ontologies Download ///

func saveAllOnto(datdir string) {
	subdir := "onto/"
	ext := ""
	ns := ""

	ns = "https://data.bioontology.org/ontologies/"
	key := "/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb"
	//	ff := "&download_format=rdf"
	//	ext = ".xml"
	//	 for _, onto := range [...]string{"BFO", "OBOREL", "MI", "ECO", "GO", } {
	//		uri := fmt.Sprintf("%s%s%s%s", ns, onto, key, ff)
	//		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, strings.ToLower(onto), ext)
	//		if err := HttpFile(uri, pth); err != nil {
	//			log.Println("main.saveAllOnto:", err)
	//		}
	//	}
	// Attn: needs to be updated prior each download!
	var subms = map[string]string{
		"OMIM": "/submissions/23", // 2023-03-07
	}
	ext = ".ttl"
	for onto, subm := range subms {
		uri := fmt.Sprintf("%s%s%s%s", ns, onto, subm, key)
		wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, strings.ToLower(onto), ext)
		//if err := HttpFile(uri, wpth); err != nil {
		if _, err := GetFile(uri, "Accept", "text", wpth); err != nil {
			log.Println("main.saveAllOnto:", err)
		}
	}

	var nss = map[string]string{
		"sio":       "http://semanticscience.org/ontology/",
		"ncbitaxon": "http://purl.obolibrary.org/obo/",
		"ro":        "http://purl.obolibrary.org/obo/",
		"bfo":       "http://purl.obolibrary.org/obo/",
		"mi":        "http://purl.obolibrary.org/obo/",
		"go-basic":  "http://purl.obolibrary.org/obo/go/",
	}
	ext = ".owl"
	for onto, ns := range nss {
		uri := fmt.Sprintf("%s%s%s", ns, onto, ext)
		wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, onto, ext)
		if _, err := HttpFile(uri, wpth); err != nil {
			// NB: GetFile() failes to download ncbitaxon.owl !
			log.Println("main.saveAllOnto:", err)
		}
	}
}

func main() {
	aP := flag.Bool("a", false, "download [a]ll") // should not normally be used
	oP := flag.Bool("o", false, "download [o]ntologies")
	OP := flag.Bool("O", false, "don't download [o]ntologies")
	mP := flag.Bool("m", false, "download ID [m]appings")
	MP := flag.Bool("M", false, "dont't download ID [m]appings")
	iP := flag.Bool("i", false, "download molecular [i]ntaraction data")
	uP := flag.Bool("u", false, "download [u]niprot data")
	sP := flag.Bool("s", false, "download [s]ignor data")
	lP := flag.Bool("l", false, "download tf[l]ink data")
	gP := flag.Bool("g", false, "download [g]ene ontology annotations")
	cP := flag.Bool("c", false, "download [c]tdb data")
	tP := flag.Bool("t", false, "download collec[t]ri data")
	flag.Parse()
	if !flag.Parsed() {
		log.Fatalln("failed to parse flags")
	}
	args := flag.Args()
	cnt := len(args)
	if cnt < 3 {
		log.Fatalln("Expecting more arguments than ", cnt)
	}
	log.Println("Started with args:", args)
	datdir := args[0]                              // path to data directory (with a trailing '/')
	rpthT := args[1]                               // path to a list of selected taxa and proteomes
	scrdir := args[2]                              // path to the script for downloading UP for one taxon
	txn2prm, err := util.MakeMap(rpthT, 1, 0, "_") // txnID -> proteomeID
	if err != nil {
		log.Fatalln("main:", err)
	}
	n := len(txn2prm)
	if n == 0 {
		log.Fatalln("main:Empty map:", rpthT)
	}
	log.Println("txn2prm:", n)

	/// Step0 ///
	// independent
	if (*aP || *oP) && !*OP {
		start := time.Now()
		saveAllOnto(datdir)
		log.Println("Done with ontos in", time.Since(start))
	}

	if *aP || *lP {
		start := time.Now()
		saveAllTflink(datdir)
		log.Println("Done with TFlink in", time.Since(start))
	}

	if *aP || *sP {
		start := time.Now()
		saveOneSignor("9606", datdir)
		log.Println("Done with Signor in", time.Since(start))
	}

	if *aP || *cP {
		start := time.Now()
		if err := saveAllCtdb(datdir); err != nil {
			panic(err)
		}
		log.Println("Done with CTDbase in", time.Since(start))
	}

	/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	if (*aP || *mP) && !*MP {
		start := time.Now()
		saveAllIdmap(datdir, txn2prm)
		log.Println("Done with idmappings in", time.Since(start))
	}

	if *aP || *iP {
		start := time.Now()
		saveAllIntact(datdir, txn2prm)
		log.Println("Done with Intact in", time.Since(start))
	}

	if *aP || *tP {
		start := time.Now()
		saveAllCtri(datdir, scrdir)
		log.Println("Done with Ctri in", time.Since(start))
	}

	if *aP || *uP {
		start := time.Now()
		uri := "http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
		wpth := datdir + "uniprot/9606.var"
		if _, err := HttpFile(uri, wpth); err != nil {
			panic(err)
		}

		saveAllUniprot(datdir, txn2prm, scrdir)
		log.Println("Done with UniProt in", time.Since(start))
	}

	if *aP || *gP {
		start := time.Now()
		uri := "http://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/proteome2taxid"
		wpth := datdir + "goa/gafpomes.tsv"
		if _, err := HttpFile(uri, wpth); err != nil {
			panic(err)
		}
		gafmap, err := util.MakeMap(wpth, 1, 2, "\t") // counting from 0
		if err != nil {
			//log.Fatalln("main:", err)
			panic(err)
		}
		n = len(gafmap)
		if n == 0 {
			msg := fmt.Sprintf("Empty map: %s", wpth)
			panic(errors.New(msg))
		}
		log.Println("gafmap:", n)
		saveAllGaf(datdir, txn2prm, gafmap)
		// saveAllGpa(datdir, txn2prm) // 10000 lines limit, do NOT use !
		log.Println("Done with Goa in", time.Since(start))
	}
}
