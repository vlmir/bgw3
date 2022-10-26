package main

import (
	"bytes"
	"errors"
	"flag"
	"fmt"
	"github.com/vlmir/bgw3/src/util"
	"io"
	"log"
	"net/http"
	"os"
	"os/exec"
	"strings"
	"time"
)

// function Rwget() recirsively identifies target files and downloads in a slecified location
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
// HttpFile() seems prone to failures in case of very large files of poor connctions

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
	} else if n == 1 {
		url = strs[0]
		cmd = exec.Command("wget", "-q", url) // Attn: all flags separately!
	} else {
		panic(errors.New(fmt.Sprintf("Expexting 1 or 2 args, got: %d", n)))
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

// HttpFile will download a url to a local file. It's efficient because it will
// write as it downloads and not load the whole file into memory.
// from: https://golangcode.com/download-a-file-from-a-url/
func HttpFile(url string, wpth string) (*bytes.Buffer, error) {
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

func saveOneUniprot(txid string, datdir string) error {
	// "https://www.uniprot.org/uniprot/?query=organism:9606&columns=id,entry%20name,genes(PREFERRED),genes(ALTERNATIVE),organism,organism-id,protein%20names,proteome,citation&format=tab" -O 9606.upt
	uri := "https://www.uniprot.org/uniprot/?query=organism:" + txid + "&columns=id,entry%20name,genes(PREFERRED),genes(ALTERNATIVE),organism,organism-id,protein%20names,proteome,citation,annotation%20score,database(GeneID)&format=tab"
	subdir := "uniprot/"
	ext := ".upt"
	wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	//if err := HttpFile(uri, wpth); err != nil {
	if _, err := GetFile(uri, "Accept", "text", wpth); err != nil {
		//if err := WgetFile(uri, wpth); err != nil {
		log.Println("saveOneUniprot(): Warning: Failed to download data for:", txid, err)
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
	uri := "http://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/" + gafpome
	subdir := "goa/"
	ext := ".gaf"
	wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	//if err := HttpFile(uri, wpth); err != nil {
	if _, err := GetFile(uri, "Accept", "text", wpth); err != nil {
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

/// Multiple Taxa Download ///

// very slow, why ??
func saveAllIdmap(datdir string, txmap util.Set2D) {
	ns := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
	for txid := range txmap {
		mask := fmt.Sprintf("_%s%s", txid, ".idmapping.gz")
		subdir := "idmapping/"
		ddir := datdir + subdir
		if err := Rwget(ns, mask, ddir); err != nil {
			log.Println("Warning: Failed to download data for:", txid, err)
		}
	}
}

func saveAllUniprot(datdir string, txmap util.Set2D) {
	for txid := range txmap {
		if err := saveOneUniprot(txid, datdir); err != nil {
			log.Println("Warning: Failed to download data for:", txid, err)
		}
	}
}

func saveAllIntact(datdir string, txmap util.Set2D) {
	for txid := range txmap {
		saveOneIntact(txid, datdir)
	}
}

func saveAllGaf(datdir string, txmap, gafmap util.Set2D) {
	for txid := range txmap {
		gpomes := gafmap[txid].Keys()
		if len(gpomes) != 1 {
			continue
		}
		saveOneGaf(txid, datdir, gpomes[0])
	}
}

func saveAllGpa(datdir string, txmap util.Set2D) {
	for txid := range txmap {
		saveOneGpa(txid, datdir)
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
		// 2021-05-2o
		"OMIM": "/submissions/20",
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
	//for _, onto := range [...]string{"sio"} {
	for onto, ns := range nss {
		uri := fmt.Sprintf("%s%s%s", ns, onto, ext)
		wpth := fmt.Sprintf("%s%s%s%s", datdir, subdir, onto, ext)
		if _, err := HttpFile(uri, wpth); err != nil {
			// NB: GetFile() failes to download ncbitaxon.owl !
			log.Println("main.saveAllOnto:", err)
		}
	}
}

// not used
//func x1tftg(pthR, pthW string, cols [5]int) error {
//	fhR, err := os.Open(pthR)
//	check(err)
//	defer fhR.Close()
//	fhW, err := os.Create(pthW)
//	check(err)
//	defer fhW.Close()
//	scanner := bufio.NewScanner(fhR)
//	for scanner.Scan() { // by default scans for '\n'
//		line := scanner.Text()
//		if len(line) == 0 {
//			continue
//		}
//		if string(line[0]) == "#" {
//			continue
//		}
//		cells := strings.Split(line, "\t")
//		n := len(cells)
//		m := 37
//		if n < m {
//			msg := "Want at least: %d cells, have: %d for: %s in %s"
//			fmt.Printf(msg, m, n, cells[0], pthR)
//			continue
//		}
//		nl := cells[0]
//		for _, i := range cols {
//			if i == 0 {
//				nl = strings.Join([]string{nl, ""}, "\t")
//			} else {
//				nl = strings.Join([]string{nl, cells[i-1]}, "\t")
//			}
//		}
//		fhW.WriteString(fmt.Sprintf("%s\n", nl))
//	}
//	return nil
//}

func main() {
	aP := flag.Bool("a", false, "download [a]ll") // should not normally be used
	oP := flag.Bool("o", false, "download [o]ntologies")
	OP := flag.Bool("O", false, "don't download [o]ntologies")
	mP := flag.Bool("m", false, "download ID [m]appings")
	MP := flag.Bool("M", false, "dont't download ID [m]appings")
	iP := flag.Bool("i", false, "download molecular [i]ntaraction data")
	uP := flag.Bool("u", false, "download [u]niprot data")
	gP := flag.Bool("g", false, "download [g]ene ontology annotations")
	tP := flag.String("t", "./taxa.tls", "selected [t]axa")
	var cnt int
	flag.Parse()
	if !flag.Parsed() {
		log.Fatalln("failed to parse flags")
	}
	args := flag.Args()
	cnt = len(args)
	if cnt < 1 {
		log.Fatalln("Expecting more arguments than ", cnt)
	}
	datdir := args[0] // path to data directory (with a trailing '/')
	//pth1 := *lP     // read list of Reference Proteomes TODO write at Step0
	pthx := *tP // path to a list of selected taxa
	log.Println("Started with args:", args)
	n := 0

	/// Step0 ///
	// independent
	if (*aP || *oP) && !*OP {
		start := time.Now()
		// Done with ontos in 54m5.856446087s - What's that ???
		saveAllOnto(datdir)
		log.Println("Done with ontos in", time.Since(start))
	}

	txnmap, err := util.MakeMap(pthx, 0, 1, ".")
	if err != nil {
		log.Fatalln("main:", err)
	}
	n = len(txnmap)
	if n == 0 {
		log.Fatalln("main:Empty map:", pthx)
	}
	log.Println("txnmap:", n)
	/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	/// Step1 ///
	if (*aP || *mP) && !*MP {
		start := time.Now()
		ns := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
		mask := "idmapping.gz"
		subdir := "idmapping/"
		ddir := datdir + subdir
		err := Rwget(ns, mask, ddir)
		if err != nil {
			panic(err)
		}
		/*
			terrrible slow this way
			~5h for 20 taxa vs ~1h for ~1500 taxa
			saveAllIdmap(datdir, txnmap)
		*/
		log.Println("Done with idmappings in", time.Since(start))
	}

	// Attn: Full stop here !!

	/// Step2 ///
	if *aP || *iP {
		start := time.Now()
		saveAllIntact(datdir, txnmap)
		log.Println("Done with IntAct in", time.Since(start))
	}

	// Attn: Full stop here !!

	/// Step3 ///
	if *aP || *uP {
		start := time.Now()
		uri := "http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
		wpth := datdir + "uniprot/9606.var"
		if _, err := HttpFile(uri, wpth); err != nil {
			panic(err)
		}

		//		saveAllUniprot(datdir, txnmap)
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
		saveAllGaf(datdir, txnmap, gafmap)
		// saveAllGpa(datdir, txnmap) // 10000 lines limit, do NOT use !
		log.Println("Done with Goa in", time.Since(start))
	}
}
