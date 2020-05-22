package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/util"
	"io"
	"log"
	"net/http"
	"os"
	"os/exec"
	"strings"
	"time"
)

func main() {
	aP := flag.Bool("a", false, "download [a]ll") // should not normally be used
	oP := flag.Bool("o", false, "download [o]ntologies")
	OP := flag.Bool("O", false, "don't download [o]ntologies")
	mP := flag.Bool("m", false, "download ID [m]appings")
	MP := flag.Bool("M", false, "dont't download ID [m]appings")
	iP := flag.Bool("i", false, "download molecular [i]ntaraction data")
	uP := flag.Bool("u", false, "download [u]niprot data")
	gP := flag.Bool("g", false, "download [g]ene ontology annotations")
	//tP := flag.Bool("t",, false, "filter [t]ranscription factor - target gene data")
	lP := flag.String("l", "./proteomes.lst", "proteome [l]ist")
	sP := flag.String("s", "./taxa.ls", "[s]elected taxa")
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
	pth0 := args[0] // path to data directory (with a trailing '/')
	pth1 := *lP     // read list of Reference Proteomes TODO write at Step0
	pthx := *sP     // path to a list of selected taxa
	log.Println("Started with args:", pth0, pth1, pthx)
	n := 0

	/// Step0 ///
	// independent
	if (*aP || *oP) && !*OP {
		start := time.Now()
		// Done with ontos in 54m5.856446087s
		getAllOnto(pth0)
		log.Println("Done with ontos in", time.Since(start))
	}

	mitmap, err := util.Makemap(pthx, 0, 1, ".")
	if err != nil {
		log.Fatalln("main:", err)
	}
	n = len(mitmap)
	if n == 0 {
		log.Fatalln("main:Empty map:", pthx)
	}
	log.Println("mitmap:", n)
	/// Step1 ///
	if (*aP || *mP) && !*MP {
		// takes ~35min
		//uri := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/*.idmapping.gz"
		/*
			var cmd *exec.Cmd
			cmd = exec.Command("pwd") // Attn: all flags separately!
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			err := cmd.Run()
		*/
		for _, txid := range mitmap.Keys() {
			uri := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/*_" + txid + ".idmapping.gz"
			if err := Wget(uri); err != nil {
				log.Println("main.main(): Failed to download idmapping data for: ", txid, err)
				continue
			}
		}
		// exec doesn't understand globbing...
		// TODO genereate pth1 here
	}

	// Attn: Full stop here !!

	/// Step2 ///
	if *aP || *iP {
		start := time.Now()
		/*
			tx2pm, err := util.Makemap(pth1, 1, 0, "_")
			if err != nil { log.Fatalln("main:", err) }
			n = len(tx2pm)
			if n == 0 { log.Fatalln("main:Empty map:", pth1) }
			log.Println("tx2pm:", n)
			getAllIntact(pth0, tx2pm)
		*/
		getAllIntact(pth0, mitmap)
		log.Println("Done with IntAct in", time.Since(start))
		// TODO generate pthx here
	}

	// Attn: Full stop here !!

	/// Step3 ///
	if *aP || *uP {
		mitmap, err := util.Makemap(pthx, 0, 1, ".")
		if err != nil {
			log.Fatalln("main:", err)
		}
		start := time.Now()
		n = len(mitmap)
		if n == 0 {
			log.Fatalln("main:Empty map:", pthx)
		}
		log.Println("mitmap:", n)

		uri := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
		expth := pth0 + "uniprot/9606.var"
		if err := Wget(uri, expth); err != nil {
			log.Println("main:Wget: Failed to download:", expth, err)
		}

		getAllUniprot(pth0, mitmap)
		log.Println("Done with UniProt in", time.Since(start))
	}

	if *aP || *gP {
		start := time.Now()
		uri := "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/proteome2taxid"
		expth := pth0 + "gafpomes.tsv"
		if err := Wget(uri, expth); err != nil {
			log.Println("main:Wget: Failed to download:", expth, err)
		}
		gafmap, err := util.Makemap(expth, 1, 2, "\t") // counting from 0
		if err != nil {
			log.Fatalln("main:", err)
		}
		n = len(gafmap)
		if n == 0 {
			log.Fatalln("main:Empty map:", expth)
		}
		log.Println("gafmap:", n)

		wgetAllGoa(pth0, mitmap, gafmap)
		//getAllGoa(pth0, mitmap, tx2pm)
		log.Println("Done with Goa in", time.Since(start))
	}
}

// from: https://siongui.github.io/2018/03/04/go-run-wget-via-shell-command/
//func Wget(url, filepath string) error {
func Wget(strs ...string) (err error) {
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
		//err = fmt.Errorf("main.Wget: Expexting 1 or 2 args, got: %d", n)
		//return err
		panic(errors.New(fmt.Sprintf("Expexting 1 or 2 args, got: %d", n)))
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

// HttpFile will download a url to a local file. It's efficient because it will
// write as it downloads and not load the whole file into memory.
// from: https://golangcode.com/download-a-file-from-a-url/
func HttpFile(url string, filepath string) error {
	// seems ~2.5 faster then wget
	// Get the data
	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	defer resp.Body.Close()
	// Create the file
	out, err := os.Create(filepath)
	if err != nil {
		return err
	}
	defer out.Close()
	// Write the body to file
	_, err = io.Copy(out, resp.Body)
	return err
}

// TODO generalize
func getFile(uri, expth string) error {
	req, err := http.NewRequest("GET", uri, nil)
	if err != nil {
		return err
	}
	req.Header.Set("Accept", "text/gpad")
	client := &http.Client{}
	resp, err := client.Do(req)
	if err != nil {
		return err
	}
	defer resp.Body.Close()
	out, err := os.Create(expth)
	if err != nil {
		return err
	}
	defer out.Close()
	_, err = io.Copy(out, resp.Body)
	return err
}

func getOneUniprot(txid string, datdir string) error {
	/// UniProt, 3702:~2.5 min
	uri := "https://www.uniprot.org/uniprot/?query=organism:" + txid + "&columns=id,entry%20name,organism,organism-id,protein%20names,proteome,citation,annotation%20score&format=tab"
	subdir := "uniprot/"
	ext := ".upt"
	expth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := HttpFile(uri, expth); err != nil {
		return err
	}
	return nil
}

func getAllUniprot(datdir string, txmap util.Set2D) {
	for txid := range txmap {
		if err := getOneUniprot(txid, datdir); err != nil {
			log.Println("Warning: Failed to download data for:", txid, err)
		}
	}
}

func getOneIntact(txid string, datdir string) error {
	uri := "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:" + txid
	subdir := "intact/"
	ext := ".mit"
	expth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := HttpFile(uri, expth); err != nil {
		log.Println("Warning: Failed to download data for:", txid, err)
		return err
	}
	return nil
}

func getAllIntact(datdir string, txmap util.Set2D) {
	for txid := range txmap {
		getOneIntact(txid, datdir)
	}
}

func getOneGoa(txid string, datdir string) error {
	uri := "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?taxonId=" + txid
	//uri = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?proteome=gcrpCan%2CgcrpIso&geneProductId=P04637&geneProductType=protein&taxonId=9606"
	subdir := "goa/"
	ext := ".gpa"
	expth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := getFile(uri, expth); err != nil {
		log.Println("Warning: Failed to download data for:", txid, err)
		return err
	}
	return nil
}

func wgetOneGoa(txid string, datdir string, gafpome string) error {
	uri := "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/" + gafpome
	subdir := "goa/"
	ext := ".gaf"
	expth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := Wget(uri, expth); err != nil {
		//err := fmt.Errorf("%s%s%s", "Warning: Failed to download data for:", txid, err)
		//return err
		log.Println("txid:", txid)
		panic(err)
	}
	return nil
}

func wgetAllGoa(datdir string, txmap util.Set2D, gafmap util.Set2D) {
	for txid := range txmap {
		gpomes := gafmap[txid].Keys()
		if len(gpomes) != 1 {
			continue
		}
		wgetOneGoa(txid, datdir, gpomes[0])
	}
}

func getAllGoa(datdir string, txmap, tx2pm util.Set2D) {
	for txid := range txmap {
		var srcs = map[string]string{
			"NCBI_TaxID": "upac",
		}
		pomes := tx2pm[txid].Keys()
		if len(pomes) != 1 {
			continue
		}
		ext := ".idmapping"
		subdir := "idmapping/"
		idmpth := fmt.Sprintf("%s%s%s%s%s%s", datdir, subdir, pomes[0], "_", txid, ext)
		//idmap, err := parse.Idmap(idmpth, srcs, 2, 0)
		idmap, err := parse.Idmap(idmpth, srcs, 2, 1, 0) // TODO test
		if err != nil {
			continue
		}
		fmt.Println(len(idmap))
		subdir = "goa/"
		ext = ".gpa"
		expth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
		out, err := os.Create(expth)
		if err != nil {
			continue
		}
		defer out.Close()
		for _, upac := range idmap[txid]["upac"].Keys() {
			qs := "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?"
			uri := fmt.Sprintf("%s%s%s%s%s", qs, "fproteome=gcrpCan%2CgcrpIso&geneProductId=", upac, "&geneProductType=protein&taxonId=", txid)
			// "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?proteome=gcrpCan%2CgcrpIso&geneProductId=P04637&geneProductType=protein&taxonId=9606"
			req, err := http.NewRequest("GET", uri, nil)
			if err != nil {
				continue
			}
			req.Header.Set("Accept", "text/gpad")
			client := &http.Client{}
			resp, err := client.Do(req)
			if err != nil {
				continue
			}
			defer resp.Body.Close()
			_, err = io.Copy(out, resp.Body)
			if err != nil {
				continue
			}
		}
	}
}

func getAllOnto(datdir string) {
	ns := "http://data.bioontology.org/ontologies/"
	key := "/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb"
	ff := "&download_format=rdf"
	subdir := "onto/"
	ext := ".xml"
	for _, onto := range [...]string{"BFO", "RO", "MI", "ECO", "GO"} {
		uri := fmt.Sprintf("%s%s%s%s", ns, onto, key, ff)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, strings.ToLower(onto), ext)
		if err := HttpFile(uri, pth); err != nil {
			log.Println("main.getAllOnto:", err)
		}
	}
	// Attn: needs to be updated !
	var subms = map[string]string{
		"OMIM":      "/submissions/16",
		"NCBITAXON": "/submissions/14",
	}
	ext = ".ttl"
	for _, onto := range [...]string{"OMIM", "NCBITAXON"} {
		uri := fmt.Sprintf("%s%s%s%s", ns, onto, subms[onto], key)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, strings.ToLower(onto), ext)
		if err := HttpFile(uri, pth); err != nil {
			log.Println("main.getAllOnto:", err)
		}
	}
	var nss = map[string]string{
		"sio": "http://semanticscience.org/ontology/",
	}
	ext = ".owl"
	for _, onto := range [...]string{"sio"} {
		uri := fmt.Sprintf("%s%s%s", nss[onto], onto, ext)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, onto, ext)
		if err := HttpFile(uri, pth); err != nil {
			log.Println("main.getAllOnto:", err)
		}
	}
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func x1tftg(pthR, pthW string, cols [5]int) error {
	fhR, err := os.Open(pthR)
	check(err)
	defer fhR.Close()
	fhW, err := os.Create(pthW)
	check(err)
	defer fhW.Close()
	scanner := bufio.NewScanner(fhR)
	for scanner.Scan() { // by default scans for '\n'
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}
		if string(line[0]) == "#" {
			continue
		}
		cells := strings.Split(line, "\t")
		n := len(cells)
		m := 37
		if n < m {
			msg := "Want at least: %d cells, have: %d for: %s in %s"
			fmt.Printf(msg, m, n, cells[0], pthR)
			continue
		}
		nl := cells[0]
		for _, i := range cols {
			if i == 0 {
				nl = strings.Join([]string{nl, ""}, "\t")
			} else {
				nl = strings.Join([]string{nl, cells[i-1]}, "\t")
			}
		}
		fhW.WriteString(fmt.Sprintf("%s\n", nl))
	}
	return nil
}
