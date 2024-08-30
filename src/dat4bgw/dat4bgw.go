package dat4bgw

import (
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/util"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
)

func gunzip(pth string) error {
	var cmd *exec.Cmd
	// Attn: all flags separately!
	cmd = exec.Command("gunzip", "-f", pth) // struct
	err := cmd.Run()
	if err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(0), pth, err)
		return errors.New(msg)
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return nil
} // gunzip

// function Rwget() recursively identifies target files and downloads in a specified location
// arg1: base URI of the source
// arg2: pattern for filtering
// download location
func Rwget(ns, mask, ddir string) error {
	// NOT used
	// recursing may take much time !
	// irreproducible set of downloaded files. TODO
	var cmd *exec.Cmd
	// Attn: all flags separately!
	cmd = exec.Command("wget", "-q", "-r", "-nd", "-A", mask, "-P", ddir, ns)
	err := cmd.Run()
	if err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(0), ns, err)
		return errors.New(msg)
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return nil
} // Rwget

// the 3 functions below have comparable performance
// Getfile is consistently fastest in the tests
// from: https://siongui.github.io/2018/03/04/go-run-wget-via-shell-command/
// NB: preserves the time stamp !!
func WgetFile(uri, pth string) (err error) {
	var cmd *exec.Cmd
	cmd = exec.Command("wget", "-q", uri, "-O", pth) // Attn: all flags separately!
	// sending to background with "-b" or "&" does not work here !
	// set pth to "-" for sending the content to STDOUT
	err = cmd.Run()
	if err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return nil
} // WgetFile

func HttpFile(uri, pth string) error {
	// NOT used
	// from: https://golangcode.com/download-a-file-from-a-uri/
	// HttpFile will download a uri to a local file. It's efficient because it will
	// write as it downloads and not load the whole file into memory.
	// HttpFile() seems prone to failures in case of very large files or poor connections

	// Get the data; resp is a ponter to a struct; resp.Body: io.ReadCloser
	//buf := new(bytes.Buffer) // *bytes.Buffer; used only if no pth provided
	// get response
	resp, err := http.Get(uri) // resp: *Reponse
	if err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	defer resp.Body.Close()
	//buf.ReadFrom(resp.Body)
	// Create the file
	wfh, err := os.Create(pth)
	if err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err) // err includes pth
		panic(errors.New(msg))
	}
	defer wfh.Close()
	// Write the body to file
	_, err = io.Copy(wfh, resp.Body)

	return nil
} // HttpFile

//func GetFile(uri, key0, val0, pth string) (*bytes.Buffer, error) {
func GetFile(uri, key0, val0, pth string) error {
	// NOT used
	// Get the data; resp is a ponter to a struct; resp.Body: io.ReadCloser
	// the field 'Header' is needed only for gpa files
	// with key0="Accept" and val0="text/gpad" works with saveOne*()

	// func NewRequest(method, uri string, body io.Reader) (*Request, error)
	// GET requests have no body
	// http.MethodGet constant can be used iso "GET"
	req, err := http.NewRequest("GET", uri, nil) // type Request struct
	if err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	req.Header.Set(key0, val0)

	client := &http.Client{} // type Client struct
	// get response
	// func (c *Client) Do(req *Request) (*Response, error)
	// Do sends an HTTP request and returns an HTTP response
	resp, err := client.Do(req) // type Response struct
	if err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	defer resp.Body.Close() // Body io.ReadCloser
	// StatusCode: https://stackoverflow.com/questions/55210301/error-handling-with-http-newrequest-in-go
	if resp.StatusCode != http.StatusOK {
		// aceptung success codes other than StausOK
		statusOK := resp.StatusCode >= 200 && resp.StatusCode < 300
		if !statusOK {
			msg := fmt.Sprintf("%s: %d", util.FN(0), resp.StatusCode)
			return errors.New(msg)
		}
	}
	if pth != "" {
		// Create the file
		wfh, err := os.Create(pth)
		if err != nil {
			msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err) // err includes pth
			panic(errors.New(msg))
		}
		defer wfh.Close()
		// Write the body to file
		_, err = io.Copy(wfh, resp.Body)
		if err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err)
			return errors.New(msg)
		}
	}
	return nil
} // GetFile

/// Single Taxon Download ///

func saveOneIdmap(txid, pome, datdir string) error {
	// Never use the 'ftp:' protocol here !!!
	ns := "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
	subdir := "idmapping/"
	ext := ".idmapping.gz"
	file := fmt.Sprintf("%s_%s%s", pome, txid, ext)
	pth := fmt.Sprintf("%s%s%s", datdir, subdir, file)
	uri := fmt.Sprintf("%s%s/%s", ns, pome, file)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	if err := gunzip(pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: pth
		return errors.New(msg)
	}
	return nil
} // saveOneIdmap

func saveOneUniprot(txid, datdir, script string) error {
	var cmd *exec.Cmd
	cmd = exec.Command("/usr/bin/python3", script, txid, datdir, "&")
	err := cmd.Run()
	if err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return nil
} // saveOneUniprot

func saveOneColtri(txid, datdir, script string) error {
	var cmd *exec.Cmd
	cmd = exec.Command(script, txid, datdir, "&")
	err := cmd.Run()
	if err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err)
		return errors.New(msg)
	}
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return nil
} // saveOneColtri

func SaveOneSignor(txid string, datdir string) error {
	// SPECIAL CASE: called in the adaptor for 9606 only
	// TODO see if all the files are limited to human data
	// TODO adding mouse and rat data ??
	subdir := "signor/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	ext := ".mi28"
	uri := "https://signor.uniroma2.it/getData.php?type=causalTab"
	pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	ext = ".tsv"
	uri = "https://signor.uniroma2.it/API/getComplexData.php"
	pth = fmt.Sprintf("%s%s%s%s", datdir, subdir, "complexes", ext)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	uri = "https://signor.uniroma2.it/API/getProteinFamilyData.php"
	pth = fmt.Sprintf("%s%s%s%s", datdir, subdir, "families", ext)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	uri = "https://signor.uniroma2.it/getPathwayData.php?relations"
	pth = fmt.Sprintf("%s%s%s%s", datdir, subdir, "pathways", ext)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	return nil
} // SaveOneSignor

func saveOneIntact(txid string, datdir string) error {
	uri := "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:" + txid
	subdir := "intact/"
	ext := ".mi25"
	pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	return nil
} // saveOneIntact

func saveOneGaf(txid string, datdir string, gafpome string) error {
	// NB: this site is NOT recognised by https !!
	// Note: txid is used only in the saved file names
	uri := "http://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/" + gafpome
	subdir := "goa/"
	ext := ".gaf"
	pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	return nil
} // saveOneGaf

func saveOneGpa(txid string, datdir string) error {
	// NOT used
	// Attn: the result is limited to 10000 annotations !!
	// pagination required
	uri := "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?taxonId=" + txid
	subdir := "goa/"
	ext := ".gpa"
	pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	// NB: GetFile MUST be used here !!
	// "text/gaf" for GAF files
	// Header is mandatory
	if err := GetFile(uri, "Accept", "text/gpad", pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	return nil
} // saveOneGpa

func saveOneTflink(uri, txid, datdir string) error {
	subdir := "tflink/"
	ext := ".tsv.gz"
	pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, txid, ext)
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	if err := gunzip(pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}
	return nil
} // saveOneTflink

/// All Taxa Download ///

func SaveAllCtdb(datdir string) error {
	subdir := "ctdb/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	ext := ".tsv.gz"
	// SPECIAL CASE: no filtering by taxon - data for all species in single files
	var uris = map[string]string{
		"chem2gene": "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz",
		"chem2dise": "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz",
		"gene2dise": "https://ctdbase.org/reports/CTD_genes_diseases.tsv.gz",
		"gene2path": "https://ctdbase.org/reports/CTD_genes_pathways.tsv.gz",
		"dise2path": "https://ctdbase.org/reports/CTD_diseases_pathways.tsv.gz",
		"chem2phen": "https://ctdbase.org/reports/CTD_pheno_term_ixns.tsv.gz",
	}
	for lbl := range uris {
		uri := uris[lbl]
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, lbl, ext)
		if err := WgetFile(uri, pth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
			return errors.New(msg)
		}
		if err := gunzip(pth); err != nil {
			msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err) // err includes: pth
			panic(errors.New(msg))
		}
	} // for lbl
	return nil
} // SaveAllCtdb

func SaveAllIdmap(datdir string, txn2prm util.Set2D) error {
	subdir := "idmapping/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	statsU := "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/STATS"
	wpth := fmt.Sprintf("%s%s%s", datdir, subdir, "stats.tsv")
	if err := WgetFile(statsU, wpth); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}

	for txid := range txn2prm {
		pome := txn2prm[txid].Keys()[0]
		if err := saveOneIdmap(txid, pome, datdir); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: txid uri||fpth
			return errors.New(msg)
		}
	}
	return nil
} // SaveAllIdmap

func SaveAllColtri(datdir, scrdir string) error {
	subdir := "coltri/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	taxa := map[string]string{
		"9606":  "human",
		"10090": "mouse",
		"10116": "rat",
	}
	script := scrdir + "download1ctri.py"
	for txid := range taxa {
		if err := saveOneColtri(txid, datdir, script); err != nil {
			// SPECIAL CASE:no URI inside err
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: txid
			return errors.New(msg)
		}
	}
	return nil
} // SaveAllColtri

func SaveAllUniprot(datdir string, txn2prm util.Set2D, scrdir string) error {
	subdir := "uniprot/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}

	uri := "http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
	pth := datdir + "uniprot/9606.var"
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
		return errors.New(msg)
	}

	script := scrdir + "download1up.py"
	for txid := range txn2prm {
		if err := saveOneUniprot(txid, datdir, script); err != nil {
			// SPECIAL CASE:no URI inside err
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: txid
			return errors.New(msg)
		}
	}
	return nil
} // SaveAllUniprot

func SaveAllIntact(datdir string, txn2prm util.Set2D) error {
	subdir := "intact/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	for txid := range txn2prm {
		if err := saveOneIntact(txid, datdir); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err uncludes: txid uri
			return errors.New(msg)
		}
	}
	return nil
}

func SaveAllGaf(datdir string, txn2prm util.Set2D) error {
	subdir := "goa/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	uri := "http://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/proteome2taxid"
	pth := datdir + "goa/gafpomes.tsv"
	// writng to pth
	if err := WgetFile(uri, pth); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err) // err includes: uri
		panic(errors.New(msg))
	}
	// reading from pth
	gafmap, err := util.MakeMap(pth, 1, 2, "\t") // counting from 0
	if err != nil {
		msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err uncludes: MakeMap pth
		panic(errors.New(msg))
	}
	for txid := range txn2prm {
		gpomes := gafmap[txid].Keys()
		if len(gpomes) != 1 {
			msg := fmt.Sprintf("%s: MultipleGoaProteomesFor %s", util.FN(0), txid)
			return errors.New(msg)
		}
		if err := saveOneGaf(txid, datdir, gpomes[0]); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err uncludes: txid uri
			return errors.New(msg)
		}
	} // txid
	return nil
} // SaveAllGaf

func SaveAllGpa(datdir string, txn2prm util.Set2D) error {
	// NOT used
	subdir := "goa/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	for txid := range txn2prm {
		if err := saveOneGpa(txid, datdir); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err uncludes: txid uri||fpth
			return errors.New(msg)
		}
	}
	return nil
} // SaveAllGpa

func SaveAllTflink(datdir string) error {
	subdir := "tflink/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}
	ns := "https://cdn.netbiol.org/tflink/download_files/TFLink_"
	taxa := map[string]string{
		"6239":   "Caenorhabditis_elegans_interactions_All_mitab_v1.0.tsv.gz",
		"7955":   "Danio_rerio_interactions_All_mitab_v1.0.tsv.gz",
		"7227":   "Drosophila_melanogaster_interactions_All_mitab_v1.0.tsv.gz",
		"9606":   "Homo_sapiens_interactions_All_mitab_v1.0.tsv.gz",
		"10090":  "Mus_musculus_interactions_All_mitab_v1.0.tsv.gz",
		"10116":  "Rattus_norvegicus_interactions_All_mitab_v1.0.tsv.gz",
		"559292": "Saccharomyces_cerevisiae_interactions_All_mitab_v1.0.tsv.gz",
	}

	for txid := range taxa {
		uri := ns + taxa[txid]
		if err := saveOneTflink(uri, txid, datdir); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err uncludes: txid uri||fpth
			return errors.New(msg)
		}
	}
	return nil
} // SaveAllTflink

/// Ontologies Download ///

func SaveAllOnto(datdir, year string) error {
	subdir := "onto/"
	if err := os.MkdirAll(filepath.Join(datdir, subdir), 0755); err != nil {
		msg := fmt.Sprintf("%s: %s: %s", util.FN(1), util.FN(0), err)
		panic(errors.New(msg))
	}

	// from BioPortal
	ontos := make(util.Set3D)
	ontos.Add("OMIM", "ext", ".ttl")
	ontos.Add("OMIM", "ns", "https://data.bioontology.org/ontologies/")
	ontos.Add("OMIM", "subm", "/submissions/"+year)
	ontos.Add("OMIM", "key", "/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb")
	for onto, vals := range ontos {
		ext := vals["ext"].Keys()[0]
		ns := vals["ns"].Keys()[0]
		key := vals["key"].Keys()[0]
		subm := vals["subm"].Keys()[0]
		uri := fmt.Sprintf("%s%s%s%s", ns, onto, subm, key)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, strings.ToLower(onto), ext)
		//if _, err := GetFile(uri, "Accept", "text", pth); err != nil { // better avoided
		if err := WgetFile(uri, pth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
			return errors.New(msg)
		}
	}

	// from OBO Foundry
	var nss = map[string]string{
		"bfo":       "http://purl.obolibrary.org/obo/",
		"go-basic":  "http://purl.obolibrary.org/obo/go/",
		"mi":        "http://purl.obolibrary.org/obo/",
		"ro":        "http://purl.obolibrary.org/obo/",
		"sio":       "http://semanticscience.org/ontology/",
		"ncbitaxon": "http://purl.obolibrary.org/obo/",
	}
	ext := ".owl"
	for onto, ns := range nss {
		uri := fmt.Sprintf("%s%s%s", ns, onto, ext)
		pth := fmt.Sprintf("%s%s%s%s", datdir, subdir, onto, ext)
		if err := WgetFile(uri, pth); err != nil {
			msg := fmt.Sprintf("%s: %s", util.FN(0), err) // err includes: uri
			return errors.New(msg)
		}
	}
	return nil
} // SaveAllOnto
