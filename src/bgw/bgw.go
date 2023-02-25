package bgw

import (
	"encoding/json"
	"github.com/vlmir/bgw3/src/util"
	"io/ioutil"
)

var CV = "3.3.0"

var Ensomes = map[string]string{
	"6239":   "ensmetazoa",
	"7227":   "ensmetazoa",
	"367110": "ensfungi",
	"559292": "ensfungi",
	"284812": "ensfungi",
	"3702":   "ensplants",
	"3055":   "ensplants",
	"39947":  "ensplants",
	"4577":   "ensplants",
	"44689":  "ensprotists",
	"36329":  "ensprotists",
}

var Upkeys = map[string]string{
	"UniProtKB-ID":      "upid",
	"Gene_Name":         "gnm",
	"Gene_Synonym":      "gsnm",
	"Ensembl":           "ensgene",
	"Ensembl_PRO":       "enspro",
	"EnsemblGenome":     "ensom",
	// "EnsemblGenome_PRO": "ensompro", // transcripts...
	"GeneID":            "ncbigene",
	"RefSeq":            "refseq",
	"UniParc":           "uniparc",
	"NCBI_TaxID":        "ncbitx",
}

var Orthokeys = map[string]string{
	// "KO": "keggortho",
	"OrthoDB": "orthodb",
}

type Column struct {
	Ind1 int
	Dlm1 string
	Ind2 int
	Dlm2 string
	Ind3 int
	Key  string
}

type Dat4bridge struct {
	Src   string
	Taxid string
	Duos  util.Set3D
	Mode  util.Set3D
	Cnts  util.Set2D
}

func (p *Dat4bridge) New() {
	d4b := *p
	d4b.Src = ""
	d4b.Taxid = ""
	d4b.Duos = make(util.Set3D)
	d4b.Mode = make(util.Set3D)
	d4b.Cnts = make(util.Set2D)
	*p = d4b
}

type Dat4rdf struct {
	Udat *util.Set3D
	Txns *util.Set3D
	Gnm  *util.Set3D
	Upac *util.Set3D
}

func (p *Dat4rdf) New() {
	d := *p
	n := make(util.Set3D)
	d.Udat = &n
	d.Txns = &n
	d.Gnm = &n
	d.Upac = &n
	*p = d
}

type Xmap struct {
	// TODO use pointers like in Dat4rdf
	Bgwg   util.Set3D
	Bgwp   util.Set3D
	Upac   util.Set3D
	Lblp   util.Set3D
	Synp   util.Set3D
	Lblg   util.Set3D
	Syng   util.Set3D
	Ensg   util.Set3D
	Ncbig  util.Set3D
	Ensp   util.Set3D
	Rfsq   util.Set3D
	Signor util.Set3D
}

func (p *Xmap) New() {
	xm := *p
	xm.Bgwg = make(util.Set3D)
	xm.Bgwp = make(util.Set3D)
	xm.Upac = make(util.Set3D)
	xm.Lblp = make(util.Set3D)
	xm.Synp = make(util.Set3D)
	xm.Lblg = make(util.Set3D)
	xm.Syng = make(util.Set3D)
	xm.Ensg = make(util.Set3D)
	xm.Ncbig = make(util.Set3D)
	xm.Ensp = make(util.Set3D)
	xm.Rfsq = make(util.Set3D)
	xm.Signor = make(util.Set3D)
	*p = xm
}

func (xmap Xmap) Unmarshal(pthj string) error {
	jdat, err := ioutil.ReadFile(pthj) // []bite
	if err != nil {
		// should be in main funcs instead
		panic(err)
	}
	if err = json.Unmarshal(jdat, &xmap); err != nil {
		panic(err)
	}
	return nil
}

func UpdatParseConf() ([]Column, []Column) {
	keys := []Column{
		{0, "; ", 0, "", 0, ""}, // UniProt canonical accessionbs from multiple proteomes
	}
	vals := []Column{
		// {0, "; ", 0, "|", 0, "upca"}, // single value
		{1, "; ", 0, "|", 0, "upid"}, // single value
		// the next two fields should be avoided!
		// {2, "; ", 0, "|", 0, "gnms"},// '; ' separated but some vals contain ';'!
		// {3, "; ", 0, "|", 0, "gsnms"},// '; ' separated but some vals contain ';'!
		{4, "; ", 0, "|", 0, "taxnm"}, // single value but 59 unique names !
		{5, "; ", 0, "|", 0, "txid"}, // single value but 59 unique names !
		{6, "!", 0, "; ", 0, "pdfns"}, // may contain '|', ';' etc, no '!'
		// the next field may contain multiple proteomes per taxon!
		// {7, "; ", -1, ": ", 0, "poms"}, // '; ' separated, 'pomeid: chrid'
		{8, "; ", -1, "; ", 0, "pubmed"}, // '; ' separated
		{9, "; ", 0, " ", 0, "score"},    // single value
	}
	return keys, vals
} // UpdatParseConf()

func TftgParseConf() ([]Column, []Column) {
	keys := []Column{
		{0, ":", 0, "--", 0, ""},
		{0, ":", 1, "--", 0, ""},
	}
	vals := []Column{
		{1, "|", 0, "|", 0, "uniprot"},
		{2, "|", 0, "|", 0, "ncbig"},
		{3, "|", -1, ";", 0, "pubmed"},
		{4, "|", 0, ";", 0, "score"},
		{5, "|", 0, ";", 0, "mode"},
	}
	return keys, vals
}

func SignorParseConf() ([]Column, []Column) {
	keys := []Column{
		{0, "|", 0, "--", 0, ""},
		{1, "|", 0, "--", 0, ""},
	}
	vals := []Column{
		{0, "|", 0, "|", 0, "Aid"},     // single value; Attn: diverse formats !!
		{1, "|", 0, "|", 0, "Bid"},     // single value; Attn: diverse formats !!
		{8, "|", 1, ":", -1, "pubmed"}, // filtering by "pubmed"
		{11, "|", 1, "\"", 0, "typeABid"},
		{11, "|", 2, "\"", 0, "typeABlbl"},
		{13, "|", 1, ":", 0, "sigid"}, // single value, unique for each line, not accepted by API
		{14, "|", 1, ":", 0, "score"},
		{20, "|", 1, "\"", 0, "typeAid"},
		{20, "|", 2, "\"", 0, "typeAlbl"},
		{21, "|", 1, "\"", 0, "typeBid"},
		{21, "|", 2, "\"", 0, "typeBlbl"},
		{44, "|", 1, "\"", 0, "reglevelid"},
		{44, "|", 2, "\"", 0, "reglevellbl"},
		{45, "|", 1, "\"", 0, "modeid"},
		{45, "|", 2, "\"", 0, "modelbl"},
	}
	return keys, vals
}

func SigPwaysParseConf() ([]Column, []Column) {
	// "|" is not used at all
	keys := []Column{
		{5, "|", 0, "--", 0, ""},  // Aid
		{10, "|", 0, "--", 0, ""}, // Bid
	}
	vals := []Column{
		{0, ";", 0, "|", 0, "pwayid"},
		{1, ";", 0, "|", 0, "pwaylbl"},
		{2, ";", 0, "|", 0, "Albl"},
		{4, ";", 0, "|", 0, "typeAlbl"},
		{5, ";", 0, "|", 0, "Aid"},
		{7, ";", 0, "|", 0, "Blbl"},
		{9, ";", 0, "|", 0, "typeBlbl"},
		{10, ";", 0, "|", 0, "Bid"},
		{12, ";", 0, "|", 0, "modelbl"},
		{16, ";", 0, "|", 0, "txid"}, // the host
		{17, ";", 0, "|", 0, "cellid"},
		{18, ";", 0, "|", 0, "tissueid"},
		{25, ";", 0, "|", 0, "pubmed"},
		{26, ";", 0, "|", 0, "isdirect"}, // "t" or "f"
		{30, ";", 0, "|", 0, "sigid"},    // interaction id, not accepted by API
		{31, ";", 0, "|", 0, "score"},
	}
	return keys, vals
}

func Keys4rgrs() util.SliceSet {
	keys := make(util.SliceSet)
	keys["Opys"] = []string{
		"reg2ptrg",
		"reg2ntrg",
		"reg2utrg",
		"reg2dtrg",
		"reg2itrg",
		"ins2cls",
		"sth2src",
		"gp2bp",
		"sth2rlm",
		"sub2cls",
		"sth2evd",
		"pcs2loc",
		"step2pway",
	}
	keys["Apys"] = []string{
		"sth2dfn",
		"sth2lbl",
		"evd2lvl",
	}
	keys["Prns"] = []string{
		"stm",
	}
	return keys
}
