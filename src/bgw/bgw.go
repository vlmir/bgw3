package bgw

import (
	"encoding/json"
	"github.com/vlmir/bgw3/src/util"
	"io/ioutil"
)

type Column struct {
	Ind1 int
	Dlm1 string
	Ind2 int
	Dlm2 string
	Ind3 int
	Key  string
}

var CV = "3.3.0"

var Upkeys = map[string]string{
	//"Gene_Name":    "gnm",
	//"Gene_Synonym": "gsnm",
	"Ensembl":       "ensgene",
	"Ensembl_PRO":   "ensprotein",
	"EnsemblGenome": "ensom",
	"GeneID":        "ncbigene",
	"RefSeq":        "refseq",
	"UniParc":       "uniparc",
}

var Orthokeys = map[string]string{
	// "KO": "keggortho",
	"OrthoDB": "orthodb",
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
	Bgwg   util.Set3D
	Bgwp   util.Set3D
	Upac   util.Set3D
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
	xm.Lblg = make(util.Set3D)
	xm.Syng = make(util.Set3D)
	xm.Ensg = make(util.Set3D)
	xm.Ncbig = make(util.Set3D)
	xm.Ensp = make(util.Set3D)
	xm.Rfsq = make(util.Set3D)
	xm.Signor = make(util.Set3D)
	*p = xm
}

// TODO to be eliminated
func NewXmap() (xmap Xmap) {
	xmap.Upac = make(util.Set3D)
	xmap.Lblg = make(util.Set3D)
	xmap.Syng = make(util.Set3D)
	xmap.Bgwp = make(util.Set3D)
	xmap.Bgwg = make(util.Set3D)
	xmap.Ensg = make(util.Set3D)
	xmap.Ncbig = make(util.Set3D)
	xmap.Ensp = make(util.Set3D)
	xmap.Rfsq = make(util.Set3D)
	return xmap
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
		{8, "|", 1, ":", -1, "pubmed"},
		{11, "|", 1, "\"", 1, "typeABid"},
		{11, "|", 2, "\"", 1, "typeABlbl"},
		{14, "|", 1, ":", 1, "score"},
		{20, "|", 1, "\"", 1, "typeAid"},
		{20, "|", 2, "\"", 1, "typeAlbl"},
		{21, "|", 1, "\"", 1, "typeBid"},
		{21, "|", 2, "\"", 1, "typeBlbl"},
		{44, "|", 1, "\"", 1, "reglevelid"},
		{44, "|", 2, "\"", 1, "reglevellbl"},
		{45, "|", 1, "\"", 1, "stmid"},
		{45, "|", 2, "\"", 1, "mode"},
	}
	return keys, vals
}

func Keys4rgrs() util.SliceSet {
	keys := make(util.SliceSet)
	keys["Opys"] = []string{
		"preg2targ",
		"nreg2targ",
		"reg2targ",
		"ins2cls",
		"sth2src",
		"gp2bp",
		"sth2rlm",
		"sub2cls",
		"sth2evd",
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
