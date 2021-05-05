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
	"KO": "keggortho",
	// "OrthoDB": "orthodb",
}

type Dat4bridge struct {
	Duos  util.Set3D
	OriAs  util.Set3D
	OriBs util.Set3D
}

func (p *Dat4bridge) New() {
	d4b := *p
	d4b.Duos = make(util.Set3D)
	d4b.OriAs = make(util.Set3D)
	d4b.OriBs = make(util.Set3D)
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
	Bgwg  util.Set3D
	Bgwp  util.Set3D
	Upac  util.Set3D
	Lblg  util.Set3D
	Syng  util.Set3D
	Ensg  util.Set3D
	Ncbig util.Set3D
	Ensp  util.Set3D
	Rfsq  util.Set3D
}

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
