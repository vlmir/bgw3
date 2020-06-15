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
	"Ensembl_PRO":  "ensp",
	"Ensembl":      "ensg",
	"GeneID":       "ncbig",
	"RefSeq":       "rfsq",
	"UniParc":      "uparc",
}

var Orthokeys = map[string]string{
	"KO":      "keggortho",
	// "OrthoDB": "orthodb",
}

type Meta struct {
	Refs  util.Set3D
	Cnfs  util.Set3D
	Signs util.Set3D
}

func NewMeta() (meta Meta) {
	meta.Refs = make(util.Set3D)
	meta.Cnfs = make(util.Set3D)
	meta.Signs = make(util.Set3D)
	return meta
}

// TODO make it work
func (meta Meta) New() {
	meta.Signs = make(util.Set3D)
	meta.Refs = make(util.Set3D)
	meta.Cnfs = make(util.Set3D)
}

type Dat4rdf struct {
	Udat *util.Set3D
	Txns *util.Set3D
	Gnm  *util.Set3D
	Upac *util.Set3D
}

type Xmap struct {
	Bgwg  util.Set3D
	Bgwp  util.Set3D
	Upac  util.Set3D
	Gsymb util.Set3D
	Ensg  util.Set3D
	Ncbig util.Set3D
	Ensp  util.Set3D
	Rfsq  util.Set3D
}

func NewXmap() (xmap Xmap) {
	xmap.Upac = make(util.Set3D)
	xmap.Gsymb = make(util.Set3D)
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
