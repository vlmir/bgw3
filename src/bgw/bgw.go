package bgw

import (
	"github.com/vlmir/bgw3/src/util" // pkg 'util'
	"encoding/json"
	"fmt"
	"io/ioutil"
)

var Orthokeys = map[string]string{
	"KO":      "keggortho",
	//"OrthoDB": "orthodb",
}

type Meta struct {
	Signs util.Set3D
	Refs  util.Set3D
	Cnfs  util.Set3D
}

func NewMeta() (meta Meta) {
	meta.Signs = make(util.Set3D)
	meta.Refs = make(util.Set3D)
	meta.Cnfs = make(util.Set3D)
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
	Upac *util.Set3D
	Upca *util.Set2D
	Gnm  *util.Set3D
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
		err = fmt.Errorf("%s%s", "bgw.Xmap.Unmarshal:outil.ReadFile: ", err)
		return err
	}
	if err = json.Unmarshal(jdat, &xmap); err != nil {
		err = fmt.Errorf("%s%s", "bgw.Xmap.Unmarshal:json.Unmarshal: ", err)
		return err
	}
	return nil
}

