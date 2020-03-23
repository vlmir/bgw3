package util

import (
	"github.com/vlmir/bgw3/src/utils" // pkg 'aux'
	"encoding/json"
	"fmt"
	"io/ioutil"
)

type Meta struct {
	Signs aux.Set3D
	Refs  aux.Set3D
	Cnfs  aux.Set3D
}

func NewMeta() (meta Meta) {
	meta.Signs = make(aux.Set3D)
	meta.Refs = make(aux.Set3D)
	meta.Cnfs = make(aux.Set3D)
	return meta
}

// TODO make it work
func (meta Meta) New() {
	meta.Signs = make(aux.Set3D)
	meta.Refs = make(aux.Set3D)
	meta.Cnfs = make(aux.Set3D)
}

type Dat4rdf struct {
	Udat *aux.Set3D
	Txns *aux.Set3D
	Upac *aux.Set3D
	Upca *aux.Set2D
	Gnm  *aux.Set3D
}

type Xmap struct {
	Bgwg  aux.Set3D
	Bgwp  aux.Set3D
	Upac  aux.Set3D
	Gsymb aux.Set3D
	Ensg  aux.Set3D
	Ncbig aux.Set3D
	Ensp  aux.Set3D
	Rfsq  aux.Set3D
}

func NewXmap() (xmap Xmap) {
	xmap.Upac = make(aux.Set3D)
	xmap.Gsymb = make(aux.Set3D)
	xmap.Bgwp = make(aux.Set3D)
	xmap.Bgwg = make(aux.Set3D)
	xmap.Ensg = make(aux.Set3D)
	xmap.Ncbig = make(aux.Set3D)
	xmap.Ensp = make(aux.Set3D)
	xmap.Rfsq = make(aux.Set3D)
	return xmap
}

func (xmap Xmap) Unmarshal(pthj string) error {
	jdat, err := ioutil.ReadFile(pthj) // []bite
	if err != nil {
		// should be in main funcs instead
		err = fmt.Errorf("%s%s", "util.Xmap.Unmarshal:outil.ReadFile: ", err)
		return err
	}
	if err = json.Unmarshal(jdat, &xmap); err != nil {
		err = fmt.Errorf("%s%s", "util.Xmap.Unmarshal:json.Unmarshal: ", err)
		return err
	}
	return nil
}

