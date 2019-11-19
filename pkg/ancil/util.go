package util

import (
	"github.com/vlmir/bgw3/pkg/utils" // pkg 'aux'
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"reflect"
)

type Meta struct {
	Signs aux.Set3D
	Refs  aux.Set3D
	Cnfs  aux.Set3D
}

/// Functions
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

func (meta Meta) Fields() (fns []string) {
	t := reflect.TypeOf(meta)
	cnt := t.NumField()
	for i := 0; i < cnt; i++ {
		tf := t.Field(i)
		n := tf.Name
		fns = append(fns, n)
	}
	return fns
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

/// Functions
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

func (xmap Xmap) Fields() {
	t := reflect.TypeOf(xmap)
	v := reflect.ValueOf(xmap)
	for i := 0; i < t.NumField(); i++ {
		tf := t.Field(i)
		vf := v.Field(i)
		n := tf.Name
		fmt.Println(i, n, vf)
	}
}

func (xmap Xmap) Counts() {
	t := reflect.TypeOf(xmap)
	for i := 0; i < t.NumField(); i++ {
		tf := t.Field(i)
		n := tf.Name
		switch n {
		case "Bgwg":
			log.Println(n, len(xmap.Bgwg))
		case "Bgwp":
			log.Println(n, len(xmap.Bgwp))
		case "Upac":
			log.Println(n, len(xmap.Upac))
		case "Gsymb":
			log.Println(n, len(xmap.Gsymb))
		case "Ensg":
			log.Println(n, len(xmap.Ensg))
		case "Ensp":
			log.Println(n, len(xmap.Ensp))
		case "Ncbig":
			log.Println(n, len(xmap.Ncbig))
		case "Rfsq":
			log.Println(n, len(xmap.Rfsq))
		}
	}
}
