package rdf

import (
	"github.com/vlmir/bgw3/src/utils" // pkg 'aux'
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"reflect"
	"strings"
)

type Zeno struct {
	Opys aux.SliceSet
	Apys aux.SliceSet
	Prns aux.SliceSet
	Uris map[string]string
}

/// Functions
func NewZeno() (zeno Zeno) {
	zeno.Opys = make(aux.SliceSet)
	zeno.Apys = make(aux.SliceSet)
	zeno.Prns = make(aux.SliceSet)
	zeno.Uris = make(map[string]string)
	return zeno
}

func FormU(u string) string {
	return strings.Join([]string{"<", u, ">"}, "")
}
func CompU(ns string, ext string) string {
	return strings.Join([]string{"<", ns, ext, ">"}, "")
}
func FormT(s string, p string, o string) string {
	return strings.Join([]string{s, p, o, ".\n"}, " ")
}
func FormL(l string) string {
	return strings.Join([]string{`"`, l, `"`}, "")
}

func Header(pitUs map[string]string, rdfmap aux.SliceSet, zeno Zeno) (string, int) {
	// pitUs returned via arg
	var sb strings.Builder
	//zvs := reflect.ValueOf(zeno) // TODO fix
	zfs := reflect.TypeOf(zeno)
	var mykeys []string
	var pdc []string
	dic := make(aux.SliceSet)
	prns := zeno.Prns // aux.SliceSet
	opys := zeno.Opys // aux.SliceSet
	apys := zeno.Apys // aux.SliceSet
	uris := zeno.Uris // map[string]string
	var rdfs string
	nln := 0
	for i := 0; i < zfs.NumField(); i++ {
		field := zfs.Field(i)
		//value := zvs.Field(i) // (type reflect.Value) must be converted to aux.SliceSet
		//fmt.Println(aux.SliceSet(value))
		key := field.Name
		mykeys = rdfmap[key] // empty slice for Uris, no complaints
		switch {
		case key == "Prns":
			pdc = opys["sub2cls"]
			dic = prns
			rdfs = "Class"
		case key == "Opys":
			pdc = opys["ppy2prn"]
			dic = opys
			rdfs = "ObjectProperty"
		case key == "Apys":
			pdc = opys["ppy2prn"]
			dic = apys
			rdfs = "AnnotationProperty"
		}
		pU := CompU(uris[pdc[0]], pdc[1])
		lblpdc := apys["sth2lbl"] // []string
		plU := CompU(uris[lblpdc[0]], lblpdc[1])
		for _, mykey := range mykeys {
			item := dic[mykey] // []string
			sU := CompU(uris[item[0]], item[1])
			pitUs[mykey] = sU
			oU := CompU(uris["rdfs"], rdfs)
			sb.WriteString(FormT(sU, pU, oU))
			nln++
			sb.WriteString(FormT(sU, plU, FormL(item[2])))
			nln++
		}
	}
	return sb.String(), nln
}

func (zeno Zeno) Unmarshal(jpth string) error {
	zdat, err := ioutil.ReadFile(jpth) // []bite
	if err != nil {
		err = fmt.Errorf("%s%s%s", "rdf.Unmarshal:ioutil.ReadFile:", jpth, err)
		return err
	}
	if len(zdat) == 0 {
		log.Fatalln("rdf.Unmarshal:No data in json:", jpth)
	}
	if err = json.Unmarshal(zdat, &zeno); err != nil {
		log.Fatalln("rdf.Unmarshal:Failed to unmarshal:", err)
		return err
	}
	return nil
}
