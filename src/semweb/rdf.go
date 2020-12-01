package rdf

import (
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/util"
	"strings"
)

var Nss = map[string]string{
	"rdf":       "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
	"rdfs":      "http://www.w3.org/2000/01/rdf-schema#",
	"owl":       "http://www.w3.org/2002/07/owl#",
	"skos":      "http://www.w3.org/2004/02/skos/core#",
	"schema":    "http://schema.org/",
	"obo":       "http://purl.obolibrary.org/obo/",
	"sio":       "http://semanticscience.org/resource/", // resolvable (2014-09-29)
	"uniprot":   "http://uniprot.org/uniprot/",          // accepts both UPID and UPAC
	"uparc":     "http://identifiers.org/uniparc/",
	"ncbig":     "http://identifiers.org/ncbigene/",
	"ncbitx":    "http://purl.bioontology.org/ontology/NCBITAXON/",
	"rfsq":      "http://identifiers.org/refseq/",
	"keggortho": "http://identifiers.org/kegg.orthology/",
	"orthodb":   "https://www.orthodb.org/",
	"omim":      "http://purl.bioontology.org/ontology/OMIM/",
	"intact":    "http://identifiers.org/intact/",
	"goa":       "http://identifiers.org/goa/",
	"pubmed":    "http://identifiers.org/pubmed/",
	"ensp":      "http://identifiers.org/ensembl/",
	"ensg":      "http://identifiers.org/ensembl/",
	"enspl":      "https://plants.ensembl.org/id/",
	"ensfu":      "https://fungi.ensembl.org/id/",
	"ensme":      "https://metazoa.ensembl.org/id/",
	"enspr":      "https://protists.ensembl.org/id/",
	"go":        "http://purl.obolibrary.org/obo/GO_",
	"bgw":       "http://rdf.biogateway.eu/",
	"gene":      "http://rdf.biogateway.eu/gene/",
	"prot":      "http://rdf.biogateway.eu/prot/",
	// "tsf_gn":    "http://rdf.biogateway.eu/prot-gene/", // for TF-TG
	// "up_up":     "http://rdf.biogateway.eu/prot-prot/",
	// "up_obo":    "http://rdf.biogateway.eu/prot-obo/",
	// "up_mim":    "http://rdf.biogateway.eu/prot-omim/",
}

//"orthodb": "https://www.orthodb.org/?query=", // accepts IDs from UP idmapping

var Opys = util.SliceSet{
	"sub2cls": {"rdfs", "subClassOf", "is subclass of"},
	"sub2ppy": {"rdfs", "subPropertyOf", "is subproperty of"},
	"sth2evd": {"sio", "SIO_000772", "has evidence"}, // PubMed only: ALL
	"sth2ori": {"schema", "evidenceOrigin", "has evidence origin"}, // DATA sources; ALL; e.g. bgwp 2 upca
	"sth2src": {"sio", "SIO_000253", "has source", "has source is a relation between an entity and another entity from which it stems from."},
	"sth2clm": {"skos", "closeMatch", "has close match"},
	"sub2set": {"obo", "BFO_0000050", "is part of"},
	"gn2txn":  {"obo", "BFO_0000052", "inheres in"},
	"gp2txn":  {"obo", "BFO_0000052", "inheres in"},
	"mbr2lst": {"schema", "memberOf", "is member of"},
	"ins2cls": {"rdf", "type", "is instance of"},
	"gn2phn":  {"obo", "RO_0002331", "involved in"},
	"gp2bp":   {"obo", "RO_0002331", "involved in", "biological_process"},
	"gp2cc":   {"obo", "BFO_0000050", "is part of", "cellular_component"},
	"gp2mf":   {"obo", "RO_0002327", "enables", "molecular_function"},
	"gn2gp":   {"sio", "SIO_010078", "encodes"},
	"tlp2tlp": {"obo", "RO_0002436", "molecularly interacts with"},
	"rgr2trg": {"obo", "RO_0002428", "involved in regulation of"},
	"sth2mtd": {"rdfs", "isDefinedBy", "is defined by"},
	"orl2orl": {"sio", "SIO_000558", "is orthologous to"},
}
var Apys = util.SliceSet{
	"sth2vrs": {"owl", "versionInfo", "current version"},
	"sth2lbl": {"skos", "prefLabel", "has name"},
	// only in entities
	"sth2dfn": {"skos", "definition", "has definition"},          // in gene prot
	"sth2syn": {"skos", "altLabel", "has synonym"},               // in gene prot
	"evd2lvl": {"schema", "evidenceLevel", "has evidence level"}, // prot tfac2gene prot2prot
	// only in bridges
	"sth2val": {"rdf", "value", "has value"}, // only tfac2gene (positive|negative)
	//"sth2cmt": {"rdfs", "comment", "has comment"}, // gene2phen?
	// to be used
	"sth2id": {"skos", "notation", "has notation"}, // not used yet
}
var Prns = util.SliceSet{
	"cls": {"rdfs", "Class", "class"},
	//"opy": {"rdfs", "ObjectProperty", "object property"},
	//"apy": {"rdfs", "AnnotationProperty", "annotation property"},
	"bag": {"rdf", "Bag", "unordered collection"},
	"stm": {"rdf", "Statement", "triple"},
	"tlp": {"sio", "SIO_010043", "protein"},
	"gn":  {"sio", "SIO_010035", "gene"},
	"chr": {"obo", "GO_0005694", "chromosome"},
	"gom": {"obo", "SO_0001026", "genome"},
}

var Uris4tftg = map[string]string{
	"extri":  "http://www.extri.org",
	"htri":   "http://www.lbbc.ibb.unesp.br/htri",
	"signor": "http://signor.uniroma2.it",
	"tfacts": "http://www.tfacts.org",
	"trrust": "http://www.grnpedia.org/trrust/",
	"cytreg": "https://www.fuxmanlab.com",
	"geredb": "http://www.thua45.cn/geredb",
}

/// Functions

func FormU(u string) string {
	util.CheckStrings(u)
	return strings.Join([]string{"<", u, ">"}, "")
}
func CompU(ns string, ext string) string {
	util.CheckStrings(ns, ext)
	return strings.Join([]string{"<", ns, ext, ">"}, "")
}
func FormT(s string, p string, o string) string {
	util.CheckStrings(s, p, o)
	return strings.Join([]string{s, p, o, ".\n"}, " ")
}
func FormL(l string) string {
	util.CheckStrings(l)
	return strings.Join([]string{`"`, l, `"`}, "")
}

// arg1: a map for filtering
func FmtURIs(rdfmap util.SliceSet) map[string]string {
	dic := make(util.SliceSet)
	fmtURIs := make(map[string]string)
	for group, urikeys := range rdfmap { // urikeys - slice of tokens
		switch {
		case group == "Prns":
			dic = Prns
		case group == "Opys":
			dic = Opys
		case group == "Apys":
			dic = Apys
		}
		for _, urikey := range urikeys {
			bits, ok := dic[urikey] // []string{nsk, uid}
			if !ok {
				msg := fmt.Sprintf("NoEntryFor: %s", urikey)
				panic(errors.New(msg))
			}
			if len(bits) < 2 {
				msg := fmt.Sprintf("Want at least 2 elements in: %v", bits)
				panic(errors.New(msg))
			}
			nsk := bits[0]
			uid := bits[1]
			ns, ok := Nss[nsk]
			if !ok {
				msg := fmt.Sprintf("NoEntryFor: %s", nsk)
				panic(errors.New(msg))
			}
			myU := CompU(ns, uid)
			util.CheckStrings(myU)
			fmtURIs[urikey] = myU
		}
	}
	n := 0
	for _, v := range rdfmap {
		n += len(v)
	}
	m := len(fmtURIs)
	if m != n {
		panic(errors.New(fmt.Sprintf("Want: %d have: %d", n, m)))
	}
	return fmtURIs
}

// arg1: a map for filtering
// return1: a string of RDF triples in the 'nt' format`
// return2: the number of lines in return1
func Capita(rdfmap util.SliceSet) (string, int) {
	dic := make(util.SliceSet)
	var pdc []string
	var rdfs string
	var plU string
	var sb strings.Builder
	nln := 0
	for group, urikeys := range rdfmap {
		switch {
		case group == "Prns":
			dic = Prns
			pdc = Opys["sub2cls"]
			rdfs = "Class"
		case group == "Opys":
			dic = Opys
			pdc = Opys["sub2ppy"]
			rdfs = "ObjectProperty"
		case group == "Apys":
			dic = Apys
			pdc = Opys["sub2ppy"]
			rdfs = "AnnotationProperty"
		}
		pU := CompU(Nss[pdc[0]], pdc[1])
		lbits := Apys["sth2lbl"]
		plU = CompU(Nss[lbits[0]], lbits[1])
		for _, urikey := range urikeys {
			util.CheckStrings(urikey)
			bits, ok := dic[urikey] // []string
			if !ok {
				msg := fmt.Sprintf("NoEntryFor: %s", urikey)
				panic(errors.New(msg))
			}
			sU := CompU(Nss[bits[0]], bits[1])
			oU := CompU(Nss["rdfs"], rdfs)
			sb.WriteString(FormT(sU, pU, oU))
			nln++
			sb.WriteString(FormT(sU, plU, FormL(bits[2])))
			nln++
		}
	}
	util.CheckStrings(sb.String())
	return sb.String(), nln
}

//func Header(uris map[string]string, rdfmap util.SliceSet) (string, int) {
//	// fmtUs returned via arg group:Uri
//	var sb strings.Builder
//	k4p := ""
//	k4o := ""
//	nln := 0
//	for group, urikeys := range(rdfmap) {
//		switch {
//		case group == "Prns":
//			k4p = "sub2cls"
//			k4o = "cls"
//		case group == "Opys":
//			k4p = "sub2ppy"
//			k4o = "opy"
//		case group == "Apys":
//			k4p = "sub2ppy"
//			k4o = "apy"
//		}
//		plU := uris["sth2lbl"]
//		pU := uris[k4p]
//		for _, urikey := range urikeys {
//			// skipping the sanity checks performed in fmtURIs()
//			sU := uris[urikey]
//			oU := uris[k4o]
//			sb.WriteString(FormT(sU, pU, oU))
//			nln++
//			l4s := urikey// TODO see how to fix the values
//			sb.WriteString(FormT(sU, plU, FormL(l4s)))
//			nln++
//		}
//	}
//	if len(sb.String()) == 0 {
//		panic(errors.New(fmt.Sprintf("Empty header")))
//	}
//	return sb.String(), nln
//}
