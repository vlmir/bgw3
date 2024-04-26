package rdf

import (
	"errors"
	"fmt"
	"github.com/vlmir/bgw3/src/util"
	"strings"
)

var Nss = map[string]string{
	"rdf":         "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
	"rdfs":        "http://www.w3.org/2000/01/rdf-schema#",
	"owl":         "http://www.w3.org/2002/07/owl#",
	"skos":        "http://www.w3.org/2004/02/skos/core#",
	"schema":      "http://schema.org/",
	"sd":          "http://www.w3.org/ns/sparql-service-description#", // Dataset, Graph, NamedGraph etc.
	"biolink":     "https://w3id.org/biolink/vocab/",                  // TODO implement
	"obo":         "http://purl.obolibrary.org/obo/",
	"sio":         "http://semanticscience.org/resource/",
	"uniprot":     "http://uniprot.org/uniprot/", // accepts both UPID and UPAC
	"uniprotkb":   "http://uniprot.org/uniprot/", // accepts both UPID and UPAC
	"uniparc":     "http://uniprot.org/uniparc/",
	"ncbigene":    "http://identifiers.org/ncbigene/",
	"ncbitx":      "http://purl.obolibrary.org/obo/NCBITaxon_",
	"refseq":      "http://identifiers.org/refseq/",
	"keggortho":   "http://identifiers.org/kegg.orthology/",
	"orthodb":     "https://www.orthodb.org/?query=", // acepts IDs from UniProt ID mappings
	"omim":        "http://purl.bioontology.org/ontology/OMIM/",
	"intact":      "http://identifiers.org/intact/",
	"goa":         "http://identifiers.org/goa/",
	"pubmed":      "http://identifiers.org/pubmed/",
	"signor":      "https://signor.uniroma2.it/relation_result.php?id=",
	"enspro":      "http://identifiers.org/ensembl/",
	"ensgene":     "http://identifiers.org/ensembl/",
	"ensplants":   "https://plants.ensembl.org/id/",
	"ensfungi":    "https://fungi.ensembl.org/id/",
	"ensmetazoa":  "https://metazoa.ensembl.org/id/",
	"ensprotists": "https://protists.ensembl.org/id/",
	"sigpway":     "https://signor.uniroma2.it/getPathwayData.php?pathway=",
	"gene":        "http://rdf.biogateway.eu/gene/",
	// "prot":        "http://rdf.biogateway.eu/prot/",
}

// Object Properties
var Opys = util.SliceSet{
	"sub2cls": {"rdfs", "subClassOf", "is subclass of"},
	"sub2ppy": {"rdfs", "subPropertyOf", "is subproperty of"},
	"sth2evd": {"sio", "SIO_000772", "has evidence"},                                                                                          // PubMed only: ALL
	"sth2ori": {"schema", "evidenceOrigin", "has evidence origin"},                                                                            // object: an external entity, currently used only by Juan
	"sth2src": {"sio", "SIO_000253", "has source", "has source is a relation between an entity and another entity from which it stems from."}, // is this the right property for linking instances to databases? TODO
	//"sth2src":   {"biolink", "provided_by", "is provided by"},
	//"sth2src":   {"biolink", "supporting_data_source", "has supporting data source"}, //  better choice, used already by Juan
	"sth2eqv": {"owl", "sameAs", "is equivalent to"},
	"sth2clm": {"skos", "closeMatch", "has close match"},
	"sth2rlm": {"skos", "relatedMatch", "has related match"},
	"sub2set": {"obo", "BFO_0000050", "is part of"},
	//"gn2txn":  {"obo", "RO_0000052", "characteristic of"},
	//"gp2txn":  {"obo", "RO_0000052", "characteristic of"},
	"be2txn": {"obo", "RO_0002162", "in taxon"}, // confusing definition
	//"be2txn":	{"biolink", "in_taxon", "is characteristic of taxon"} // declared equivalent to RO_0002162
	"mbr2lst": {"schema", "memberOf", "is member of"},
	"ins2cls": {"rdf", "type", "is instance of"},
	"stm2sbj": {"rdf", "subject", "has subject"},
	"stm2obj": {"rdf", "object", "has object"},
	"stm2pdc": {"rdf", "predicate", "has predicate"},
	"gn2phn":  {"obo", "RO_0002331", "involved in"},
	"gp2bp":   {"obo", "RO_0002331", "involved in", "biological_process"},
	"gp2cc":   {"obo", "BFO_0000050", "is part of", "cellular_component"},
	"gp2mf":   {"obo", "RO_0002327", "enables", "molecular_function"},
	"gn2gp":   {"sio", "SIO_010078", "encodes"},
	//"gn2gp":     {"biolink", "has_gene_product", "has gene product"}, // TODO
	"tlp2tlp":   {"obo", "RO_0002436", "molecularly interacts with"},
	"reg2utrg":  {"obo", "RO_0002428", "involved in regulation of"},
	"reg2ptrg":  {"obo", "RO_0002429", "involved in positive regulation of"},
	"reg2ntrg":  {"obo", "RO_0002430", "involved in negative regulation of"},
	"sth2mtd":   {"rdfs", "isDefinedBy", "is defined by"},
	"orl2orl":   {"sio", "SIO_000558", "is orthologous to"},
	"mi2bp":     {"obo", "RO_0002331", "involved in"},
	"mi2loc":    {"obo", "BFO_0000066", "occurs in"}, // sigpways only
	"step2pway": {"obo", "BFO_0000050", "is part of"},
	"reg2dtrg":  {"obo", "RO_0002578", "directly regulates"},
	"reg2itrg":  {"obo", "RO_0012012", "indirectly regulates"},
	// not used yet SIO_000208"
	// "evd4sth": {"sio", "SIO_000208", "is supporting evidence for"},
}

// Datatype Properties
var Dpys = util.SliceSet{
	"gi2^": {"obo", "GENO_0000894", "start position"},
	"gi2#": {"obo", "GENO_0000895", "end position"},
}

// Annotation Properties
var Apys = util.SliceSet{
	"sth2vrs": {"owl", "versionInfo", "current version"},
	"sth2lbl": {"skos", "prefLabel", "has name"},
	// only in entities
	"sth2dfn": {"skos", "definition", "has definition"},          // in gene prot
	"sth2syn": {"skos", "altLabel", "has synonym"},               // in gene prot
	"evd2lvl": {"schema", "evidenceLevel", "has evidence level"}, // prot tfac2gene prot2prot
	// only in bridges
	// "sth2val": {"rdf", "value", "has value"}, // only tfac2gene (positive|negative)
	//"sth2cmt": {"rdfs", "comment", "has comment"}, // gene2phen?
	// to be used
	"sth2id": {"skos", "notation", "has notation"}, // TODO
}

// Parental classes
var Prns = util.SliceSet{
	"cls": {"rdfs", "Class", "class"},
	//"opy": {"rdfs", "ObjectProperty", "object property"},
	//"apy": {"rdfs", "AnnotationProperty", "annotation property"},
	// "bag": {"rdf", "Bag", "unordered collection"},
	"stm": {"rdf", "Statement", "Triple"},
	"tlp": {"sio", "SIO_010043", "Protein"},
	"gn":  {"sio", "SIO_010035", "Gene"},
	// "chr": {"obo", "GO_0005694", "Chromosome"},
	"chr": {"obo", "SO_0000340", "Chromosome"},
	// "gom": {"obo", "SO_0001026", "Genome"},
}

var Uris4rgrtrg = map[string]string{
	"signor": "http://signor.uniroma2.it",
}

var Uris4tftg = map[string]string{
	//	"extri":  "http://www.extri.org",
	//	"htri":   "http://www.lbbc.ibb.unesp.br/htri",
	//	"intact": "https://www.ebi.ac.uk/intact/interactomes",
	//	"ntnu":   "NTNU",
	//	"tfacts": "http://www.tfacts.org",
	//	"trrust": "http://www.grnpedia.org/trrust",
	//	"cytreg": "https://www.fuxmanlab.com",
	//	"geredb": "http://www.thua45.cn/geredb",
	//	"signor": "http://signor.uniroma2.it",
	"tflink": "https://tflink.net",
	"coltri": "https://github.com/saezlab/CollecTRI",
}

/// Functions

// he 6 functions below panic if util.TrimString() returns an error
func FormL(s0 string) string {
	// used only in FormT() background
	p0 := &s0
	err := util.TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("rdf.FormL(%s): %s", "s0", err)
		panic(errors.New(msg))
	}
	return strings.Join([]string{`"`, *p0, `"`}, "")
}

func FormU(s0 string) string {
	p0 := &s0
	err := util.TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("rdf.FormU(%s): %s", "s0", err)
		panic(errors.New(msg))
	}
	return strings.Join([]string{"<", *p0, ">"}, "")
}

func CompU(s0, s1 string) string {
	p0 := &s0
	err := util.TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("rdf.CompU(%s, _): %s", "s0", err)
		panic(errors.New(msg))
	}
	p1 := &s1
	err = util.TrimString(p1)
	if err != nil {
		msg := fmt.Sprintf("rdf.CompU(_, %s): %s", "s1", err)
		panic(errors.New(msg))
	}
	return strings.Join([]string{"<", *p0, *p1, ">"}, "")
}

func FormT(s0, s1, s2 string) string {
	p0 := &s0
	err := util.TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("rdf.FormT(%s, _, _): %s", "s0", err)
		panic(errors.New(msg))
	}
	p1 := &s1
	err = util.TrimString(p1)
	if err != nil {
		msg := fmt.Sprintf("rdf.FormT(_, %s, _): %s", "s1", err)
		panic(errors.New(msg))
	}
	p2 := &s2
	err = util.TrimString(p2)
	if err != nil {
		msg := fmt.Sprintf("rdf.FormT(_, _, %s): %s", "s2", err)
		panic(errors.New(msg))
	}
	return strings.Join([]string{*p0, *p1, *p2, ".\n"}, " ")
}

// arg1: a map for filtering
func FmtURIs(rdfmap util.SliceSet) map[string]string {
	set := make(util.SliceSet)
	fmtURIs := make(map[string]string)
	for group, urikeys := range rdfmap { // urikeys - slice of tokens
		if len(urikeys) == 0 {
			msg := fmt.Sprintf("rdf.FmtURIs(): NoKeys")
			panic(errors.New(msg))
		}
		switch {
		case group == "Prns":
			set = Prns
		case group == "Opys":
			set = Opys
		case group == "Apys":
			set = Apys
		}
		for _, urikey := range urikeys {
			bits, ok := set[urikey] // []string{nsk, uid}
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
} // FmtURIs

// arg1: a map for filtering
// return1: a string of RDF triples in the 'nt' format`
// return2: the number of lines in return1
func Capita(rdfmap util.SliceSet) (string, int) {
	ss := make(util.SliceSet)
	var pdc []string
	var top string
	var plU string
	var sb strings.Builder
	nln := 0
	groups := rdfmap.Keys() // sorted, non-empty
	for _, group := range groups {
		urikeys := rdfmap[group]
		switch {
		case group == "Prns":
			ss = Prns
			pdc = Opys["ins2cls"]
			top = "Class"
		case group == "Opys":
			ss = Opys
			pdc = Opys["ins2cls"]
			top = "ObjectProperty"
		case group == "Apys":
			ss = Apys
			pdc = Opys["ins2cls"]
			top = "AnnotationProperty"
		}
		pU := CompU(Nss[pdc[0]], pdc[1])
		lbits := Apys["sth2lbl"]
		plU = CompU(Nss[lbits[0]], lbits[1])
		for _, urikey := range urikeys {
			bits, ok := ss[urikey] // []string
			if !ok {
				msg := fmt.Sprintf("NoEntryFor: %s", urikey)
				panic(errors.New(msg))
			}
			sU := CompU(Nss[bits[0]], bits[1])
			oU := CompU(Nss["owl"], top)
			sb.WriteString(FormT(sU, pU, oU))
			nln++
			sb.WriteString(FormT(sU, plU, FormL(bits[2])))
			nln++
		}
	}
	return sb.String(), nln
} // Capita
