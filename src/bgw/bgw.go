package bgw

import (
	"encoding/json"
	"github.com/vlmir/bgw3/src/util"
	"io/ioutil"
)

var CV = "3.3.0"

var Coltri = map[string]string{
	"9606":  "human",
	"10090": "mouse",
	"10116": "rat",
}

var Tflink = map[string]string{
	"6239":   "Caenorhabditis_elegans_interactions_All_mitab_v1.0.tsv.gz",
	"7955":   "Danio_rerio_interactions_All_mitab_v1.0.tsv.gz",
	"7227":   "Drosophila_melanogaster_interactions_All_mitab_v1.0.tsv.gz",
	"9606":   "Homo_sapiens_interactions_All_mitab_v1.0.tsv.gz",
	"10090":  "Mus_musculus_interactions_All_mitab_v1.0.tsv.gz",
	"10116":  "Rattus_norvegicus_interactions_All_mitab_v1.0.tsv.gz",
	"559292": "Saccharomyces_cerevisiae_interactions_All_mitab_v1.0.tsv.gz",
}

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
	"UniProtKB-ID":  "upid",
	"Gene_Name":     "gnm",
	"Gene_Synonym":  "gsnm",
	"Ensembl":       "ensgene",
	"Ensembl_PRO":   "enspro",
	"EnsemblGenome": "ensom",
	// "EnsemblGenome_PRO": "ensompro", // transcripts...
	"GeneID":     "ncbigene",
	"RefSeq":     "refseq",
	"UniParc":    "uniparc",
	"NCBI_TaxID": "ncbitx",
}

var Orthokeys = map[string]string{
	// "KO": "keggortho",
	"OrthoDB": "orthodb",
}

type Column struct {
	// TODO swap Ind2 Dlm2 ?
	Ind1 int
	Dlm1 string
	Ind2 int
	Dlm2 string
	Ind3 int
	Key  string
}

type SrcConf struct {
	Keys []Column
	Vals []Column
}

type Dat4bridge struct {
	Out   map[string]string
	Src   string
	Taxid string
	Duos  util.Set3D
	Mode  util.Set3D
	Cnts  util.Set2D
}

func (p *Dat4bridge) New() {
	d4b := *p
	d4b.Out = map[string]string{}
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

var UpdatConf = SrcConf{
	Keys: []Column{
		{0, "; ", 0, "", 0, ""}, // UniProt canonical accessionbs from multiple proteomes
	},
	Vals: []Column{
		// {0, "; ", 0, "|", 0, "upca"}, // single value
		{1, "; ", 0, "|", 0, "upid"}, // single value
		// gnms: '|' separated for CHLRE else '; ' with no '|';
		{2, "; ", -1, "|", 0, "gnms"}, // '; ' separated but some vals contain ';'!
		// gsnms: '|' separated for CHLRE; ' ' separated for DROME, '|' do occrur;
		{2, "; ", -1, "|", 0, "gsnms"}, // '; ' separated but some vals contain ';'!
		{4, "; ", 0, "|", 0, "taxnm"},  // single value but 59 unique names !
		{5, "; ", 0, "|", 0, "txid"},   // single value but 59 unique names !
		{6, "!", 0, "; ", 0, "pdfns"},  // may contain '|', ';' etc, no '!'
		// the next field may contain multiple proteomes per taxon!
		// {7, "; ", -1, ": ", 0, "poms"}, // '; ' separated, 'pomeid: chrid'
		{8, "; ", -1, "; ", 0, "pubmed"}, // '; ' separated
		{9, "; ", 0, " ", 0, "score"},    // single value
	},
}

//func TftgParseConf() ([]Column, []Column) {
//	// not used anymore
//	keys := []Column{
//		{0, ":", 0, "--", 0, ""},
//		{0, ":", 1, "--", 0, ""},
//	}
//	vals := []Column{
//		{1, "|", 0, "|", 0, "uniprot"},
//		{2, "|", 0, "|", 0, "ncbig"},
//		{3, "|", -1, ";", 0, "pubmed"},
//		{4, "|", 0, ";", 0, "score"},
//		{5, "|", 0, ";", 0, "mode"},
//	}
//	return keys, vals
//} // TftgParseConf

var ColtriConf = SrcConf{
	Keys: []Column{
		{0, ";", 0, "--", 0, ""},
		{1, ";", 0, "--", 0, ""},
	},
	Vals: []Column{
		{0, ";", 0, ":", 0, "Aupca"},       // no iso-form present 23-11-17
		{1, ";", 0, ":", 0, "Bupca"},       // no iso-form present 23-11-17
		{2, ";", 0, ":", 0, "Aglbl"},       // single values
		{3, ";", 0, ":", 0, "Bglbl"},       // single values
		{5, ";", 0, ":", 0, "pos"},         // single values, True|False
		{6, ";", 0, ":", 0, "neg"},         // single values, True|False
		{11, ";", 1, ":", -1, "CollecTRI"}, // all filtered
		{15, ";", 0, ";", 0, "score"},      // number of refs
		{16, ";", 0, ":", 0, "pubmed"},     // all refs
	},
}

var TflinkConf = SrcConf{
	Keys: []Column{
		{4, "\"", 1, "--", 0, ""},
		{5, "\"", 1, "--", 0, ""},
	},
	Vals: []Column{
		{0, "|", 1, ":", 0, "uniprot"}, // prot id; single value
		{3, "|", 1, ":", 0, "ncbig"},   // gene id; single value
		{6, "|", 1, "\"", 1, "mtdid"},  // multiple values
		{8, "|", 1, ":", -1, "pubmed"}, // filtering by "pubmed"
		// {11, "|", 1, "\"", 0, "typeABid"}, // single value, all MI:2232 (molecular association)
		// {11, "|", 3, "\"", 0, "typeABlbl"},
		// {13, "|", 1, ":", 0, "oregannoid"}, // not for all, multiple values
		// {20, "|", 1, "\"", 0, "typeAid"}, // single value, all MI:0326 (protein)
		// {20, "|", 3, "\"", 0, "typeAlbl"},
		// {21, "|", 1, "\"", 0, "typeBid"}, // single value, all MI:0250 (gene)
		// {21, "|", 3, "\"", 0, "typeBlbl"},
		{25, "|", 1, "\"", 1, "mode"},  // multiple values
		{27, "|", 1, "\"", 1, "score"}, // single values, TODO adjust ??
	},
}

var IntactConf = SrcConf{
	Keys: []Column{
		{0, ":", 1, "--", -1, "uniprotkb"}, // sorting and filtering by "uniprotkb"
		{1, ":", 1, "--", -1, "uniprotkb"}, // the last join string is defining
	},
	Vals: []Column{
		{0, "|", 1, ":", 0, "uniprotkb"}, // filtering by "uniprotkb"
		{1, "|", 1, ":", 0, "uniprotkb"}, // filtering by "uniprotkb"
		{6, "|", 1, "\"", 1, "mtd"},      // psi-mi ids
		{8, "|", 1, ":", -1, "pubmed"},
		{11, "|", 1, "\"", 0, "typeABid"},       // psi-mi id; 9606: single values
		{11, "|", 2, "\"", 0, "typeABlbl"},      // not used
		{13, "|", 1, ":", -1, "intact"},         // intact:EBI-20559053|imex:IM-26397-1
		{14, "|", 1, ":", -1, "intact-miscore"}, // author score:D|intact-miscore:0.37
	},
}

var SignorConf = SrcConf{
	Keys: []Column{
		// must be this way for complexes and families
		{0, "|", 0, "--", 0, ""},
		{1, "|", 0, "--", 0, ""},
	},
	Vals: []Column{
		{0, "|", 0, "|", 0, "Aid"},     // single value; Attn: diverse formats !!
		{1, "|", 0, "|", 0, "Bid"},     // single value; Attn: diverse formats !!
		{6, "|", 1, "\"", 1, "mtd"},    // psi-mi ids
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
	},
}

var SigPwaysConf = SrcConf{
	// "|" is not used at all
	Keys: []Column{
		{5, "|", 0, "--", 0, ""},  // Aid
		{10, "|", 0, "--", 0, ""}, // Bid
	},
	Vals: []Column{
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
	},
}

var SigMapConf = SrcConf{
	Keys:  []Column{
		{0, ";", 0, "", 0, ""},
	},
	Vals: []Column{
		{2, ";", -1, ";", 0, "ids"}, // no iso-form present 23-12-04
	},
}
