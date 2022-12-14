package export

import (
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/util"
	"testing"
)

// TODO use Test_Tfac2gene as a paradigm

func Test_GeneProt(t *testing.T) {
	// TODO implement properly without duplicating the tests
	type tt struct {
		arg1 bgw.Dat4rdf
		arg2 string
		arg3 string
		arg4 string
		val1 int
		val2 int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	var idmkeys = map[string]string{
		"Ensembl_PRO":   "ensprotein",
		"Ensembl":       "ensgene",
		"EnsemblGenome": "ensom",
		"GeneID":        "ncbigene",
		"RefSeq":        "refseq",
		"UniParc":       "uniparc",
	}
	txn2prm := make(util.Set2D)
	txn2prm.Add("9606", "UP000005640")
	txn2prm.Add("7227", "UP000000803")
	upacs, _ := parse.UpIdMap(pth+"idmapping/UP000005640_9606.idmapping", idmkeys)
	arg01, _ := parse.UpTab(pth+"uniprot/9606.upt", upacs, txn2prm)
	arg01.Upac = &upacs
	arg02 := wpth + "gene/9606.nt"
	arg03 := wpth + "prot/9606.nt"
	arg04 := wpth + "xmap/9606.json"
	tts := []tt{
		//		{arg01, arg02, arg03, arg04, 35, 76},
		{arg01, arg02, arg03, arg04, 50, 68},
	}
	for i, tt := range tts {
		n, m, _ := GeneProt(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
		if m != tt.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val2,
				"\n\thave", m,
			)
		}
	}
	upacs, _ = parse.UpIdMap(pth+"idmapping/UP000000803_7227.idmapping", idmkeys)
	arg11, _ := parse.UpTab(pth+"uniprot/7227.upt", upacs, txn2prm)
	arg11.Upac = &upacs
	arg12 := wpth + "gene/7227.nt"
	arg13 := wpth + "prot/7227.nt"
	arg14 := wpth + "xmap/7227.json"
	tts = []tt{
		//		{arg11, arg12, arg13, arg14, 14, 33},
		{arg11, arg12, arg13, arg14, 21, 25},
	}
	for i, tt := range tts {
		n, m, _ := GeneProt(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
		if m != tt.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val2,
				"\n\thave", m,
			)
		}
	}
}


func TestSigPways(t *testing.T) {
	type tt struct {
		arg1 *bgw.Dat4bridge
		arg2 *bgw.Xmap
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	// parsing
	keys, vals := bgw.SigPwaysParseConf()
	var d4b0 bgw.Dat4bridge
	d4b0.New()
	_ = parse.Tab2struct(pth+"signor/sigpcrels.tsv", keys, vals, &d4b0)
	xmap := bgw.NewXmap()
	xmap.Upac.Add("O15393", "bgwp", "9606/O15393")
	xmap.Upac.Add("P14210", "bgwp", "9606/P14210")
	xmap.Upac.Add("P08581", "bgwp", "9606/P08581")
	xmap.Upac.Add("P10275", "bgwp", "9606/P10275")
	xmap.Upac.Add("P31749", "bgwp", "9606/P31749") // PF24
	xmap.Upac.Add("P31751", "bgwp", "9606/P31751") // PF24
	/// exporting
	srcs := []string{"signor",}

	// for Signor complexes and protein families
	sigmap := make(util.Set3D)
	subdir := "signor/"
	dpth := pth + subdir
	ss0 := []string{dpth + "complexes.map", dpth + "families.map"}
	parse.Sig2up(sigmap, ss0)
	xmap.Signor = sigmap

	tts := []tt{
		{&d4b0, &xmap, pth + "OUT/export/", 5}, // 3 if filtered by the host
	}
	pdck := "reg2utrg"
	for i, tt := range tts {
		(*tt.arg1).Src = srcs[i]
		(*tt.arg1).Taxid = "9606"
		SigPways(tt.arg1, tt.arg2, tt.arg3)
		cnts := (*tt.arg1).Cnts
		if cnts[pdck][srcs[i]] != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", cnts[pdck][srcs[i]],
			)
		}
	}
} // TestSigPways

func TestRgr2trg(t *testing.T) {
	type tt struct {
		arg1 *bgw.Dat4bridge
		arg2 *bgw.Xmap
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	// parsing
	keys, vals := bgw.SignorParseConf()
	var d4b0 bgw.Dat4bridge
	d4b0.New()
	_ = parse.Tab2struct(pth+"signor/9606.mi28", keys, vals, &d4b0)
	xmap := bgw.NewXmap() //TODO see which entries are really necessary
	xmap.Upac.Add("P27361", "bgwp", "9606/P27361")
	xmap.Upac.Add("P48431", "bgwp", "9606/P48431")
	xmap.Upac.Add("Q9BTC0", "bgwp", "9606/Q9BTC0")
	xmap.Upac.Add("P08648", "bgwp", "9606/P08648")
	xmap.Bgwp.Add("9606/P08648", "bgwg", "9606/GENEX")
	xmap.Upac.Add("P10275", "bgwp", "9606/P10275")
	xmap.Upac.Add("P19838", "bgwp", "9606/P19838")
	xmap.Upac.Add("Q04206", "bgwp", "9606/Q04206")
	xmap.Upac.Add("P24385", "bgwp", "9606/P24385")
	xmap.Upac.Add("P04637", "bgwp", "9606/P04637")
	xmap.Upac.Add("Q01081", "bgwp", "9606/Q01081")
	/// exporting
	srcs := []string{"signor",}

	// for Signor complexes and protein families
	sigmap := make(util.Set3D)
	subdir := "signor/"
	dpth := pth + subdir
	ss0 := []string{dpth + "complexes.map", dpth + "families.map"}
	parse.Sig2up(sigmap, ss0)
	xmap.Signor = sigmap

	tts := []tt{
		{&d4b0, &xmap, pth + "OUT/export/", 2},
	}
	pdck := "reg2ptrg"
	for i, tt := range tts {
		(*tt.arg1).Src = srcs[i]
		(*tt.arg1).Taxid = "9606"
		Rgr2trg(tt.arg1, tt.arg2, tt.arg3)
		cnts := (*tt.arg1).Cnts
		if cnts[pdck][srcs[i]] != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", cnts[pdck][srcs[i]],
			)
		}
	}
} // TestRgr2trg

func Test_Tfac2gene(t *testing.T) {
	type tt struct {
		arg1 *bgw.Dat4bridge
		arg2 *bgw.Xmap
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	// parsing
	var d4b0 bgw.Dat4bridge
	d4b0.New()
	keys, vals := bgw.TftgParseConf()
	_ = parse.Tab2struct(pth+"static/tfacts/9606.f2g", keys, vals, &d4b0)
	var d4b1 bgw.Dat4bridge
	d4b1.New()
	keys, vals = bgw.TftgParseConf()
	_ = parse.Tab2struct(pth+"static/ntnu/9606.f2g", keys, vals, &d4b1)
	xmap := bgw.NewXmap()
	xmap.Upac.Add("P01100", "bgwp", "9606/P01100")
	xmap.Upac.Add("P04637", "bgwp", "9606/P04637")
	xmap.Ncbig.Add("4322", "bgwg", "9606/MMP13")
	xmap.Ncbig.Add("7157", "bgwg", "9606/TP53")
	xmap.Bgwg.Add("9606/MMP13", "bgwp", "9606/P01100")
	xmap.Bgwg.Add("9606/TP53", "bgwp", "9606/P04637")
	/// exporting
	srcs := []string{"tfacts", "ntnu"}

	tts := []tt{
		{&d4b0, &xmap, pth + "OUT/export/", 2},
		{&d4b1, &xmap, pth + "OUT/export/", 1},
	}
	pdck := "reg2ptrg"
	for i, tt := range tts {
		(*tt.arg1).Src = srcs[i]
		(*tt.arg1).Taxid = "9606"
		Tfac2gene(tt.arg1, tt.arg2, tt.arg3)
		cnts := (*tt.arg1).Cnts
		if cnts[pdck][srcs[i]] != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", cnts[pdck][srcs[i]],
			)
		}
	}
}

func Test_UpVar(t *testing.T) {
	type tt struct {
		arg1 util.Set3D
		arg3 util.Set3D
		arg4 string
		val1 int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	arg3 := make(util.Set3D)
	arg3.Add("TP53", "bgwg", "9606/TP53")
	arg1, _ := parse.UpVar(pth + "uniprot/P04637.var")
	arg4 := wpth + "gene2phen/9606.nt"
	t2s := []tt{
		{arg1, arg3, arg4, 11},
	}
	for i, tt := range t2s {
		n, err := Gene2phen(tt.arg1, tt.arg3, tt.arg4)
		if err != nil {
			return
		}
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ", len(tt.arg1), tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
	}
}

func Test_MiTab(t *testing.T) {
	type tt struct {
		arg1 util.Set3D
		arg2 util.Set3D
		arg3 string
		val1 int
	}

	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	arg2 := make(util.Set3D)
	arg2.Add("P04637", "bgwp", "9606/P04637")
	arg1, _ := parse.MiTab(pth+"intact/9606.mit", arg2)
	arg3 := wpth + "prot2prot/9606.nt"
	tts := []tt{
		{arg1, arg2, arg3, 73},
	}
	for i, tt := range tts {
		n, err := Prot2prot(tt.arg1, tt.arg2, tt.arg3)
		if err != nil {
			return
		}
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ", len(tt.arg1), tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
	}
}

func Test_Prot2go(t *testing.T) {
	type tt struct {
		arg1 util.Set3D
		arg2 util.Set3D
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	arg2 := make(util.Set3D)
	arg2.Add("P04637", "bgwp", "9606/P04637")
	bps, ccs, mfs, _ := parse.Gaf(pth+"goa/9606.gaf", arg2)
	out := [3]string{"prot2bp/9606.nt", "prot2cc/9606.nt", "prot2mf/9606.nt"}
	tts := []tt{
		{bps, arg2, wpth + out[0], 1961},
		{ccs, arg2, wpth + out[1], 260},
		{mfs, arg2, wpth + out[2], 853},
	}
	for i, tt := range tts {
		n, err := Prot2go(tt.arg1, tt.arg2, tt.arg3)
		if err != nil {
			return
		}
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ", len(tt.arg1), tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
	}
}

func Test_Ortho(t *testing.T) {
	type tt struct {
		arg1 util.Set3D
		arg2 util.Set3D
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	arg1 := make(util.Set3D)
	arg2 := make(util.Set3D)
	arg1.Add("uniprot!P02340--uniprot!P04637", "KO", "K04451")
	arg1.Add("uniprot!P02340--uniprot!P04637", "OrthoDB", "257530at2759")
	// arg2.Add("P04637", "bgwp", "9606/P04637#UPI000002ED67")
	arg2.Add("P04637", "bgwp", "9606/P04637")
	// arg2.Add("P02340", "bgwp", "10090/P02340#PI00000002B3")
	arg2.Add("P02340", "bgwp", "10090/P02340")
	arg3 := wpth + "ortho/10090-9606.nt"
	tts := []tt{
		{arg1, arg2, arg3, 12},
	}
	for i, tt := range tts {
		n, err := Ortho(tt.arg1, tt.arg2, tt.arg3)
		if err != nil {
			return
		}
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ", len(tt.arg1), tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
	}
}
