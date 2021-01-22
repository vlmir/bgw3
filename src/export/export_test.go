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
	xpth := pth + "output/"
	var idmkeys = map[string]string{
		"Ensembl_PRO":   "ensp",
		"Ensembl":       "ensg",
		"EnsemblGenome": "ensom",
		"GeneID":        "ncbig",
		"RefSeq":        "rfsq",
		"UniParc":       "uparc",
	}
	txn2prm := make(util.Set2D)
	txn2prm.Add("9606", "UP000005640")
	txn2prm.Add("7227", "UP000000803")
	upacs, _ := parse.UpIdMap(pth+"idmapping/UP000005640_9606.idmapping", idmkeys)
	arg01, _ := parse.UpTab(pth+"uniprot/9606.upt", upacs, txn2prm)
	arg01.Upac = &upacs
	arg02 := xpth + "gene/export0.nt"
	arg03 := xpth + "prot/export0.nt"
	arg04 := xpth + "xmap/export0.json"
	tts := []tt{
		//		{arg01, arg02, arg03, arg04, 35, 76},
		{arg01, arg02, arg03, arg04, 34, 52},
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
	arg12 := xpth + "gene/export1.nt"
	arg13 := xpth + "prot/export1.nt"
	arg14 := xpth + "xmap/export1.json"
	tts = []tt{
		//		{arg11, arg12, arg13, arg14, 14, 33},
		{arg11, arg12, arg13, arg14, 19, 34},
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

func Test_Tfac2gene(t *testing.T) {
	type tt struct {
		arg1 util.Set4D
		arg2 bgw.Xmap
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	xpth := pth + "output/"
	arg01 := make(util.Set4D)
	arg01.Add("TP53--TP53", "uniprot", "P04637", "tfacts")
	arg01.Add("TP53--TP53", "ncbig", "7157", "tfacts")
	arg01.Add("TP53--TP53", "tfacts", "pubmed", "P04637")
	arg01.Add("TP53--TP53", "tfacts", "pubmed", "04637")
	arg01.Add("TP53--TP53", "tfacts", "confidence", "High")
	arg01.Add("TP53--TP53", "tfacts", "mode", "UP")
	arg01.Add("TP53--TP53", "tfacts", "mode", "Unknown")

	arg02 := bgw.NewXmap()
	arg02.Upac.Add("P04637", "bgwp", "9606/P04637")
	arg02.Lblg.Add("TP53", "bgwg", "9606/TP53")
	arg02.Lblg.Add("7157", "bgwg", "9606/TP53")
	arg03 := xpth + "tfac2gene/export.nt"
	tts := []tt{
		{arg01, arg02, arg03, 20},
	}
	for i, tt := range tts {
		n, _ := Tfac2gene(tt.arg1, tt.arg2, tt.arg3)
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ", len(tt.arg1), tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", n,
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
	xpth := pth + "output/"
	arg3 := make(util.Set3D)
	arg3.Add("TP53", "bgwg", "9606/TP53")
	arg1, _ := parse.UpVar(pth+"test.var", arg3)
	arg4 := xpth + "gene2phen/export.nt"
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
	xpth := pth + "output/"
	arg2 := make(util.Set3D)
	arg2.Add("P04637", "bgwp", "9606/P04637")
	arg1, _ := parse.MiTab(pth+"test.mit", arg2)
	arg3 := xpth + "prot2prot/export.nt"
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
	xpth := pth + "output/"
	arg2 := make(util.Set3D)
	arg2.Add("P04637", "bgwp", "9606/P04637")
	bps, ccs, mfs, _ := parse.Gaf(pth+"test.gaf", arg2)
	out := [3]string{"prot2bp/export.nt", "prot2cc/export.nt", "prot2mf/export.nt"}
	tts := []tt{
		{bps, arg2, xpth + out[0], 1961},
		{ccs, arg2, xpth + out[1], 260},
		{mfs, arg2, xpth + out[2], 853},
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
	xpth := pth + "output/"
	arg1 := make(util.Set3D)
	arg2 := make(util.Set3D)
	arg1.Add("uniprot!P02340--uniprot!P04637", "KO", "K04451")
	arg1.Add("uniprot!P02340--uniprot!P04637", "OrthoDB", "257530at2759")
	// arg2.Add("P04637", "bgwp", "9606/P04637#UPI000002ED67")
	arg2.Add("P04637", "bgwp", "9606/P04637")
	// arg2.Add("P02340", "bgwp", "10090/P02340#PI00000002B3")
	arg2.Add("P02340", "bgwp", "10090/P02340")
	arg3 := xpth + "ortho/export.nt"
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
