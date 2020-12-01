package export

import (
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/util"
	"testing"
)

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
		"Ensembl_PRO": "ensp",
		"Ensembl":     "ensg",
	"EnsemblGenome":     "ensom",
		"GeneID":      "ncbig",
		"RefSeq":      "rfsq",
		"UniParc":     "uparc",
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
	t1s := []tt{
		{arg01, arg02, arg03, arg04, 35, 76},
	}
	for i, tt := range t1s {
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
	t1s = []tt{
		{arg11, arg12, arg13, arg14, 14, 33},
	}
	for i, tt := range t1s {
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
		arg1 map[string]util.Set3D
		arg2 util.Set3D
		arg3 util.Set3D
		arg4 string
		val1 int
	}
	pth := "../../tdata/"
	xpth := pth + "output/"
	arg1 := make(map[string]util.Set3D)
	set1 := make(util.Set3D)
	set1.Add("TP53--TP53", "uniprot", "P04637")
	set1.Add("TP53--TP53", "ncbig", "7157")
	set1.Add("TP53--TP53", "pubmed", "P04637")
	set1.Add("TP53--TP53", "pubmed", "04637")
	set1.Add("TP53--TP53", "confidence", "High")
	set1.Add("TP53--TP53", "mode", "UP")
	set1.Add("TP53--TP53", "mode", "Unknown")
	set1.Add("tfacts", "uri", "http://www.tfacts.org")
	arg1["tfacts"] = set1
	set2 := make(util.Set3D)
	set2.Add("TP53--TP53", "uniprot", "P04637")
	set2.Add("TP53--TP53", "ncbig", "7157")
	set2.Add("TP53--TP53", "pubmed", "04637")
	set2.Add("test", "uri", "http://www.test.org")
	arg1["test"] = set2

	arg2 := make(util.Set3D)
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	arg3 := make(util.Set3D)
	arg3.Add("TP53", "bgwg", "9606/chr-17/TP53")
	arg3.Add("7157", "bgwg", "9606/chr-17/TP53")
		xmap := bgw.NewXmap()
		xmap.Upac = arg2
		xmap.Lblg = arg3
	arg4 := xpth + "tfac2gene/export.nt"
	t3s := []tt{
		{arg1, arg2, arg3, arg4, 24},
	}
	for i, tt := range t3s {
		//n, _ := Tfac2gene(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		n, _ := Tfac2gene(tt.arg1, xmap, tt.arg4)
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
	arg3.Add("TP53", "bgwg", "9606/chr-17/TP53")
	arg1, _ := parse.UpVar(pth+"test.var", arg3)
	arg4 := xpth + "gene2phen/export.nt"
	t2s := []tt{
		{arg1, arg3, arg4, 10},
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
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	arg1, _ := parse.MiTab(pth+"test.mit", arg2)
	arg3 := xpth + "prot2prot/export.nt"
	tts := []tt{
		{arg1, arg2, arg3, 72},
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
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	bps, ccs, mfs, _ := parse.Gaf(pth+"test.gaf", arg2)
	out := [3]string{"prot2bp/export.nt", "prot2cc/export.nt", "prot2mf/export.nt"}
	tts := []tt{
		{bps, arg2, xpth + out[0], 1811},
		{ccs, arg2, xpth + out[1], 243},
		{mfs, arg2, xpth + out[2], 814},
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
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	arg2.Add("P02340", "bgwp", "10090/chr-11/Tp53/UPI00000002B3")
	arg3 := xpth + "ortho/export.nt"
	tts := []tt{
		{arg1, arg2, arg3, 11},
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
