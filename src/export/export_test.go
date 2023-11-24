package export

import (
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/parse"
	"github.com/vlmir/bgw3/src/util"
	"testing"
)

// TODO use Test_Tfac2gene as a paradigm

// the 2 tests below MUST come first !
func Test_Gene(t *testing.T) {
	// TODO implement properly without duplicating the tests
	type tt struct {
		arg0 string // path to UP data
		arg1 string
		arg2 string
		arg3 *bgw.Xmap
		val1 int // number of gene names added
		val2 int // number of gene synonyms added
	}
	var xmap bgw.Xmap
	xmap.New()
	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	arg00 := pth + "uniprot/9606.upt"
	arg01 := pth + "idmapping/UP000005640_9606.idmapping"
	arg02 := wpth + "gene/9606.nt"
	arg03 := &xmap
	arg10 := pth + "uniprot/7227.upt"
	arg11 := pth + "idmapping/UP000000803_7227.idmapping"
	arg12 := wpth + "gene/7227.nt"
	arg13 := &xmap
	tts := []tt{
		// {arg00, arg01, arg02, arg03, 2, 1},
		// {arg00, arg01, arg02, arg03, 4, 1},
		{arg00, arg01, arg02, arg03, 4, 3},
		// {arg00, arg11, arg12, arg13, 3, 2},
		// {arg10, arg11, arg12, arg13, 6, 2},
		{arg10, arg11, arg12, arg13, 6, 5}, // cumulative
	}
	for i, tt := range tts {
		_ = Gene(tt.arg0, tt.arg1, tt.arg2, tt.arg3)
		n := len(xmap.Lblg.Keys())
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
		m := len(xmap.Syng.Keys())
		if m != tt.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val2,
				"\n\thave", m,
			)
		}
	}
} // Test_Gene

func Test_Prot(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 string
		arg4 *bgw.Xmap
		val1 int // number of prot names added
		val2 int // number of prot synonyms added
	}
	var xmap bgw.Xmap
	xmap.New()
	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	arg01 := pth + "uniprot/9606.upt"
	arg02 := pth + "idmapping/UP000005640_9606.idmapping"
	arg03 := wpth + "prot/9606.nt"
	arg04 := &xmap
	arg11 := pth + "uniprot/7227.upt"
	arg12 := pth + "idmapping/UP000000803_7227.idmapping"
	arg13 := wpth + "prot/7227.nt"
	arg14 := &xmap
	tts := []tt{
		{arg01, arg02, arg03, arg04, 2, 5},
		{arg11, arg12, arg13, arg14, 3, 8}, // cumulative !
	}
	for i, tt := range tts {
		_ = Prot(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		n := len(xmap.Lblp.Keys())
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
		m := len(xmap.Synp.Keys())
		if m != tt.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val2,
				"\n\thave", m,
			)
		}
	}
} // Test_Prot

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
	_ = parse.Tab2struct(pth+"signor/sigpcrels.tsv", keys, vals, &d4b0, "\t")
	var xmap bgw.Xmap
	xmap.New()
	xmap.Upac.Add("O15393", "bgwp", "O15393")
	xmap.Upac.Add("P14210", "bgwp", "P14210")
	xmap.Upac.Add("P08581", "bgwp", "P08581")
	xmap.Upac.Add("P10275", "bgwp", "P10275")
	xmap.Upac.Add("P31749", "bgwp", "P31749") // PF24
	xmap.Upac.Add("P31751", "bgwp", "P31751") // PF24
	/// exporting
	srcs := []string{"signor"}

	// for Signor complexes and protein families
	sigmap := make(util.Set3D)
	subdir := "signor/"
	dpth := pth + subdir
	ss0 := []string{dpth + "complexes.map", dpth + "families.map"}
	parse.Sig2up(sigmap, ss0)
	xmap.Signor = sigmap

	tts := []tt{
		// {&d4b0, &xmap, pth + "OUT/export/", 6}, // 3 if filtered by the host
		{&d4b0, &xmap, pth + "OUT/export/", 5},
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
	_ = parse.Tab2struct(pth+"signor/9606.mi28", keys, vals, &d4b0, "\t")
	var xmap bgw.Xmap
	xmap.New()
	xmap.Upac.Add("P27361", "bgwp", "P27361")
	xmap.Upac.Add("P48431", "bgwp", "P48431") // for testing isoforms
	xmap.Upac.Add("Q9BTC0", "bgwp", "Q9BTC0")
	xmap.Upac.Add("P08648", "bgwp", "P08648")
	xmap.Bgwp.Add("P08648", "bgwg", "9606/GENEX")
	xmap.Upac.Add("P10275", "bgwp", "P10275")
	xmap.Upac.Add("P04637", "bgwp", "P04637")
	xmap.Upac.Add("Q01081", "bgwp", "Q01081")
	/// exporting
	srcs := []string{"signor"}

	// for Signor complexes and protein families
	sigmap := make(util.Set3D)
	subdir := "signor/"
	dpth := pth + subdir
	ss0 := []string{dpth + "complexes.map", dpth + "families.map"}
	parse.Sig2up(sigmap, ss0)
	xmap.Signor = sigmap

	tts := []tt{
		{&d4b0, &xmap, pth + "OUT/export/", 4},
	}
	pdck := "reg2utrg"
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

func Test_Prot2prot(t *testing.T) {
	type tt struct {
		arg1 *bgw.Dat4bridge
		arg2 *bgw.Xmap
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	// parsing
	keys, vals := bgw.IntactParseConf()
	var d4b0 bgw.Dat4bridge
	d4b0.New()
	_ = parse.Tab2struct(pth+"intact/9606.mi25", keys, vals, &d4b0, "\t")
	var d4b1 bgw.Dat4bridge
	d4b1.New()
	_ = parse.Tab2struct(pth+"intact/367110.mi25", keys, vals, &d4b1, "\t")

	var xmap bgw.Xmap
	xmap.New()
	xmap.Upac.Add("P04637", "bgwp", "P04637")
	xmap.Upac.Add("Q8X225", "bgwp", "Q8X225") // for 367110
	xmap.Upac.Add("P07041", "bgwp", "P07041") // for 367110
	xmap.Upac.Add("Q01371", "bgwp", "Q01371") // for 367110
	/// exporting
	srcs := []string{"intact", "intact"}
	taxa := []string{"9606", "367110"}
	tts := []tt{
		{&d4b0, &xmap, pth + "OUT/export/", 1},
		{&d4b1, &xmap, pth + "OUT/export/", 2},
	}
	pdck := "tlp2tlp"
	for i, tt := range tts {
		(*tt.arg1).Src = srcs[i]
		(*tt.arg1).Taxid = taxa[i]
		Prot2prot(tt.arg1, tt.arg2, tt.arg3)
		cnts := (*tt.arg1).Cnts
		if cnts[pdck][srcs[i]] != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", cnts[pdck][srcs[i]],
			)
		}
	}
} // Prot2prot

// output: 3 RDF files, each one specific for a particular predicate
// URIs conditioned on the predicate
// multiple subjects and objects are possible
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
	_ = parse.Tab2struct(pth+"static/tfacts/9606.f2g", keys, vals, &d4b0, "\t")
	var d4b1 bgw.Dat4bridge
	d4b1.New()
	keys, vals = bgw.TftgParseConf()
	_ = parse.Tab2struct(pth+"static/ntnu/9606.f2g", keys, vals, &d4b1, "\t")
	var d4b2 bgw.Dat4bridge
	d4b2.New()
	keys, vals = bgw.TflinkParseConf()
	_ = parse.Tab2struct(pth+"tflink/9606.tsv", keys, vals, &d4b2, "\t")
	var xmap bgw.Xmap
	xmap.New()
	xmap.Upac.Add("P01100", "bgwp", "P01100")
	xmap.Upac.Add("P04637", "bgwp", "P04637")
	xmap.Upac.Add("P37231", "bgwp", "P37231") // for tflink
	xmap.Ncbig.Add("4322", "bgwg", "9606/MMP13")
	xmap.Ncbig.Add("7157", "bgwg", "9606/TP53")
	xmap.Ncbig.Add("5915", "bgwg", "9606/RARB") // for tflink
	xmap.Bgwg.Add("9606/MMP13", "bgwp", "P01100")
	xmap.Bgwg.Add("9606/TP53", "bgwp", "P04637")
	xmap.Bgwg.Add("9606/RARB", "bgwp", "P10826") // for tflink
	/// exporting
	srcs := []string{"tfacts", "ntnu", "tflink"}

	tts := []tt{
		{&d4b0, &xmap, pth + "OUT/export/", 2},
		{&d4b1, &xmap, pth + "OUT/export/", 1},
		{&d4b2, &xmap, pth + "OUT/export/", 1}, // sic!
	}
	pdck := "reg2utrg"
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
	arg2.Add("P04637", "bgwp", "P04637")
	// no isoforms in gaf files
	bps, ccs, mfs, _ := parse.Gaf(pth+"goa/9606.gaf", arg2)
	out := [3]string{"prot2bp/9606.nt", "prot2cc/9606.nt", "prot2mf/9606.nt"}
	tts := []tt{
		{bps, arg2, wpth + out[0], 13},
		{ccs, arg2, wpth + out[1], 14},
		{mfs, arg2, wpth + out[2], 13},
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
		arg3 string
		val1 int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/export/"
	arg1 := make(util.Set3D)
	// arg1.Add("uniprot!P02340--uniprot!P04637", "KO", "K04451")
	arg1.Add("uniprot!P02340--uniprot!P04637", "OrthoDB", "257530at2759") // no isoforms
	arg3 := wpth + "ortho/10090-9606.nt"
	tts := []tt{
		{arg1, arg3, 14},
	}
	for i, tt := range tts {
		n, err := Ortho(tt.arg1, tt.arg3)
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
