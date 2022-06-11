package parse

import (
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/util"
	"log"
	"testing"
	"time"
)

func Test_Sig2up(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 util.Set3D
		arg2 []string
		val1 error
	}

	mp := make(util.Set3D)
	mps := [...]util.Set3D{mp}
	pth := "../../tdata/"
	subdir := "signor/"
	dpth := pth + subdir
	ss0 := []string{dpth + "sig-c.head3", dpth + "sig-pf.head3"}
	sss := [...][]string{ss0}

	tts := []tt{
		{mps[0], sss[0], nil},
	}
	for i, tt := range tts {
		err := Sig2up(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with Sig2up() in", time.Since(mystart))
}

func Test_Tab2struct(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 []bgw.Column
		arg3 []bgw.Column
		arg4 *bgw.Dat4bridge
		val  int
	}
	keys0 := []bgw.Column{
		{0, ":", 1, "--", 0, ""},
		{1, ":", 1, "--", 0, ""},
	}
	vals0 := []bgw.Column{
		{6, "|", 1, "\"", 1, "mtd"},
		{8, "|", 1, ":", -1, "pubmed"},
	}
	keys1, vals1 := bgw.TftgParseConf()
	keys2, vals2 := bgw.SignorParseConf()
	var d4b bgw.Dat4bridge
	d4b.New()
	arg4 := &d4b
	pth := "../../tdata/"
	tts := []tt{
		{pth + "intact/test.mit", keys0, vals0, arg4, 2},
		{pth + "tfacts/test.f2g", keys1, vals1, arg4, 5},
		{pth + "signor/9606.mi28", keys2, vals2, arg4, 12},
	}
	keys := []string{
		"P04637--P04637",
		"AP1--SPP1",
		"uniprotkb:Q9BTC0--uniprotkb:P08648",
	}
	for i, tt := range tts {
		_ = Tab2struct(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		if len(d4b.Duos[keys[i]]) != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", len(d4b.Duos[keys[i]]),
			)
		}
	}
}

func Test_Tab2set3D(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 []bgw.Column
		arg3 []bgw.Column
		val  int
	}
	arg2_1 := []bgw.Column{
		{0, ":", 1, "--", 0, ""},
		{1, ":", 1, "--", 0, ""},
	}
	arg3_1 := []bgw.Column{
		{6, "|", 1, "\"", 1, "mtd"},
		{8, "|", 1, ":", -1, "pubmed"},
	}
	arg2_2, arg3_2 := bgw.TftgParseConf()
	pth := "../../tdata/"
	tts := []tt{
		{pth + "intact/test.mit", arg2_1, arg3_1, 2},
		{pth + "tfacts/test.f2g", arg2_2, arg3_2, 5},
	}
	keys := []string{
		"P04637--P04637",
		"AP1--SPP1",
		"uniprotkb:Q9BTC0--uniprotkb:P08648",
	}
	for i, tt := range tts {
		out, _ := Tab2set3D(tt.arg1, tt.arg2, tt.arg3)
		if len(out[keys[i]]) != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", len(out[keys[i]]),
			)
		}
	}
}

func Test_UpIdMap(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 map[string]string
		val1 int
	}
	pth := "../../tdata/idmapping/"
	idms := []tt{
		{pth + "UP000005640_9606.idmapping", map[string]string{"UniParc": "test"}, 5},
		{pth + "UP000000803_7227.idmapping", map[string]string{"EnsemblGenome": "test"}, 1},
	}
	for i, tt := range idms {
		out, _ := UpIdMap(tt.arg1, tt.arg2)
		if len(out) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", len(out),
			)
		}
	}
}

func Test_UpTab(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 util.Set3D
		arg3 util.Set2D
		val1 int
		val2 int
		val3 int
	}
	pth := "../../tdata/"
	txn2prm := make(util.Set2D)
	txn2prm.Add("9606", "UP000005640")
	txn2prm.Add("7227", "UP000000803")
	upt := make(util.Set3D)
	upts := []tt{
		{pth + "uniprot/9606.upt", upt, txn2prm, 2, 1, 2},
		{pth + "uniprot/7227.upt", upt, txn2prm, 1, 1, 1},
	}
	for i, tt := range upts {
		tt.arg2.Add("P04637", "upac", "P04637-2")
		tt.arg2.Add("FOOFOO", "upac", "FOOFOO-2")
		tt.arg2.Add("P00528", "upac", "P00528")
		out, _ := UpTab(tt.arg1, tt.arg2, txn2prm)
		set1 := *out.Udat
		set2 := *out.Txns
		set3 := *out.Gnm
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
		if len(set2) != tt.val2 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val2,
				"\n\thave", len(set2),
			)
		}
		if len(set3) != tt.val3 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val3,
				"\n\thave", len(set3),
			)
		}
	}
}

func Test_UpVar(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 util.Set3D
		val1 int
	}
	pth := "../../tdata/"
	upvar := make(util.Set3D)
	upvars := []tt{
		{pth + "uniprot/P04637.var", upvar, 1},
	}
	for i, tt := range upvars {
		tt.arg2.Add("TP53", "test", "t")
		set1, _ := UpVar(tt.arg1, tt.arg2)
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func Test_MiTab(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 util.Set3D
		val1 int
	}
	pth := "../../tdata/"
	mit := make(util.Set3D)
	mits := []tt{
		{pth + "intact/test.mit", mit, 1},
	}
	for i, tt := range mits {
		tt.arg2.Add("P04637", "test", "t")
		set1, _ := MiTab(tt.arg1, tt.arg2)
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func Test_Gaf(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 util.Set3D
		val1 int
		val2 int
		val3 int
	}
	pth := "../../tdata/"
	var set = make(util.Set3D)
	var gafs = []tt{
		{pth + "goa/test.gaf", set, 150, 17, 39},
	}
	for i, tt := range gafs {
		tt.arg2.Add("P04637", "test", "t")
		set1, set2, set3, _ := Gaf(tt.arg1, tt.arg2)
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
		if len(set2) != tt.val2 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val2,
				"\n\thave", len(set2),
			)
		}
		if len(set3) != tt.val3 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val3,
				"\n\thave", len(set3),
			)
		}
	}
}

func Test_OrthoDuo(t *testing.T) {
	idmkeys := bgw.Orthokeys
	arg1 := "../../tdata/"
	arg3 := "9606"
	arg2 := "10090"
	var arg4 [5]util.Set2D // txn2prm
	arg5 := idmkeys
	var val1 [5]int
	n := 1 // number of tests
	for i := 0; i < n; i++ {
		arg4[i] = make(util.Set2D)
	}
	arg4[0].Add("9606", "ortho")
	arg4[0].Add("10090", "ortho")
	val1[0] = 2
	for i := 0; i < n; i++ {
		out, _ := OrthoDuo(arg1, arg2, arg3, arg4[i], arg5)
		if len(out) != val1[i] {
			t.Error(
				"For test", i+1, ": ", arg1, arg2, arg3, arg4[i],
				"\n\twant", val1[i],
				"\n\thave", len(out),
			)
		}
	}
}
