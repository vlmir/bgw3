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
	ss0 := []string{dpth + "complexes.map", dpth + "families.map"}
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
	d4bs := make([]bgw.Dat4bridge, 0, 5)
	keys0, vals0 := bgw.IntactParseConf()
	keys1, vals1 := bgw.TftgParseConf()
	keys2, vals2 := bgw.SignorParseConf()
	keys3, vals3 := bgw.TflinkParseConf()
	var d4b0 bgw.Dat4bridge
	var d4b1 bgw.Dat4bridge
	var d4b2 bgw.Dat4bridge
	var d4b3 bgw.Dat4bridge
	d4b0.New()
	d4b1.New()
	d4b2.New()
	d4b3.New()
	d4bs = append(d4bs, d4b0, d4b1, d4b2, d4b3)
	pth := "../../tdata/"
	tts := []tt{
		{pth + "intact/9606.mi25", keys0, vals0, &d4bs[0], 7}, // 8? TODO
		{pth + "static/tfacts/9606.f2g", keys1, vals1, &d4bs[1], 5},
		{pth + "signor/9606.mi28", keys2, vals2, &d4bs[2], 16}, // 15->16 - added mtd field
		{pth + "tflink/9606.tsv", keys3, vals3, &d4bs[3], 6},
	}
	keys := []string{
		"P04637--P04637",
		"AP1--SPP1",
		"uniprotkb:Q9BTC0--uniprotkb:P08648",
		//		"uniprot:Q9H9S0--geneid:22943",
		"NANOG--DKK1",
	}
	for i, tt := range tts {
		_ = Tab2struct(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		if len(d4bs[i].Duos[keys[i]]) != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", len(d4bs[i].Duos[keys[i]]),
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
		{pth + "intact/9606.mi25", arg2_1, arg3_1, 2},
		{pth + "static/tfacts/9606.f2g", arg2_2, arg3_2, 5},
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
		{pth + "goa/9606.gaf", set, 150, 17, 39},
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
	arg1 := "../../tdata/"
	arg3 := "9606"
	arg2 := "10090"
	var arg4 [5]util.Set2D // txn2prm
	var val1 [5]int
	n := 1 // number of tests
	for i := 0; i < n; i++ {
		arg4[i] = make(util.Set2D)
	}
	arg4[0].Add("9606", "ortho")
	arg4[0].Add("10090", "ortho")
	val1[0] = 2
	for i := 0; i < n; i++ {
		out, _ := OrthoDuo(arg1, arg2, arg3, arg4[i])
		if len(out) != val1[i] {
			t.Error(
				"For test", i+1, ": ", arg1, arg2, arg3, arg4[i],
				"\n\twant", val1[i],
				"\n\thave", len(out),
			)
		}
	}
}
