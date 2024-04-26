package parse

import (
	"github.com/vlmir/bgw3/src/bgw"
	"github.com/vlmir/bgw3/src/util"
	"testing"
)

func Test_Tab2struct(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 []bgw.Column
		arg3 []bgw.Column
		arg4 *bgw.Dat4bridge
		arg5 string
		val  int // the number of fields associated with a key
	}
	d4bs := make([]bgw.Dat4bridge, 0, 5)
	keys0, vals0 := bgw.IntactParseConf()
	keys1, vals1 := bgw.TftgParseConf()
	keys2, vals2 := bgw.SignorParseConf()
	keys3, vals3 := bgw.TflinkParseConf()
	keys4, vals4 := bgw.ColtriParseConf()
	var d4b0 bgw.Dat4bridge
	var d4b1 bgw.Dat4bridge
	var d4b2 bgw.Dat4bridge
	var d4b3 bgw.Dat4bridge
	var d4b4 bgw.Dat4bridge
	d4b0.New()
	d4b1.New()
	d4b2.New()
	d4b3.New()
	d4b4.New()
	d4bs = append(d4bs, d4b0, d4b1, d4b2, d4b3, d4b4)
	pth := "../../tdata/"
	tts := []tt{
		{pth + "intact/9606.mi25", keys0, vals0, &d4bs[0], "\t", 7}, // 8? TODO
		{pth + "static/tfacts/9606.f2g", keys1, vals1, &d4bs[1], "\t", 5},
		{pth + "signor/9606.mi28", keys2, vals2, &d4bs[2], "\t", 16}, // 15->16 - added mtd field
		{pth + "tflink/9606.tsv", keys3, vals3, &d4bs[3], "\t", 6},
		{pth + "coltri/9606.csv", keys4, vals4, &d4bs[4], ",", 9},
	}
	keys := []string{
		"P04637--P04637",
		"AP1--SPP1",
		"uniprotkb:Q9BTC0--uniprotkb:P08648",
		"NANOG--DKK1",
		"Q16254--Q01094",
	}
	for i, tt := range tts {
		_ = Tab2struct(tt.arg1, tt.arg2, tt.arg3, tt.arg4, tt.arg5)
		if len(d4bs[i].Duos[keys[i]]) != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", len(d4bs[i].Duos[keys[i]]),
			)
		}
	}
} // Tab2struct

func Test_Tab2set3D(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 []bgw.Column
		arg3 []bgw.Column
		val  int
	}
	arg2_1, arg3_1 := bgw.IntactParseConf()
	arg2_2, arg3_2 := bgw.TftgParseConf()
	arg2_3, arg3_3 := bgw.SigMapParseConf()
	arg2_4, arg3_4 := bgw.SigMapParseConf()
	pth := "../../tdata/"
	tts := []tt{
		{pth + "intact/9606.mi25", arg2_1, arg3_1, 7},
		{pth + "static/tfacts/9606.f2g", arg2_2, arg3_2, 5},
		{pth + "signor/complexes.tsv", arg2_3, arg3_3, 1},
		{pth + "signor/families.tsv", arg2_4, arg3_4, 1},
	}
	keys := []string{
		"P04637--P04637",
		"AP1--SPP1",
		"SIGNOR-C1",
		"SIGNOR-PF1",
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
} // Tab2set3D

func Test_UpVar(t *testing.T) {
	type tt struct {
		arg1 string
		val1 int
	}
	pth := "../../tdata/"
	upvars := []tt{
		{pth + "uniprot/P04637.var", 1},
	}
	for i, tt := range upvars {
		set1, _ := UpVar(tt.arg1)
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func Test_Gaf(t *testing.T) {
	type tt struct {
		arg1 string
		val1 int
		val2 int
		val3 int
	}
	pth := "../../tdata/"
	var gafs = []tt{
		// {pth + "goa/9606.gaf", set, 150, 17, 39},
		{pth + "goa/9606.gaf", 1, 1, 1},
	}
	for i, tt := range gafs {
		set1, set2, set3, _ := Gaf(tt.arg1)
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
		if len(set2) != tt.val2 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val2,
				"\n\thave", len(set2),
			)
		}
		if len(set3) != tt.val3 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val3,
				"\n\thave", len(set3),
			)
		}
	}
}

func Test_OrthoDuo(t *testing.T) {
	// TODO re-implement
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
