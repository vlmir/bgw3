package parse

import (
	"github.com/vlmir/bgw3/src/bgw"
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
		{pth + "intact/9606.mi25", bgw.IntactConf.Keys, bgw.IntactConf.Vals, &d4bs[0], "\t", 7},  // 8? TODO
		{pth + "signor/9606.mi28", bgw.SignorConf.Keys, bgw.SignorConf.Vals, &d4bs[1], "\t", 16}, // 15->16 - added mtd field
		{pth + "tflink/9606.tsv", bgw.TflinkConf.Keys, bgw.TflinkConf.Vals, &d4bs[2], "\t", 6},
		{pth + "coltri/9606.csv", bgw.ColtriConf.Keys, bgw.ColtriConf.Vals, &d4bs[3], ",", 9},
		{pth + "atregnet/3702.tsv", bgw.AtregnetConf.Keys, bgw.AtregnetConf.Vals, &d4bs[4], "\t", 5},
	}
	keys := []string{
		"P04637--P04637",
		"uniprotkb:Q9BTC0--uniprotkb:P08648",
		"NANOG--DKK1",
		"Q16254--Q01094",
		"AT3G22170--AT1G02400", // Q9LIE5 - 3702/GA2OX6
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
	pth := "../../tdata/"
	tts := []tt{
		{pth + "intact/9606.mi25", bgw.IntactConf.Keys, bgw.IntactConf.Vals, 7},
		{pth + "signor/complexes.tsv", bgw.SigMapConf.Keys, bgw.SigMapConf.Vals, 1},
		{pth + "signor/families.tsv", bgw.SigMapConf.Keys, bgw.SigMapConf.Vals, 1},
	}
	keys := []string{
		"P04637--P04637",
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
