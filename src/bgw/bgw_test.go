package bgw

import (
	"testing"
)

func TestDat4bridgeNew(t *testing.T) {
	type tt struct {
		val1 int
		val2 int
	}
	var d4b Dat4bridge
	d4b.New()
	d4b.Duos.Add("abc", "def", "ghi")
	d4b.Duos.Add("abc", "klm", "nop")
	tts := []tt{
		{2, 1},
	}
	for i, tt := range tts {
		val1 := len(d4b.Duos["abc"])
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
		val2 := len(d4b.Duos["abc"]["def"])
		if val2 != tt.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val2,
				"\n\thave", val2,
			)
		}
	}
}

func TestDat4rdfNew(t *testing.T) {
	type tt struct {
		val1 int
	}
	var d Dat4rdf
	d.New()
	(*d.Udat).Add("abc", "def", "ghi")
	(*d.Udat).Add("abc", "klm", "nop")
	tts := []tt{
		{2},
	}
	for i, tt := range tts {
		val1 := len((*d.Udat)["abc"])
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func TestXmapNew(t *testing.T) {
	type tt struct {
		val1 int
	}
	var x Xmap
	x.New()
	x.Upac.Add("abc", "def", "ghi")
	x.Upac.Add("abc", "klm", "nop")
	tts := []tt{
		{2},
	}
	for i, tt := range tts {
		val1 := len(x.Upac["abc"])
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Unmarshal(t *testing.T) {
	type tt struct {
		arg1 string
		val1 int
	}
	pth := "../../tdata/"
	tts := []tt{
		// {pth + "OUT/export/xmap/9606.json", 7},
		{pth + "OUT/rdf4bgw/xmap/9606.json", 3}, // rdf4bgw.Geneprot() new
	}
	for i, tt := range tts {
		var xmap Xmap
		xmap.New()
		xmap.Unmarshal(tt.arg1)
		val1 := len(xmap.Bgwp)
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}
