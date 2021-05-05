package bgw

import (
	"testing"
)

func TestDat4bridgeNew(t *testing.T) {
	var d4b Dat4bridge
	d4b.New()
	d4b.Duos.Add("abc", "def", "ghi")
	d4b.Duos.Add("abc", "klm", "nop")
	type tt struct {
		val1 int
		val2 int
	}
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
	var d Dat4rdf
	d.New()
	(*d.Udat).Add("abc", "def", "ghi")
	(*d.Udat).Add("abc", "klm", "nop")
	type tt struct {
		val1 int
	}
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

func Test_Unmarshal(t *testing.T) {

	pth := "../../tdata/"
	type tt struct {
		arg1 string
		val1 int
	}
	tts := []tt{
		{pth + "xmap.json", 9},
	}
	for i, tt := range tts {
		xmap := NewXmap()
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
