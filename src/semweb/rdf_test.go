package rdf

import (
	"github.com/vlmir/bgw3/src/util"
	"testing"
)

func Test_FmtURIs(t *testing.T) {
	type tt struct {
		arg1 util.SliceSet
		val1 int
	}
	arg1 := make(util.SliceSet)
	arg1["Opys"] = []string{
		"sub2cls",
		"ins2cls",
	}
	arg1["Apys"] = []string{
		"sth2lbl",
	}
	arg1["Prns"] = []string{
		"cls",
	}
	tts := []tt{
		{arg1, 4},
	}
	for i, tt := range tts {
		val1 := len(FmtURIs(tt.arg1))
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Capita(t *testing.T) {
	type tt struct {
		arg1 util.SliceSet
		val1 int
	}
	arg1 := make(util.SliceSet)
	arg1["Opys"] = []string{
		"sub2cls",
	}
	arg1["Apys"] = []string{
		"sth2lbl",
	}
	arg1["Prns"] = []string{
		"stm",
	}
	tts := []tt{
		{arg1, 6},
	}
	for i, tt := range tts {
		_, n := Capita(tt.arg1)
		val1 := n
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}
