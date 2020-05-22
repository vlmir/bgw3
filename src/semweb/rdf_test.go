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
		"ppy2prn",
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

//func Test_Header(t *testing.T) {
//	type tt struct {
//		arg1 map[string]string
//		arg2 util.SliceSet
//		val1 int
//	}
//	arg2 := make(util.SliceSet)
//	arg2["Opys"] = []string{
//		"sub2cls",
//		"ppy2prn",
//	}
//	arg2["Apys"] = []string{
//		"sth2lbl",
//	}
//	arg2["Prns"] = []string{
//		"cls",
//		"opy",
//		"apy",
//	}
//	arg1 := FmtURIs(arg2)
//	tts := []tt{
//		{arg1, arg2, 12},
//	}
//	for i, tt := range tts {
//		_, n := Header(tt.arg1, tt.arg2)
//		val1 := n
//		if val1 != tt.val1 {
//			t.Error(
//				"For test", i+1, ": ",
//				"\n\twant", tt.val1,
//				"\n\thave", val1,
//			)
//		}
//	}
//}
