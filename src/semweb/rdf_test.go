package rdf
import (
	"testing"
	"github.com/vlmir/bgw3/src/utils" // pkg 'aux'"
)

func TestZenoUnmarshal(t *testing.T) {
	pth := "../../tdata/"
	type tt struct {
		arg1 string
		val1 int
	}
	tts := []tt{
		{pth + "zeno.json", 33},
	}
	for i, tt := range tts {
		zeno := NewZeno()
		zeno.Unmarshal(tt.arg1)
		val1 := len(zeno.Uris)
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func TestHeader(t *testing.T) {
	pth := "../../tdata/"
	type tt struct {
		arg1 map[string]string
		arg2 aux.SliceSet
		arg3 Zeno
		val1 int
		val2 int
	}
	uris := make(map[string]string)
	zeno := NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	set1 := make(aux.SliceSet)
	set1["Opys"] = []string{
		"sub2cls",
	}
	set1["Apys"] = []string{
		"sth2lbl",
	}
	set1["Prns"] = []string{
		"stm",
	}
	tts := []tt{
		{uris, set1, zeno, 6, 3},
	}
	for i, tt := range tts {
		_, n := Header(tt.arg1, tt.arg2, tt.arg3)
		val1 := n
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
		val2 := len(tt.arg1)
		if val2 != tt.val2 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val2,
				"\n\thave", val2,
			)
		}
	}
}

