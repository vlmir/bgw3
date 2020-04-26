package parse

import (
	"github.com/vlmir/bgw3/src/bgw" // pkg 'bgw'
	"github.com/vlmir/bgw3/src/util" // pkg 'util'
	"testing"
)

func Test_Idmap(t *testing.T) {
type tt struct {
	arg1 string
	arg2 map[string]string
	arg3 int
	arg4 int
	arg5 int
	val  int
}
	pth := "../../tdata/"
	t1s := []tt{
		{pth + "test.idm", map[string]string{"NCBI_TaxID": "test"}, 2, 1, 0, 1},
		{pth + "test.idm", map[string]string{"NCBI_TaxID": "test"}, 0, 1, 2, 1},
		{pth + "test.idm", map[string]string{"UniParc": "test"}, 0, 1, 2, 9},
	}
	for i, tt := range t1s {
		idm, _ := Idmap(tt.arg1, tt.arg2, tt.arg3, tt.arg4, tt.arg5)
		if len(idm) != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3, tt.arg4, tt.arg5,
				"\n\twant", tt.val,
				"\n\thave", len(idm),
			)
		}
	}
}

func Test_Upidmap(t *testing.T) {
type tt struct {
	arg1 string
	arg2 map[string]string
	val1 int
	val2 int
	val3 int
}
	pth := "../../tdata/"
	idms := []tt{
		{pth + "test.idm", map[string]string{"UniParc": "test"}, 1, 9, 0},
	}
	for i, tt := range idms {
		set1, set2, set3, _ := Upidmap(tt.arg1, tt.arg2)
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

func Test_Updat(t *testing.T) {
type tt struct {
	arg1 string
	arg2 util.Set2D
	val1 int
	val2 int
}
	pth := "../../tdata/"
	upt := make(util.Set2D)
	upts := []tt{
		{pth + "test.upt", upt, 1, 1},
	}
	for i, tt := range upts {
		tt.arg2.Add("P04637", "P04637-2")
		set1, set2, _ := Updat(tt.arg1, tt.arg2)
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
	}
}

func Test_Upvar(t *testing.T) {
type tt struct {
	arg1 string
	arg2 util.Set3D
	val1 int
}
	pth := "../../tdata/"
	upvar := make(util.Set3D)
	upvars := []tt{
		{pth + "test.var", upvar, 1},
	}
	for i, tt := range upvars {
		tt.arg2.Add("TP53", "test", "t")
		set1 := Upvar(tt.arg1, tt.arg2)
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func Test_Mitab(t *testing.T) {
type tt struct {
	arg1 string
	arg2 util.Set3D
	val1 int
}
	pth := "../../tdata/"
	mit := make(util.Set3D)
	mits := []tt{
		{pth + "test.mit", mit, 1},
	}
	for i, tt := range mits {
		tt.arg2.Add("P04637", "test", "t")
		set1 := Mitab(tt.arg1, tt.arg2)
		if len(set1) != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func Test_Tftg(t *testing.T) {
type tt struct {
	arg1 string
	arg2 util.Set3D
	arg3 util.Set3D
	val1 int
}
	pth := "../../tdata/"
	s1 := make(util.Set3D)
	s2 := make(util.Set3D)
	var f2gs = []tt{
		{pth + "test.f2g", s1, s2, 1},
	}
	for i, tt := range f2gs {
		tt.arg2.Add("P04637", "test", "TP53")
		tt.arg3.Add("TP53", "test", "P04637")
		set1, _ := Tftg(tt.arg1, tt.arg2, tt.arg3)
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
		{pth + "test.gaf", set, 150, 17, 39},
	}
	for i, tt := range gafs {
		tt.arg2.Add("P04637", "test", "t")
		set1, set2, set3 := Gaf(tt.arg1, tt.arg2)
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

func Test_Orthoduo(t *testing.T){
	idmkeys := bgw.Orthokeys
	arg1 := "../../tdata/"
	arg3 := "9606"
	arg2 := "10090"
	var arg4 [5]util.Set2D // tx2pm
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
		out, _ := Orthoduo(arg1, arg2, arg3, arg4[i], arg5)
		if len(out) != val1[i] {
			t.Error(
				"For test", i+1, ": ", arg1, arg2, arg3, arg4[i],
				"\n\twant", val1[i],
				"\n\thave", len(out),
			)
		}
	}
}
