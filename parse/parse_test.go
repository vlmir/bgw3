package parse

import (
	"aux"
	"testing"
)

type t1 struct {
	arg1 string
	arg2 map[string]string
	arg3 int
	arg4 int
	arg5 int
	val  int
}
type t2 struct {
	arg1 string
	arg2 map[string]string
	val1 int
	val2 int
	val3 int
}
type t3 struct {
	arg1 string
	arg2 aux.Set2D
	val1 int
	val2 int
}
type t4 struct {
	arg1 string
	arg2 aux.Set3D
	val1 int
}
type t5 struct {
	arg1 string
	arg2 aux.Set3D
	arg3 aux.Set3D
	val1 int
}
type t6 struct {
	arg1 string
	arg2 aux.Set3D
	val1 int
	val2 int
	val3 int
}

func TestIdmap(t *testing.T) {
	pth := "../tdata/"
	t1s := []t1{
		{pth + "test.idm", map[string]string{"NCBI_TaxID": "test"}, 2, 1, 0, 1},
		{pth + "test.idm", map[string]string{"NCBI_TaxID": "test"}, 0, 1, 2, 1},
		{pth + "test.idm", map[string]string{"UniParc": "test"}, 0, 1, 2, 9},
	}
	for i, t1 := range t1s {
		idm, _ := Idmap(t1.arg1, t1.arg2, t1.arg3, t1.arg4, t1.arg5)
		if len(idm) != t1.val {
			t.Error(
				"For test", i+1, ": ", t1.arg1, t1.arg2, t1.arg3, t1.arg4, t1.arg5,
				"\n\twant", t1.val,
				"\n\thave", len(idm),
			)
		}
	}
}

func TestUpidmap(t *testing.T) {
	pth := "../tdata/"
	idms := []t2{
		{pth + "test.idm", map[string]string{"UniParc": "test"}, 1, 9, 0},
	}
	for i, t2 := range idms {
		set1, set2, set3, _ := Upidmap(t2.arg1, t2.arg2)
		if len(set1) != t2.val1 {
			t.Error(
				"For test", i+1, ": ", t2.arg1, t2.arg2,
				"\n\twant", t2.val1,
				"\n\thave", len(set1),
			)
		}
		if len(set2) != t2.val2 {
			t.Error(
				"For test", i+1, ": ", t2.arg1, t2.arg2,
				"\n\twant", t2.val2,
				"\n\thave", len(set2),
			)
		}
		if len(set3) != t2.val3 {
			t.Error(
				"For test", i+1, ": ", t2.arg1, t2.arg2,
				"\n\twant", t2.val3,
				"\n\thave", len(set3),
			)
		}
	}
}

func TestUpdat(t *testing.T) {
	pth := "../tdata/"
	upt := make(aux.Set2D)
	upts := []t3{
		{pth + "test.upt", upt, 1, 1},
	}
	for i, t3 := range upts {
		t3.arg2.Add("P04637", "P04637-2")
		set1, set2, _ := Updat(t3.arg1, t3.arg2)
		if len(set1) != t3.val1 {
			t.Error(
				"For test", i+1, ": ", t3.arg1, t3.arg2,
				"\n\twant", t3.val1,
				"\n\thave", len(set1),
			)
		}
		if len(set2) != t3.val2 {
			t.Error(
				"For test", i+1, ": ", t3.arg1, t3.arg2,
				"\n\twant", t3.val2,
				"\n\thave", len(set2),
			)
		}
	}
}

func TestUpvar(t *testing.T) {
	pth := "../tdata/"
	upvar := make(aux.Set3D)
	upvars := []t4{
		{pth + "test.var", upvar, 1},
	}
	for i, t4 := range upvars {
		t4.arg2.Add("TP53", "test", "t")
		set1 := Upvar(t4.arg1, t4.arg2)
		if len(set1) != t4.val1 {
			t.Error(
				"For test", i+1, ": ", t4.arg1, t4.arg2,
				"\n\twant", t4.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func TestMitab(t *testing.T) {
	pth := "../tdata/"
	mit := make(aux.Set3D)
	mits := []t4{
		{pth + "test.mit", mit, 1},
	}
	for i, t4 := range mits {
		t4.arg2.Add("P04637", "test", "t")
		set1 := Mitab(t4.arg1, t4.arg2)
		if len(set1) != t4.val1 {
			t.Error(
				"For test", i+1, ": ", t4.arg1, t4.arg2,
				"\n\twant", t4.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func TestTftg(t *testing.T) {
	pth := "../tdata/"
	s1 := make(aux.Set3D)
	s2 := make(aux.Set3D)
	var f2gs = []t5{
		{pth + "test.f2g", s1, s2, 1},
	}
	for i, t5 := range f2gs {
		t5.arg2.Add("P04637", "test", "TP53")
		t5.arg3.Add("TP53", "test", "P04637")
		set1, _ := Tftg(t5.arg1, t5.arg2, t5.arg3)
		if len(set1) != t5.val1 {
			t.Error(
				"For test", i+1, ": ", t5.arg1, t5.arg2,
				"\n\twant", t5.val1,
				"\n\thave", len(set1),
			)
		}
	}
}

func TestGaf(t *testing.T) {
	pth := "../tdata/"
	var set = make(aux.Set3D)
	var gafs = []t6{
		{pth + "test.gaf", set, 150, 17, 39},
	}
	for i, t6 := range gafs {
		t6.arg2.Add("P04637", "test", "t")
		set1, set2, set3 := Gaf(t6.arg1, t6.arg2)
		if len(set1) != t6.val1 {
			t.Error(
				"For test", i+1, ": ", t6.arg1, t6.arg2,
				"\n\twant", t6.val1,
				"\n\thave", len(set1),
			)
		}
		if len(set2) != t6.val2 {
			t.Error(
				"For test", i+1, ": ", t6.arg1, t6.arg2,
				"\n\twant", t6.val2,
				"\n\thave", len(set2),
			)
		}
		if len(set3) != t6.val3 {
			t.Error(
				"For test", i+1, ": ", t6.arg1, t6.arg2,
				"\n\twant", t6.val3,
				"\n\thave", len(set3),
			)
		}
	}
}
