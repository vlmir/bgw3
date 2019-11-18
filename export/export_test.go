package export

import (
	"aux"
	"bgw/parse"
	"bgw/semweb"
	"bgw/util"
	"testing"
)

/// for: GeneProt()
type t1 struct {
	arg1 util.Dat4rdf
	arg2 string
	arg3 string
	arg4 string
	arg5 rdf.Zeno
	val1 int
	val2 int
}

/// for: Upvar()
type t2 struct {
	arg1 aux.Set3D
	arg2 aux.Set3D
	arg3 aux.Set3D
	arg4 string
	arg5 rdf.Zeno
	val1 int
}

/// for: Tftg()
type t3 struct {
	arg1 aux.Set3D
	arg2 util.Meta
	arg3 aux.Set3D
	arg4 aux.Set3D
	arg5 string
	arg6 rdf.Zeno
	val1 int
}

/// for: Mitab(), Gaf()
type t4 struct {
	arg1 aux.Set3D
	arg2 aux.Set3D
	arg3 string
	arg4 rdf.Zeno
	val1 int
}

func TestOrtho(t *testing.T) {
	pth := "../tdata/"
	var idmkeys = map[string]string{
		"KO":      "kego",
		"OrthoDB": "orthodb",
	}
	arg1, _ := parse.Idmap(pth+"test.idm", idmkeys, 1, 2, 0)
	arg2 := make(aux.Set3D)
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	t4s := []t4{
		{arg1, arg2, "/dev/null", zeno, 18},
	}
	for i, t4 := range t4s {
		n, err := Ortho(t4.arg1, t4.arg2, t4.arg3, t4.arg4)
		if err != nil {
			return
		}
		if n != t4.val1 {
			t.Error(
				"For test", i+1, ": ", len(t4.arg1), t4.arg2, t4.arg3,
				"\n\twant", t4.val1,
				"\n\thave", n,
			)
		}
	}
}

func TestUpdat(t *testing.T) {
	pth := "../tdata/"
	var idmkeys = map[string]string{
		"Gene_Name":    "gnm",
		"Gene_Synonym": "gsnm",
		//"Gene_ORFName": "gsnm", // no impact on triple number TODO
		"Ensembl_PRO": "ensp",
		"Ensembl":     "ensg",
		"GeneID":      "ncbig",
		"RefSeq":      "rfsq",
		"UniParc":     "uparc",
	}
	upcas, upacs, gnms, _ := parse.Upidmap(pth+"test.idm", idmkeys)
	updat, txns, _ := parse.Updat(pth+"test.upt", upcas)
	var arg1 util.Dat4rdf
	arg1.Udat = &updat
	arg1.Txns = &txns
	arg1.Upac = &upacs
	arg1.Upca = &upcas
	arg1.Gnm = &gnms
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	t1s := []t1{
		{arg1, "/dev/null", "/dev/null", "/dev/null", zeno, 57, 1888},
	}
	for i, t1 := range t1s {
		n, m, _ := GeneProt(t1.arg1, t1.arg2, t1.arg3, t1.arg4, t1.arg5)
		if n != t1.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", t1.val1,
				"\n\thave", n,
			)
		}
		if m != t1.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", t1.val2,
				"\n\thave", m,
			)
		}
	}
}

func TestTftg(t *testing.T) {
	pth := "../tdata/"
	arg3 := make(aux.Set3D)
	arg3.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	arg4 := make(aux.Set3D)
	arg4.Add("TP53", "bgwg", "9606/chr-17/TP53")
	arg1, arg2 := parse.Tftg(pth+"test.f2g", arg3, arg4)
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	t3s := []t3{
		{arg1, arg2, arg3, arg4, "/dev/null", zeno, 350},
	}
	for i, t3 := range t3s {
		n, _ := Tftg(t3.arg1, t3.arg2, t3.arg3, t3.arg4, t3.arg5, t3.arg6)
		if n != t3.val1 {
			t.Error(
				"For test", i+1, ": ", len(t3.arg1), t3.arg3, t3.arg4,
				"\n\twant", t3.val1,
				"\n\thave", n,
			)
		}
	}
}

func TestUpvar(t *testing.T) {
	pth := "../tdata/"
	arg2 := make(aux.Set3D)
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	arg3 := make(aux.Set3D)
	arg3.Add("TP53", "bgwg", "9606/chr-17/TP53")
	arg1 := parse.Upvar(pth+"test.var", arg3)
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	t2s := []t2{
		{arg1, arg2, arg3, "/dev/null", zeno, 28},
	}
	for i, t2 := range t2s {
		n, err := Upvar(t2.arg1, t2.arg2, t2.arg3, t2.arg4, t2.arg5)
		if err != nil {
			return
		}
		if n != t2.val1 {
			t.Error(
				"For test", i+1, ": ", len(t2.arg1), t2.arg2, t2.arg3,
				"\n\twant", t2.val1,
				"\n\thave", n,
			)
		}
	}
}

func TestMitab(t *testing.T) {
	pth := "../tdata/"
	arg2 := make(aux.Set3D)
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	arg1 := parse.Mitab(pth+"test.mit", arg2)
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	t4s := []t4{
		{arg1, arg2, "/dev/null", zeno, 95},
	}
	for i, t4 := range t4s {
		n, err := Mitab(t4.arg1, t4.arg2, t4.arg3, t4.arg4)
		if err != nil {
			return
		}
		if n != t4.val1 {
			t.Error(
				"For test", i+1, ": ", len(t4.arg1), t4.arg2, t4.arg3,
				"\n\twant", t4.val1,
				"\n\thave", n,
			)
		}
	}
}

func TestGaf(t *testing.T) {
	pth := "../tdata/"
	arg2 := make(aux.Set3D)
	arg2.Add("P04637", "bgwp", "9606/chr-17/TP53/UPI000002ED67")
	bps, ccs, mfs := parse.Gaf(pth+"test.gaf", arg2)
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	t4s := []t4{
		{bps, arg2, "/dev/null", zeno, 1836},
		{ccs, arg2, "/dev/null", zeno, 268},
		{mfs, arg2, "/dev/null", zeno, 839},
	}
	for i, t4 := range t4s {
		n, err := Goa(t4.arg1, t4.arg2, t4.arg3, t4.arg4)
		if err != nil {
			return
		}
		if n != t4.val1 {
			t.Error(
				"For test", i+1, ": ", len(t4.arg1), t4.arg2, t4.arg3,
				"\n\twant", t4.val1,
				"\n\thave", n,
			)
		}
	}
}
