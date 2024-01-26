// package main
package rdf4bgw

import (
	"github.com/vlmir/bgw3/src/util" // pkg 'util'
	"testing"
)

// must be the first one !
func Test_Geneprot(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		val1 error
	}

	pth := "../../tdata/"
	wpth := pth + "OUT/rdf4bgw/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "UP000005640")
	//arg3[0].Add("7227", "UP000000803")
	tts := []tt{
		// {pth, wpth, arg3[0], 35, 76},
		// {pth, wpth, arg3[0], 50, 68},
		{pth, wpth, arg3[0], nil},
	}

	for i, tt := range tts {
		err := Geneprot(tt.arg1, tt.arg2, tt.arg3)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
		//		if np != tt.val2 {
		//			t.Error(
		//				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
		//				"\n\twant", tt.val2,
		//				"\n\thave", np,
		//			)
		//		}
	}
}

func Test_Prot2go(t *testing.T) {
	// TODO
}

func Test_Gene2phen(t *testing.T) {
	// TODO
}

func Test_Reg2pway(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		val  int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/rdf4bgw/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "testprome")
	tts := []tt{
		{pth, wpth, arg3[0], 1}, // p53-bcl2 interaction
	}
	pdck := "reg2ntrg"
	srck := "signor"
	for i, tt := range tts {
		cnts, _ := Reg2pway(tt.arg1, tt.arg2, tt.arg3) // process ALL files in srck dir
		if cnts[pdck][srck] != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", cnts[pdck][srck],
			)
		}
	}
}

func Test_Rgr2trg(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		val  int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/rdf4bgw/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "testprome")
	tts := []tt{
		{pth, wpth, arg3[0], 1}, // p53-bcl2
	}
	pdck := "reg2ntrg"
	srck := "signor"
	for i, tt := range tts {
		cnts, _ := Rgr2trg(tt.arg1, tt.arg2, tt.arg3) // process ALL files in srck dir
		if cnts[pdck][srck] != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", cnts[pdck][srck],
			)
		}
	}
}

func Test_tfac2gene(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		val  int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/rdf4bgw/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "testprome")
	tts := []tt{
		// {pth, wpth, arg3[0], 2},
		// {pth, wpth, arg3[0], 4}, // new tests
		{pth, wpth, arg3[0], 1}, // new tests
	}
	pdck := "reg2ptrg"
	// srck := "tfacts"
	srck := "coltri"
	// srck := "tflink" // no interactions after filtering by xmap
	for i, tt := range tts {
		cnts, _ := tfac2gene(tt.arg1, tt.arg2, tt.arg3) // process ALL files in srck dir
		if cnts[pdck][srck] != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", cnts[pdck][srck],
			)
		}
	}
}

func Test_Prot2prot(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		val  int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/rdf4bgw/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "testprome")
	tts := []tt{
		{pth, wpth, arg3[0], 1},
	}
	pdck := "tlp2tlp"
	srck := "intact"
	for i, tt := range tts {
		cnts, _ := Prot2prot(tt.arg1, tt.arg2, tt.arg3) // process ALL files in srck dir
		if cnts[pdck][srck] != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", cnts[pdck][srck],
			)
		}
	}
}

func Test_ortho(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		val1 int
	}
	pth := "../../tdata/"
	wpth := pth + "OUT/rdf4bgw/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "ortho")
	arg3[0].Add("10090", "ortho")
	tts := []tt{
		{pth, wpth, arg3[0], 28},
	}
	for i, tt := range tts {
		n, _ := ortho(tt.arg1, tt.arg2, tt.arg3)
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
	}
}
