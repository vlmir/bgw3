package main

import (
	"github.com/vlmir/bgw3/src/util" // pkg 'util'
	"testing"
)

// must be the first one !
func Test_geneprot(t *testing.T) {
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
		err := geneprot(tt.arg1, tt.arg2, tt.arg3)
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

func Test_gene2phen(t *testing.T) {
	// TODO
}

func Test_rgr2trg(t *testing.T) {
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
		{pth, wpth, arg3[0], 0}, // used to be '1' TODO
	}
	pdck := "reg2ntrg" // sic! positive interactions eliminated due to xmap
	src := "signor"
	for i, tt := range tts {
		cnts, _ := rgr2trg(tt.arg1, tt.arg2, tt.arg3) // process ALL files in src dir
		if cnts[pdck][src] != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", cnts[pdck][src],
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
		{pth, wpth, arg3[0], 3}, // new tests
	}
	pdck := "reg2ptrg"
	src := "tfacts"
	// src := "tflink" // no interactions after filtering by xmap
	for i, tt := range tts {
		cnts, _ := tfac2gene(tt.arg1, tt.arg2, tt.arg3) // process ALL files in src dir
		if cnts[pdck][src] != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", cnts[pdck][src],
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
