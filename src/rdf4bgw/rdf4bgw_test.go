package main

import (
	"github.com/vlmir/bgw3/src/util" // pkg 'util'
	"testing"
)

func Test_geneprot(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		val1 int
		val2 int
	}

	pth := "../../tdata/"
	wpth := pth + "OUT/rdf4bgw/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "UP000005640")
	//arg3[0].Add("7227", "UP000000803")
	tts := []tt{
		//		{pth, wpth, arg3[0], 35, 76},
		{pth, wpth, arg3[0], 50, 68},
	}

	for i, tt := range tts {
		ng, np, _ := geneprot(tt.arg1, tt.arg2, tt.arg3)
		if ng != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", ng,
			)
		}
		if np != tt.val2 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val2,
				"\n\thave", np,
			)
		}
	}
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
		{pth, wpth, arg3[0], 2}, // Rgr2trg(): 4, Tfac2gene(): 2
	}
	pdck := "preg2targ"
	src := "tfacts"
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
		{pth, wpth, arg3[0], 12},
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
