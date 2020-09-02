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
	xpth := pth + "output/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "UP000005640")
	tts := []tt{
		{pth, xpth, arg3[0], 35, 76},
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

func Test_tfac2gene(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		arg4 map[string]string
		val  int
	}

	pth := "../../tdata/"
	xpth := pth + "output/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "1")
	var arg4 [5]map[string]string
	arg4[0] = map[string]string{
		"tfacts": "http://www.tfacts.org",
		"htri":   "http://www.htri.org",
	}
	tts := []tt{
		{pth, xpth, arg3[0], arg4[0], 41},
	}

	for i, tt := range tts {
		n, _ := tfac2gene(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		if n != tt.val {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val,
				"\n\thave", n,
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
	xpth := pth + "output/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "ortho")
	arg3[0].Add("10090", "ortho")
	tts := []tt{
		{pth, xpth, arg3[0], 11},
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
