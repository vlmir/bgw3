package main

import (
	"fmt"
	"github.com/vlmir/bgw3/src/util" // pkg 'util'
	"github.com/vlmir/bgw3/src/semweb"
	"testing"
)

func Test_orthoduo(t *testing.T){
	type tt struct {
		arg1 string
		arg2 string
		arg3 util.Set2D
		arg4 util.Set2D
		arg5 rdf.Zeno
		val1 int
	}
	pth := "../../tdata/"
	var arg3 [5]util.Set2D // txmap
	arg3[0] = make(util.Set2D)
	arg3[0].Add("9606", "1")
	arg3[0].Add("10090", "1")
	var arg4 [5]util.Set2D // tx2pm
	arg4[0] = make(util.Set2D)
	arg4[0].Add("9606", "ortho")
	arg4[0].Add("10090", "ortho")
	zeno := rdf.NewZeno()
	zeno.Unmarshal(pth + "zeno.json")
	tts := []tt{
		{pth, pth, arg3[0], arg4[0], zeno, 11},
	}
	for i, tt := range tts {
		n, err := orthoduo(tt.arg1, tt.arg2, tt.arg3, tt.arg4, tt.arg5)
		fmt.Println("err:", err)
		if n != tt.val1 {
			t.Error(
				"For test", i+1, ": ", len(tt.arg1), tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", n,
			)
		}
	}
}

