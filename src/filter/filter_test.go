package main

import (
	"testing"
)

type tt struct {
	arg1 string
	arg2 string
	arg3 [5]int
	val1 error
}

func Test_subset(t *testing.T) {
	var fields = map[string][5]int{
		"htri": {2, 3, 9, 10, 0},
	}
	pth := "../../tdata/"
	tts := []tt{
		{pth + "htri.ori", pth + "output/filter/htri.f2g", fields["htri"], nil},
	}
	for i, tt := range tts {
		err := subset(tt.arg1, tt.arg2, tt.arg3)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
}
