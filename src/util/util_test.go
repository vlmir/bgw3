package util

import (
	"fmt"
	"testing"
)

func Test_TrimString(t *testing.T) {
	s := "\n\n"
	if err := TrimString(&s); err != nil {
		fmt.Println(err)
	}
}

func Test_Set1D_Add(t *testing.T) {
	type tt struct {
		arg1 string
		val1 int
	}
	tts := []tt{
		{"key0", 1},
	}
	for i, tt := range tts {
		set := make(Set1D)
		set.Add(tt.arg1)
		val1 := set[tt.arg1]
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Set1D_Keys(t *testing.T) {
	type tt struct {
		arg1 string
		val1 int
	}
	tts := []tt{
		{"key0", 1},
	}
	for i, tt := range tts {
		set := make(Set1D)
		set[tt.arg1] = 1
		val1 := len(set.Keys())
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Set2D_Add(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		val1 int
	}
	tts := []tt{
		{"key0", "key1", 1},
	}
	for i, tt := range tts {
		set := make(Set2D)
		set.Add(tt.arg1, tt.arg2)
		val1 := set[tt.arg1][tt.arg2]
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Set2D_Keys(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		val1 int
	}
	tts := []tt{
		{"key0", "key1", 1},
	}
	for i, tt := range tts {
		set := make(Set2D)
		set.Add(tt.arg1, tt.arg2)
		val1 := len(set.Keys())
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Set3D_Add(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 string
		val1 int
	}
	tts := []tt{
		{"key0", "key1", "key3", 1},
	}
	for i, tt := range tts {
		set := make(Set3D)
		set.Add(tt.arg1, tt.arg2, tt.arg3)
		val1 := set[tt.arg1][tt.arg2][tt.arg3]
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Set4D_Add(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 string
		arg4 string
		val1 int
	}
	tts := []tt{
		{"key0", "key1", "key3", "key4", 1},
	}
	for i, tt := range tts {
		set := make(Set4D)
		set.Add(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		val1 := set[tt.arg1][tt.arg2][tt.arg3][tt.arg4]
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_Set3D_Keys(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		arg3 string
		val1 int
	}
	tts := []tt{
		{"key0", "key1", "key3", 1},
	}
	for i, tt := range tts {
		set := make(Set3D)
		set.Add(tt.arg1, tt.arg2, tt.arg3)
		val1 := len(set.Keys())
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_SliceSet_Add(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 string
		val1 int
	}
	tts := []tt{
		{"key", "str", 1},
	}
	for i, tt := range tts {
		set := make(SliceSet)
		set.Add(tt.arg1, tt.arg2)
		val1 := len(set[tt.arg1])
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func Test_X1type(t *testing.T) {
	type tt struct {
		arg1 []string
		arg2 string
		arg3 string
		val1 int
	}
	s := []string{"one|a", "two|b", "one|c", "two|d"}
	tts := []tt{
		{s, "one", "|", 2},
	}
	for i, tt := range tts {
		r := X1type(tt.arg1, tt.arg2, tt.arg3)
		val1 := len(r)
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

//func Test_FilterByValues(t *testing.T) {
//	type tt struct {
//		arg1 string
//		arg2 map[string]string
//		arg3 int
//		arg4 int
//		arg5 int
//		val  int
//	}
//	pth := "../../tdata/idmapping/"
//	t1s := []tt{
//		{pth + "UP000005640_9606.idmapping", map[string]string{"NCBI_TaxID": "test"}, 2, 1, 0, 1},
//		{pth + "UP000005640_9606.idmapping", map[string]string{"NCBI_TaxID": "test"}, 0, 1, 2, 1},
//		{pth + "UP000005640_9606.idmapping", map[string]string{"UniParc": "test"}, 0, 1, 2, 5},
//	}
//	for i, tt := range t1s {
//		idm, _ := FilterByValues(tt.arg1, tt.arg2, tt.arg3, tt.arg4, tt.arg5)
//		if len(idm) != tt.val {
//			t.Error(
//				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3, tt.arg4, tt.arg5,
//				"\n\twant", tt.val,
//				"\n\thave", len(idm),
//			)
//		}
//	}
//}

func Test_MakeMap(t *testing.T) {
	pth := "../../tdata/idmapping/"
	type tt struct {
		arg1 string
		arg2 int
		arg3 int
		arg4 string
		val1 int
	}
	tts := []tt{
		{pth + "UP000005640_9606.idmapping", 0, 1, "\t", 8},
	}
	for i, tt := range tts {
		r, _ := MakeMap(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		val1 := len(r)
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ":[", tt.arg1, tt.arg2, tt.arg3, tt.arg4,
				"]\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}
