package util
import (
	"testing"
)

func TestSet1DKeys(t *testing.T) {
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
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func TestSet2DAdd(t *testing.T) {
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
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func TestSet2DKeys(t *testing.T) {
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
				"For test", i+1, ": ", tt.arg1, tt.arg2,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func TestSet3DAdd(t *testing.T) {
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
				"For test", i+1, ": ", tt.arg1, tt.arg2, tt.arg3,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

func TestSet3DKeys(t *testing.T) {
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

func TestSliceSetAdd(t *testing.T) {
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

func TestX1type(t *testing.T) {
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

func TestMakemap(t *testing.T) {
	pth := "../../tdata/"
	type tt struct {
		arg1 string
		arg2 int
		arg3 int
		arg4 string
		val1 int
	}
	tts := []tt{
		{pth + "test.idm", 0, 1, "\t", 9},
	}
	for i, tt := range tts {
		r, _ := Makemap(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
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

