package bgw
import (
	"testing"
)

func Test_Unmarshal(t *testing.T) {
	pth := "../../tdata/"
	type tt struct {
		arg1 string
		val1 int
	}
	tts := []tt{
		{pth + "xmap.json", 9},
	}
	for i, tt := range tts {
		xmap := NewXmap()
		xmap.Unmarshal(tt.arg1)
		val1 := len(xmap.Bgwp)
		if val1 != tt.val1 {
			t.Error(
				"For test", i+1, ": ", tt.arg1,
				"\n\twant", tt.val1,
				"\n\thave", val1,
			)
		}
	}
}

