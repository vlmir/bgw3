package util

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"sort"
	"strings"
)

// TODO trim all strings ? better to be done in the parsers

type Set1D map[string]int
type Set2D map[string]Set1D
type Set3D map[string]Set2D
type Set4D map[string]Set3D
type SliceSet map[string][]string

// FN retrieves the names of the running||calling functions
func FN(skip int) string {
	pc, _, _, ok := runtime.Caller(skip + 1)
	if !ok {
		return ""
	}
	f := runtime.FuncForPC(pc)
	if f == nil {
		return ""
	}
	return f.Name()
}

func TrimString(p *string) error {
	*p = strings.TrimSpace(*p)
	if *p == "" {
		msg := fmt.Sprintf("util.TrimString(): EmptyString")
		return errors.New(msg)
	}
	return nil
}

func (m Set1D) Add(key0 string) {
	p0 := &key0
	err := TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("util.(Set1D)Add(%s): %s", "key0", err)
		panic(errors.New(msg))
	}
	_, ok := m[*p0]
	if !ok {
		m[*p0]++
	}
}
func (m Set2D) Add(key0, key1 string) {
	p0 := &key0
	err := TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("util.(Set2D)Add(%s, _): %s", "key0", err)
		panic(errors.New(msg))
	}
	p1 := &key1
	err = TrimString(p1)
	if err != nil {
		msg := fmt.Sprintf("util.(Set2D)Add(_, %s): %s", "key1", err)
		panic(errors.New(msg))
	}
	mm, ok := m[*p0]
	if !ok {
		mm = make(map[string]int)
		m[*p0] = mm
	}
	mm[*p1]++
}
func (m Set3D) Add(key0, key1, key2 string) {
	p0 := &key0
	err := TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("util.(Set3D)Add(%s, _, _): %s", "key0", err)
		panic(errors.New(msg))
	}
	p1 := &key1
	err = TrimString(p1)
	if err != nil {
		msg := fmt.Sprintf("util.(Set3D)Add(_, %s, _): %s", "key1", err)
		panic(errors.New(msg))
	}
	p2 := &key2
	err = TrimString(p2)
	if err != nil {
		msg := fmt.Sprintf("util.(Set3D)Add(_, _, %s): %s", "key2", err)
		panic(errors.New(msg))
	}
	mm, ok := m[*p0]
	if !ok {
		mm = make(Set2D)
		m[*p0] = mm
	}
	mmm, ok := mm[*p1]
	if !ok {
		mmm = make(Set1D)
		mm[*p1] = mmm
	}
	mmm[*p2]++
}
func (m Set4D) Add(key0, key1, key2, key3 string) {
	p0 := &key0
	err := TrimString(p0)
	if err != nil {
		msg := fmt.Sprintf("util.(Set4D)Add(%s, _, _, _): %s", "key0", err)
		panic(errors.New(msg))
	}
	p1 := &key1
	err = TrimString(p1)
	if err != nil {
		msg := fmt.Sprintf("util.(Set4D)Add(_, %s, _, _): %s", "key1", err)
		panic(errors.New(msg))
	}
	p2 := &key2
	err = TrimString(p2)
	if err != nil {
		msg := fmt.Sprintf("util.(Set4D)Add(_, _, %s, _): %s", "key2", err)
		panic(errors.New(msg))
	}
	p3 := &key3
	err = TrimString(p3)
	if err != nil {
		msg := fmt.Sprintf("util.(Set4D)Add(_, _, _, %s): %s", "key3", err)
		panic(errors.New(msg))
	}
	mm, ok := m[*p0]
	if !ok {
		mm = make(Set3D)
		m[*p0] = mm
	}
	mmm, ok := mm[*p1]
	if !ok {
		mmm = make(Set2D)
		mm[*p1] = mmm
	}
	mmmm, ok := mmm[*p2]
	if !ok {
		mmmm = make(Set1D)
		mmm[*p2] = mmmm
	}
	mmmm[*p3]++
}
func (s SliceSet) Add(key, value string) {
	_, ok := s[key]
	if !ok {
		s[key] = make([]string, 0, 20)
	}
	s[key] = append(s[key], value)
}

func (m Set1D) Keys() []string {
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}
func (m Set2D) Keys() []string {
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}
func (m Set3D) Keys() []string {
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}
func (m Set4D) Keys() []string {
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}
func (m SliceSet) Keys() []string {
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// not used
func (s SliceSet) Peek(key string) (string, bool) {
	slice, ok := s[key]
	if !ok || len(slice) == 0 {
		return "", false
	}
	return s[key][0], true
}

/*
input: slice of srings containig pairs of values, all delimited by the same char
*/
func X1type(s []string, key string, dlm string) []string {
	var out []string
	for _, cell := range s {
		bits := strings.Split(cell, dlm)
		if len(bits) != 2 {
			continue
		}
		if bits[0] == key {
			out = append(out, bits[1])
		}
	}
	return out
}

/*
Input: slice of strings containing double quoted OBO IDs
not used
*/
func X1quoted(s []string, key string, dlm string) []string {
	var out []string
	for _, cell := range s {
		bits := strings.Split(cell, "\"")
		if len(bits) <= 2 {
			continue
		}
		shreds := strings.Split(bits[0], dlm)
		if shreds[0] == key {
			id := bits[1]
			out = append(out, id)
		}
	}
	return out
}

// MakeMap generates a map from a file delimited by an arbitrary string
func MakeMap(pth string, keyind int, valind int, dlm string) (Set2D, error) {
	set := make(Set2D)
	fh, err := os.Open(pth)
	if err != nil {
		return set, err
	}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	for scanner.Scan() { // by default scans for '\n'
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}
		if string(line[0]) == "#" {
			continue
		}
		cells := strings.Split(line, dlm)
		if keyind >= len(cells) {
			continue
		}
		if valind >= len(cells) {
			continue
		}
		key := cells[keyind]
		val := cells[valind]
		set.Add(key, val)
	}
	return set, nil
}

func StripParQuots(s string) string {
	l := len(s)
	s = s[1 : l-1]
	return s
}
func Index(vs []string, t string) int {
	for i, v := range vs {
		if v == t {
			return i
		}
	}
	return -1
}

func Includes(vs []string, t string) bool {
	return Index(vs, t) >= 0
}

func Shared(slA, slB []string) []string {
	s := make([]string, 0)
	for _, strA := range slA {
		for _, strB := range slB {
			if strA == strB {
				s = append(s, strA)
			}
		}
	}
	return s
}

// write a test
func IsDigital(s string) bool {
	// https://programming-idioms.org/idiom/137/check-if-string-contains-only-digits/1739/go
	b := true
	for _, c := range s {
		if c < '0' || c > '9' {
			b = false
			break
		}
	}
	return b
}

func Gzip(pth string) error {
	// TODO tests
	var cmd *exec.Cmd
	// Attn: all flags separately!
	cmd = exec.Command("gzip", "-f", pth) // struct
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run() // returns error
}
