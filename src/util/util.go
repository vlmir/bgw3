package util

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"sort"
	"strings"
)

// TODO trim all strings ? better to be done in the parsers

type Set1D map[string]int
type Set2D map[string]Set1D
type Set3D map[string]Set2D
type Set4D map[string]Set3D
type SliceSet map[string][]string

// TODO return error instead of panicing ?
func CheckStrings(s ...string) error {
	// just checks, does ot trim
	for i, v := range s {
		if strings.TrimSpace(v) == "" {
			// panic(errors.New(fmt.Sprintf("%v: index: %d: EmptyString", s, i)))
			return errors.New(fmt.Sprintf("%v: index: %d: EmptyString", s, i))
		}
	}
	return nil
}

func (m Set1D) Add(key0 string) {
	CheckStrings(key0)
	_, ok := m[key0]
	if !ok {
		m[key0]++
	}
}

func (m Set1D) Keys() []string {
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

func (m Set2D) Add(key0, key1 string) {
	CheckStrings(key0, key1)
	mm, ok := m[key0]
	if !ok {
		mm = make(map[string]int)
		m[key0] = mm
	}
	mm[key1]++
}

func (m SliceSet) Keys() []string {
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

func (m Set3D) Add(key0, key1, key2 string) {
	CheckStrings(key0, key1, key2)
	mm, ok := m[key0]
	if !ok {
		mm = make(Set2D)
		m[key0] = mm
	}
	mmm, ok := mm[key1]
	if !ok {
		mmm = make(Set1D)
		mm[key1] = mmm
	}
	mmm[key2]++
}

func (m Set4D) Add(key0, key1, key2, key3 string) {
	CheckStrings(key0, key1, key2, key3)
	mm, ok := m[key0]
	if !ok {
		mm = make(Set3D)
		m[key0] = mm
	}
	mmm, ok := mm[key1]
	if !ok {
		mmm = make(Set2D)
		mm[key1] = mmm
	}
	mmmm, ok := mmm[key2]
	if !ok {
		mmmm = make(Set1D)
		mmm[key2] = mmmm
	}
	mmmm[key3]++
}

func (m Set4D) Keys() []string {
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

func (m Set1D) Check() {
	if len(m) == 0 {
		panic(errors.New("Empty map!"))
	}
}

func (m Set2D) Check() {
	if len(m) == 0 {
		panic(errors.New("Empty map!"))
	}
}

func (m Set3D) Check() {
	if len(m) == 0 {
		panic(errors.New("Empty map!"))
	}
}

func (m Set4D) Check() {
	if len(m) == 0 {
		panic(errors.New("Empty map!"))
	}
}

func (s SliceSet) Add(key, value string) {
	_, ok := s[key]
	if !ok {
		s[key] = make([]string, 0, 20)
	}
	s[key] = append(s[key], value)
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

// FilterByValues filters a tab-delimited file by values in a specified field
//func FilterByValues(rpth string, srcs map[string]string, ind1, ind2, ind3 int) (Set3D, error) {
//	// used only in parse.orthosolo() and dat4bgw
//	// looks like a problem with memory leakage TODO fis
//	out := make(Set3D)
//	fh, err := os.Open(rpth)
//	CheckE(err)
//	defer fh.Close()
//	scanner := bufio.NewScanner(fh)
//	indmax := 0
//	inds := [3]int{ind1, ind2, ind3}
//	for _, v := range inds {
//		if v > indmax {
//			indmax = v
//		}
//	}
//	for scanner.Scan() { // by default scans for '\n'
//		cells := strings.Split(scanner.Text(), "\t")
//		if len(cells) <= indmax {
//			continue
//		}
//		_, ok := srcs[cells[ind2]] // filtering
//		if !ok {
//			continue
//		} // filtering by srcs
//		key1 := strings.Replace(cells[ind1], "\"", "''", -1) // present e.g. in 44689
//		key2 := strings.Replace(cells[ind2], "\"", "''", -1) // present e.g. in 44689
//		key3 := strings.Replace(cells[ind3], "\"", "''", -1) // present e.g. in 44689
//		out.Add(key1, key2, key3)
//	}
//	return out, nil
//}

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

func CheckE(e error) {
	if e != nil {
		panic(e)
	}
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
