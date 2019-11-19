package aux

import (
	"bufio"
	"os"
	"strings"
)

type Set1D map[string]int

//type Set2D map[string]map[string]int
type Set2D map[string]Set1D
type Set3D map[string]Set2D
type SliceSet map[string][]string

func (m Set2D) Add(key0, key1 string) {
	mm, ok := m[key0]
	if !ok {
		mm = make(map[string]int)
		m[key0] = mm
	}
	mm[key1]++
}

func (m Set3D) Add(key0, key1, key2 string) {
	mm, ok := m[key0]
	if !ok {
		//mm = make(map[string]map[string]int)
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

func Makemap(pth string, keyind int, valind int, dlm string) (Set2D, error) {
	set := make(Set2D)
	fh, err := os.Open(pth)
	if err != nil {
		return set, err
	}
	defer fh.Close()
	scanner := bufio.NewScanner(fh)
	for scanner.Scan() { // by default scans for '\n'
		cells := strings.Split(scanner.Text(), dlm)
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

func (s SliceSet) Add(key, value string) {
	_, ok := s[key]
	if !ok {
		s[key] = make([]string, 0, 20)
	}
	s[key] = append(s[key], value)
}

//not used
func (s SliceSet) Peek(key string) (string, bool) {
	slice, ok := s[key]
	if !ok || len(slice) == 0 {
		return "", false
	}
	return s[key][0], true
}

func (m Set1D) Keys() []string {
	// extracting map keys
	//keys := make([]string, 0, len(m))
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	return keys
}

func (m Set2D) Keys() []string {
	// extracting map keys
	//keys := make([]string, 0, len(m))
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	return keys
}

func (m Set3D) Keys() []string {
	// extracting map keys
	//keys := make([]string, 0, len(m))
	var keys []string
	for k := range m {
		keys = append(keys, k)
	}
	return keys
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
