package main

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"strings"
)

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func subset(pthR, pthW string, inds [5]int) error {
	fhR, err := os.Open(pthR)
	check(err)
	defer fhR.Close()
	fhW, err := os.Create(pthW)
	check(err)
	defer fhW.Close()
	scanner := bufio.NewScanner(fhR)
	for scanner.Scan() { // by default scans for '\n'
		line := scanner.Text()
		if len(line) == 0 {
			continue
		}
		if string(line[0]) == "#" {
			continue
		}
		cells := strings.Split(line, "\t")
		n := len(cells)
		m := 37
		if n < m {
			msg := "Want at least: %d cells, have: %d for: %s in %s"
			fmt.Printf(msg, m, n, cells[0], pthR)
			continue
		}
		subline := cells[0] // pair ID, e.g. MYC:TERT
		// adding sorce specific fields
		for _, i := range inds {
			if i == 0 { // the value not provided by the source
				subline = strings.Join([]string{subline, ""}, "\t")
			} else {
				subline = strings.Join([]string{subline, cells[i-1]}, "\t")
			}
		}
		if len(subline) == 0 {
			msg := fmt.Sprintf("%s", "EmptyString")
			panic(errors.New(msg))
		}
		fhW.WriteString(fmt.Sprintf("%s\n", subline))
	}
	return nil
}

// columns
// 1: TF IDs, 2: TG IDs 3: PUBMED IDs 4: evidence scores 5: modes
// Note: currently a single TG in the whole data set
var src2ind = map[string][5]int{
	"cytreg": {2, 3, 33, 0, 32},
	"extri":  {2, 3, 5, 4, 0},
	"geredb": {2, 3, 37, 0, 36},
	//"goa":    {2, 3, 0, 0, 21}, // no PUBMED refs
	"htri":   {2, 3, 9, 10, 0},
	"intact": {2, 3, 23, 0, 0},
	"ntnu": {2, 3, 40, 0, 39},
	"signor": {2, 3, 28, 0, 27},
	"tfacts": {2, 3, 18, 19, 15},
	"trrust": {2, 3, 13, 0, 12},
}

func main() {
	datdir := os.Args[1] // location of source specific files (with a trailing slash)
	for src, inds := range src2ind {
		pthR := fmt.Sprintf("%s%s/9606.ori", datdir, src)
		pthW := fmt.Sprintf("%s%s/9606.f2g", datdir, src)
		err := subset(pthR, pthW, inds)
		check(err)
	}
}
