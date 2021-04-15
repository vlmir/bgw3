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
		subline := cells[0]
		for _, i := range inds {
			if i == 0 {
				subline = strings.Join([]string{subline, ""}, "\t")
			} else {
				subline = strings.Join([]string{subline, cells[i-1]}, "\t")
			}
		}
		if len(subline) == 0 {
			msg := fmt.Sprintf("%s", "EmptyString")
			panic(errors.New(msg))}
		fhW.WriteString(fmt.Sprintf("%s\n", subline))
	}
	return nil
}

var src2ind = map[string][5]int{
	"cytreg": {2, 3, 33, 0, 32},
	"extri":  {2, 3, 5, 4, 0},
	"geredb": {2, 3, 37, 0, 36},
	"goa":    {2, 3, 0, 0, 21},
	"htri":   {2, 3, 9, 10, 0},
	"intact": {2, 3, 23, 0, 0},
	"signor": {2, 3, 28, 0, 27},
	"tfacts": {2, 3, 18, 19, 15},
	"trrust": {2, 3, 13, 0, 12},
}

func main() {
	datdir := "/home/mironov/data/bgw/tftg/9606/"
	for src, inds := range src2ind {
		pthR := fmt.Sprintf("%s%s.tsv", datdir, src)
		pthW := fmt.Sprintf("%s%s.f2g", datdir, src)
		err := subset(pthR, pthW, inds)
		check(err)
	}
}
