package main

import (
	"log"
	"os"
	"path/filepath"
	"testing"
	"time"
)

func Test_Init(t *testing.T) {
	pth := "../../tdata/"
	dpth := filepath.Join(pth, "dat4bgw/")
	if err := os.RemoveAll(dpth); err != nil {
		log.Println(err)
	}
	// NB: the leading '0' is necessary !
	if err := os.MkdirAll(filepath.Join(dpth, "uniprot/"), 0755); err != nil {
		log.Println(err)
	}
	if err := os.MkdirAll(filepath.Join(dpth, "goa/"), 0755); err != nil {
		log.Println(err)
	}
	if err := os.MkdirAll(filepath.Join(dpth, "intact/"), 0755); err != nil {
		log.Println(err)
	}
}

/*
 */
func Test_getOneUniprot(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		val1 error
	}
	pth := "../../tdata/"
	txids := [...]string{"1343077", "2018007"}
	datdirs := [...]string{pth + "dat4bgw/", pth + "dat4bgw/"}
	tts := []tt{
		{txids[0], datdirs[0], nil},
		//{txids[1], datdirs[1], nil},
	}
	for i, tt := range tts {
		err := saveOneUniprot(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with saveOneUniprot in", time.Since(mystart))
}

func Test_getOneIntact(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		val1 error
	}
	pth := "../../tdata/"
	txids := [...]string{"1672772", "9031"}
	datdirs := [...]string{pth + "dat4bgw/", pth + "dat4bgw/"}
	tts := []tt{
		{txids[0], datdirs[0], nil},
		{txids[1], datdirs[1], nil},
	}
	for i, tt := range tts {
		err := saveOneIntact(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with saveOneIntact in", time.Since(mystart))
}

func Test_getOneGaf(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		arg3 string
		val1 error
	}
	pth := "../../tdata/"
	txids := [...]string{"1672772"}
	datdirs := [...]string{pth + "dat4bgw/", pth + "dat4bgw/"}
	pomes := [...]string{"493.P_falciparum.goa"}
	tts := []tt{
		{txids[0], datdirs[0], pomes[0], nil},
	}
	for i, tt := range tts {
		err := saveOneGaf(tt.arg1, tt.arg2, tt.arg3)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with saveOneGaf in", time.Since(mystart))
}
func Test_getOneGpa(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		val1 error
	}
	pth := "../../tdata/"
	txids := [...]string{"1672772"}
	datdirs := [...]string{pth + "dat4bgw/", pth + "dat4bgw/"}
	tts := []tt{
		{txids[0], datdirs[0], nil},
	}
	for i, tt := range tts {
		err := saveOneGpa(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with saveOneGpa in", time.Since(mystart))
}

func Test_HttpFile(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		val1 error
		val2 int
	}

	//pth := "../../tdata/"
	//dpth := pth + "dat4bgw/"
	ns01 := "http://purl.obolibrary.org/obo/"
	uri01 := ns01 + "bfo.owl"
	ns02 := "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:"
	uri02 := ns02 + "1672772"
	uris := [...]string{uri01, uri02}
	//pths :=[...]string{dpth+"HttpFile.01", dpth+"HttpFile.02"}
	pths := [...]string{"", ""}

	tts := []tt{
		{uris[0], pths[0], nil, 157931},
		{uris[1], pths[1], nil, 0},
	}
	for i, tt := range tts {
		buf, err := HttpFile(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
		if len(buf.String()) != tt.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val2,
				"\n\thave", buf,
			)
		}
	}
	log.Println("Done with HttpFile in", time.Since(mystart))
}

func Test_GetFile(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		arg3 string
		arg4 string
		val1 error
		val2 int
	}

	//pth := "../../tdata/"
	//dpth := pth + "dat4bgw/"
	//uri = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?proteome=gcrpCan%2CgcrpIso&geneProductId=P04637&geneProductType=protein&taxonId=9606"
	ns01 := "http://purl.obolibrary.org/obo/"
	uri01 := ns01 + "bfo.owl"
	ns02 := "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:"
	uri02 := ns02 + "1672772"
	uris := [...]string{uri01, uri02}
	hact01 := "Accept"
	htype01 := "text"
	hact02 := "Accept"
	htype02 := "text"
	hacts := [...]string{hact01, hact02}
	htypes := [...]string{htype01, htype02}
	pths := [...]string{"", ""}
	//pths :=[...]string{dpth+"GetFile.01", dpth+"GetFile.02"}

	tts := []tt{
		{uris[0], hacts[0], htypes[0], pths[0], nil, 157931},
		{uris[1], hacts[1], htypes[1], pths[1], nil, 0},
	}
	for i, tt := range tts {
		buf, err := GetFile(tt.arg1, tt.arg2, tt.arg3, tt.arg4)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
		if len(buf.String()) != tt.val2 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val2,
				"\n\thave", buf,
			)
		}
	}
	log.Println("Done with GetFile in", time.Since(mystart))
}

func Test_WgetFile(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		val1 error
	}

	pth := "../../tdata/"
	dpth := pth + "dat4bgw/"
	ns01 := "http://purl.obolibrary.org/obo/"
	uri01 := ns01 + "bfo.owl"
	ns02 := "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:"
	uri02 := ns02 + "1672772"
	uris := [...]string{uri01, uri02}
	pths := [...]string{dpth + "WgetFile.01", dpth + "WgetFile.02"}

	tts := []tt{
		{uris[0], pths[0], nil},
		{uris[1], pths[1], nil},
	}
	for i, tt := range tts {
		err := WgetFile(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with WgetFile in", time.Since(mystart))
}

func Test_Rwget(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		arg3 string
		val1 error
	}

	pth := "../../tdata/"
	subdir := "dat4bgw/"
	dpth := pth + subdir
	ns01 := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping"
	nss := [...]string{ns01}
	mask01 := "README"
	masks := [...]string{mask01}

	tts := []tt{
		{nss[0], masks[0], dpth, nil},
	}
	for i, tt := range tts {
		err := Rwget(tt.arg1, tt.arg2, tt.arg3)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with Rwget in", time.Since(mystart))
}

func Test_Cleanup(t *testing.T) {
	pth := "../../tdata/"
	dpth := filepath.Join(pth, "dat4bgw/")
	if err := os.RemoveAll(dpth); err != nil {
		log.Println(err)
	}
}

