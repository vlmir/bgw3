package main

import (
	// "github.com/vlmir/bgw3/src/util"
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
	if err := os.MkdirAll(filepath.Join(dpth, "signor/"), 0755); err != nil {
		log.Println(err)
	}
}

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

/*
// takes too much time
func Test_getOneSignor(t *testing.T) {
	mystart := time.Now()
	type tt struct {
		arg1 string
		arg2 string
		val1 error
	}
	pth := "../../tdata/"
	txids := [...]string{"9606"}
	datdirs := [...]string{pth + "dat4bgw/"}
	tts := []tt{
		{txids[0], datdirs[0], nil},
	}
	for i, tt := range tts {
		err := saveOneSignor(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
	log.Println("Done with saveOneSignor in", time.Since(mystart))
}
*/

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

/*
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
	// subdir := "idmapping/"
	dpth := pth + subdir
//	ns01 := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping"
 ns01 := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
	ns02 := "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
	nss := [...]string{ns01, ns02}
	mask01 := "README"
	mask02 := "STATS.gz"
	masks := [...]string{mask01, mask02}

	tts := []tt{
		{nss[0], masks[0], dpth, nil},
		// {nss[1], masks[1], dpth, nil},
	}
	for i, tt := range tts {
		// saves nothing and fails: have exit status 8
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
*/
/*
mironov@nt-biogateway ~/r/b/t/dat4bgw> ll
total 172K
drwxr-xr-x 2 mironov mironov 4.0K Mar  9 01:32 goa/
drwxr-xr-x 2 mironov mironov 4.0K Mar  9 01:31 intact/
drwxr-xr-x 2 mironov mironov 4.0K Mar  9 01:31 signor/
drwxr-xr-x 2 mironov mironov 4.0K Mar  9 01:31 uniprot/
-rw-rw-r-- 1 mironov mironov 155K Mar  9 01:32 WgetFile.01
-rw-rw-r-- 1 mironov mironov    0 Mar  9 01:32 WgetFile.02

mironov@nt-biogateway ~/r/b/t/dat4bgw> ll goa/
total 6.2M
-rw-rw-r-- 1 mironov mironov 6.2M Mar  9 01:32 1672772.gaf
-rw-rw-r-- 1 mironov mironov  941 Mar  9 01:32 1672772.gpa
mironov@nt-biogateway ~/r/b/t/dat4bgw> ll intact/
total 352K
-rw-rw-r-- 1 mironov mironov    0 Mar  9 01:31 1672772.mi25
-rw-rw-r-- 1 mironov mironov 352K Mar  9 01:32 9031.mi25
mironov@nt-biogateway ~/r/b/t/dat4bgw> ll signor/
total 0
mironov@nt-biogateway ~/r/b/t/dat4bgw> ll uniprot/
total 0
-rw-rw-r-- 1 mironov mironov 0 Mar  9 01:31 1343077.upt
*/
func Test_Cleanup(t *testing.T) {
	pth := "../../tdata/"
	dpth := filepath.Join(pth, "dat4bgw/")
	if err := os.RemoveAll(dpth); err != nil {
		log.Println(err)
	}
}

/*
// 'panic: test timed out after 10m0s' after saving one of the 2 files
func Test_saveAllIdmap(t *testing.T) {
	type tt struct {
		arg1 string
		arg2 util.Set2D
		val1 error
	}
	txmap := make(util.Set2D)
	txmap.Add("185453", "test")
	txmap.Add("7739", "test")
	pth := "../../tdata/"
	tts := []tt{
		{pth, txmap, nil},
	}
	for i, tt := range tts {
		err := saveAllIdmap(tt.arg1, tt.arg2)
		if err != tt.val1 {
			t.Error(
				"For test", i+1, ": ",
				"\n\twant", tt.val1,
				"\n\thave", err,
			)
		}
	}
}
*/
