package main

import (
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strings"
)

func replaceN(pam string) []string {
	res := make([]string, 0)

	if strings.Contains(pam, "N") {
		for _, i := range []string{"A", "T", "C", "G"} {
			res = append(res, strings.Replace(pam, "N", i, -1))
		}
	} else {
		res = append(res, pam)
	}
	return res
}

// PoolLibrary handle all the spacer or targets
type PoolLibrary struct {
	Key      string                       `json:"key"`
	Scaffold string                       `json:"scaffold"`
	Spacer   map[string]map[string]string `json:"spacer"` // spacer and before
	Before   map[string]map[string]string `json:"before"` // before and list of pam
	Target   map[string]string            `json:"target"`
	Behind   map[string]map[string]string `json:"behind"` // before+pam+spacer with behinds
	Reverse  string
}

func LoadPoolLibrary(path string) map[string]*PoolLibrary {
	// open files
	f, err := os.Open(path)
	if err != nil {
		sugar.Fatal(err)
	}
	defer func() { _ = f.Close() }()

	// prepare readers and headers
	var gzReader *gzip.Reader
	var reader *csv.Reader
	if strings.HasSuffix(path, ".gz") {
		gzReader, err = gzip.NewReader(f)
		if err != nil {
			sugar.Fatal(err)
		}
		defer func() { _ = gzReader.Close() }()
		reader = csv.NewReader(gzReader)
	} else {
		reader = csv.NewReader(f)
	}

	headers, err := reader.Read()
	if err != nil {
		sugar.Fatal(err)
	}

	var libraries = make(map[string]*PoolLibrary)

	// optionally, resize scanner's capacity for lines over 64K, see next example
	for {
		temp, err := reader.Read()
		if err == io.EOF {
			break
		} else if err != nil {
			sugar.Fatal(err)
		}
		row := make(map[string]string)

		// 读完一行
		for idx, header := range headers {
			row[strings.ToLower(header)] = temp[idx]
		}

		if id_, ok := row["id"]; ok {
			row["guide"] = id_
		}

		key := strings.Split(row["guide"], "_")[0]
		key = strings.Split(key, "-")[0]

		_, ok := libraries[key]
		if !ok {
			libraries[key] = &PoolLibrary{
				Key:      key,
				Scaffold: row["scaffold"],
				Before:   make(map[string]map[string]string),
				Spacer:   make(map[string]map[string]string),
				Behind:   make(map[string]map[string]string),
				Target:   make(map[string]string),
				Reverse:  row["primer_R"],
			}
		}

		for _, pam := range replaceN(row["pam"]) {

			if before, ok := libraries[key].Before[row["before"]]; ok {
				before[pam] = ""
				libraries[key].Before[row["before"]] = before
			} else {
				libraries[key].Before[row["before"]] = map[string]string{pam: ""}
			}

			for _, spacer := range replaceN(row["spacer"]) {
				libraries[key].Spacer[spacer] = map[string]string{row["before"]: ""}
				libraries[key].Target[spacer] = row["target"]

				combined := fmt.Sprintf("%v%v%v", row["before"], pam, spacer)

				if behind, ok := libraries[key].Behind[combined]; ok {
					behind[row["behind"]] = ""
					libraries[key].Behind[combined] = behind
				} else {
					libraries[key].Behind[combined] = map[string]string{row["behind"]: ""}
				}

			}
		}
	}

	return libraries
}

// DesignSize 代表设计的质粒各元件的大小，及特定元件之间的距离
type DesignSize struct {
	TargetLen      int // gDNA target的长度
	PrimerLen      int // Reverse Primer的长度
	BehindLen      int // behind target length
	BeforeLen      int // before target length
	PAMLen         int // PAM的长度
	SpacerLen      int // spacer的长度
	SpacerDistance int // spacer与before target的间距
	PrimerRLen     int
	UMILen         int
}
