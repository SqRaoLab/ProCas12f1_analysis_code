package main

import (
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"sort"

	"github.com/schollz/progressbar/v3"
)

func createMergeProgressbar(files []string) *progressbar.ProgressBar {

	fileSize := int64(0)
	for _, file := range files {
		if stat, err := os.Stat(file); !os.IsNotExist(err) {
			fileSize += stat.Size()
		}
	}

	return progressBar(fileSize, true)
}

// 后续的文件融合策略
func merge(pattern, outPath string) {
	tempMatches, _ := filepath.Glob(pattern)

	matches := make([]string, 0)
	for _, i := range tempMatches {
		if i == outPath {
			continue
		}
		matches = append(matches, i)
	}

	sort.Strings(matches)

	if len(matches) < 2 {
		_ = os.Rename(matches[0], outPath)
	} else {
		bar := createMergeProgressbar(matches)

		outFile, err := os.Create(outPath)
		if err != nil {
			sugar.Fatal(err)
		}
		gz := gzip.NewWriter(outFile)

		for _, path := range matches {
			inFile, err := os.Open(path)
			if err != nil {
				sugar.Fatal(err)
			}

			pbarReader := &progressReader{inFile, bar}

			gzReader, err := gzip.NewReader(pbarReader)
			if err != nil {
				sugar.Fatal(err)
			}

			_, _ = io.Copy(gz, gzReader)

			_ = gzReader.Close()
			_ = inFile.Close()
		}
		_ = gz.Close()
		_ = outFile.Close()
		_ = bar.Finish()
	}
}
