package main

import (
	"encoding/json"
	"log"
	"os"
	"path/filepath"
	"sort"
	"sync"

	"github.com/voxelbrain/goptions"
)

func loadInputFiles(path string) (map[string]int64, map[string][]string) {
	fqPatterns := map[string]string{}
	content, err := os.ReadFile(path)
	if err != nil {
		sugar.Fatal(err)
	}
	err = json.Unmarshal(content, &fqPatterns)
	if err != nil {
		sugar.Fatal(err)
	}

	// 收集所有文件对
	// var allPairs = make(map[string][]string)
	fileSize := make(map[string]int64)
	allPairs := make(map[string][]string)
	for group, pattern := range fqPatterns {
		matches, _ := filepath.Glob(pattern)
		sort.Strings(matches)

		if len(matches) < 1 {
			sugar.Fatalf("%v not found", pattern)
		}

		allPairs[group] = matches

		fileSize[group] = 0
		for _, f := range matches {
			if stat, ok := os.Stat(f); !os.IsNotExist(ok) {
				fileSize[group] += stat.Size()
			}
		}
	}
	return fileSize, allPairs
}

func main() {
	options := defaultParams()
	goptions.ParseAndFail(options)
	setLogger(filepath.Join(options.Output, "splitFq.log"))

	if options.Version {
		sugar.Info("Version: 0.2.6")
		os.Exit(0)
	}

	// 创建输出目录
	if err := os.MkdirAll(options.Output, 0755); err != nil {
		log.Fatal(err)
	}

	// 加载 barcode 规则
	rules, samples, err := loadBarcodes(options.Barcode)
	if err != nil {
		sugar.Fatal("Load barcodes error: ", err)
	}

	if options.Debug {
		sugar.Infof("%v", rules)
	}

	// 定义 FASTQ 路径模式
	fileSize, allPairs := loadInputFiles(options.Pool)

	// 统计总共要读取的文件大小
	total := int64(0)
	for key := range rules {
		fs, ok := fileSize[key]
		if ok {
			total += fs
		}
	}

	if options.Debug {
		sugar.Infof("%v", allPairs)
	}

	// read processes
	readChan := make(chan *ReadPair, 50000)
	var workerWg sync.WaitGroup
	var readerWg sync.WaitGroup
	var writerWg sync.WaitGroup

	if options.Merge {
		sugar.Infof("running in merge mode")
	}

	// create writers for different files
	writers := make(map[string]chan *ReadPair, len(samples))
	for sample := range samples {
		writerWg.Add(1)
		writers[sample] = make(chan *ReadPair, 50000)
		go writer(writers[sample], options.Output, &writerWg, options.Merge)
	}

	// create reader and processor
	bar := progressBar(total, true)
	for key, rule := range rules {
		fs, ok := allPairs[key]
		if ok {
			readerWg.Add(1)
			go readFastqInChannel(fs, readChan, &readerWg, bar)

			for i := 0; i < 6; i++ {
				workerWg.Add(1)
				go processFilePair(readChan, writers, rule, &workerWg, options.MinMeanQual, options.MinQ30Percent, options.Trim)
			}
		}
	}

	readerWg.Wait()
	close(readChan)

	workerWg.Wait()
	for _, channel := range writers {
		close(channel)
	}

	writerWg.Wait()
	_ = bar.Finish()
}
