package main

import (
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"os"
	"strings"
	"sync"

	"github.com/voxelbrain/goptions"
	"go.uber.org/zap"
)

var (
	sugar      *zap.SugaredLogger
	designSize *DesignSize
	poolLibMap map[string]*PoolLibrary
)

// reverseComplement 返回 DNA 序列的反向互补链
func reverseComplement(seq string) string {
	// 定义互补映射
	complement := map[rune]rune{
		'A': 'T',
		'T': 'A',
		'G': 'C',
		'C': 'G',
		'a': 't', // 可选：支持小写
		't': 'a',
		'g': 'c',
		'c': 'g',
	}

	var builder strings.Builder
	builder.Grow(len(seq)) // 预分配空间，提升性能

	// 从后往前遍历，实现反向 + 互补
	for i := len(seq) - 1; i >= 0; i-- {
		r := rune(seq[i])
		if comp, exists := complement[r]; exists {
			builder.WriteRune(comp)
		} else {
			// 可选：处理非法字符，如 'N' 或其他
			builder.WriteRune(r) // 或 builder.WriteRune('N') 表示未知
		}
	}

	return builder.String()
}

// processPairs 处理读段对
func processPairs(pairChan <-chan *ReadPair, resultChan chan<- []string, key string, umi string) {
	for pair := range pairChan {
		for _, res := range designSize.decodeString(pair.R1+reverseComplement(pair.R2), umi, poolLibMap[key]) {
			if umi != "" && res.UMI == "" {
				continue
			}
			row := []string{key, pair.ID}
			row = append(row, res.String(umi != "")...)
			resultChan <- row
			break
		}
	}
}

// --- 主函数 ---

func main() {
	options := defaultParams()
	goptions.ParseAndFail(options)
	setLogger(options.Debug)

	if options.Version {
		fmt.Println("Version: 0.2.2")
		os.Exit(0)
	}

	if options.R1 == "" || options.R2 == "" || options.Cas == "" {
		sugar.Errorf("-1, -2, -c all reqquired")
		goptions.PrintHelp()
		os.Exit(0)
	}

	designSize = &DesignSize{
		TargetLen:      options.TargetLen,
		PAMLen:         options.PAMLen,
		SpacerLen:      options.SpacerLen,
		SpacerDistance: options.SpacerDistance,
		BeforeLen:      options.BeforeLen,
		BehindLen:      options.BehindLen,
		PrimerLen:      options.PrimerLen,
		UMILen:         options.UMILen,
	}

	// 定义条形码集
	poolLibMap = LoadPoolLibrary(options.Pool)
	for key := range poolLibMap {
		sugar.Debugf("Library: %v loaded", key)
	}

	if _, ok := poolLibMap[options.Cas]; !ok {
		sugar.Fatalf("%s not found in library, please chech the input cas protein name", options.Cas)
	}

	//sugar.Infof("%v", poolLibMap[options.Cas])

	// 创建channel
	r1Chan := make(chan []string, 100) // 缓冲channel
	r2Chan := make(chan []string, 100)
	pairChan := make(chan *ReadPair, 1000)           // 更大的缓冲区以平衡生产者和消费者
	resultChan := make(chan []string, options.Procs) // 缓冲区大小等于worker数

	var wg sync.WaitGroup

	// 启动读取goroutine
	wg.Add(2)
	if options.ReadCmd != "" {
		go readFromCmd(options.R1, options.ReadCmd, r1Chan, &wg, true)
		go readFromCmd(options.R2, options.ReadCmd, r2Chan, &wg, false)
	} else {
		go readFastqInChannel(options.R1, r1Chan, &wg, true)
		go readFastqInChannel(options.R2, r2Chan, &wg, false)
	}

	// 启动配对goroutine
	go func() {
		defer close(pairChan)

		r1Cache := make(map[string]string)
		r2Cache := make(map[string]string)

		for {
			r1Record, ok1 := <-r1Chan
			r2Record, ok2 := <-r2Chan

			if !ok1 || !ok2 {
				break // 如果任一channel关闭，则停止
			}

			// 简单的ID匹配检查 (实际应解析@后的ID)
			// 这里假设顺序严格对应，且ID行格式为 @...
			if len(r1Record) > 0 && len(r2Record) > 0 {
				id1 := strings.SplitN(r1Record[0], " ", 2)[0] // 取第一个空格前的部分
				id2 := strings.SplitN(r2Record[0], " ", 2)[0]

				if strings.HasSuffix(id1, "/1") || strings.HasSuffix(id1, "/2") {
					id1 = strings.Split(id1, "/")[0]
					id2 = strings.Split(id2, "/")[0]
				}

				pair := &ReadPair{}

				if id1 == id2 {
					pair.ID = id1
					pair.R1 = r1Record[1]
					pair.R2 = r2Record[1]
				} else if _, ok := r1Cache[id2]; ok {
					pair.ID = id2
					pair.R1 = r1Cache[id2]
					pair.R2 = r2Record[1]
					delete(r1Cache, id2)
				} else if _, ok := r2Cache[id1]; ok {
					pair.ID = id1
					pair.R1 = r1Record[1]
					pair.R2 = r2Cache[id1]
					delete(r2Cache, id1)
				} else {
					r1Cache[id1] = r1Record[1]
					r2Cache[id2] = r2Record[1]
					continue
				}
				//sugar.Debugf("Pair: %v", pair)
				pairChan <- pair
			}
		}
	}()

	// 启动worker goroutine
	var workerWg sync.WaitGroup
	for i := 0; i < options.Procs; i++ {
		workerWg.Add(1)
		go func(id int) {
			defer workerWg.Done()
			processPairs(pairChan, resultChan, options.Cas, options.UMIAnchor)
		}(i)
	}

	// 等待读取完成
	go func() {
		wg.Wait()
	}()

	// 等待所有worker完成
	go func() {
		workerWg.Wait()
		close(resultChan) // 关闭结果channel
	}()

	// --- 写入结果 ---
	outFile, err := os.Create(options.Output)
	if err != nil {
		sugar.Fatalf("Error creating output file %s: %v", options.Output, err)
	}

	var writer *csv.Writer
	var gzFile *gzip.Writer
	if strings.HasSuffix(options.Output, ".gz") {
		gzFile = gzip.NewWriter(outFile)

		// 创建 TSV 写入器（使用 tab 作为分隔符）
		writer = csv.NewWriter(gzFile)
	} else {
		// 创建 TSV 写入器（使用 tab 作为分隔符）
		writer = csv.NewWriter(outFile)
	}

	// 写入表头 pam, spacer, behind, before
	temp := &SeqComponent{}
	header := []string{"CasProtein", "ReadID"}
	header = append(header, temp.Header(options.UMIAnchor != "")...)
	if err := writer.Write(header); err != nil {
		sugar.Fatalf("Error writing header: %v", err)
	}

	// do not flush manually, to prevent EOF error
	count := 0
	for res := range resultChan {
		_ = writer.Write(res)

		count++
		if count%10000 == 0 {
			writer.Flush()
		}
	}

	if err != nil {
		sugar.Infof("Error final flush to output file: %v", err)
	}

	// close file to prevent EOF error
	writer.Flush()
	if gzFile != nil {
		_ = gzFile.Close()
	}
	_ = outFile.Close()

	sugar.Infof("Results written to %s", options.Output)
}
